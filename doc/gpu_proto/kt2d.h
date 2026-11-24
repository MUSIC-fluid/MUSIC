// Portable kernels for the ideal, boost-invariant (2+1)D KT update of MUSIC
// (face_reconstruction 0 and 1, Heun RK2, the velocity root solve).
//
// One source for three compilers:
//   C++ (host, OpenMP)  : the reference, in double or float
//   CUDA                : __device__ functions, called from gpu/proto/kt2d_cuda.cu
//   Metal Shading Lang. : prepended with the MSL prelude by kt2d_metal.mm
//
// It mirrors src/advance.cpp (computeFaceFlux, computeFaceFluxPrimitive,
// MakeDeltaQI, FirstRKStepT, MaxSpeed) and src/reconst.cpp
// (ReconstIt_velocity_Newton, solve_v_from_guess, solve_v_Hybrid) for
// u^eta = 0, n_Q = n_S = 0. MUSIC's u^tau refinement for v > v_critical
// is not ported (the double results agree with MUSIC to 5e-9 without it).
// Errors are returned as flags instead of messages and exit().

#ifndef GPU_PROTO_KT2D_H_
#define GPU_PROTO_KT2D_H_

#if defined(__METAL_VERSION__)
#include <metal_stdlib>
using namespace metal;
#define KT_DEV inline
#define KT_GLOBAL device
#define KT_CONST constant
#define KT_THREAD thread
#elif defined(__CUDACC__)
#define KT_DEV __device__ inline
#define KT_GLOBAL
#define KT_CONST
#define KT_THREAD
#else
#include <cmath>
#define KT_DEV inline
#define KT_GLOBAL
#define KT_CONST
#define KT_THREAD
#endif

namespace kt {

// error flags, OR-ed atomically into one word
#if defined(__METAL_VERSION__)
#define KT_FLAGS device atomic_int *
inline void flagOr(KT_FLAGS f, int v) {
    atomic_fetch_or_explicit(f, v, memory_order_relaxed);
}
#elif defined(__CUDACC__)
#define KT_FLAGS int *
__device__ inline void flagOr(KT_FLAGS f, int v) { atomicOr(f, v); }
#else
#define KT_FLAGS int *
inline void flagOr(KT_FLAGS f, int v) { __atomic_fetch_or(f, v, __ATOMIC_RELAXED); }
#endif

#if !defined(__METAL_VERSION__) && !defined(__CUDACC__)
using std::fabs;
using std::fmax;
using std::fmin;
using std::isfinite;
using std::sqrt;
#endif

// equation of state: 0 = conformal P = e/3, 1 = table P(e) on a uniform grid
// in ln e (linear interpolation in e between the nodes)
template <typename Real>
struct Params {
    int nx, ny;
    int mode;   // face_reconstruction: 0 = limit tau T^{tau mu}, 1 = e, n, u^i
    int eosType;
    int eosN;   // table nodes
    int rkFlag;
    int storeQ;  // 1: the state also holds q = tau T^{tau mu}, tau J^tau
                 // (components 5-8), and the update starts from it
    Real dx, dy;
    Real dtau;
    Real tau;   // tau_n of the step
    Real theta;  // minmod theta
    Real absErr, relErr;
    Real eosLnE0, eosDlnE;  // ln e of node 0 and the node spacing
};

template <typename Real>
struct EosState {
    Real p, dpde, cs2;
};

template <typename Real>
KT_DEV EosState<Real> eosEval(
    KT_CONST const Params<Real> &P, KT_GLOBAL const Real *table, Real e) {
    EosState<Real> s;
    if (P.eosType == 0) {
        s.dpde = Real(1) / Real(3);
        s.p = s.dpde * e;
        s.cs2 = s.dpde;
        return s;
    }
    // table[2 i] = e_i, table[2 i + 1] = P_i
    const Real x = (log(fmax(e, Real(1e-30))) - P.eosLnE0) / P.eosDlnE;
    int i = int(x);
    i = (i < 0 ? 0 : (i > P.eosN - 2 ? P.eosN - 2 : i));
    const Real e0 = table[2 * i], p0 = table[2 * i + 1];
    const Real e1 = table[2 * i + 2], p1 = table[2 * i + 3];
    s.dpde = (p1 - p0) / (e1 - e0);
    s.p = p0 + s.dpde * (e - e0);
    if (e < e0) s.p = p0 * e / e0;  // below the table: P proportional to e
    s.cs2 = fmin(fmax(s.dpde, Real(0.01)), Real(1) / Real(3));
    return s;
}

template <typename Real>
struct Cell {
    Real e, nb, u0, ux, uy;
};

template <typename Real>
KT_DEV Real minmodDx(Real up1, Real u, Real um1, Real theta) {
    const Real diffup = (up1 - u) * theta;
    const Real diffdown = (u - um1) * theta;
    const Real diffmid = (up1 - um1) * Real(0.5);
    if (diffup == Real(0)) return Real(0);
    return diffup
           * fmax(Real(0), fmin(Real(1), fmin(diffdown / diffup,
                                              diffmid / diffup)));
}

// tau T^{tau mu}, tau J^tau components 0 (tt), 1 (tx), 2 (ty), 3 (J)
// and the flux components T^{d mu}, J^d along direction d = 1 (x), 2 (y)
template <typename Real>
KT_DEV Real TJ(const KT_THREAD Cell<Real> &c, int alpha, int nu, Real p) {
    const Real unu = (nu == 0 ? c.u0 : (nu == 1 ? c.ux : c.uy));
    if (alpha == 3) return c.nb * unu;
    const Real umu = (alpha == 0 ? c.u0 : (alpha == 1 ? c.ux : c.uy));
    Real g = Real(0);
    if (alpha == nu) g = (alpha == 0 ? Real(-1) : Real(1));
    return (c.e + p) * umu * unu + p * g;
}

template <typename Real>
KT_DEV Real maxSpeed(Real tau, int direc, const KT_THREAD Cell<Real> &c,
                     const KT_THREAD EosState<Real> &s, KT_THREAD int &flag) {
    const Real utau = c.u0;
    const Real utau2 = utau * utau;
    const Real ux = fabs(direc == 1 ? c.ux : c.uy);
    const Real ut2mux2 = utau2 - ux * ux;
    const Real vs2 = s.cs2;
    const Real numSqrt = (ut2mux2 - (ut2mux2 - Real(1)) * vs2) * vs2;
    Real num;
    if (numSqrt >= Real(0)) {
        num = utau * ux * (Real(1) - vs2) + sqrt(numSqrt);
    } else {
        const Real h = s.p + c.e;
        num = sqrt(fmax(Real(0), -(h * s.dpde * h
                                   * (s.dpde * (Real(-1) + ut2mux2) - ut2mux2))))
              - h * (Real(-1) + s.dpde) * utau * ux;
        if (s.dpde >= Real(0.001)) flag |= 1;
    }
    const Real den = utau2 * (Real(1) - vs2) + vs2;
    Real f = num / fmax(den, Real(1e-30));
    if (f > Real(1)) {
        flag |= 2;
        f = Real(1);
    }
    if (f < ux / utau) f = ux / utau;
    return f;
}

// ----- velocity solve: src/reconst.cpp ------------------------------------
template <typename Real>
KT_DEV Real velocityF(KT_CONST const Params<Real> &P,
                      KT_GLOBAL const Real *table, Real v, Real T00, Real M,
                      Real J0) {
    const Real e = T00 - v * M;
    const EosState<Real> s = eosEval(P, table, e);
    return v - M / (T00 + s.p);
}

template <typename Real>
KT_DEV int solveVFromGuess(KT_CONST const Params<Real> &P,
                           KT_GLOBAL const Real *table, Real vGuess, Real T00,
                           Real M, Real J0, KT_THREAD Real &vSol) {
    if (!(vGuess >= Real(0) && vGuess < Real(1))) return 0;
    Real v = vGuess;
    Real fv = velocityF(P, table, v, T00, M, J0);
    Real dfdv = Real(1) - Real(0.2) * v * v;
    if (!isfinite(fv)) return 0;
    if (fabs(fv) < P.absErr) {
        vSol = v;
        return 1;
    }
    for (int iter = 0; iter < 8; iter++) {
        if (!isfinite(dfdv) || dfdv == Real(0)) return 0;
        const Real dv = fv / dfdv;
        const Real vNext = v - dv;
        if (!(vNext >= Real(0) && vNext < Real(1))) return 0;
        const Real absE = fabs(dv);
        if (absE <= P.absErr || absE / (vNext + P.absErr) <= P.relErr) {
            vSol = vNext;
            return 1;
        }
        const Real fNext = velocityF(P, table, vNext, T00, M, J0);
        if (!isfinite(fNext)) return 0;
        dfdv = (fNext - fv) / (vNext - v + P.absErr);
        v = vNext;
        fv = fNext;
    }
    return 0;
}

template <typename Real>
KT_DEV void velocityFdf(KT_CONST const Params<Real> &P,
                        KT_GLOBAL const Real *table, Real v, Real T00, Real M,
                        Real J0, KT_THREAD Real &fv, KT_THREAD Real &dfdv) {
    fv = velocityF(P, table, v, T00, M, J0);
    const Real h = fmax(Real(1e-3) * v, Real(1e-8));
    const Real v1 = (v > h) ? v - h : v + h;
    const Real fv1 = velocityF(P, table, v1, T00, M, J0);
    dfdv = (fv - fv1) / (v - v1);
}

template <typename Real>
KT_DEV int solveVHybrid(KT_CONST const Params<Real> &P,
                        KT_GLOBAL const Real *table, Real vGuess, Real T00,
                        Real M, Real J0, KT_THREAD Real &vSol) {
    Real vl = fmax(Real(0), vGuess - Real(0.05));
    Real vh = fmin(Real(1), vGuess + Real(0.05));
    Real fl = velocityF(P, table, vl, T00, M, J0);
    Real fh = velocityF(P, table, vh, T00, M, J0);
    if (fl * fh > Real(0)) {
        vl = Real(0);
        vh = Real(1);
        fl = velocityF(P, table, vl, T00, M, J0);
        fh = velocityF(P, table, vh, T00, M, J0);
        if (fl * fh > Real(0)) return 0;
    }
    if (fabs(fl) < P.absErr) {
        vSol = vl;
        return 1;
    }
    if (fabs(fh) < P.absErr) {
        vSol = vh;
        return 1;
    }
    Real dvPrev = vh - vl;
    Real dvCurr = dvPrev;
    Real vr = (vh + vl) / Real(2);
    Real fv = Real(0), dfdv = Real(0);
    if (dvPrev > Real(0.1) * vr) {
        velocityFdf(P, table, vr, T00, M, J0, fv, dfdv);
    } else {
        fv = velocityF(P, table, vr, T00, M, J0);
        dfdv = (fh - fl) / (vh - vl);
    }
    Real absE = Real(10), relE = Real(10);
    int iter = 0;
    do {
        iter++;
        const Real vPrev = vr;
        const Real fPrev = fv;
        if (((vr - vh) * dfdv - fv) * ((vr - vl) * dfdv - fv) > Real(0)
            || (fabs(Real(2) * fv) > fabs(dvPrev * dfdv))) {
            dvPrev = dvCurr;
            dvCurr = (vh - vl) / Real(2);
            vr = vl + dvCurr;
        } else {
            dvPrev = dvCurr;
            dvCurr = fv / dfdv;
            vr = vr - dvCurr;
        }
        absE = fabs(dvCurr);
        relE = absE / (vr + P.absErr);
        if (fabs(vr - vPrev) > Real(0.1) * vr) {
            velocityFdf(P, table, vr, T00, M, J0, fv, dfdv);
        } else {
            fv = velocityF(P, table, vr, T00, M, J0);
            dfdv = (fv - fPrev) / (vr - vPrev + P.absErr);
        }
        if (fv * fl < Real(0)) {
            vh = vr;
            fh = fv;
        } else {
            vl = vr;
            fl = fv;
        }
        if (iter > 100) return 0;
    } while (absE > P.absErr && relE > P.relErr);
    vSol = vr;
    return 1;
}

// q = tau T^{tau mu}, tau J^tau at time tau. guess: the cell on the same side
// (its state is the fallback). flag bits: 4 revert, 8 regulate, 16 capped
template <typename Real>
KT_DEV Cell<Real> reconstruct(KT_CONST const Params<Real> &P,
                              KT_GLOBAL const Real *table, Real tau,
                              const KT_THREAD Real *q,
                              const KT_THREAD Cell<Real> &guess, Real pGuess,
                              KT_THREAD int &flag) {
    Real T00 = q[0] / tau, Mx = q[1] / tau, My = q[2] / tau;
    const Real J0 = q[3] / tau;
    Real M = sqrt(Mx * Mx + My * My);
    Cell<Real> c;
    if (T00 >= P.absErr && M >= T00) {  // superluminal: cap |M| at 0.999 T00
        const Real scale = Real(0.999) * T00 / M;
        Mx *= scale;
        My *= scale;
        M = Real(0.999) * T00;
        flag |= 16;
    }
    if (T00 < P.absErr) {  // regulate: keep n_B = J^tau
        c.e = fmax(P.absErr, T00);
        c.nb = J0;
        c.u0 = Real(1);
        c.ux = Real(0);
        c.uy = Real(0);
        flag |= 8;
        return c;
    }
    Real vGuess = sqrt(Real(1) - Real(1) / (guess.u0 * guess.u0 + P.absErr));
    if (pGuess > Real(0)) vGuess = M / (T00 + pGuess);
    Real v = Real(0);
    int ok = solveVFromGuess(P, table, vGuess, T00, M, J0, v);
    if (!ok) ok = solveVHybrid(P, table, vGuess, T00, M, J0, v);
    if (!ok) {
        flag |= 4;
        return guess;
    }
    const Real u0 = Real(1) / (sqrt(Real(1) - v * v) + v * P.absErr);
    if (fabs(u0 - guess.u0) / guess.u0 > Real(100)) {
        flag |= 4;
        return guess;
    }
    c.e = T00 - v * M;
    c.nb = J0 / u0;
    const EosState<Real> s = eosEval(P, table, c.e);
    const Real fac = u0 / (T00 + s.p);
    c.u0 = u0;
    c.ux = Mx * fac;
    c.uy = My * fac;
    const Real u2 = c.ux * c.ux + c.uy * c.uy;
    // u0^2 - 1 = (v u0)^2: written without the cancellation of u0^2 - 1,
    // which in float loses all digits at small v (MUSIC: u0 * u0 - 1)
    const Real vu0sq = v * v * u0 * u0;
    if (fabs(vu0sq - u2) > P.absErr) {
        const Real sf = sqrt(vu0sq / (u2 + P.absErr));
        c.ux *= sf;
        c.uy *= sf;
    }
    return c;
}

// ----- state access ---------------------------------------------------------
template <typename Real>
KT_DEV Cell<Real> loadCell(KT_GLOBAL const Real *s, int n, int ix, int iy,
                           int nx, int ny) {
    ix = (ix < 0 ? 0 : (ix > nx - 1 ? nx - 1 : ix));  // zero-gradient edges
    iy = (iy < 0 ? 0 : (iy > ny - 1 ? ny - 1 : iy));
    const int i = ix * ny + iy;
    Cell<Real> c;
    c.e = s[i];
    c.nb = s[n + i];
    c.u0 = s[2 * n + i];
    c.ux = s[3 * n + i];
    c.uy = s[4 * n + i];
    return c;
}

template <typename Real>
KT_DEV void storeCell(KT_GLOBAL Real *s, int n, int i,
                      const KT_THREAD Cell<Real> &c) {
    s[i] = c.e;
    s[n + i] = c.nb;
    s[2 * n + i] = c.u0;
    s[3 * n + i] = c.ux;
    s[4 * n + i] = c.uy;
}

// ----- KT flux through one face --------------------------------------------
// direction 1: face (fx, fy) between cells fx - 1 and fx along x, fx in
// [0, nx]; direction 2 likewise along y. flux holds 4 components per face,
// stored as flux[alpha * nFaces + face]
template <typename Real>
KT_DEV void faceFlux(KT_CONST const Params<Real> &P,
                     KT_GLOBAL const Real *table, KT_GLOBAL const Real *s,
                     int direction, int fx, int fy, KT_GLOBAL Real *flux,
                     KT_FLAGS flags) {
    const int nx = P.nx, ny = P.ny, n = nx * ny;
    const Real tau = P.tau + Real(P.rkFlag) * P.dtau;
    const int dx = (direction == 1), dy = (direction == 2);
    const Cell<Real> cLL = loadCell(s, n, fx - 2 * dx, fy - 2 * dy, nx, ny);
    const Cell<Real> cL = loadCell(s, n, fx - dx, fy - dy, nx, ny);
    const Cell<Real> cR = loadCell(s, n, fx, fy, nx, ny);
    const Cell<Real> cRR = loadCell(s, n, fx + dx, fy + dy, nx, ny);
    int flag = 0;
    Cell<Real> gL, gR;
    Real qL[4], qR[4];
    EosState<Real> sL, sR;
    if (P.mode == 1) {
        const Real th = P.theta;
        gL.e = cL.e + Real(0.5) * minmodDx(cR.e, cL.e, cLL.e, th);
        gR.e = cR.e - Real(0.5) * minmodDx(cRR.e, cR.e, cL.e, th);
        gL.nb = cL.nb + Real(0.5) * minmodDx(cR.nb, cL.nb, cLL.nb, th);
        gR.nb = cR.nb - Real(0.5) * minmodDx(cRR.nb, cR.nb, cL.nb, th);
        gL.ux = cL.ux + Real(0.5) * minmodDx(cR.ux, cL.ux, cLL.ux, th);
        gR.ux = cR.ux - Real(0.5) * minmodDx(cRR.ux, cR.ux, cL.ux, th);
        gL.uy = cL.uy + Real(0.5) * minmodDx(cR.uy, cL.uy, cLL.uy, th);
        gR.uy = cR.uy - Real(0.5) * minmodDx(cRR.uy, cR.uy, cL.uy, th);
        gL.u0 = sqrt(Real(1) + gL.ux * gL.ux + gL.uy * gL.uy);
        gR.u0 = sqrt(Real(1) + gR.ux * gR.ux + gR.uy * gR.uy);
        sL = eosEval(P, table, gL.e);
        sR = eosEval(P, table, gR.e);
        for (int a = 0; a < 4; a++) {
            qL[a] = tau * TJ(gL, a, 0, sL.p);
            qR[a] = tau * TJ(gR, a, 0, sR.p);
        }
    } else {
        const Real pLL = eosEval(P, table, cLL.e).p;
        const Real pL = eosEval(P, table, cL.e).p;
        const Real pR = eosEval(P, table, cR.e).p;
        const Real pRR = eosEval(P, table, cRR.e).p;
        for (int a = 0; a < 4; a++) {
            const Real g0 = tau * TJ(cLL, a, 0, pLL);
            const Real g1 = tau * TJ(cL, a, 0, pL);
            const Real g2 = tau * TJ(cR, a, 0, pR);
            const Real g3 = tau * TJ(cRR, a, 0, pRR);
            qL[a] = g1 + Real(0.5) * minmodDx(g2, g1, g0, P.theta);
            qR[a] = g2 - Real(0.5) * minmodDx(g3, g2, g1, P.theta);
        }
        gL = reconstruct(P, table, tau, qL, cL, pL, flag);
        gR = reconstruct(P, table, tau, qR, cR, pR, flag);
        sL = eosEval(P, table, gL.e);
        sR = eosEval(P, table, gR.e);
    }
    const Real a = fmax(maxSpeed(tau, direction, gL, sL, flag),
                        maxSpeed(tau, direction, gR, sR, flag));
    const int nfx = nx + dx, nfy = ny + dy;
    const int face = fx * nfy + fy;
    const int nFaces = nfx * nfy;
    for (int al = 0; al < 4; al++) {
        const Real FL = tau * TJ(gL, al, direction, sL.p);
        const Real FR = tau * TJ(gR, al, direction, sR.p);
        flux[al * nFaces + face] =
            Real(0.5) * ((FL + FR) - a * (qR[al] - qL[al]));
    }
    if (flag) flagOr(flags, flag);
}

// ----- the cell update of one RK stage (FirstRKStepT) ---------------------
// sCurr: the stage's state (tau_n or the predictor at tau_{n+1}),
// sPrev: the state at tau_n (used in the second stage), sNext: output
template <typename Real>
KT_DEV void cellUpdate(KT_CONST const Params<Real> &P,
                       KT_GLOBAL const Real *table, KT_GLOBAL const Real *sCurr,
                       KT_GLOBAL const Real *sPrev, KT_GLOBAL const Real *fluxX,
                       KT_GLOBAL const Real *fluxY, KT_GLOBAL Real *sNext,
                       int ix, int iy, KT_FLAGS flags) {
    const int nx = P.nx, ny = P.ny, n = nx * ny;
    const Real tauRk = P.tau + Real(P.rkFlag) * P.dtau;
    const Cell<Real> c = loadCell(sCurr, n, ix, iy, nx, ny);
    const Real p = eosEval(P, table, c.e).p;
    Real qi[4];
    const int i0 = ix * ny + iy;
    for (int a = 0; a < 4; a++)
        qi[a] = (P.storeQ ? sCurr[(5 + a) * n + i0] : TJ(c, a, 0, p) * tauRk);
    const int nFx = (nx + 1) * ny, nFy = nx * (ny + 1);
    for (int a = 0; a < 4; a++) {
        const Real Fxm = fluxX[a * nFx + ix * ny + iy];
        const Real Fxp = fluxX[a * nFx + (ix + 1) * ny + iy];
        const Real Fym = fluxY[a * nFy + ix * (ny + 1) + iy];
        const Real Fyp = fluxY[a * nFy + ix * (ny + 1) + iy + 1];
        qi[a] += ((Fxm - Fxp) / P.dx + (Fym - Fyp) / P.dy) * P.dtau;
    }
    // the eta faces of a boost-invariant grid: -P (tau^2 T^{eta eta})
    qi[0] -= p * P.dtau;
    if (P.rkFlag == 1) {
        const Cell<Real> c0 = loadCell(sPrev, n, ix, iy, nx, ny);
        const Real p0 = eosEval(P, table, c0.e).p;
        for (int a = 0; a < 4; a++) {
            const Real q0 =
                (P.storeQ ? sPrev[(5 + a) * n + i0] : TJ(c0, a, 0, p0) * P.tau);
            qi[a] = Real(0.5) * (qi[a] + q0);
        }
    }
    int flag = 0;
    const Cell<Real> out =
        reconstruct(P, table, P.tau + P.dtau, qi, c, p, flag);
    storeCell(sNext, n, i0, out);
    if (P.storeQ) {
        // a capped, regulated or reverted cell keeps the q of its state
        const Real tauN = P.tau + P.dtau;
        const Real pOut = eosEval(P, table, out.e).p;
        for (int a = 0; a < 4; a++)
            sNext[(5 + a) * n + i0] = (flag ? TJ(out, a, 0, pOut) * tauN : qi[a]);
    }
    if (flag) flagOr(flags, flag);
}

}  // namespace kt

#endif  // GPU_PROTO_KT2D_H_
