// Copyright 2026 Chun Shen
//
// The ideal part of one RK stage on the device: the KT face fluxes
// (Advance::computeFaceFlux, computeFaceFluxPrimitive, MaxSpeed), the cell
// update (Advance::FirstRKStepT with MakeDeltaQI, without sources and
// viscous terms) and the root solve of Reconst (ReconstIt_shell,
// ReconstIt_velocity_Newton, solve_v_from_guess, solve_v_Hybrid,
// solve_u0_Hybrid). Each function follows its host version operation by
// operation, so that in double the device reproduces MUSIC; where MUSIC
// prints a message and exits, the kernels set a flag (device_defs.h).
//
// Data layout: a state slice holds kNComp components of n = nx ny neta
// cells, component c of cell (ix, iy, ieta) at
//   c n + ix + nx (iy + ny ieta);
// a neighbour past the grid edge is the edge cell (MUSIC's ghost cells).
// The faces along direction d (1 = x, 2 = y, 3 = eta) are indexed like the
// cells on a grid with one more point along d; face f lies between cells
// f - 1 and f, and component alpha of its flux is at alpha nFaces + face.

#ifndef SRC_DEVICE_DEVICE_IDEAL_H_
#define SRC_DEVICE_DEVICE_IDEAL_H_

#ifndef MUSIC_DEVICE_BUNDLED
#include "device_defs.h"
#include "device_eos.h"
#endif

namespace music_device {

template <typename Real>
struct StepParams {
    int nx, ny, neta;
    int mode;    // face_reconstruction
    int rkFlag;  // 0: first stage, 1: second
    int storeQ;  // the update starts from the stored tau T^{tau mu}
    int maxIter;
    Real tau;  // tau_n of the step
    Real dtau, dx, dy, deta;
    Real theta;  // minmod theta
    Real absErr, relErr, vCritical;
    Real coshDeta, sinhDeta;  // the eta terms of MakeDeltaQI
    int nComp;  // components of a state slice: kNComp, or kNCompVisc
    // ----- the viscous update (device_visc.h), viscosity_flag = 1 -----
    int shear, bulk;    // turn_on_shear, turn_on_bulk
    int secondOrder;    // include_second_order_terms
    int vorticity;      // include_vorticity_terms
    int duIdeal;        // Du_from_ideal_EoM
    int resum;          // FlagResumTransportCoeff
    int muBShear;       // muB_dependent_shear_to_s
    int bulkRelaxType;  // bulk_relaxation_type
    int questRevert;    // Initial_profile is neither 0 nor 1
    int shearT, bulkT;  // T_dependent_shear_to_s, T_dependent_bulk_to_s
    int nWComp;         // Diss::nWComponents
    // DATA.gmunu, read at run time as on the host (see device_visc.h)
    Real g[4][4];
    Real duIdealEDilute, resumAlpha, questRevertStrength;
    Real shearRelaxFactor, bulkRelaxFactor;
    Real shearToS, shear2Min, shear2Slope, shear2Curv;
    Real shear3TKink, shear3LowSlope, shear3HighSlope, shear3AtKink;
    Real shearMuBf0p2, shearMuBf0p4, shearMuBAlpha, shearMuBSlope,
        shearMuBScale;
    Real bulk2Norm, bulk2Width, bulk2Peak;
    Real bulk3Max, bulk3Width, bulk3TPeak, bulk3Lambda;
    Real bulk10Max, bulk10MaxMuB0p2, bulk10MaxMuB0p4, bulk10Tpeak,
        bulk10TpeakMuBCurv, bulk10WidthLow, bulk10WidthHigh;
    // the stage adds the source terms qi_source (kNTJ components of the
    // cells, computed by the host as in Advance::FirstRKStepT)
    int hasSource;
};

template <typename Real>
struct Cell {
    Real e, rhob, rhoq, rhos;
    Real u[4];
};

template <typename Real>
struct EosState {
    EosGrad<Real> g;
    bool valid;
};

MD_FN int cellIndex(int nx, int ny, int neta, int ix, int iy, int ieta) {
    ix = dmax(0, dmin(nx - 1, ix));
    iy = dmax(0, dmin(ny - 1, iy));
    ieta = dmax(0, dmin(neta - 1, ieta));
    return ix + nx * (iy + ny * ieta);
}

template <typename Real>
MD_FN Cell<Real> loadCell(MD_GLOBAL const Real *s, int n, int i) {
    Cell<Real> c;
    c.e = s[i];
    c.rhob = s[n + i];
    c.rhoq = s[2 * n + i];
    c.rhos = s[3 * n + i];
    for (int k = 0; k < 4; k++) c.u[k] = s[(4 + k) * n + i];
    return c;
}

template <typename Real>
MD_FN void storeCell(
    MD_GLOBAL Real *s, int n, int i, MD_THREAD const Cell<Real> &c) {
    s[i] = c.e;
    s[n + i] = c.rhob;
    s[2 * n + i] = c.rhoq;
    s[3 * n + i] = c.rhos;
    for (int k = 0; k < 4; k++) s[(4 + k) * n + i] = c.u[k];
}

//! Advance::get_TJb
template <typename Real>
MD_FN Real TJb(MD_THREAD const Cell<Real> &c, int mu, int nu, Real p) {
    const Real u_nu = c.u[nu];
    if (mu == 4) return c.rhob * u_nu;
    if (mu == 5) return c.rhoq * u_nu;
    if (mu == 6) return c.rhos * u_nu;
    const Real e = c.e;
    Real gfac = Real(0);
    Real u_mu = Real(0);
    if (mu == nu) {
        u_mu = u_nu;
        gfac = Real(1);
        if (mu == 0) gfac = Real(-1);
    } else {
        u_mu = c.u[mu];
    }
    return (e + p) * u_mu * u_nu + p * gfac;
}

//! Minmod::minmod_dx
template <typename Real>
MD_FN Real minmodDx(Real up1, Real u, Real um1, Real theta) {
    const Real diffup = (up1 - u) * theta;
    const Real diffdown = (u - um1) * theta;
    const Real diffmid = (up1 - um1) * Real(0.5);
    if (diffup == Real(0)) return Real(0);
    return diffup
           * dmax(Real(0),
                  dmin(Real(1), dmin(diffdown / diffup, diffmid / diffup)));
}

// ----- the root solve: Reconst ---------------------------------------------

template <typename Real>
struct Solver {
    // pointers: MSL has no reference members
    MD_CONST const StepParams<Real> *P_;
    MD_CONST const EosDesc<Real> *d_;
    MD_GLOBAL const Real *tab;
    // J^0 of the electric charge and the strangeness of the cell solved
    // (set in velocityNewton): the 4D EoS depends on n_Q, n_S
    Real J0Q_, J0S_;

    //! the pressure at e, n_B = rhob and n_{Q,S} = J0{Q,S} f, as the
    //! solves of Reconst compute them
    MD_FN Real pressure(Real e, Real rhob, Real f) const {
        return eosPressure(*d_, tab, e, rhob, J0Q_ * f, J0S_ * f);
    }

    // Reconst::reconst_velocity_f
    MD_FN Real velocityF(Real v, Real T00, Real M, Real J0B) const {
        const Real epsilon = T00 - v * M;
        const Real temp = sqrt(Real(1) - v * v);
        const Real rhob = J0B * temp;
        const Real p = pressure(epsilon, rhob, temp);
        return v - M / (T00 + p);
    }

    // Reconst::reconst_velocity_fdf
    MD_FN void velocityFdf(
        Real v, Real T00, Real M, Real J0B, MD_THREAD Real &fv,
        MD_THREAD Real &dfdv) const {
        fv = velocityF(v, T00, M, J0B);
        const Real h = dmax(Real(1e-3) * v, Real(1e-8));
        const Real v1 = (v > h) ? v - h : v + h;
        const Real fv1 = velocityF(v1, T00, M, J0B);
        dfdv = (fv - fv1) / (v - v1);
    }

    // Reconst::solve_v_from_guess
    MD_FN int solveVFromGuess(
        Real vGuess, Real T00, Real M, Real J0B, MD_THREAD Real &vSol) const {
        if (!(vGuess >= Real(0) && vGuess < Real(1))) return 0;
        Real v = vGuess;
        Real fv = velocityF(v, T00, M, J0B);
        Real dfdv = Real(1) - Real(0.2) * v * v;
        if (!isfinite(fv)) return 0;
        if (fabs(fv) < P_->absErr) {
            vSol = v;
            return 1;
        }
        for (int iter = 0; iter < 8; iter++) {
            if (!isfinite(dfdv) || dfdv == Real(0)) return 0;
            const Real dv = fv / dfdv;
            const Real vNext = v - dv;
            if (!(vNext >= Real(0) && vNext < Real(1))) return 0;
            const Real absE = fabs(dv);
            const Real relE = absE / (vNext + P_->absErr);
            if (absE <= P_->absErr || relE <= P_->relErr) {
                vSol = vNext;
                return 1;
            }
            const Real fNext = velocityF(vNext, T00, M, J0B);
            if (!isfinite(fNext)) return 0;
            dfdv = (fNext - fv) / (vNext - v + P_->absErr);
            v = vNext;
            fv = fNext;
        }
        return 0;
    }

    // Reconst::solve_v_Hybrid
    MD_FN int solveVHybrid(
        Real vGuess, Real T00, Real M, Real J0B, MD_THREAD Real &vSol,
        MD_THREAD int &flag) const {
        Real vl = dmax(Real(0), vGuess - Real(0.05));
        Real vh = dmin(Real(1), vGuess + Real(0.05));
        Real fl = velocityF(vl, T00, M, J0B);
        Real fh = velocityF(vh, T00, M, J0B);
        if (fl * fh > Real(0)) {
            vl = Real(0);
            vh = Real(1);
            fl = velocityF(vl, T00, M, J0B);
            fh = velocityF(vh, T00, M, J0B);
            if (fl * fh > Real(0)) {
                flag |= kFlagHybridFail;
                return 0;
            }
        }
        if (fabs(fl) < P_->absErr) {
            vSol = vl;
            return 1;
        }
        if (fabs(fh) < P_->absErr) {
            vSol = vh;
            return 1;
        }
        Real dvPrev = vh - vl;
        Real dvCurr = dvPrev;
        Real vr = (vh + vl) / Real(2);
        Real fv = Real(0);
        Real dfdv = Real(0);
        if (dvPrev > Real(0.1) * vr) {
            velocityFdf(vr, T00, M, J0B, fv, dfdv);
        } else {
            fv = velocityF(vr, T00, M, J0B);
            dfdv = (fh - fl) / (vh - vl);
        }
        Real absE = Real(10);
        Real relE = Real(10);
        int iter = 0;
        int status = 1;
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
            relE = absE / (vr + P_->absErr);
            if (fabs(vr - vPrev) > Real(0.1) * vr) {
                velocityFdf(vr, T00, M, J0B, fv, dfdv);
            } else {
                fv = velocityF(vr, T00, M, J0B);
                dfdv = (fv - fPrev) / (vr - vPrev + P_->absErr);
            }
            if (fv * fl < Real(0)) {
                vh = vr;
                fh = fv;
            } else {
                vl = vr;
                fl = fv;
            }
            if (iter > P_->maxIter) {
                status = 0;
                break;
            }
        } while (absE > P_->absErr && relE > P_->relErr);
        vSol = vr;
        return status;
    }

    // Reconst::reconst_u0_f
    MD_FN Real u0F(Real u0, Real T00, Real K00, Real M, Real J0B) const {
        const Real u0_inv = Real(1) / u0;
        const Real v = sqrt(Real(1) - u0_inv * u0_inv);
        const Real epsilon = T00 - v * M;
        const Real rhob = J0B * u0_inv;
        const Real p = pressure(epsilon, rhob, u0_inv);
        const Real temp1 = Real(1) - K00 / ((T00 + p) * (T00 + p));
        return u0 - Real(1) / sqrt(temp1);
    }

    // Reconst::solve_u0_Hybrid
    MD_FN int solveU0Hybrid(
        Real u0Guess, Real T00, Real K00, Real M, Real J0B,
        MD_THREAD Real &u0Sol) const {
        int status = 1;
        Real l = Real(1);
        Real h = Real(1e5);
        if (u0Guess >= Real(1)) {
            l = dmax(l, Real(0.99) * u0Guess);
            h = Real(1.01) * u0Guess;
        }
        Real fl = u0F(l, T00, K00, M, J0B);
        Real fh = u0F(h, T00, K00, M, J0B);
        if (fabs(fl) < P_->absErr) {
            u0Sol = l;
            return 1;
        }
        if (fabs(fh) < P_->absErr) {
            u0Sol = h;
            return 1;
        }
        if (fl * fh > Real(0)) return 0;
        Real duPrev = h - l;
        Real duCurr = duPrev;
        Real root = (h + l) / Real(2);
        Real fu = u0F(root, T00, K00, M, J0B);
        Real df = (fh - fl) / (h - l);
        Real absE = Real(10);
        Real relE = Real(10);
        int iter = 0;
        do {
            iter++;
            const Real uPrev = root;
            const Real fPrev = fu;
            if (((root - h) * df - fu) * ((root - l) * df - fu) > Real(0)
                || (fabs(Real(2) * fu) > fabs(duPrev * df))) {
                duPrev = duCurr;
                duCurr = (h - l) / Real(2);
                root = l + duCurr;
            } else {
                duPrev = duCurr;
                duCurr = fu / df;
                root = root - duCurr;
            }
            if (root != root) {
                status = 0;
                break;
            }
            fu = u0F(root, T00, K00, M, J0B);
            df = (fu - fPrev) / (root - uPrev);
            absE = duCurr;
            relE = duCurr / root;
            if (fu * fl < Real(0)) {
                h = root;
                fh = fu;
            } else {
                l = root;
                fl = fu;
            }
            if (iter > P_->maxIter) {
                status = 0;
                break;
            }
        } while (fabs(absE) > P_->absErr && fabs(relE) > P_->relErr);
        u0Sol = root;
        return status;
    }

    // Reconst::ReconstIt_velocity_Newton: 1 solved, -1 revert, -2 regulate
    MD_FN int velocityNewton(
        MD_THREAD Cell<Real> &g, MD_THREAD const Real *q,
        MD_THREAD const Cell<Real> &gpt, MD_THREAD EosState<Real> *eosOut,
        Real pGuess, MD_THREAD int &flag) const {
        if (eosOut) eosOut->valid = false;
        const Real K00 = q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
        const Real M = sqrt(K00);
        const Real T00 = q[0];
        const Real J0B = q[4], J0Q = q[5], J0S = q[6];
        Solver<Real> s = *this;
        s.J0Q_ = J0Q;
        s.J0S_ = J0S;
        if (T00 < P_->absErr) return -2;
        if (T00 < M) return -1;
        Real vGuess = sqrt(Real(1) - Real(1) / (gpt.u[0] * gpt.u[0] + P_->absErr));
        if (pGuess > Real(0)) vGuess = M / (T00 + pGuess);
        Real v = Real(0);
        int st = s.solveVFromGuess(vGuess, T00, M, J0B, v);
        if (st == 0) st = s.solveVHybrid(vGuess, T00, M, J0B, v, flag);
        if (st == 0) return -1;
        Real u[4];
        u[0] = Real(1) / (sqrt(Real(1) - v * v) + v * P_->absErr);
        Real epsilon = T00 - v * sqrt(K00);
        Real rhob = J0B / u[0];
        Real rhoq = J0Q / u[0];
        Real rhos = J0S / u[0];
        if (v > P_->vCritical && epsilon > Real(1e-6)) {
            Real u0Sol = u[0];
            if (s.solveU0Hybrid(u[0], T00, K00, M, J0B, u0Sol) == 1) {
                u[0] = u0Sol;
                epsilon =
                    T00 - sqrt((Real(1) - Real(1) / (u0Sol * u0Sol)) * K00);
                rhob = J0B / u0Sol;
                rhoq = J0Q / u0Sol;
                rhos = J0S / u0Sol;
            }
        }
        if (fabs(u[0] - gpt.u[0]) / gpt.u[0] > Real(100)) return -1;
        g.e = epsilon;
        g.rhob = rhob;
        g.rhoq = rhoq;
        g.rhos = rhos;
        Real p;
        if (eosOut) {
            eosOut->g = eosPressureGrad(*d_, tab, epsilon, rhob, rhoq, rhos);
            eosOut->valid = true;
            p = eosOut->g.p;
        } else {
            p = eosPressure(*d_, tab, epsilon, rhob, rhoq, rhos);
        }
        const Real fac = u[0] / (T00 + p);
        u[1] = q[1] * fac;
        u[2] = q[2] * fac;
        u[3] = q[3] * fac;
        const Real u2 = u[1] * u[1] + u[2] * u[2] + u[3] * u[3];
        if (fabs(u[0] * u[0] - u2 - Real(1)) > P_->absErr) {
            const Real sf = sqrt((u[0] * u[0] - Real(1)) / (u2 + P_->absErr));
            u[1] *= sf;
            u[2] *= sf;
            u[3] *= sf;
        }
        for (int k = 0; k < 4; k++) g.u[k] = u[k];
        return 1;
    }

    // Reconst::ReconstIt_shell
    MD_FN Cell<Real> reconstruct(
        Real tau, MD_THREAD const Real *tauq, MD_THREAD const Cell<Real> &gpt,
        MD_THREAD EosState<Real> *eosOut, Real pGuess,
        MD_THREAD int &flag) const {
        if (eosOut) eosOut->valid = false;
        Real q[7];
        for (int i = 0; i < 7; i++) q[i] = tauq[i] / tau;
        const Real M = sqrt(q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
        if (q[0] >= P_->absErr && M >= q[0]) {
            flag |= kFlagMomentumCap;
            const Real mm = Real(0.999) * q[0];
            const Real Mmax2 = mm * mm;
            const Real Mperp2 = q[1] * q[1] + q[2] * q[2];
            const Real Meta2 = q[3] * q[3];
            if (Meta2 < Mmax2 && Mperp2 > Real(0)) {
                const Real scale = sqrt((Mmax2 - Meta2) / Mperp2);
                q[1] *= scale;
                q[2] *= scale;
            } else {
                const Real scale = Real(0.999) * q[0] / M;
                for (int i = 1; i < 4; i++) q[i] *= scale;
            }
        }
        Cell<Real> g = gpt;
        const int st = velocityNewton(g, q, gpt, eosOut, pGuess, flag);
        if (st == -1) {
            g = gpt;
            flag |= kFlagRevert;
        } else if (st == -2) {
            g.e = dmax(P_->absErr, q[0]);
            g.rhob = q[4];
            g.rhoq = q[5];
            g.rhos = q[6];
            g.u[0] = Real(1);
            g.u[1] = Real(0);
            g.u[2] = Real(0);
            g.u[3] = Real(0);
            flag |= kFlagRegulate;
        }
        return g;
    }
};

// ----- KT fluxes -----------------------------------------------------------

//! Advance::MaxSpeed
template <typename Real>
MD_FN Real maxSpeed(
    Real tau, int direc, MD_THREAD const Cell<Real> &c,
    MD_THREAD const EosGrad<Real> &s, MD_THREAD int &flag) {
    const Real utau = c.u[0];
    const Real utau2 = utau * utau;
    const Real ux = fabs(c.u[direc]);
    const Real ut2mux2 = utau2 - ux * ux;
    const Real eps = c.e;
    const Real vs2 = s.cs2;
    const Real dpde = s.dpde;
    const Real numSqrt = (ut2mux2 - (ut2mux2 - Real(1)) * vs2) * vs2;
    Real num;
    if (numSqrt >= Real(0)) {
        num = utau * ux * (Real(1) - vs2) + sqrt(numSqrt);
    } else {
        const Real h = s.p + eps;
        if (dpde < Real(0.001)) {
            num = (sqrt(-(h * dpde * h * (dpde * (Real(-1) + ut2mux2) - ut2mux2)))
                   - h * (Real(-1) + dpde) * utau * ux);
        } else {
            flag |= kFlagMaxSpeedFatal;
            num = Real(0);
        }
    }
    const Real den = utau2 * (Real(1) - vs2) + vs2;
    Real f = num / dmax(den, Real(1e-16));
    if (f < Real(0)) {
        flag |= kFlagMaxSpeedFatal;
    } else if (f < ux / utau) {
        if (num != Real(0)) {
            if (fabs(f - ux / utau) < Real(0.0001)) {
                f = ux / utau;
            } else {
                flag |= kFlagMaxSpeedFatal;
            }
        }
    } else if (f > Real(1)) {
        flag |= kFlagMaxSpeedFatal;
        f = Real(1);
    }
    if (direc == 3) f *= Real(1) / tau;
    return f;
}

//! the pressure of every cell of a slice (Advance::computePressure)
template <typename Real>
MD_FN void kernelPressure(
    MD_CONST const StepParams<Real> &P, MD_CONST const EosDesc<Real> &d,
    MD_GLOBAL const Real *tab, MD_GLOBAL const Real *s, MD_GLOBAL Real *p,
    int i) {
    const int n = P.nx * P.ny * P.neta;
    p[i] = eosPressure(d, tab, s[i], s[n + i], s[2 * n + i], s[3 * n + i]);
}

MD_FN int faceCount(int nx, int ny, int neta, int direction) {
    return (nx + (direction == 1)) * (ny + (direction == 2))
           * (neta + (direction == 3));
}

//! Advance::computeFaceFlux for face number f along direction
template <typename Real>
MD_FN void kernelFace(
    MD_CONST const StepParams<Real> &P, MD_CONST const EosDesc<Real> &d,
    MD_GLOBAL const Real *tab, MD_GLOBAL const Real *s,
    MD_GLOBAL const Real *pressure, int direction, int f, MD_GLOBAL Real *flux,
    MD_FLAGS flags) {
    const int nx = P.nx, ny = P.ny, neta = P.neta;
    const int n = nx * ny * neta;
    const int dx = (direction == 1), dy = (direction == 2),
              deta = (direction == 3);
    const int nfx = nx + dx, nfy = ny + dy;
    const int fx = f % nfx;
    const int fy = (f / nfx) % nfy;
    const int feta = f / (nfx * nfy);
    const Real tau = P.tau + Real(P.rkFlag) * P.dtau;
    const Real tauFac = (direction == 3 ? Real(1) : tau);
    const int ILL =
        cellIndex(nx, ny, neta, fx - 2 * dx, fy - 2 * dy, feta - 2 * deta);
    const int IL = cellIndex(nx, ny, neta, fx - dx, fy - dy, feta - deta);
    const int IR = cellIndex(nx, ny, neta, fx, fy, feta);
    const int IRR = cellIndex(nx, ny, neta, fx + dx, fy + dy, feta + deta);
    const Cell<Real> cLL = loadCell(s, n, ILL);
    const Cell<Real> cL = loadCell(s, n, IL);
    const Cell<Real> cR = loadCell(s, n, IR);
    const Cell<Real> cRR = loadCell(s, n, IRR);
    const Solver<Real> solver = {&P, &d, tab};
    int flag = 0;
    Cell<Real> gL, gR;
    EosState<Real> eL, eR;
    Real qL[7], qR[7];
    if (P.mode == 1) {
        // Advance::computeFaceFluxPrimitive
        const Real th = P.theta;
        gL.e = cL.e + Real(0.5) * minmodDx(cR.e, cL.e, cLL.e, th);
        gR.e = cR.e - Real(0.5) * minmodDx(cRR.e, cR.e, cL.e, th);
        gL.rhob = cL.rhob + Real(0.5) * minmodDx(cR.rhob, cL.rhob, cLL.rhob, th);
        gR.rhob = cR.rhob - Real(0.5) * minmodDx(cRR.rhob, cR.rhob, cL.rhob, th);
        gL.rhoq = cL.rhoq + Real(0.5) * minmodDx(cR.rhoq, cL.rhoq, cLL.rhoq, th);
        gR.rhoq = cR.rhoq - Real(0.5) * minmodDx(cRR.rhoq, cR.rhoq, cL.rhoq, th);
        gL.rhos = cL.rhos + Real(0.5) * minmodDx(cR.rhos, cL.rhos, cLL.rhos, th);
        gR.rhos = cR.rhos - Real(0.5) * minmodDx(cRR.rhos, cR.rhos, cL.rhos, th);
        Real u2L = Real(0), u2R = Real(0);
        for (int i = 1; i < 4; i++) {
            gL.u[i] = cL.u[i] + Real(0.5) * minmodDx(cR.u[i], cL.u[i], cLL.u[i], th);
            gR.u[i] = cR.u[i] - Real(0.5) * minmodDx(cRR.u[i], cR.u[i], cL.u[i], th);
            u2L += gL.u[i] * gL.u[i];
            u2R += gR.u[i] * gR.u[i];
        }
        gL.u[0] = sqrt(Real(1) + u2L);
        gR.u[0] = sqrt(Real(1) + u2R);
        eL.g = eosPressureGrad(d, tab, gL.e, gL.rhob, gL.rhoq, gL.rhos);
        eR.g = eosPressureGrad(d, tab, gR.e, gR.rhob, gR.rhoq, gR.rhos);
        for (int a = 0; a < 7; a++) {
            qL[a] = tau * TJb(gL, a, 0, eL.g.p);
            qR[a] = tau * TJb(gR, a, 0, eR.g.p);
        }
    } else {
        const Real pLL = pressure[ILL], pL = pressure[IL];
        const Real pR = pressure[IR], pRR = pressure[IRR];
        for (int a = 0; a < 7; a++) {
            const Real gLL = tau * TJb(cLL, a, 0, pLL);
            const Real gL_ = tau * TJb(cL, a, 0, pL);
            const Real gR_ = tau * TJb(cR, a, 0, pR);
            const Real gRR = tau * TJb(cRR, a, 0, pRR);
            const Real fL = Real(0.5) * minmodDx(gR_, gL_, gLL, P.theta);
            const Real fR = Real(-0.5) * minmodDx(gRR, gR_, gL_, P.theta);
            qL[a] = gL_ + fL;
            qR[a] = gR_ + fR;
        }
        gL = solver.reconstruct(tau, qL, cL, &eL, pL, flag);
        gR = solver.reconstruct(tau, qR, cR, &eR, pR, flag);
        if (!eL.valid)
            eL.g = eosPressureGrad(d, tab, gL.e, gL.rhob, gL.rhoq, gL.rhos);
        if (!eR.valid)
            eR.g = eosPressureGrad(d, tab, gR.e, gR.rhob, gR.rhoq, gR.rhos);
    }
    const Real aL = maxSpeed(tau, direction, gL, eL.g, flag);
    const Real aR = maxSpeed(tau, direction, gR, eR.g, flag);
    const Real a = dmax(aL, aR);
    const int nFaces = faceCount(nx, ny, neta, direction);
    for (int al = 0; al < 7; al++) {
        const Real FL = tauFac * TJb(gL, al, direction, eL.g.p);
        const Real FR = tauFac * TJb(gR, al, direction, eR.g.p);
        flux[al * nFaces + f] = Real(0.5) * ((FL + FR) - a * (qR[al] - qL[al]));
    }
    if (flag) flagOr(flags, flag);
}

//! the first part of Advance::FirstRKStepT for cell i: tau T^{tau mu} of
//! the stage's slice (or the stored one) plus MakeDeltaQI
template <typename Real>
MD_FN void cellQiFlux(
    MD_CONST const StepParams<Real> &P, MD_GLOBAL const Real *sCurr,
    MD_GLOBAL const Real *pCurr, MD_GLOBAL const Real *fluxX,
    MD_GLOBAL const Real *fluxY, MD_GLOBAL const Real *fluxEta, int i,
    MD_THREAD Real *qi) {
    const int nx = P.nx, ny = P.ny, neta = P.neta;
    const int n = nx * ny * neta;
    const int ix = i % nx;
    const int iy = (i / nx) % ny;
    const int ieta = i / (nx * ny);
    const Real tau_rk = P.tau + Real(P.rkFlag) * P.dtau;
    const Cell<Real> cellCurr = loadCell(sCurr, n, i);
    const Real pressure = pCurr[i];
    for (int a = 0; a < 7; a++) {
        qi[a] = (P.storeQ ? sCurr[(kNPrim + a) * n + i]
                          : TJb(cellCurr, a, 0, pressure) * tau_rk);
    }
    // MakeDeltaQI
    Real rhs[7];
    Real Tm[4], Tp[4];
    for (int a = 0; a < 7; a++) rhs[a] = Real(0);
    for (int k = 0; k < 4; k++) {
        Tm[k] = Real(0);
        Tp[k] = Real(0);
    }
    for (int direction = 1; direction <= 3; direction++) {
        const int dx = (direction == 1), dy = (direction == 2),
                  deta = (direction == 3);
        const int nfx = nx + dx, nfy = ny + dy;
        const int nFaces = faceCount(nx, ny, neta, direction);
        const int Imh = ix + nfx * (iy + nfy * ieta);
        const int Iph = (ix + dx) + nfx * ((iy + dy) + nfy * (ieta + deta));
        MD_GLOBAL const Real *F =
            (direction == 1 ? fluxX : (direction == 2 ? fluxY : fluxEta));
        const Real delta = (direction == 1 ? P.dx : (direction == 2 ? P.dy : P.deta));
        for (int a = 0; a < 7; a++) {
            const Real Fiph = F[a * nFaces + Iph];
            const Real Fimh = F[a * nFaces + Imh];
            if (direction == 3 && (a == 0 || a == 3)) {
                Tm[a] = Fimh;
                Tp[a] = Fiph;
            } else {
                const Real DFmmp = (Fimh - Fiph) / delta;
                rhs[a] += DFmmp * P.dtau;
            }
        }
    }
    rhs[0] += (((Tm[0] - Tp[0]) * P.coshDeta - (Tm[3] + Tp[3]) * P.sinhDeta)
               * P.dtau);
    rhs[3] += (((Tm[3] - Tp[3]) * P.coshDeta - (Tm[0] + Tp[0]) * P.sinhDeta)
               * P.dtau);
    for (int a = 0; a < 7; a++) qi[a] += rhs[a];
}

//! the RK average of Advance::FirstRKStepT: qi += q(tau_n), qi /= 2 in the
//! second stage (sPrev, pPrev: the slice at tau_n)
template <typename Real>
MD_FN void cellRkAverage(
    MD_CONST const StepParams<Real> &P, MD_GLOBAL const Real *sPrev,
    MD_GLOBAL const Real *pPrev, int i, MD_THREAD Real *qi) {
    if (P.rkFlag != 1) return;
    const int n = P.nx * P.ny * P.neta;
    const Cell<Real> cellPrev = loadCell(sPrev, n, i);
    const Real pressurePrev = pPrev[i];
    // one statement as in FirstRKStepT: compilers fuse the multiply-add
    // within a statement (clang) or across them (gcc), and the device
    // reproduces the host only with the same grouping
    for (int a = 0; a < 7; a++) {
        if (P.storeQ) {
            qi[a] += sPrev[(kNPrim + a) * n + i];
        } else {
            qi[a] += Real(P.rkFlag) * TJb(cellPrev, a, 0, pressurePrev) * P.tau;
        }
        qi[a] *= Real(1) / (Real(1) + Real(P.rkFlag));
    }
}

//! the source part of Advance::FirstRKStepT: qi += qi_source dtau
template <typename Real>
MD_FN void cellAddSource(
    MD_CONST const StepParams<Real> &P, MD_GLOBAL const Real *src, int i,
    MD_THREAD Real *qi) {
    if (!P.hasSource) return;
    const int n = P.nx * P.ny * P.neta;
    for (int a = 0; a < 7; a++) qi[a] += src[a * n + i] * P.dtau;
}

//! Advance::FirstRKStepT (no viscous terms) with MakeDeltaQI and
//! reconstructAndStore, for cell i
template <typename Real>
MD_FN void kernelCell(
    MD_CONST const StepParams<Real> &P, MD_CONST const EosDesc<Real> &d,
    MD_GLOBAL const Real *tab, MD_GLOBAL const Real *sCurr,
    MD_GLOBAL const Real *sPrev, MD_GLOBAL const Real *pCurr,
    MD_GLOBAL const Real *pPrev, MD_GLOBAL const Real *fluxX,
    MD_GLOBAL const Real *fluxY, MD_GLOBAL const Real *fluxEta,
    MD_GLOBAL const Real *src, MD_GLOBAL Real *sNext, int i, MD_FLAGS flags) {
    const int n = P.nx * P.ny * P.neta;
    const Cell<Real> cellCurr = loadCell(sCurr, n, i);
    const Real pressure = pCurr[i];
    Real qi[7];
    cellQiFlux(P, sCurr, pCurr, fluxX, fluxY, fluxEta, i, qi);
    cellAddSource(P, src, i, qi);
    cellRkAverage(P, sPrev, pPrev, i, qi);
    const Solver<Real> solver = {&P, &d, tab};
    int flag = 0;
    const Cell<Real> out = solver.reconstruct(
        P.tau + P.dtau, qi, cellCurr, (MD_THREAD EosState<Real> *)0, pressure,
        flag);
    storeCell(sNext, n, i, out);
    if (P.storeQ) {
        // a capped, regulated or reverted cell keeps the q of its state
        const bool changed =
            (flag & (kFlagMomentumCap | kFlagRegulate | kFlagRevert)) != 0;
        const Real tauNext = P.tau + P.dtau;
        const Real pOut =
            (changed ? eosPressure(d, tab, out.e, out.rhob, out.rhoq, out.rhos)
                     : Real(0));
        for (int a = 0; a < 7; a++) {
            sNext[(kNPrim + a) * n + i] =
                (changed ? TJb(out, a, 0, pOut) * tauNext : qi[a]);
        }
    }
    if (flag) flagOr(flags, flag);
}

}  // namespace music_device

#endif  // SRC_DEVICE_DEVICE_IDEAL_H_
