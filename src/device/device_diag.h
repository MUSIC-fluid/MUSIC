// Copyright 2026 Chun Shen
//
// Diagnostics on the device, so that the host Fields need not be filled in
// every step: kernelConservationRow is the sum over one (eta, x) row of
// Cell_info::check_conservation_law (the 7 sums and the 5 maxima), for the
// current slice and the previous one. The host adds the rows in the order
// of check_conservation_law (in double, backend 1 reproduces it bit for
// bit).

#ifndef SRC_DEVICE_DEVICE_DIAG_H_
#define SRC_DEVICE_DEVICE_DIAG_H_

#if !defined(__METAL_VERSION__)
#include "device_visc.h"
#endif

namespace music_device {

//! values per row: N_B, N_Q, N_S, T^{tau t}, T^{tau x}, T^{tau y},
//! T^{tau z}, then the maxima of e, n_B, n_Q, n_S, T
enum { kConsSums = 7, kConsVals = 12 };

//! row = ix + nx ieta; P.tau: the tau of sCurr; sPrev: the previous slice
//! (read for W and Pi only with the viscous layout)
template <typename Real>
MD_FN void kernelConservationRow(
    MD_CONST const StepParams<Real> &P, MD_CONST const EosDesc<Real> &d,
    MD_GLOBAL const Real *tab, MD_GLOBAL const Real *sCurr,
    MD_GLOBAL const Real *sPrev, MD_GLOBAL Real *out, int row) {
    const int nx = P.nx, ny = P.ny;
    const int n = nx * ny * P.neta;
    const int ix = row % nx;
    const int ieta = row / nx;
    const bool viscous = (MD_RUN(P, nComp) == kNCompVisc);
    const Real tau = P.tau;
    const Real eta_s = P.etaMin + P.deta * Real(ieta);
    const Real cosh_eta = cosh(eta_s);
    const Real sinh_eta = sinh(eta_s);
    Real rN_B = Real(0.0), rN_Q = Real(0.0), rN_S = Real(0.0);
    Real rT_t = Real(0.0), rT_x = Real(0.0), rT_y = Real(0.0),
         rT_z = Real(0.0);
    Real eps_max = Real(0.0);
    Real rhob_max = Real(0.0);
    Real rhoq_max = Real(-10.0);
    Real rhos_max = Real(-10.0);
    Real T_max = Real(0.0);
    for (int iy = 0; iy < ny; iy++) {
        const int Idx = ix + nx * (iy + ny * ieta);
        const Real wCurr = Real(0.5) * (Real(1.) + P.dtau / tau);
        const Real wPrev = Real(0.5);
        const Real e_local = sCurr[Idx];
        const Real rhob = sCurr[n + Idx];
        const Real rhoq = sCurr[2 * n + Idx];
        const Real rhos = sCurr[3 * n + Idx];
        const Real pressure = eosPressure(d, tab, e_local, rhob, rhoq, rhos);
        Real PiCurr = Real(0.), PiPrev = Real(0.);
        if (viscous) {
            PiCurr = sCurr[kOffPi * n + Idx] + sCurr[kOffPiChem * n + Idx];
            PiPrev = sPrev[kOffPi * n + Idx] + sPrev[kOffPiChem * n + Idx];
        }
        const Real u0 = sCurr[4 * n + Idx];
        const Real u0Prev = sPrev[4 * n + Idx];
        Real T[4];
        for (int a = 0; a < 4; a++) {
            const Real g = (a == 0 ? Real(-1.) : Real(0.));
            const Real ua = sCurr[(4 + a) * n + Idx];
            const Real WCurr = (viscous ? sCurr[(kOffW + a) * n + Idx] : Real(0.));
            const Real WPrev = (viscous ? sPrev[(kOffW + a) * n + Idx] : Real(0.));
            const Real VCurr = WCurr + PiCurr * (g + u0 * ua);
            const Real VPrev =
                WPrev + PiPrev * (g + u0Prev * sPrev[(4 + a) * n + Idx]);
            T[a] = (e_local + pressure) * u0 * ua
                   + (a == 0 ? -pressure : Real(0.)) + wCurr * VCurr
                   + wPrev * VPrev;
        }
        const Real q0Curr = (viscous ? sCurr[(kOffW + 10) * n + Idx] : Real(0.));
        const Real q0Prev = (viscous ? sPrev[(kOffW + 10) * n + Idx] : Real(0.));
        rN_B += rhob * u0 + wCurr * q0Curr + wPrev * q0Prev;
        rN_Q += rhoq * u0;
        rN_S += rhos * u0;
        const Real T_tau_tau = T[0];
        const Real T_tau_eta = T[3];
        rT_t += T_tau_tau * cosh_eta + T_tau_eta * sinh_eta;
        rT_x += T[1];
        rT_y += T[2];
        rT_z += T_tau_tau * sinh_eta + T_tau_eta * cosh_eta;
        Real Tl, muB;
        eosTMuB(d, tab, e_local, rhob, rhoq, rhos, Tl, muB);
        eps_max = dmax(eps_max, e_local);
        rhob_max = dmax(rhob_max, rhob);
        T_max = dmax(T_max, Tl);
        rhoq_max = dmax(rhoq_max, rhoq);
        rhos_max = dmax(rhos_max, rhos);
    }
    MD_GLOBAL Real *r = out + kConsVals * row;
    r[0] = rN_B;
    r[1] = rN_Q;
    r[2] = rN_S;
    r[3] = rT_t;
    r[4] = rT_x;
    r[5] = rT_y;
    r[6] = rT_z;
    r[7] = eps_max;
    r[8] = rhob_max;
    r[9] = rhoq_max;
    r[10] = rhos_max;
    r[11] = T_max;
}

}  // namespace music_device

#endif  // SRC_DEVICE_DEVICE_DIAG_H_
