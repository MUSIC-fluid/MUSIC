// Copyright 2026 Chun Shen
//
// The viscous part of one RK stage on the device (shear and bulk; no baryon
// diffusion, no chemical bulk, no FlagCausalityMethod = 1):
//   kernelWFace: Diss::computeWDivergences, the values of one face (KT
//       fluxes of u^i W, face averages of tau W^{i alpha} and of the bulk
//       part)
//   kernelCellVisc: the cell loop of Advance::AdvanceIt with
//       viscosity_flag = 1 for one cell: FirstRKStepT with d_m W^{m alpha},
//       U_derivative (MakedU, theta, Du, sigma, omega), Diss::Make_uWRHS,
//       the thermal variables, FirstRKStepW (Make_uWSource, Make_uPiSource,
//       QuestRevert, QuestRevertResummedTransCoeff) and
//       stage2ViscousFixedPoint.
// As in device_ideal.h, every function follows its host version statement
// by statement: clang fuses a*b + c within a statement, and the device
// reproduces the host bit for bit (in double, with clang) only with the
// same statements. For the same reason DATA.gmunu is read from StepParams
// at run time, as the host reads it from DATA, rather than as constants the
// compiler could fold into a different contraction.
//
// Slices of viscous runs have kNCompVisc components: those of device_ideal.h,
// then W^{mu nu} (Fields::Wmunu_ order) at kOffW, Pi at kOffPi and Pi_chem at
// kOffPiChem.

#ifndef SRC_DEVICE_DEVICE_VISC_H_
#define SRC_DEVICE_DEVICE_VISC_H_

#ifndef MUSIC_DEVICE_BUNDLED
#include "device_ideal.h"
#endif

namespace music_device {

//! Util::map_2d_idx_to_1d
MD_FN int wIdx(int a, int b) {
    if (a == 4) return 10 + b;
    if (a > b) {
        const int t = a;
        a = b;
        b = t;
    }
    if (a == 0) return b;
    if (a == 1) return 3 + b;
    if (a == 2) return 5 + b;
    return 9;
}

//! Util::map_1d_idx_to_2d for idx_1d = 4 ... 8
MD_FN void wIdx2d(int idx1d, MD_THREAD int &a, MD_THREAD int &b) {
    a = (idx1d < 7 ? 1 : 2);
    b = (idx1d < 7 ? idx1d - 3 : idx1d - 5);
}

//! the slice component of the W flux component idx (Diss piIdxArr:
//! 4, 5, 6, 7, 8, 9 = Pi, 14 = Pi_chem, 11, 12, 13)
MD_FN int wFluxComp(int idx) {
    if (idx < 5) return kOffW + 4 + idx;
    if (idx == 5) return kOffPi;
    if (idx == 6) return kOffPiChem;
    return kOffW + idx + 4;
}

//! the viscous state of a cell
template <typename Real>
struct VCell {
    Cell<Real> c;
    Real W[14];
    Real pi, piChem;
};

template <typename Real>
MD_FN VCell<Real> loadVCell(MD_GLOBAL const Real *s, int n, int i) {
    VCell<Real> v;
    v.c = loadCell(s, n, i);
    for (int k = 0; k < 14; k++) v.W[k] = s[(kOffW + k) * n + i];
    v.pi = s[kOffPi * n + i];
    v.piChem = s[kOffPiChem * n + i];
    return v;
}

// ----- Diss::computeWDivergences: one face ---------------------------------

//! out[k nFaces + f], k < kNWFaceVals: the flux of W component k < nWComp,
//! then (at kNWFlux + alpha) the face averages
template <typename Real>
MD_FN void kernelWFace(
    MD_CONST const StepParams<Real> &P, MD_GLOBAL const Real *s,
    int direction, int f, MD_GLOBAL Real *out) {
    const int nx = P.nx, ny = P.ny, neta = P.neta;
    const int n = nx * ny * neta;
    const int dx = (direction == 1), dy = (direction == 2),
              deta = (direction == 3);
    const int nfx = nx + dx, nfy = ny + dy;
    const int fx = f % nfx;
    const int fy = (f / nfx) % nfy;
    const int feta = f / (nfx * nfy);
    const int nFaces = faceCount(nx, ny, neta, direction);
    const Real tau = P.tau + Real(P.rkFlag) * P.dtau;
    const Real tauFac = (direction == 3 ? Real(1.0) : tau);
    const bool withBulk = (MD_RUN(P, bulk) == 1 || MD_RUN(P, bulkChem) == 1);
    const int IL = cellIndex(nx, ny, neta, fx - dx, fy - dy, feta - deta);
    const int IR = cellIndex(nx, ny, neta, fx, fy, feta);
    const int ILL =
        cellIndex(nx, ny, neta, fx - 2 * dx, fy - 2 * dy, feta - 2 * deta);
    const int IRR = cellIndex(nx, ny, neta, fx + dx, fy + dy, feta + deta);
    MD_GLOBAL const Real *u0 = s + 4 * n;
    MD_GLOBAL const Real *ud = s + (4 + direction) * n;
    const Real aL = fabs(ud[IL]) / u0[IL];
    const Real aR = fabs(ud[IR]) / u0[IR];
    const Real ax = dmax(aL, aR);
    for (int idx = 0; idx < MD_RUN(P, nWComp); idx++) {
        MD_GLOBAL const Real *W = s + wFluxComp(idx) * n;
        const Real fLL = W[ILL] * ud[ILL];
        const Real fL = W[IL] * ud[IL];
        const Real fR = W[IR] * ud[IR];
        const Real fRR = W[IRR] * ud[IRR];
        const Real gLL = W[ILL] * u0[ILL];
        const Real gL = W[IL] * u0[IL];
        const Real gR = W[IR] * u0[IR];
        const Real gRR = W[IRR] * u0[IRR];
        const Real uWR = fR - Real(0.5) * minmodDx(fRR, fR, fL, P.theta);
        const Real uWL = fL + Real(0.5) * minmodDx(fR, fL, fLL, P.theta);
        const Real WR = gR - Real(0.5) * minmodDx(gRR, gR, gL, P.theta);
        const Real WL = gL + Real(0.5) * minmodDx(gR, gL, gLL, P.theta);
        out[idx * nFaces + f] = ((uWR + uWL) - ax * (WR - WL)) * Real(0.5);
    }
    for (int alpha = 0; alpha < 5; alpha++) {
        MD_GLOBAL const Real *W = s + (kOffW + wIdx(alpha, direction)) * n;
        const Real sL = W[IL] * tauFac;
        const Real sR = W[IR] * tauFac;
        out[(kNWFlux + alpha) * nFaces + f] = (sR + sL) * Real(0.5);
    }
    if (!withBulk) return;
    const Real PiL = s[kOffPi * n + IL] + s[kOffPiChem * n + IL];
    const Real PiR = s[kOffPi * n + IR] + s[kOffPiChem * n + IR];
    for (int alpha = 0; alpha < 4; alpha++) {
        MD_GLOBAL const Real *ua = s + (4 + alpha) * n;
        const Real gfac1 = (alpha == direction ? Real(1.0) : Real(0.0));
        const Real bL = PiL * (gfac1 + ua[IL] * ud[IL]) * tauFac;
        const Real bR = PiR * (gfac1 + ua[IR] * ud[IR]) * tauFac;
        out[(kNWFlux + 5 + alpha) * nFaces + f] = (bR + bL) * Real(0.5);
    }
}

// ----- TransportCoeffs -------------------------------------------------------

template <typename Real>
struct Transport {
    MD_CONST const StepParams<Real> *P_;

    MD_FN Real etaOverSDefault(Real T) const {
        const Real Ttr = Real(0.18) / Real(0.19733);
        const Real Tfrac = T / Ttr;
        Real eta_over_s;
        if (Tfrac < Real(1.)) {
            eta_over_s =
                (P_->shearToS + Real(0.0594) * (Real(1.) - Tfrac)
                 + Real(0.544) * (Real(1.) - Tfrac * Tfrac));
        } else {
            eta_over_s =
                (P_->shearToS + Real(0.288) * (Tfrac - Real(1.))
                 + Real(0.0818) * (Tfrac * Tfrac - Real(1.)));
        }
        return (eta_over_s);
    }

    MD_FN Real etaOverSDuke(Real T_in_fm) const {
        Real T_in_GeV = T_in_fm * Real(0.19733);
        Real Ttr_in_GeV = Real(0.154);
        Real Tfrac = T_in_GeV / Ttr_in_GeV;
        Real eta_over_s =
            (P_->shear2Min
             + (P_->shear2Slope) * (T_in_GeV - Ttr_in_GeV)
                   * pow(Tfrac, P_->shear2Curv));
        return eta_over_s;
    }

    MD_FN Real etaOverSSims(Real T_in_fm) const {
        Real T_in_GeV = T_in_fm * Real(0.19733);
        Real T_kink_in_GeV = P_->shear3TKink;
        Real low_T_slope = P_->shear3LowSlope;
        Real high_T_slope = P_->shear3HighSlope;
        Real eta_over_s_at_kink = P_->shear3AtKink;
        const Real eta_over_s_min = Real(1e-3);
        Real eta_over_s;
        if (T_in_GeV < T_kink_in_GeV) {
            eta_over_s = (eta_over_s_at_kink
                          + low_T_slope * (T_in_GeV - T_kink_in_GeV));
        } else {
            eta_over_s = (eta_over_s_at_kink
                          + high_T_slope * (T_in_GeV - T_kink_in_GeV));
        }
        eta_over_s = dmax(eta_over_s, eta_over_s_min);
        return (eta_over_s);
    }

    MD_FN Real shearProfileT(Real T_in_fm) const {
        const Real T_in_GeV = T_in_fm * Real(0.19733);
        const Real Tc = Real(0.165);
        Real f_T = Real(1.0);
        if (T_in_GeV < Tc) {
            const Real Tslope = Real(1.2);
            const Real Tlow = Real(0.1);
            f_T += Tslope * (Tc - T_in_GeV) / (Tc - Tlow);
        } else {
            const Real Tslope2 = Real(0.0);
            const Real Thigh = Real(0.4);
            f_T += Tslope2 * (T_in_GeV - Tc) / (Thigh - Tc);
        }
        return (f_T);
    }

    MD_FN Real shearMuBPiecewise(Real muB_in_fm) const {
        const Real muB_in_GeV = fabs(muB_in_fm) * Real(0.19733);
        const Real f0 = Real(1.0);
        const Real f1 = P_->shearMuBf0p2;
        const Real f2 = P_->shearMuBf0p4;
        Real f_muB = f0;
        if (muB_in_GeV < Real(0.2)) {
            f_muB = f0 + (f1 - f0) / Real(0.2) * muB_in_GeV;
        } else if (muB_in_GeV < Real(0.4)) {
            f_muB = f1 + (f2 - f1) / Real(0.2) * (muB_in_GeV - Real(0.2));
        } else {
            f_muB = f2;
        }
        return (f_muB);
    }

    MD_FN Real shearMuBProfile(Real muB_in_fm) const {
        const Real muB_in_GeV = muB_in_fm * Real(0.19733);
        const Real alpha = P_->shearMuBAlpha;
        const Real muB_slope = P_->shearMuBSlope;
        const Real muB_scale = P_->shearMuBScale;
        Real f_muB =
            Real(1.) + muB_slope * pow(muB_in_GeV / muB_scale, alpha);
        return (f_muB);
    }

    //! TransportCoeffs::get_eta_over_s (not T_dependent_shear_to_s = 23)
    MD_FN Real etaOverS(Real T, Real muB) const {
        Real eta_over_s = P_->shearToS;
        const int shear_T_ = P_->shearT;
        if (shear_T_ == 0) {
            eta_over_s = P_->shearToS;
        } else if (shear_T_ == 3) {
            eta_over_s = etaOverSSims(T);
        } else if (shear_T_ == 1) {
            eta_over_s = etaOverSDefault(T);
        } else if (shear_T_ == 2) {
            eta_over_s = etaOverSDuke(T);
        } else if (shear_T_ == 11) {
            eta_over_s *= shearProfileT(T);
        } else {
            eta_over_s = P_->shearToS;
        }
        if (P_->muBShear == 7) {
            eta_over_s *= shearMuBPiecewise(muB);
        } else if (P_->muBShear == 10) {
            eta_over_s *= shearMuBProfile(muB);
        }
        eta_over_s = dmin(Real(0.5), dmax(Real(0.), eta_over_s));
        return eta_over_s;
    }

    MD_FN Real zetaAsymGaussian(
        Real T_in_fm, Real B_norm, Real B_width1, Real B_width2,
        Real Tpeak) const {
        const Real T_in_GeV = T_in_fm * Real(0.19733);
        Real Tdiff = T_in_GeV - Tpeak;
        if (Tdiff > Real(0.)) {
            Tdiff = Tdiff / B_width2;
        } else {
            Tdiff = Tdiff / B_width1;
        }
        Real bulk = B_norm * exp(-Tdiff * Tdiff);
        return (bulk);
    }

    MD_FN Real zetaDuke(Real T_in_fm) const {
        const Real A = P_->bulk2Norm;
        const Real G = P_->bulk2Width;
        const Real Tpeak_in_GeV = P_->bulk2Peak;
        const Real T_in_GeV = T_in_fm * Real(0.19733);
        const Real diff_ratio = (T_in_GeV - Tpeak_in_GeV) / G;
        return A / (Real(1) + diff_ratio * diff_ratio);
    }

    MD_FN Real zetaSims(Real T_in_fm) const {
        const Real T_in_GeV = T_in_fm * Real(0.19733);
        const Real max = P_->bulk3Max;
        const Real width = P_->bulk3Width;
        const Real T_peak_in_GeV = P_->bulk3TPeak;
        const Real lambda = P_->bulk3Lambda;
        const Real diff = T_in_GeV - T_peak_in_GeV;
        const Real sign = Real((diff > Real(0)) - (diff < Real(0)));
        const Real diff_ratio = (diff) / (width * (lambda * sign + Real(1)));
        return max / (Real(1) + diff_ratio * diff_ratio);
    }

    MD_FN Real zetaDefault(Real T_in_fm) const {
        const Real T_in_GeV = T_in_fm * Real(0.19733);
        Real Ttr = Real(0.18);
        Real dummy = T_in_GeV / Ttr;
        Real A1 = Real(-13.77), A2 = Real(27.55), A3 = Real(13.45);
        Real lambda1 = Real(0.9), lambda2 = Real(0.25), lambda3 = Real(0.9),
             lambda4 = Real(0.22);
        Real sigma1 = Real(0.025), sigma2 = Real(0.13), sigma3 = Real(0.0025),
             sigma4 = Real(0.022);
        Real bulk = A1 * dummy * dummy + A2 * dummy - A3;
        if (T_in_GeV < Real(0.995) * Ttr) {
            bulk =
                (lambda3 * exp((dummy - Real(1)) / sigma3)
                 + lambda4 * exp((dummy - Real(1)) / sigma4) + Real(0.03));
        }
        if (T_in_GeV > Real(1.05) * Ttr) {
            bulk =
                (lambda1 * exp(-(dummy - Real(1)) / sigma1)
                 + lambda2 * exp(-(dummy - Real(1)) / sigma2) + Real(0.001));
        }
        return (bulk);
    }

    MD_FN Real zetaBigBroadP(Real T_in_fm) const {
        const Real T_in_GeV = T_in_fm * Real(0.19733);
        const Real B_norm = Real(0.24);
        const Real B_width = Real(1.5);
        const Real Tpeak = Real(0.165);
        const Real Ttilde = (T_in_GeV / Tpeak - Real(1.)) / B_width;
        Real bulk = B_norm / (Ttilde * Ttilde + Real(1.));
        if (T_in_GeV < Tpeak) {
            const Real Tdiff = (T_in_GeV - Tpeak) / Real(0.01);
            bulk = B_norm * exp(-Tdiff * Tdiff);
        }
        return (bulk);
    }

    //! TransportCoeffs::get_zeta_over_s (not T_dependent_bulk_to_s = 23)
    MD_FN Real zetaOverS(Real T, Real mu_B) const {
        Real muB_in_GeV = mu_B * Real(0.19733);
        Real zeta_over_s = Real(0.);
        const int bulk_T_ = P_->bulkT;
        if (bulk_T_ == 10) {
            Real peak_norm = P_->bulk10Max;
            if (muB_in_GeV < Real(0.2)) {
                peak_norm =
                    (P_->bulk10Max
                     + (P_->bulk10MaxMuB0p2 - P_->bulk10Max) / Real(0.2)
                           * muB_in_GeV);
            } else if (muB_in_GeV < Real(0.4)) {
                peak_norm =
                    (P_->bulk10MaxMuB0p2
                     + (P_->bulk10MaxMuB0p4 - P_->bulk10MaxMuB0p2) / Real(0.2)
                           * (muB_in_GeV - Real(0.2)));
            } else {
                peak_norm = P_->bulk10MaxMuB0p4;
            }
            const Real Tpeak =
                (P_->bulk10Tpeak
                 + P_->bulk10TpeakMuBCurv * muB_in_GeV * muB_in_GeV);
            const Real B_width1 = P_->bulk10WidthLow;
            const Real B_width2 = P_->bulk10WidthHigh;
            zeta_over_s =
                zetaAsymGaussian(T, peak_norm, B_width1, B_width2, Tpeak);
        } else if (bulk_T_ == 8) {
            const Real peak_norm = Real(0.13);
            const Real B_width1 = Real(0.01);
            const Real B_width2 = Real(0.12);
            const Real Tpeak = Real(0.160);
            zeta_over_s =
                zetaAsymGaussian(T, peak_norm, B_width1, B_width2, Tpeak);
        } else if (bulk_T_ == 1) {
            zeta_over_s = zetaDefault(T);
        } else if (bulk_T_ == 2) {
            zeta_over_s = zetaDuke(T);
        } else if (bulk_T_ == 3) {
            zeta_over_s = zetaSims(T);
        } else if (bulk_T_ == 7) {
            zeta_over_s = zetaBigBroadP(T);
        } else if (bulk_T_ == 9) {
            const Real peak_norm = Real(0.175);
            const Real B_width1 = Real(0.01);
            const Real B_width2 = Real(0.12);
            const Real Tpeak = Real(0.160);
            zeta_over_s =
                zetaAsymGaussian(T, peak_norm, B_width1, B_width2, Tpeak);
        }
        zeta_over_s = dmax(Real(0.), zeta_over_s);
        return zeta_over_s;
    }

    //! TransportCoeffs::getResummedCorrFactor
    MD_FN Real resummedCorrFactor(Real Rshear, Real Rbulk) const {
        Real resummedCorrection = Real(1e-4);
        Real r_combined = P_->resumAlpha * (Rshear + Rbulk);
        if (r_combined < Real(1.0)) {
            Real atanhFactor = atanh(r_combined);
            resummedCorrection = (dmax(
                resummedCorrection,
                Real(1.) / (Real(1.) + atanhFactor * atanhFactor)));
        }
        return (resummedCorrection);
    }
};

//! Util::UnpackVecToMatrix of the 10 W^{mu nu} (or sigma^{mu nu})
template <typename Real>
MD_FN void unpackSym(MD_THREAD const Real *v, MD_THREAD Real (*m)[4]) {
    m[0][0] = v[0];
    m[0][1] = v[1];
    m[0][2] = v[2];
    m[0][3] = v[3];
    m[1][0] = v[1];
    m[1][1] = v[4];
    m[1][2] = v[5];
    m[1][3] = v[6];
    m[2][0] = v[2];
    m[2][1] = v[5];
    m[2][2] = v[7];
    m[2][3] = v[8];
    m[3][0] = v[3];
    m[3][1] = v[6];
    m[3][2] = v[8];
    m[3][3] = v[9];
}

//! Util::UnpackVecToMatrix of a VorticityVec
template <typename Real>
MD_FN void unpackAnti(MD_THREAD const Real *v, MD_THREAD Real (*m)[4]) {
    m[0][0] = Real(0.0);
    m[0][1] = v[0];
    m[0][2] = v[1];
    m[0][3] = v[2];
    m[1][0] = -m[0][1];
    m[1][1] = Real(0.0);
    m[1][2] = v[3];
    m[1][3] = v[4];
    m[2][0] = -m[0][2];
    m[2][1] = -m[1][2];
    m[2][2] = Real(0.0);
    m[2][3] = v[5];
    m[3][0] = -m[0][3];
    m[3][1] = -m[1][3];
    m[3][2] = -m[2][3];
    m[3][3] = Real(0.0);
}

//! the gradients of one cell (U_derivative), shared by the stage-2 fixed
//! point
template <typename Real>
struct Gradients {
    Real theta;
    Real a[5];       // DumuVec
    Real sigma[10];  // VelocityShearVec
    Real omega[6];   // VorticityVec
    Real wRhs[10];   // Diss::Make_uWRHS
    Real dMuBT[4];   // DmuMuBoverTVec: partial^mu (mu_B / T)
};

template <typename Real, int EosType = -1>
struct ViscSolver {
    MD_CONST const StepParams<Real> *P_;
    MD_CONST const EosDesc<Real> *d_;
    MD_GLOBAL const Real *tab;

    //! U_derivative::idealEoMWeight
    MD_FN Real idealEoMWeight(Real eps) const {
        if (P_->duIdeal != 1 || eps < Real(1e-8)) return Real(0.);
        const Real eDilute = P_->duIdealEDilute / Real(0.19733);
        const Real blendRange = Real(5.);
        if (eps >= blendRange * eDilute) return Real(0.);
        if (eps <= eDilute) return Real(1.);
        return Real(1.) - log(eps / eDilute) / log(blendRange);
    }

    //! U_derivative::MakedU, calculate_expansion_rate, calculate_Du_supmu,
    //! calculate_velocity_shear_tensor,
    //! calculate_kinetic_vorticity_with_spatial_projector, then
    //! Diss::Make_uWRHS
    //! aCurr, aPrev: mu_B / T of sCurr and sPrev (kernelMuBOverT), read
    //! only with baryon diffusion
    MD_FN void gradients(
        Real tau, MD_GLOBAL const Real *sCurr, MD_GLOBAL const Real *sPrev,
        MD_GLOBAL const Real *pCurr, MD_GLOBAL const Real *aCurr,
        MD_GLOBAL const Real *aPrev, int i, MD_THREAD const Real *wFluxDiv,
        MD_THREAD Gradients<Real> &G) const {
        MD_CONST const StepParams<Real> &P = *P_;
        const int nx = P.nx, ny = P.ny, neta = P.neta;
        const int n = nx * ny * neta;
        const int ix = i % nx;
        const int iy = (i / nx) % ny;
        const int ieta = i / (nx * ny);
        Real dUsup[5][4];
        for (int m = 0; m < 5; m++)
            for (int k = 0; k < 4; k++) dUsup[m][k] = Real(0.0);
        Real dPsup[4] = {Real(0.0), Real(0.0), Real(0.0), Real(0.0)};
        // MakeDSpatial
        {
            const Real delta[4] = {Real(0.0), P.dx, P.dy, P.deta * tau};
            for (int dir = 1; dir <= 3; dir++) {
                const int ox = (dir == 1), oy = (dir == 2), oe = (dir == 3);
                const int Ic = cellIndex(nx, ny, neta, ix, iy, ieta);
                const int Ip1 =
                    cellIndex(nx, ny, neta, ix + ox, iy + oy, ieta + oe);
                const int Im1 =
                    cellIndex(nx, ny, neta, ix - ox, iy - oy, ieta - oe);
                for (int m = 1; m < 4; m++) {
                    MD_GLOBAL const Real *um = sCurr + (4 + m) * n;
                    const Real f = um[Ic];
                    const Real fp1 = um[Ip1];
                    const Real fm1 = um[Im1];
                    dUsup[m][dir] =
                        (minmodDx(fp1, f, fm1, P.theta) / delta[dir]);
                }
                if (MD_RUN(P, duIdeal) == 1) {
                    const Real Pc = pCurr[Ic];
                    const Real Pp1 = pCurr[Ip1];
                    const Real Pm1 = pCurr[Im1];
                    dPsup[dir] = (minmodDx(Pp1, Pc, Pm1, P.theta) / delta[dir]);
                }
                if (MD_RUN(P, diff) == 1) {
                    // partial_n (mu_B / T)
                    const Real f = aCurr[Ic];
                    const Real fp1 = aCurr[Ip1];
                    const Real fm1 = aCurr[Im1];
                    dUsup[4][dir] = (minmodDx(fp1, f, fm1, P.theta) / delta[dir]);
                }
            }
            for (int nn = 1; nn < 4; nn++) {
                Real f = Real(0.0);
                for (int m = 1; m < 4; m++) {
                    f += dUsup[m][nn] * (sCurr[(4 + m) * n + i]);
                }
                f /= sCurr[4 * n + i];
                dUsup[0][nn] = f;
            }
        }
        // MakeDTau
        Real u[4];
        for (int m = 0; m < 4; m++) u[m] = sCurr[(4 + m) * n + i];
        {
            for (int m = 0; m < 4; m++) {
                Real uCurr = sCurr[(4 + m) * n + i];
                Real uPrev = sPrev[(4 + m) * n + i];
                Real f = (uCurr - uPrev) / P.dtau;
                dUsup[m][0] = -f;
            }
            Real f = Real(0.0);
            for (int m = 1; m < 4; m++) {
                f += dUsup[m][0] * sCurr[(4 + m) * n + i];
            }
            dUsup[0][0] = f / sCurr[4 * n + i];
            if (MD_RUN(P, diff) == 1) {
                // the backward difference of mu_B / T
                const Real tildemu = aCurr[i];
                const Real tildemu_prev = aPrev[i];
                f = (tildemu - tildemu_prev) / (P.dtau);
                dUsup[4][0] = -f;
            }
        }
        // MakeDTauIdealEoM
        if (MD_RUN(P, duIdeal) == 1) {
            const Real eps = sCurr[i];
            const Real weight = idealEoMWeight(eps);
            if (weight > Real(0.)) {
                const EosGrad<Real> eg = eosPressureGrad<Real, EosType>(
                    *d_, tab, eps, sCurr[n + i], sCurr[2 * n + i],
                    sCurr[3 * n + i]);
                const Real Pp = eg.p;
                const Real cs2 = eg.cs2;
                const Real enthalpy = dmax(eps + Pp, tinyNormal<Real>());
                const Real K = enthalpy * cs2;
                Real spatialDiv = Real(0.);
                for (int nn = 1; nn < 4; nn++) spatialDiv += dUsup[nn][nn];
                Real dtu[4];
                Real th = Real(0.);
                Real theta0 = Real(0.), slope = Real(0.);
                // c_s^2 v^2 >= 1: no solution, keep the backward difference
                bool causal = true;
                for (int pass = 0; pass < 3; pass++) {
                    if (pass == 1) th = Real(1.);
                    if (pass == 2) {
                        if (!(slope < Real(1.))) {
                            causal = false;
                            break;
                        }
                        th = theta0 / (Real(1.) - slope);
                    }
                    Real udtu = Real(0.);
                    for (int ii = 1; ii < 4; ii++) {
                        Real conv = Real(0.);
                        for (int nn = 1; nn < 4; nn++)
                            conv += u[nn] * dUsup[ii][nn];
                        const Real a = -(dPsup[ii] - u[ii] * K * th) / enthalpy;
                        const Real christoffel =
                            (ii == 3 ? u[0] * u[3] / tau : Real(0.));
                        dtu[ii] = (a - conv - christoffel) / u[0];
                        udtu += u[ii] * dtu[ii];
                    }
                    dtu[0] = udtu / u[0];
                    const Real thetaOut = dtu[0] + spatialDiv + u[0] / tau;
                    if (pass == 0) theta0 = thetaOut;
                    if (pass == 1) slope = thetaOut - theta0;
                }
                bool finite = causal;
                for (int m = 0; m < 4; m++) {
                    if (!isfinite(dtu[m])) finite = false;
                }
                if (finite) {
                    for (int m = 0; m < 4; m++) {
                        dUsup[m][0] =
                            (Real(1.) - weight) * dUsup[m][0] - weight * dtu[m];
                    }
                }
            }
        }
        // calculate_expansion_rate
        {
            Real partial_mu_u_supmu = Real(0.0);
            for (int mu = 0; mu < 4; mu++) {
                Real gfac = (mu == 0 ? Real(-1.0) : Real(1.0));
                partial_mu_u_supmu += dUsup[mu][mu] * gfac;
            }
            G.theta = partial_mu_u_supmu + u[0] / tau;
        }
        // calculate_Du_supmu
        for (int mu = 0; mu <= 4; mu++) {
            Real u_supnu_partial_nu_u_supmu = Real(0.0);
            for (int nu = 0; nu < 4; nu++) {
                Real tfac = (nu == 0 ? Real(-1.0) : Real(1.0));
                u_supnu_partial_nu_u_supmu += (tfac * u[nu] * dUsup[mu][nu]);
            }
            G.a[mu] = u_supnu_partial_nu_u_supmu;
        }
        // calculate_velocity_shear_tensor
        {
            Real gfac = Real(0.0);
            Real sigma_local[4][4];
            const Real theta_u_local = G.theta;
            for (int a = 1; a < 4; a++) {
                for (int b = a; b < 4; b++) {
                    if (b == a) {
                        gfac = Real(1.0);
                    } else {
                        gfac = Real(0.0);
                    }
                    sigma_local[a][b] =
                        ((dUsup[a][b] + dUsup[b][a]) / Real(2.)
                         - (gfac + u[a] * u[b]) * theta_u_local / Real(3.)
                         + u[0] / tau * P.g[a][3] * P.g[b][3]
                         + u[3] * u[0] / tau / Real(2.)
                               * (P.g[a][3] * u[b] + P.g[b][3] * u[a])
                         + (u[a] * G.a[b] + u[b] * G.a[a]) / Real(2.));
                    sigma_local[b][a] = sigma_local[a][b];
                }
            }
            sigma_local[3][3] =
                ((Real(2.)
                      * (u[1] * u[2] * sigma_local[1][2]
                         + u[1] * u[3] * sigma_local[1][3]
                         + u[2] * u[3] * sigma_local[2][3])
                  - (u[0] * u[0] - u[1] * u[1]) * sigma_local[1][1]
                  - (u[0] * u[0] - u[2] * u[2]) * sigma_local[2][2])
                 / dmax(u[0] * u[0] - u[3] * u[3], Real(0.5)));
            for (int a = 1; a < 4; a++) {
                Real temp = Real(0.0);
                for (int b = 1; b < 4; b++) {
                    temp += sigma_local[a][b] * u[b];
                }
                sigma_local[0][a] = temp / u[0];
            }
            Real temp = Real(0.0);
            for (int a = 1; a < 4; a++) {
                temp += sigma_local[0][a] * u[a];
            }
            sigma_local[0][0] = temp / u[0];
            G.sigma[0] = sigma_local[0][0];
            G.sigma[1] = sigma_local[0][1];
            G.sigma[2] = sigma_local[0][2];
            G.sigma[3] = sigma_local[0][3];
            G.sigma[4] = sigma_local[1][1];
            G.sigma[5] = sigma_local[1][2];
            G.sigma[6] = sigma_local[1][3];
            G.sigma[7] = sigma_local[2][2];
            G.sigma[8] = sigma_local[2][3];
            G.sigma[9] = sigma_local[3][3];
        }
        // calculate_kinetic_vorticity_with_spatial_projector
        if (MD_RUN(P, vorticity)) {
            Real omega_local[4][4];
            for (int mu = 0; mu < 4; mu++) {
                omega_local[mu][mu] = Real(0.0);
            }
            for (int mu = 1; mu < 4; mu++) {
                for (int nu = mu + 1; nu < 4; nu++) {
                    omega_local[mu][nu] =
                        Real(0.5)
                        * ((dUsup[nu][mu] - dUsup[mu][nu])
                           + (u[mu] * G.a[nu] - u[nu] * G.a[mu])
                           + u[3] * u[0] / tau
                                 * (u[mu] * P.g[nu][3] - u[nu] * P.g[mu][3]));
                    omega_local[nu][mu] = -omega_local[mu][nu];
                }
            }
            for (int mu = 1; mu < 4; mu++) {
                Real temp = Real(0.0);
                for (int nu = 1; nu < 4; nu++) {
                    temp += omega_local[mu][nu] * u[nu];
                }
                omega_local[0][mu] = temp / u[0];
            }
            G.omega[0] = omega_local[0][1];
            G.omega[1] = omega_local[0][2];
            G.omega[2] = omega_local[0][3];
            G.omega[3] = omega_local[1][2];
            G.omega[4] = omega_local[1][3];
            G.omega[5] = omega_local[2][3];
        } else {
            for (int k = 0; k < 6; k++) G.omega[k] = Real(0.);
        }
        // Diss::Make_uWRHS
        {
            const VCell<Real> grid_pt = loadVCell(sCurr, n, i);
            Real Wmunu_local[4][4];
            unpackSym(grid_pt.W, Wmunu_local);
            for (int k = 0; k < 10; k++) G.wRhs[k] = Real(0.);
            // P.nWComp at run time, not MD_RUN: with this bound a compile-
            // time constant, the Metal compiler (M4, macOS 26) gives the
            // update with FlagResetCausality = 1 an error of 5e-5 per step
            // (gpu_backend_viscous, "causality"); the same loop with the
            // bound from P, and the kernels on the CPU in float, agree with
            // the CPU code to rounding
            for (int idx = 0; idx < P.nWComp; idx++) {
                G.wRhs[idx] += wFluxDiv[idx];
            }
            int mu = 0;
            int nu = 0;
            for (int idx_1d = 4; idx_1d < 9; idx_1d++) {
                wIdx2d(idx_1d, mu, nu);
                Real tempf =
                    (-(P.g[3][mu]) * (Wmunu_local[0][nu])
                     - (P.g[3][nu]) * (Wmunu_local[0][mu])
                     + (P.g[0][mu]) * (Wmunu_local[3][nu])
                     + (P.g[0][nu]) * (Wmunu_local[3][mu])
                     + (Wmunu_local[3][nu]) * (grid_pt.c.u[mu])
                           * (grid_pt.c.u[0])
                     + (Wmunu_local[3][mu]) * (grid_pt.c.u[nu])
                           * (grid_pt.c.u[0])
                     - (Wmunu_local[0][nu]) * (grid_pt.c.u[mu])
                           * (grid_pt.c.u[3])
                     - (Wmunu_local[0][mu]) * (grid_pt.c.u[nu])
                           * (grid_pt.c.u[3]))
                    * (grid_pt.c.u[3] / tau);
                for (int ic = 0; ic < 4; ic++) {
                    const Real ic_fac = (ic == 0 ? Real(-1.0) : Real(1.0));
                    tempf +=
                        ((Wmunu_local[ic][nu]) * (grid_pt.c.u[mu]) * (G.a[ic])
                             * ic_fac
                         + (Wmunu_local[ic][mu]) * (grid_pt.c.u[nu])
                               * (G.a[ic]) * ic_fac);
                }
                G.wRhs[idx_1d - 4] +=
                    (tempf * (P.dtau)
                     + (-(grid_pt.c.u[0] * Wmunu_local[mu][nu]) / tau
                        + (G.theta * Wmunu_local[mu][nu]))
                           * (P.dtau));
            }
            if (MD_RUN(P, bulk) == 1) {
                G.wRhs[5] -= (grid_pt.pi) * (grid_pt.c.u[0]) / tau * P.dtau;
                G.wRhs[5] += (grid_pt.pi) * G.theta * P.dtau;
            }
            if (MD_RUN(P, bulkChem) == 1) {
                G.wRhs[6] -=
                    (grid_pt.piChem) * (grid_pt.c.u[0]) / tau * P.dtau;
                G.wRhs[6] += (grid_pt.piChem) * G.theta * P.dtau;
            }
        }
        for (int mu = 0; mu < 4; mu++) G.dMuBT[mu] = dUsup[4][mu];
    }

    //! Diss::D_muB_over_T_ideal (the EoS at n_Q = n_S = 0, as on the host)
    MD_FN Real dMuBOverTIdeal(
        Real epsilon, Real rhob, Real pressure, Real theta) const {
        const Real small_eps = Real(1e-16);
        if (epsilon < small_eps) return Real(0.);
        if (sizeof(Real) == 4 && d_->type == 0) {
            // the two terms below cancel exactly for the conformal ideal
            // gas, where mu_B / T depends on n / s alone and the ideal flow
            // conserves n / s: with float's rounding in each 1e-4
            // difference only noise remained. One difference along
            // (e + P, n), the same derivative, avoids the cancellation.
            // (The tables of the other EoSs are piecewise linear, and a
            // step this large crosses their kinks.)
            const Real h = Real(1e-3);
            const Real dEp = h * (epsilon + pressure);
            const Real dNp = h * rhob;
            const Real Dalpha =
                -theta
                * (alphaB(epsilon + dEp, rhob + dNp)
                   - alphaB(epsilon - dEp, rhob - dNp))
                / (Real(2.) * h);
            return isfinite(Dalpha) ? Dalpha : Real(0.);
        }
        const Real de = Real(1e-4) * epsilon;
        const Real dalpha_de =
            (alphaB(epsilon + de, rhob) - alphaB(epsilon - de, rhob))
            / (Real(2.) * de);
        Real dalpha_dn = Real(0.);
        if (fabs(rhob) > small_eps) {
            const Real dn = Real(1e-4) * fabs(rhob);
            dalpha_dn = (alphaB(epsilon, rhob + dn) - alphaB(epsilon, rhob - dn))
                        / (Real(2.) * dn);
        }
        const Real Dalpha =
            -theta * ((epsilon + pressure) * dalpha_de + rhob * dalpha_dn);
        return isfinite(Dalpha) ? Dalpha : Real(0.);
    }
    MD_FN Real alphaB(Real e, Real nb) const {
        Real T, muB;
        eosTMuB<Real, EosType>(*d_, tab, e, nb, Real(0.), Real(0.), T, muB);
        return (muB / dmax(T, Real(1e-16)));
    }

    //! Diss::Make_uqSource
    MD_FN void uqSource(
        Real tau, MD_THREAD const VCell<Real> &grid_pt,
        MD_THREAD const Gradients<Real> &G, MD_THREAD const Thermal<Real> &th,
        MD_THREAD Real *sourceTerms) const {
        MD_CONST const StepParams<Real> &P = *P_;
        const Real small_eps = Real(1e-16);
        Real epsilon = th.e;
        Real rhob = th.rhob;
        Real pressure = th.p;
        Real T = th.T;
        Real kappa_coefficient = P.kappaCoeff;
        Real tau_rho = kappa_coefficient / dmax(T, small_eps);
        tau_rho = dmin(Real(10.), dmax(Real(3.) * P.dtau, tau_rho));
        Real mub = th.muB;
        Real alpha = mub / dmax(T, small_eps);
        Real denorm_safe = copysign(
            dmax(fabs(Real(3.) * T * tanh(alpha)), small_eps),
            Real(3.) * T * tanh(alpha));
        Real kappa = kappa_coefficient
                     * (rhob / denorm_safe
                        - rhob * rhob / dmax(epsilon + pressure, small_eps));
        if (MD_RUN(P, profile1)) {
            Real denorm_safe1 = copysign(dmax(fabs(mub), small_eps), mub);
            kappa = kappa_coefficient * (rhob / denorm_safe1);
        }
        kappa = dmax(Real(0.), kappa);
        Real q[4];
        for (int ii = 0; ii < 4; ii++) q[ii] = grid_pt.W[10 + ii];
        Real transport_coeff = Real(1.0) * tau_rho;
        Real transport_coeff_2 = Real(3.) / Real(5.) * tau_rho;
        Real sigma[4][4];
        unpackSym(G.sigma, sigma);
        const Real DmuBoverT =
            dMuBOverTIdeal(epsilon, rhob, pressure, G.theta);
        for (int nu = 1; nu < 4; nu++) {
            Real NS = kappa * (G.dMuBT[nu] + grid_pt.c.u[nu] * DmuBoverT);
            Real Nonlinear1 = -transport_coeff * q[nu] * G.theta;
            Real temptemp = Real(0.0);
            for (int ii = 0; ii < 4; ii++) {
                temptemp += q[ii] * sigma[ii][nu] * P.g[ii][ii];
            }
            Real Nonlinear2 = -transport_coeff_2 * temptemp;
            Real Nonlinear3 = Real(0.0);
            if (MD_RUN(P, vorticity)) {
                Real transport_coeff_3 = Real(1.0) * tau_rho;
                Real omega[4][4];
                unpackAnti(G.omega, omega);
                Real temp3 = Real(0.0);
                for (int ii = 0; ii < 4; ii++) {
                    temp3 += q[ii] * omega[ii][nu] * P.g[ii][ii];
                }
                Nonlinear3 = -transport_coeff_3 * temp3;
            }
            sourceTerms[nu - 1] =
                ((-q[nu] - NS + Nonlinear1 + Nonlinear2 + Nonlinear3) / tau_rho);
            if (MD_RUN(P, profile1)) {
                sourceTerms[nu - 1] = (-q[nu] - NS) / tau_rho;
            }
            sourceTerms[nu - 1] += (G.theta - grid_pt.c.u[0] / tau) * q[nu];
            Real tempf =
                ((P.g[nu][0] + grid_pt.c.u[nu] * grid_pt.c.u[0])
                     * grid_pt.c.u[3] * q[3] / tau
                 - (P.g[nu][3] + grid_pt.c.u[nu] * grid_pt.c.u[3])
                       * grid_pt.c.u[3] * q[0] / tau);
            sourceTerms[nu - 1] += tempf;
            tempf = Real(0.0);
            for (int ii = 0; ii < 4; ii++) {
                const Real gmn = (ii == 0 ? Real(-1.) : Real(1.));
                tempf += q[ii] * gmn * G.a[ii];
            }
            sourceTerms[nu - 1] += (grid_pt.c.u[nu]) * tempf;
        }
    }

    //! Diss::Make_uPiChemSource (without the messages)
    MD_FN Real uPiChemSource(
        MD_THREAD const VCell<Real> &grid_pt,
        MD_THREAD const Gradients<Real> &G,
        MD_THREAD const Thermal<Real> &th) const {
        MD_CONST const StepParams<Real> &P = *P_;
        const Real small_eps = Real(1e-16);
        const Real epsilon = th.e;
        const Real pressure = th.p;
        const Real cs2 = th.cs2;
        const Real temperature = th.T;
        const Real Iqcd = (epsilon - Real(3.0) * pressure);
        const Real Iqcd_safe = copysign(dmax(fabs(Iqcd), small_eps), Iqcd);
        const Real y_q_inner = Real(1.0) - Real(3.0) * grid_pt.piChem / Iqcd_safe;
        const Real y_q_sqrt = dmin(Real(1.0), fabs(y_q_inner));
        const Real chem_rate_C = dmax(P.chemRateC, small_eps);
        Real tau_Pi_chem =
            Real(2.0)
            / dmax(
                chem_rate_C * temperature * (Real(1.0) + y_q_sqrt)
                    * (sqrt(Real(2.0)) + y_q_sqrt),
                small_eps);
        tau_Pi_chem =
            dmin(Real(10.0), dmax(Real(3.0) * P.dtau, tau_Pi_chem));
        const Real dI3de = (Real(1.0) - Real(3.0) * cs2);
        const Real zeta_raw = tau_Pi_chem * (Real(1.0) / Real(3.0))
                              * (Real(1.0) - y_q_sqrt) * dI3de
                              * (epsilon + pressure);
        const Real zeta = dmax(Real(0.0), zeta_raw);
        const Real NS_term = -zeta * G.theta;
        const Real relax_term = -(grid_pt.piChem);
        return (NS_term + relax_term) / tau_Pi_chem;
    }

    //! Diss::computeInverseReynoldsNumbers
    MD_FN void inverseReynolds(
        Real enthalpy, MD_THREAD const VCell<Real> &grid_pt,
        MD_THREAD Real &R_shear, MD_THREAD Real &R_bulk) const {
        Real pi_00 = grid_pt.W[0];
        Real pi_01 = grid_pt.W[1];
        Real pi_02 = grid_pt.W[2];
        Real pi_03 = grid_pt.W[3];
        Real pi_11 = grid_pt.W[4];
        Real pi_12 = grid_pt.W[5];
        Real pi_13 = grid_pt.W[6];
        Real pi_22 = grid_pt.W[7];
        Real pi_23 = grid_pt.W[8];
        Real pi_33 = grid_pt.W[9];
        Real pisize =
            (pi_00 * pi_00 + pi_11 * pi_11 + pi_22 * pi_22 + pi_33 * pi_33
             - Real(2.) * (pi_01 * pi_01 + pi_02 * pi_02 + pi_03 * pi_03)
             + Real(2.) * (pi_12 * pi_12 + pi_13 * pi_13 + pi_23 * pi_23));
        const Real h = dmax(enthalpy, tinyNormal<Real>());
        R_shear = sqrt(pisize) / h;
        R_bulk = fabs(grid_pt.pi) / h;
    }

    //! Diss::Make_uWSource
    MD_FN void uWSource(
        MD_THREAD const VCell<Real> &grid_pt,
        MD_THREAD const Gradients<Real> &G, MD_THREAD const Thermal<Real> &th,
        MD_THREAD Real *sourceTerms) const {
        MD_CONST const StepParams<Real> &P = *P_;
        const Transport<Real> tc = {P_};
        Real sigma[4][4];
        Real Wmunu[4][4];
        unpackSym(G.sigma, sigma);
        unpackSym(grid_pt.W, Wmunu);
        const Real theta_local = G.theta;
        Real epsilon = th.e;
        Real T = th.T;
        Real muB = th.muB;
        Real shear_to_s = tc.etaOverS(T, muB);
        bool include_WWterm = false;
        bool include_Wsigma_term = false;
        if (MD_RUN(P, secondOrder) == 1) {
            include_WWterm = true;
            include_Wsigma_term = true;
        }
        Real pressure = th.p;
        Real shear = Real(0.);
        if (MD_RUN(P, muBShear) == 0) {
            Real entropy = th.s;
            shear = shear_to_s * entropy;
        } else {
            shear = shear_to_s * (epsilon + pressure)
                    / dmax(T, Real(1e-16));
        }
        Real tau_pi =
            (P.shearRelaxFactor * shear
             / dmax(epsilon + pressure, Real(1e-16)));
        tau_pi = dmin(Real(10.), dmax(Real(3.) * P.dtau, tau_pi));
        Real transport_coefficient =
            (Real(9. / 70.) * P.shearRelaxFactor
             / dmax(epsilon + pressure, Real(1e-16)) * Real(4. / 5.));
        Real transport_coefficient2 = (Real(4. / 3.) * tau_pi);
        Real transport_coefficient3 = (Real(10. / 7.) * tau_pi);
        Real transport_coefficient_b = (Real(6. / 5.) * tau_pi);
        Real transport_coefficient2_b = Real(0.);
        Real resummedCorrection = Real(1);
        if (MD_RUN(P, resum)) {
            Real R_shear = Real(0.);
            Real R_bulk = Real(0.);
            inverseReynolds(epsilon + pressure, grid_pt, R_shear, R_bulk);
            resummedCorrection = (tc.resummedCorrFactor(R_shear, R_bulk));
            Real resummedCorrSq = resummedCorrection * resummedCorrection;
            shear *= resummedCorrection;
            transport_coefficient *= resummedCorrection;
            transport_coefficient2 *= resummedCorrSq;
            transport_coefficient3 *= resummedCorrSq;
            transport_coefficient_b *= resummedCorrSq;
            transport_coefficient2_b *= resummedCorrection;
        }
        int mu = 0;
        int nu = 0;
        for (int idx_1d = 4; idx_1d < 9; idx_1d++) {
            wIdx2d(idx_1d, mu, nu);
            Real tempf =
                (-(Real(1.0) + transport_coefficient2 * theta_local)
                 * Wmunu[mu][nu]);
            Real NS_term = Real(-2.) * shear * sigma[mu][nu];
            Real Coupling_to_Bulk = Real(0.0);
            if (MD_RUN(P, secondOrder) == 1) {
                Real Bulk_Sigma = grid_pt.pi * sigma[mu][nu];
                Real Bulk_W = grid_pt.pi * Wmunu[mu][nu];
                Real Bulk_Sigma_term = Bulk_Sigma * transport_coefficient_b;
                Real Bulk_W_term = Bulk_W * transport_coefficient2_b;
                Coupling_to_Bulk = -Bulk_Sigma_term + Bulk_W_term;
            }
            Real Vorticity_term = Real(0.0);
            if (MD_RUN(P, vorticity)) {
                Real transport_coefficient4 = Real(2.) * tau_pi;
                if (MD_RUN(P, resum)) {
                    transport_coefficient4 *= resummedCorrection;
                }
                Real omega[4][4];
                unpackAnti(G.omega, omega);
                Real term1_Vorticity =
                    (-Wmunu[mu][0] * omega[nu][0] - Wmunu[nu][0] * omega[mu][0]
                     + Wmunu[mu][1] * omega[nu][1] + Wmunu[nu][1] * omega[mu][1]
                     + Wmunu[mu][2] * omega[nu][2] + Wmunu[nu][2] * omega[mu][2]
                     + Wmunu[mu][3] * omega[nu][3]
                     + Wmunu[nu][3] * omega[mu][3])
                    / Real(2.);
                Vorticity_term = transport_coefficient4 * term1_Vorticity;
            }
            Real Wsigma_term = Real(0.0);
            if (include_Wsigma_term) {
                Real Wsigma =
                    (Wmunu[0][0] * sigma[0][0] + Wmunu[1][1] * sigma[1][1]
                     + Wmunu[2][2] * sigma[2][2] + Wmunu[3][3] * sigma[3][3]
                     - Real(2.)
                           * (Wmunu[0][1] * sigma[0][1]
                              + Wmunu[0][2] * sigma[0][2]
                              + Wmunu[0][3] * sigma[0][3])
                     + Real(2.)
                           * (Wmunu[1][2] * sigma[1][2]
                              + Wmunu[1][3] * sigma[1][3]
                              + Wmunu[2][3] * sigma[2][3]));
                Real term1_Wsigma =
                    (-Wmunu[mu][0] * sigma[nu][0] - Wmunu[nu][0] * sigma[mu][0]
                     + Wmunu[mu][1] * sigma[nu][1] + Wmunu[nu][1] * sigma[mu][1]
                     + Wmunu[mu][2] * sigma[nu][2] + Wmunu[nu][2] * sigma[mu][2]
                     + Wmunu[mu][3] * sigma[nu][3]
                     + Wmunu[nu][3] * sigma[mu][3])
                    / Real(2.);
                Real term2_Wsigma =
                    (-(Real(1.) / Real(3.))
                     * (P.g[mu][nu] + grid_pt.c.u[mu] * grid_pt.c.u[nu])
                     * Wsigma);
                term1_Wsigma = transport_coefficient3 * term1_Wsigma;
                term2_Wsigma = transport_coefficient3 * term2_Wsigma;
                Wsigma_term = -term1_Wsigma - term2_Wsigma;
            }
            Real WW_term = Real(0.0);
            if (include_WWterm) {
                Real Wsquare =
                    (Wmunu[0][0] * Wmunu[0][0] + Wmunu[1][1] * Wmunu[1][1]
                     + Wmunu[2][2] * Wmunu[2][2] + Wmunu[3][3] * Wmunu[3][3]
                     - Real(2.)
                           * (Wmunu[0][1] * Wmunu[0][1]
                              + Wmunu[0][2] * Wmunu[0][2]
                              + Wmunu[0][3] * Wmunu[0][3])
                     + Real(2.)
                           * (Wmunu[1][2] * Wmunu[1][2]
                              + Wmunu[1][3] * Wmunu[1][3]
                              + Wmunu[2][3] * Wmunu[2][3]));
                Real term1_WW =
                    (-Wmunu[mu][0] * Wmunu[nu][0] + Wmunu[mu][1] * Wmunu[nu][1]
                     + Wmunu[mu][2] * Wmunu[nu][2]
                     + Wmunu[mu][3] * Wmunu[nu][3]);
                Real term2_WW =
                    (-(Real(1.) / Real(3.))
                     * (P.g[mu][nu] + grid_pt.c.u[mu] * grid_pt.c.u[nu])
                     * Wsquare);
                term1_WW = term1_WW * transport_coefficient;
                term2_WW = term2_WW * transport_coefficient;
                WW_term = -term1_WW - term2_WW;
            }
            sourceTerms[idx_1d - 4] =
                ((NS_term + tempf + Vorticity_term + Wsigma_term + WW_term
                  + Coupling_to_Bulk)
                 / tau_pi);
        }
    }

    //! Diss::Make_uPiSource
    MD_FN Real uPiSource(
        MD_THREAD const VCell<Real> &grid_pt,
        MD_THREAD const Gradients<Real> &G,
        MD_THREAD const Thermal<Real> &th) const {
        MD_CONST const StepParams<Real> &P = *P_;
        const Transport<Real> tc = {P_};
        const Real theta_local = G.theta;
        Real tempf;
        Real bulk;
        Real Bulk_Relax_time = Real(0.2);
        Real transport_coeff1, transport_coeff2;
        Real transport_coeff1_s, transport_coeff2_s;
        Real NS_term, BB_term;
        Real Final_Answer;
        bool include_BBterm = false;
        bool include_coupling_to_shear = false;
        if (MD_RUN(P, secondOrder) == 1) {
            include_BBterm = true;
            include_coupling_to_shear = true;
        }
        Real epsilon = th.e;
        Real temperature = th.T;
        Real mu_B = th.muB;
        Real cs2 = th.cs2;
        Real pressure = th.p;
        bulk = tc.zetaOverS(temperature, mu_B);
        bulk = bulk * (epsilon + pressure) / dmax(temperature, Real(1e-16));
        if (MD_RUN(P, bulkRelaxType) == 0) {
            Real csfactor = dmax(Real(1.) / Real(3.) - cs2, Real(1e-16));
            Bulk_Relax_time =
                (P.bulkRelaxFactor * bulk / (csfactor * csfactor)
                 / dmax(epsilon + pressure, Real(1e-16)));
        } else if (MD_RUN(P, bulkRelaxType) == 1) {
            Real csfactor = dmax(Real(1.) / Real(3.) - cs2, Real(1e-16));
            Bulk_Relax_time =
                (P.bulkRelaxFactor * bulk / csfactor
                 / dmax(epsilon + pressure, Real(1e-16)));
        } else if (MD_RUN(P, bulkRelaxType) == 2) {
            Bulk_Relax_time =
                (P.bulkRelaxFactor * bulk
                 / dmax(epsilon + pressure, Real(1e-16)));
        }
        Bulk_Relax_time =
            (dmin(Real(10.), dmax(Real(3.) * P.dtau, Bulk_Relax_time)));
        transport_coeff1 = (Real(2. / 3.) * Bulk_Relax_time);
        transport_coeff2 = (Real(0.) * Bulk_Relax_time);
        transport_coeff1_s =
            (Real(8. / 5.) * (Real(1.) / Real(3.) - cs2) * Bulk_Relax_time);
        if (MD_RUN(P, bulkRelaxType) == 2) {
            transport_coeff1_s =
                (Real(8. / 5.) * Real(1.) / Real(3.) * Bulk_Relax_time);
        }
        transport_coeff2_s = Real(0.);
        if (MD_RUN(P, resum)) {
            Real resummedCorrection = Real(1);
            Real R_shear = Real(0.);
            Real R_bulk = Real(0.);
            inverseReynolds(epsilon + pressure, grid_pt, R_shear, R_bulk);
            resummedCorrection = (tc.resummedCorrFactor(R_shear, R_bulk));
            bulk *= resummedCorrection;
            transport_coeff1 *= resummedCorrection;
            transport_coeff2 *= resummedCorrection;
            transport_coeff1_s *= resummedCorrection;
            transport_coeff2_s *= resummedCorrection;
        }
        NS_term = -bulk * theta_local;
        tempf = (-(grid_pt.pi) - transport_coeff1 * theta_local * (grid_pt.pi));
        if (include_BBterm == 1) {
            BB_term = transport_coeff2 * (grid_pt.pi) * (grid_pt.pi);
        } else {
            BB_term = Real(0.0);
        }
        Real Wsigma, WW, Shear_Sigma_term, Shear_Shear_term, Coupling_to_Shear;
        if (include_coupling_to_shear) {
            Real sigma[4][4];
            Real Wmunu[4][4];
            unpackSym(G.sigma, sigma);
            unpackSym(grid_pt.W, Wmunu);
            Wsigma =
                (Wmunu[0][0] * sigma[0][0] + Wmunu[1][1] * sigma[1][1]
                 + Wmunu[2][2] * sigma[2][2] + Wmunu[3][3] * sigma[3][3]
                 - Real(2.)
                       * (Wmunu[0][1] * sigma[0][1] + Wmunu[0][2] * sigma[0][2]
                          + Wmunu[0][3] * sigma[0][3])
                 + Real(2.)
                       * (Wmunu[1][2] * sigma[1][2] + Wmunu[1][3] * sigma[1][3]
                          + Wmunu[2][3] * sigma[2][3]));
            WW =
                (Wmunu[0][0] * Wmunu[0][0] + Wmunu[1][1] * Wmunu[1][1]
                 + Wmunu[2][2] * Wmunu[2][2] + Wmunu[3][3] * Wmunu[3][3]
                 - Real(2.)
                       * (Wmunu[0][1] * Wmunu[0][1] + Wmunu[0][2] * Wmunu[0][2]
                          + Wmunu[0][3] * Wmunu[0][3])
                 + Real(2.)
                       * (Wmunu[1][2] * Wmunu[1][2] + Wmunu[1][3] * Wmunu[1][3]
                          + Wmunu[2][3] * Wmunu[2][3]));
            Shear_Sigma_term = Wsigma * transport_coeff1_s;
            Shear_Shear_term = WW * transport_coeff2_s;
            Coupling_to_Shear = -Shear_Sigma_term + Shear_Shear_term;
        } else {
            Coupling_to_Shear = Real(0.0);
        }
        Final_Answer = NS_term + tempf + BB_term + Coupling_to_Shear;
        return Final_Answer / (Bulk_Relax_time);
    }

    //! the factor of Advance::QuestRevert and QuestRevert_qmu
    MD_FN Real questFactor(Real e_local) const {
        MD_CONST const StepParams<Real> &P = *P_;
        Real eps_scale = Real(0.1);
        Real xi = Real(0.05);
        Real factor;
        if (sizeof(Real) == 8) {
            factor = Real(10.) * P.questRevertStrength
                     * (Real(1.) / (exp(-(e_local - eps_scale) / xi) + Real(1.))
                        - Real(1.) / (exp(eps_scale / xi) + Real(1.)));
        } else {
            // the same difference of two logistic functions without the
            // cancellation of the form above, which in dilute cells (both
            // terms ~ 0.119, their difference ~ 2 e) loses log10(0.06 / e)
            // digits: all of float's at e ~ 1e-6
            const Real em = exp(-(e_local - eps_scale) / xi);
            const Real e0 = exp(eps_scale / xi);
            factor = Real(10.) * P.questRevertStrength
                     * (-e0 * dexpm1(-e_local / xi))
                     / ((em + Real(1.)) * (e0 + Real(1.)));
        }
        return (factor);
    }

    //! Advance::QuestRevert_qmu (without the messages)
    MD_FN void questRevertQmu(MD_THREAD VCell<Real> &grid_pt) const {
        Real factor = questFactor(grid_pt.c.e);
        Real q_mu_local[4];
        for (int ii = 0; ii < 4; ii++) q_mu_local[ii] = grid_pt.W[10 + ii];
        Real q_size = Real(0.0);
        for (int ii = 0; ii < 4; ii++) {
            Real gfac = (ii == 0 ? Real(-1.0) : Real(1.0));
            q_size += gfac * q_mu_local[ii] * q_mu_local[ii];
        }
        if (q_size < Real(0.0)) {
            for (int ii = 0; ii < 4; ii++) grid_pt.W[10 + ii] = Real(0.0);
        }
        Real rhob_local = grid_pt.c.rhob;
        Real rho_q =
            sqrt(q_size / dmax(rhob_local * rhob_local, tinyNormal<Real>()))
            / dmax(factor, tinyNormal<Real>());
        Real rho_q_max = Real(0.1);
        if (rho_q > rho_q_max) {
            for (int ii = 0; ii < 4; ii++) {
                grid_pt.W[10 + ii] = (rho_q_max / rho_q) * q_mu_local[ii];
            }
        }
    }

    //! Advance::QuestRevert (without the messages)
    MD_FN void questRevert(MD_THREAD VCell<Real> &grid_pt) const {
        Real e_local = grid_pt.c.e;
        Real rhob = grid_pt.c.rhob;
        Real factor = questFactor(e_local);
        Real factor_bulk = factor;
        Real pi_00 = grid_pt.W[0];
        Real pi_01 = grid_pt.W[1];
        Real pi_02 = grid_pt.W[2];
        Real pi_03 = grid_pt.W[3];
        Real pi_11 = grid_pt.W[4];
        Real pi_12 = grid_pt.W[5];
        Real pi_13 = grid_pt.W[6];
        Real pi_22 = grid_pt.W[7];
        Real pi_23 = grid_pt.W[8];
        Real pi_33 = grid_pt.W[9];
        Real pisize =
            (pi_00 * pi_00 + pi_11 * pi_11 + pi_22 * pi_22 + pi_33 * pi_33
             - Real(2.) * (pi_01 * pi_01 + pi_02 * pi_02 + pi_03 * pi_03)
             + Real(2.) * (pi_12 * pi_12 + pi_13 * pi_13 + pi_23 * pi_23));
        Real pi_local = grid_pt.pi;
        Real bulksize = Real(3.) * pi_local * pi_local;
        Real p_local = eosPressure<Real, EosType>(
            *d_, tab, e_local, rhob, grid_pt.c.rhoq, grid_pt.c.rhos);
        Real eq_size = dmax(
            e_local * e_local + Real(3.) * p_local * p_local,
            tinyNormal<Real>());
        factor = dmax(factor, tinyNormal<Real>());
        factor_bulk = factor;
        Real rho_shear = sqrt(pisize / eq_size) / factor;
        Real rho_bulk = sqrt(bulksize / eq_size) / factor_bulk;
        Real rho_shear_max = Real(0.1);
        if (isnan(rho_shear)) {
            for (int mu = 0; mu < 10; mu++) {
                grid_pt.W[mu] = Real(0.0);
            }
        } else if (rho_shear > rho_shear_max) {
            for (int mu = 0; mu < 10; mu++) {
                grid_pt.W[mu] = (rho_shear_max / rho_shear) * grid_pt.W[mu];
            }
        }
        Real rho_bulk_max = Real(0.1);
        if (rho_bulk > rho_bulk_max) {
            grid_pt.pi = (rho_bulk_max / rho_bulk) * grid_pt.pi;
        }
        Real rho_bulkChem =
            (sqrt(Real(3.) * grid_pt.piChem * grid_pt.piChem / eq_size)
             / factor_bulk);
        Real rho_bulkChem_max = Real(0.1);
        if (rho_bulkChem > rho_bulkChem_max) {
            grid_pt.piChem = (rho_bulkChem_max / rho_bulkChem) * grid_pt.piChem;
        }
    }

    //! Advance::solveEigenvaluesWmunu: Lambdas = (min, -min - max, max) of
    //! the real eigenvalues of A = W^mu_nu, in closed form: A u = 0
    //! (transversality), so the eigenvalues are 0 and the three roots of
    //! the cubic with the power sums tr A, tr A^2, tr A^3. They are real: on
    //! the space orthogonal to u, A is the spatial W of the local rest
    //! frame, a symmetric matrix. The host code does the same operations.
    MD_FN void wEigenvalues(
        MD_THREAD const VCell<Real> &g, MD_THREAD Real *Lambdas) const {
        Real A[4][4];
        A[0][0] = -g.W[0];
        A[0][1] = g.W[1];
        A[0][2] = g.W[2];
        A[0][3] = g.W[3];
        A[1][0] = -g.W[1];
        A[1][1] = g.W[4];
        A[1][2] = g.W[5];
        A[1][3] = g.W[6];
        A[2][0] = -g.W[2];
        A[2][1] = g.W[5];
        A[2][2] = g.W[7];
        A[2][3] = g.W[8];
        A[3][0] = -g.W[3];
        A[3][1] = g.W[6];
        A[3][2] = g.W[8];
        A[3][3] = g.W[9];
        Real t1 = Real(0.), t2 = Real(0.), t3 = Real(0.);
        for (int i = 0; i < 4; i++) {
            t1 += A[i][i];
            for (int j = 0; j < 4; j++) {
                t2 += A[i][j] * A[j][i];
                Real a2 = Real(0.);
                for (int k = 0; k < 4; k++) a2 += A[i][k] * A[k][j];
                t3 += a2 * A[j][i];
            }
        }
        // roots m + y_k: sum y = 0, sum y^2 = 6 q, prod y = 2 r
        const Real m = t1 / Real(3.);
        const Real q = (t2 - Real(3.) * m * m) / Real(6.);
        const Real sumY3 = t3 - Real(3.) * m * t2 + Real(6.) * m * m * m;
        const Real r = sumY3 / Real(6.);
        Real roots[3] = {m, m, m};
        if (q > Real(0.)) {
            const Real sq = sqrt(q);
            // q * sq underflows to 0 only for (numerically) equal roots
            const Real qsq = q * sq;
            const Real c =
                (qsq > Real(0.) ? dmax(Real(-1.), dmin(Real(1.), r / qsq))
                                : Real(0.));
            const Real phi = acos(c) / Real(3.);
            const Real twoPiOver3 = Real(2.0943951023931954923);
            roots[0] = m + Real(2.) * sq * cos(phi);
            roots[1] = m + Real(2.) * sq * cos(phi + twoPiOver3);
            roots[2] = m + Real(2.) * sq * cos(phi + Real(2.) * twoPiOver3);
        }
        Real lmin = Real(0.), lmax = Real(0.);
        for (int k = 0; k < 3; k++) {
            lmin = dmin(lmin, roots[k]);
            lmax = dmax(lmax, roots[k]);
        }
        Lambdas[0] = lmin;
        Lambdas[1] = -lmin - lmax;
        Lambdas[2] = lmax;
    }

    //! Advance::nCausalityConstraints
    MD_FN void causalityConstraints(
        MD_THREAD VCell<Real> &grid_pt, MD_THREAD const Real *Lambdas) const {
        Real eps = grid_pt.c.e;
        Real rhob = grid_pt.c.rhob;
        Real cs2 = eosCs2<Real, EosType>(
            *d_, tab, eps, rhob, grid_pt.c.rhoq, grid_pt.c.rhos);
        Real Pr = eosPressure<Real, EosType>(
            *d_, tab, eps, rhob, grid_pt.c.rhoq, grid_pt.c.rhos);
        // TransportCoeffs: the second-order coefficients
        const Real tau_pipi = Real(10.) / Real(7.);
        const Real delta_pipi = Real(4.) / Real(3.);
        const Real lambda_pibulkPi = Real(6.) / Real(5.);
        const Real lambda_bulkPipi = Real(8.) / Real(5.);
        const Real delta_bulkPibulkPi = Real(2.) / Real(3.);
        Real transportPart_n13 =
            Real(2.) * Real(1.) / P_->shearRelaxFactor;
        Real viscousPart1_n13 = lambda_pibulkPi;
        Real viscousPart2_n13 = -Real(1.) / Real(2.) * tau_pipi;
        Real transportPart_n56 =
            cs2 + Real(4.) / Real(3.) * Real(1.) / P_->shearRelaxFactor
            + Real(1.) / P_->bulkRelaxFactor * (Real(1.) / Real(3.) - cs2)
                  * (Real(1.) / Real(3.) - cs2);
        Real viscousPart1_n56 =
            Real(2.) / Real(3.) * lambda_pibulkPi + delta_bulkPibulkPi + cs2;
        Real viscousPart2_n56 =
            delta_pipi + Real(1.) / Real(3.) * tau_pipi
            + lambda_bulkPipi * (Real(1.) / Real(3.) - cs2) + cs2;

        const Real enthalpy = dmax(eps + Pr, tinyNormal<Real>());
        Real n1 = transportPart_n13
                  + viscousPart1_n13 * grid_pt.pi / enthalpy
                  + viscousPart2_n13 * fabs(Lambdas[0]) / enthalpy;
        Real n3 = transportPart_n13
                  + viscousPart1_n13 * grid_pt.pi / enthalpy
                  + viscousPart2_n13 * Lambdas[2] / enthalpy;
        Real n5 = transportPart_n56
                  + viscousPart1_n56 * grid_pt.pi / enthalpy
                  + viscousPart2_n56 * Lambdas[0] / enthalpy;
        Real n6 = (Real(1.) - transportPart_n56)
                  + (Real(1.) - viscousPart1_n56) * grid_pt.pi / enthalpy
                  + (Real(1.) - viscousPart2_n56) * Lambdas[2] / enthalpy;

        Real minAlp = Real(1.);
        Real alp = Real(1.);
        if (n1 < Real(0.)) {
            alp = -transportPart_n13 / (n1 - transportPart_n13);
            minAlp = dmax(Real(0.), dmin(minAlp, alp));
        }
        if (n3 < Real(0.)) {
            alp = -transportPart_n13 / (n3 - transportPart_n13);
            minAlp = dmax(Real(0.), dmin(minAlp, alp));
        }
        if (n5 < Real(0.)) {
            alp = -transportPart_n56 / (n5 - transportPart_n56);
            minAlp = dmax(Real(0.), dmin(minAlp, alp));
        }
        if (n6 < Real(0.)) {
            const Real den6 = n6 - (Real(1.) - transportPart_n56);
            alp = (den6 < Real(0.) ? -(Real(1.) - transportPart_n56) / den6
                                   : Real(0.));
            minAlp = dmax(Real(0.), dmin(minAlp, alp));
        }
        for (int k = 0; k < 14; k++) grid_pt.W[k] = grid_pt.W[k] * minAlp;
        grid_pt.pi = grid_pt.pi * minAlp;
    }

    //! Advance::QuestRevertResummedTransCoeff
    MD_FN void questRevertResummed(MD_THREAD VCell<Real> &grid_pt) const {
        Real e_local = grid_pt.c.e;
        Real rhob = grid_pt.c.rhob;
        Real p_local = eosPressure<Real, EosType>(
            *d_, tab, e_local, rhob, grid_pt.c.rhoq, grid_pt.c.rhos);
        Real pi_00 = grid_pt.W[0];
        Real pi_01 = grid_pt.W[1];
        Real pi_02 = grid_pt.W[2];
        Real pi_03 = grid_pt.W[3];
        Real pi_11 = grid_pt.W[4];
        Real pi_12 = grid_pt.W[5];
        Real pi_13 = grid_pt.W[6];
        Real pi_22 = grid_pt.W[7];
        Real pi_23 = grid_pt.W[8];
        Real pi_33 = grid_pt.W[9];
        Real pisize =
            (pi_00 * pi_00 + pi_11 * pi_11 + pi_22 * pi_22 + pi_33 * pi_33
             - Real(2.) * (pi_01 * pi_01 + pi_02 * pi_02 + pi_03 * pi_03)
             + Real(2.) * (pi_12 * pi_12 + pi_13 * pi_13 + pi_23 * pi_23));
        Real pi_local = grid_pt.pi;
        Real Rshear =
            sqrt(pisize) / (dmax(e_local + p_local, Real(0.)) + Real(1e-16));
        Real Rbulk =
            fabs(pi_local) / (dmax(e_local + p_local, Real(0.)) + Real(1e-16));
        Real r_combined = Real(1.5) * (Rshear + Rbulk) + Real(1e-16);
        Real renorm = dmin(Real(1.), (Real(1.) - Real(1e-6)) / r_combined);
        for (int mu = 0; mu < 10; mu++) {
            grid_pt.W[mu] = renorm * grid_pt.W[mu];
        }
        grid_pt.pi = renorm * grid_pt.pi;
    }

    //! Advance::FirstRKStepW: grid_f holds the primitives of the next
    //! slice; its W and Pi are overwritten
    MD_FN void firstRKStepW(
        MD_THREAD const VCell<Real> &grid_c, MD_THREAD const VCell<Real> &prev,
        MD_THREAD const Gradients<Real> &G, MD_THREAD const Thermal<Real> &th,
        MD_THREAD VCell<Real> &grid_f) const {
        MD_CONST const StepParams<Real> &P = *P_;
        const Real rk_flag = Real(P.rkFlag);
        Real tempf;
        Real u0Prev = prev.c.u[0];
        if (MD_RUN(P, shear) == 1) {
            Real sourceTerms[5];
            uWSource(grid_c, G, th, sourceTerms);
            for (int idx_1d = 4; idx_1d < 9; idx_1d++) {
                Real WmunuPrev = prev.W[idx_1d];
                tempf =
                    ((Real(1.) - rk_flag) * (grid_c.W[idx_1d] * grid_c.c.u[0])
                     + rk_flag * (WmunuPrev * u0Prev));
                tempf += sourceTerms[idx_1d - 4] * (P.dtau);
                tempf += G.wRhs[idx_1d - 4];
                tempf += rk_flag * ((grid_c.W[idx_1d]) * (grid_c.c.u[0]));
                tempf *= Real(1.) / (Real(1.) + rk_flag);
                grid_f.W[idx_1d] = tempf / (grid_f.c.u[0]);
            }
        } else {
            for (int idx_1d = 4; idx_1d < 9; idx_1d++) {
                grid_f.W[idx_1d] = Real(0.0);
            }
        }
        if (MD_RUN(P, bulk) == 1) {
            Real piBulkPrev = prev.pi;
            tempf =
                ((Real(1.) - rk_flag) * (grid_c.pi * grid_c.c.u[0])
                 + rk_flag * piBulkPrev * u0Prev);
            Real temps = uPiSource(grid_c, G, th);
            tempf += temps * (P.dtau);
            tempf += G.wRhs[5];
            tempf += rk_flag * ((grid_c.pi) * (grid_c.c.u[0]));
            tempf *= Real(1.) / (Real(1.) + rk_flag);
            grid_f.pi = tempf / (grid_f.c.u[0]);
        } else {
            grid_f.pi = Real(0.0);
            grid_f.piChem = Real(0.0);
        }
        if (MD_RUN(P, bulkChem) == 1) {
            Real piBulkChemPrev_local = prev.piChem;
            tempf =
                ((Real(1.) - rk_flag) * (grid_c.piChem * grid_c.c.u[0])
                 + rk_flag * piBulkChemPrev_local * u0Prev);
            Real temps_chem = uPiChemSource(grid_c, G, th);
            tempf += temps_chem * (P.dtau);
            tempf += G.wRhs[6];
            tempf += rk_flag * ((grid_c.piChem) * (grid_c.c.u[0]));
            tempf *= Real(1.) / (Real(1.) + rk_flag);
            grid_f.piChem = tempf / (grid_f.c.u[0]);
        }
        if (MD_RUN(P, diff) == 1) {
            const Real tau_now = P.tau + rk_flag * P.dtau;
            Real sourceTerms[3] = {Real(0.), Real(0.), Real(0.)};
            uqSource(tau_now, grid_c, G, th, sourceTerms);
            for (int idx_1d = 11; idx_1d < 14; idx_1d++) {
                Real WmunuPrev = prev.W[idx_1d];
                tempf =
                    ((Real(1.) - rk_flag) * (grid_c.W[idx_1d] * grid_c.c.u[0])
                     + rk_flag * WmunuPrev * u0Prev);
                tempf += sourceTerms[idx_1d - 11] * (P.dtau);
                tempf += G.wRhs[idx_1d - 4];
                tempf += rk_flag * (grid_c.W[idx_1d] * grid_c.c.u[0]);
                tempf *= Real(1.) / (Real(1.) + rk_flag);
                grid_f.W[idx_1d] = tempf / (grid_f.c.u[0]);
            }
        } else {
            for (int idx_1d = 10; idx_1d < 14; idx_1d++) {
                grid_f.W[idx_1d] = Real(0.0);
            }
        }
        Real u0sq = grid_f.c.u[0] * grid_f.c.u[0];
        grid_f.W[9] =
            ((Real(2.)
                  * (grid_f.c.u[1] * grid_f.c.u[2] * grid_f.W[5]
                     + grid_f.c.u[1] * grid_f.c.u[3] * grid_f.W[6]
                     + grid_f.c.u[2] * grid_f.c.u[3] * grid_f.W[8])
              - (u0sq - grid_f.c.u[1] * grid_f.c.u[1]) * grid_f.W[4]
              - (u0sq - grid_f.c.u[2] * grid_f.c.u[2]) * grid_f.W[7])
             / dmax(u0sq - grid_f.c.u[3] * grid_f.c.u[3], Real(0.5)));
        for (int mu = 1; mu < 4; mu++) {
            tempf = Real(0.0);
            for (int nu = 1; nu < 4; nu++) {
                int idx_1d = wIdx(mu, nu);
                tempf += grid_f.W[idx_1d] * grid_f.c.u[nu];
            }
            grid_f.W[mu] = tempf / (grid_f.c.u[0]);
        }
        tempf = Real(0.0);
        for (int nu = 1; nu < 4; nu++) tempf += grid_f.W[nu] * grid_f.c.u[nu];
        grid_f.W[0] = tempf / (grid_f.c.u[0]);
        // q^tau from the transversality
        tempf = Real(0.0);
        for (int idx_1d = 11; idx_1d < 14; idx_1d++) {
            tempf += grid_f.W[idx_1d] * grid_f.c.u[idx_1d - 10];
        }
        grid_f.W[10] = Real(MD_RUN(P, diff)) * tempf / (grid_f.c.u[0]);
        if (MD_RUN(P, questRevert)) {
            questRevert(grid_f);
            if (MD_RUN(P, diff) == 1) questRevertQmu(grid_f);
            if (MD_RUN(P, resum)) questRevertResummed(grid_f);
            if (MD_RUN(P, causality)) {
                Real Lambdas[3];
                wEigenvalues(grid_f, Lambdas);
                causalityConstraints(grid_f, Lambdas);
            }
        }
    }

    //! Advance::viscousTimeComponent
    MD_FN Real viscousTimeComponent(
        MD_THREAD const VCell<Real> &v, int alpha) const {
        if (alpha == 4) return v.W[10];
        Real val = v.W[alpha];
        if (P_->bulk == 1 || P_->bulkChem == 1) {
            const Real g = (alpha == 0 ? Real(-1.) : Real(0.));
            const Real Pi = v.pi + v.piChem;
            val += Pi * (g + v.c.u[alpha] * v.c.u[0]);
        }
        return val;
    }
};

//! mu_B / T of cell i of slice s, as eos.get_muB / eos.get_temperature
//! (Advance::computeMuBOverT): the derivatives of the baryon diffusion
//! read it from a, so that the neighbours do not call the EoS again
template <typename Real>
MD_FN void kernelMuBOverT(
    MD_CONST const StepParams<Real> &P, MD_CONST const EosDesc<Real> &d,
    MD_GLOBAL const Real *tab, MD_GLOBAL const Real *s, MD_GLOBAL Real *a,
    int i) {
    const int n = P.nx * P.ny * P.neta;
    Real T, muB;
    eosTMuB<Real>(d, tab, s[i], s[n + i], s[2 * n + i], s[3 * n + i], T, muB);
    a[i] = muB / dmax(T, Real(1e-16));
}

//! the cell loop of Advance::AdvanceIt (viscosity_flag = 1) for cell i.
//! sCurr: the stage's slice; sPrev: arenaFieldsPrev (tau_{n-1} in the
//! first stage, tau_n in the second); pCurr, pPrev their pressures; aCurr,
//! aPrev their mu_B / T (kernelMuBOverT; read only with diffusion); flux*:
//! the ideal KT fluxes; wFace*: kernelWFace of the three directions; src:
//! the source terms (if P.hasSource). kStage: -1 reads the stage from
//! P.rkFlag (host, CUDA); 0 or 1 compiles one stage only (Metal, one pipeline
//! per stage), which removes the fixed-point code from the first stage's
//! kernel and its register cost (9% of that kernel with EoS 20). It gates
//! only the fixed-point branch; every arithmetic use of P.rkFlag stays a
//! run-time value, so the arithmetic is the same
template <typename Real, int kStage = -1, int EosType = -1>
MD_FN void kernelCellVisc(
    MD_CONST const StepParams<Real> &P, MD_CONST const EosDesc<Real> &d,
    MD_GLOBAL const Real *tab, MD_GLOBAL const Real *sCurr,
    MD_GLOBAL const Real *sPrev, MD_GLOBAL const Real *pCurr,
    MD_GLOBAL const Real *pPrev, MD_GLOBAL const Real *aCurr,
    MD_GLOBAL const Real *aPrev, MD_GLOBAL const Real *fluxX,
    MD_GLOBAL const Real *fluxY, MD_GLOBAL const Real *fluxEta,
    MD_GLOBAL const Real *wFaceX, MD_GLOBAL const Real *wFaceY,
    MD_GLOBAL const Real *wFaceEta, MD_GLOBAL const Real *src,
    MD_GLOBAL Real *sNext, int i, MD_FLAGS flags) {
    const int nx = P.nx, ny = P.ny, neta = P.neta;
    const int n = nx * ny * neta;
    const int ix = i % nx;
    const int iy = (i / nx) % ny;
    const int ieta = i / (nx * ny);
    const Real tau = P.tau + Real(P.rkFlag) * P.dtau;
    const bool withBulk = (MD_RUN(P, bulk) == 1 || MD_RUN(P, bulkChem) == 1);
    const VCell<Real> curr = loadVCell(sCurr, n, i);
    const VCell<Real> prev = loadVCell(sPrev, n, i);

    // ----- Diss::computeWDivergences for this cell ----------------------
    Real wSource[9];
    Real wFluxDiv[kNWFlux];
    for (int k = 0; k < kNWFlux; k++) wFluxDiv[k] = Real(0.);
    for (int alpha = 0; alpha < 5; alpha++) {
        const int idx_1d_alpha0 = wIdx(alpha, 0);
        Real dWdtau =
            ((curr.W[idx_1d_alpha0] - prev.W[idx_1d_alpha0]) / P.dtau);
        Real dPidtau = Real(0.0);
        Real Pi_alpha0 = Real(0.0);
        if (alpha < 4 && withBulk) {
            Real gfac = (alpha == 0 ? Real(-1.0) : Real(0.0));
            Real PiCurr = curr.pi + curr.piChem;
            Real PiPrev = prev.pi + prev.piChem;
            Pi_alpha0 = PiCurr * (gfac + curr.c.u[alpha] * curr.c.u[0]);
            dPidtau =
                (Pi_alpha0 - PiPrev * (gfac + prev.c.u[alpha] * prev.c.u[0]))
                / P.dtau;
        }
        Real sf = tau * dWdtau + curr.W[idx_1d_alpha0];
        Real bf = tau * dPidtau + Pi_alpha0;
        Real dd = Real(0.);
        dd += sf + bf;
        wSource[alpha] = dd;
    }
    for (int direction = 1; direction <= 3; direction++) {
        const int dx = (direction == 1), dy = (direction == 2),
                  deta = (direction == 3);
        const int nfx = nx + dx, nfy = ny + dy;
        const int nFaces = faceCount(nx, ny, neta, direction);
        const int Imh = ix + nfx * (iy + nfy * ieta);
        const int Iph = (ix + dx) + nfx * ((iy + dy) + nfy * (ieta + deta));
        MD_GLOBAL const Real *F =
            (direction == 1 ? wFaceX : (direction == 2 ? wFaceY : wFaceEta));
        const Real deltaFlux =
            (direction == 1 ? P.dx : (direction == 2 ? P.dy : P.deta * tau));
        const Real deltaSrc =
            (direction == 1 ? P.dx : (direction == 2 ? P.dy : P.deta));
        for (int idx = 0; idx < MD_RUN(P, nWComp); idx++) {
            const Real HW =
                (F[idx * nFaces + Iph] - F[idx * nFaces + Imh]) / deltaFlux;
            if (direction == 1) wFluxDiv[idx] = Real(0.);
            wFluxDiv[idx] += (-HW) * P.dtau;
        }
        for (int alpha = 0; alpha < 5; alpha++) {
            const Real W_p = F[(kNWFlux + alpha) * nFaces + Iph];
            const Real W_m = F[(kNWFlux + alpha) * nFaces + Imh];
            Real Pi_m = Real(0);
            Real Pi_p = Real(0);
            if (withBulk && alpha < 4) {
                Pi_p = F[(kNWFlux + 5 + alpha) * nFaces + Iph];
                Pi_m = F[(kNWFlux + 5 + alpha) * nFaces + Imh];
            }
            if (direction == 3 && (alpha == 0 || alpha == 3)) {
                const int k = (alpha == 0 ? 0 : 1);
                wSource[5 + k] = W_p + Pi_p;
                wSource[7 + k] = W_m + Pi_m;
            } else {
                wSource[alpha] += (W_p - W_m + Pi_p - Pi_m) / deltaSrc;
            }
        }
    }
    // Diss::MakeWSourceFromDivergence
    Real dwmn[7];
    for (int k = 0; k < 7; k++) dwmn[k] = Real(0.);
    for (int alpha = 0; alpha < 5; alpha++) dwmn[alpha] = wSource[alpha];
    {
        const Real W_eta_p0 = wSource[5];
        const Real W_eta_p3 = wSource[6];
        const Real W_eta_m0 = wSource[7];
        const Real W_eta_m3 = wSource[8];
        const Real cosh_deta = P.coshDeta;
        const Real sinh_deta = P.sinhDeta;
        dwmn[0] +=
            ((W_eta_p0 - W_eta_m0) * cosh_deta
             + (W_eta_p3 + W_eta_m3) * sinh_deta);
        dwmn[3] +=
            ((W_eta_p3 - W_eta_m3) * cosh_deta
             + (W_eta_p0 + W_eta_m0) * sinh_deta);
    }

    // ----- FirstRKStepT ---------------------------------------------------
    Real qi[7];
    cellQiFlux(P, sCurr, pCurr, fluxX, fluxY, fluxEta, i, qi);
    for (int alpha = 0; alpha < 7; alpha++) {
        qi[alpha] -= dwmn[alpha] * (P.dtau);
    }
    cellAddSource(P, src, i, qi);
    cellRkAverage(P, sPrev, pPrev, i, qi);
    const Solver<Real, EosType> solver = {&P, &d, tab};
    const Real pressure = pCurr[i];
    const Real tauNext = P.tau + P.dtau;
    int flagAll = 0;
    int flag = 0;
    VCell<Real> next;
    next.c = solver.reconstruct(
        tauNext, qi, curr.c, (MD_THREAD EosState<Real> *)0, pressure, flag);
    flagAll |= flag;
    next.piChem = (withBulk ? sNext[kOffPiChem * n + i] : Real(0.));
    for (int k = 0; k < 14; k++) next.W[k] = Real(0.);
    next.pi = Real(0.);
    Real qiUsed[7];
    for (int k = 0; k < 7; k++) qiUsed[k] = qi[k];
    int flagUsed = flag;

    // ----- the viscous update ---------------------------------------------
    // a frozen cell (Advance::frozenViscousUpdate): tau W, tau Pi and
    // tau Pi_chem are kept in the first stage, and the second keeps the
    // predictor
    if (isFrozenCell(P, sCurr, ix, iy, ieta)) {
        const Real s = P.tau / tauNext;
        for (int k = 0; k < 14; k++)
            next.W[k] = (P.rkFlag == 0 ? s * curr.W[k] : curr.W[k]);
        next.pi = (P.rkFlag == 0 ? s * curr.pi : curr.pi);
        next.piChem = (P.rkFlag == 0 ? s * curr.piChem : curr.piChem);
    } else {
        const ViscSolver<Real, EosType> vs = {&P, &d, tab};
        Gradients<Real> G;
        vs.gradients(tau, sCurr, sPrev, pCurr, aCurr, aPrev, i, wFluxDiv, G);
        const Thermal<Real> th = eosThermal<Real, EosType>(
            d, tab, curr.c.e, curr.c.rhob, curr.c.rhoq, curr.c.rhos);
        vs.firstRKStepW(curr, prev, G, th, next);

        if (kStage < 0 ? P.rkFlag == 1 : kStage == 1) {
            // Advance::stage2ViscousFixedPoint
            const Real tau_np2 = P.tau + Real(2.) * P.dtau;
            const Real scaleT = dmax(fabs(qi[0]), tinyNormal<Real>());
            const Real scaleJ = fabs(qi[4]) + Real(1e-12) * scaleT;
            const Real tolerance = Real(1e-8);
            const int maxIterations = 16;
            const Real bigNumber = realMax<Real>();
            Real shift[5] = {Real(0.), Real(0.), Real(0.), Real(0.), Real(0.)};
            Real best[5] = {Real(0.), Real(0.), Real(0.), Real(0.), Real(0.)};
            Real bestMiss = bigNumber;
            Real prevMiss = bigNumber;
            Real lastMiss = Real(0.);
            Real relaxation = Real(1.);
            bool accelerate = true;
            bool havePrev = false;
            Real prevShift[5] = {Real(0.), Real(0.), Real(0.), Real(0.),
                                 Real(0.)};
            Real prevResidual[5] = {Real(0.), Real(0.), Real(0.), Real(0.),
                                    Real(0.)};
            bool converged = false;
            for (int iter = 0; iter <= maxIterations; iter++) {
                Real target[5];
                Real miss = Real(0.);
                bool conv = true;
                for (int alpha = 0; alpha < 5; alpha++) {
                    target[alpha] =
                        Real(0.5) * tau_np2
                        * (vs.viscousTimeComponent(next, alpha)
                           - vs.viscousTimeComponent(curr, alpha));
                    const Real scale = (alpha < 4 ? scaleT : scaleJ);
                    const Real dd = fabs(target[alpha] - shift[alpha]);
                    miss += dd / scale;
                    if (dd > tolerance * scale) conv = false;
                }
                lastMiss = miss;
                if (conv) {
                    converged = true;
                    break;
                }
                if (miss < bestMiss) {
                    bestMiss = miss;
                    for (int k = 0; k < 5; k++) best[k] = shift[k];
                }
                if (iter == maxIterations) break;
                if (iter > 0 && miss > prevMiss) {
                    relaxation *= Real(0.5);
                    accelerate = false;
                }
                prevMiss = miss;
                Real residual[5];
                for (int alpha = 0; alpha < 5; alpha++) {
                    residual[alpha] = target[alpha] - shift[alpha];
                }
                Real nextShift[5];
                for (int k = 0; k < 5; k++) nextShift[k] = shift[k];
                bool mixed = false;
                if (accelerate && havePrev) {
                    Real num = Real(0.);
                    Real den = Real(0.);
                    for (int alpha = 0; alpha < 5; alpha++) {
                        const Real w = Real(1.) / (alpha < 4 ? scaleT : scaleJ);
                        const Real dr =
                            (residual[alpha] - prevResidual[alpha]) * w;
                        num += dr * residual[alpha] * w;
                        den += dr * dr;
                    }
                    if (den > Real(0.)) {
                        const Real gamma = num / den;
                        for (int alpha = 0; alpha < 5; alpha++) {
                            nextShift[alpha] =
                                shift[alpha] + residual[alpha]
                                - gamma
                                      * ((shift[alpha] - prevShift[alpha])
                                         + (residual[alpha]
                                            - prevResidual[alpha]));
                        }
                        mixed = true;
                    }
                }
                if (!mixed) {
                    for (int alpha = 0; alpha < 5; alpha++) {
                        nextShift[alpha] =
                            shift[alpha] + relaxation * residual[alpha];
                    }
                }
                for (int k = 0; k < 5; k++) {
                    prevShift[k] = shift[k];
                    prevResidual[k] = residual[k];
                }
                havePrev = true;
                for (int k = 0; k < 5; k++) shift[k] = nextShift[k];
                // applyShift
                for (int k = 0; k < 7; k++) qiUsed[k] = qi[k];
                for (int alpha = 0; alpha < 5; alpha++)
                    qiUsed[alpha] -= shift[alpha];
                flag = 0;
                next.c = solver.reconstruct(
                    tauNext, qiUsed, curr.c, (MD_THREAD EosState<Real> *)0,
                    pressure, flag);
                flagAll |= flag;
                flagUsed = flag;
                vs.firstRKStepW(curr, prev, G, th, next);
            }
            if (!converged && lastMiss > bestMiss) {
                for (int k = 0; k < 7; k++) qiUsed[k] = qi[k];
                for (int alpha = 0; alpha < 5; alpha++)
                    qiUsed[alpha] -= best[alpha];
                flag = 0;
                next.c = solver.reconstruct(
                    tauNext, qiUsed, curr.c, (MD_THREAD EosState<Real> *)0,
                    pressure, flag);
                flagAll |= flag;
                flagUsed = flag;
                vs.firstRKStepW(curr, prev, G, th, next);
            }
        }

    }

    storeCell(sNext, n, i, next.c);
    for (int k = 0; k < 14; k++) sNext[(kOffW + k) * n + i] = next.W[k];
    sNext[kOffPi * n + i] = next.pi;
    sNext[kOffPiChem * n + i] = next.piChem;
    if (MD_RUN(P, storeQ)) {
        const bool changed =
            (flagUsed & (kFlagMomentumCap | kFlagRegulate | kFlagRevert)) != 0;
        const Real pOut =
            (changed ? eosPressure<Real, EosType>(
                           d, tab, next.c.e, next.c.rhob, next.c.rhoq,
                           next.c.rhos)
                     : Real(0));
        for (int a = 0; a < 7; a++) {
            sNext[(kNPrim + a) * n + i] =
                (changed ? TJb(next.c, a, 0, pOut) * tauNext : qiUsed[a]);
        }
    }
    if (flagAll) flagOr(flags, flagAll);
}

}  // namespace music_device

#endif  // SRC_DEVICE_DEVICE_VISC_H_
