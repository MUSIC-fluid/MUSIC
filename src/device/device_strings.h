// Copyright 2026 Chun Shen
//
// The strings source (HydroSourceStrings, Initial_profile 13) on the
// device: kernelStringSource evaluates qi_source of Advance::FirstRKStepT
// for every cell of a stage from the u of the stage's input slice, so the
// host neither evaluates the source nor copies u and the source across.
//
// The host packs the strings of the current tau frame (the four lists of
// HydroSourceStrings::prepare_list_for_current_tau_frame) once per step
// into one buffer (see DeviceHydro::packStrings): a header (kStrHead
// values, enum StrHead), then kStrRec values per string (enum StrRec), in
// the order of HydroSourceStrings::QCD_strings_list. Each record has one
// flag per list, so every list is walked in the host's order. The
// quantities that depend on the string alone (the normalisations, the skip
// distances) are computed on the host with the host's expressions.
//
// Every function copies its host version statement by statement (in
// double, backend 1 reproduces the host source bit for bit).

#ifndef SRC_DEVICE_DEVICE_STRINGS_H_
#define SRC_DEVICE_DEVICE_STRINGS_H_

#if !defined(__METAL_VERSION__)
#include "device_ideal.h"
#endif

namespace music_device {

enum StrHead {
    kStrNStrings,     // strings in the buffer
    kStrNEnergy,      // QCD_strings_list_current_tau.size()
    kStrNRemnant,     // QCD_strings_remnant_list_current_tau.size()
    kStrNBaryon,      // QCD_strings_baryon_list_current_tau.size()
    kStrNElectric,    // QCD_strings_electric_list_current_tau.size()
    kStrPreEqFlow,    // preEqFlowFactor_
    kStrQuench,       // parton_quench_factor
    kStrUseBQ,        // DATA.use_BQ_ratios
    kStrRhob,         // DATA.turn_on_rhob
    kStrQS,           // DATA.turn_on_QS
    kStrDtau,         // DATA.delta_tau
    kStrXSize,        // DATA.x_size, y_size, eta_size
    kStrYSize,
    kStrEtaSize,
    kStrDx,           // DATA.delta_x, delta_y, delta_eta
    kStrDy,
    kStrDeta,
    kStrHbarc,        // Util::hbarc
    kStrHead
};

enum StrRec {
    // membership of the four lists
    kRInEnergy,
    kRInRemnant,
    kRInBaryon,
    kRInElectric,
    // QCD_string
    kRNorm,
    kRTilt,
    kRERemL,
    kRERemR,
    kRTauForm,
    kRSigmaX,
    kRSigmaEta,
    kRTau0,
    kREta0,
    kRXl,
    kRYl,
    kRXr,
    kRYr,
    kRTauEndL,
    kRTauEndR,
    kREtaL,
    kREtaR,
    kRYL,
    kRYR,
    kRRemL,
    kRRemR,
    kRTauBL,
    kRTauBR,
    kREtaBL,
    kREtaBR,
    kRYBL,
    kRYBR,
    kRFracL,
    kRFracR,
    kRTauQL,
    kRTauQR,
    kRQL,
    kRQR,
    kREtaQL,
    kREtaQR,
    // per string, as in the loops of HydroSourceStrings
    kRPrefPrep,    // energy: 1 / (2 pi preFlowNorm)
    kRPrefEtas,    // 1 / (sqrt(2 pi) sigma_eta)
    kRPrefactors,  // remnants: prefactor_prep * prefactor_etas
    kRPrefPrepB,   // n_B, n_Q: 1 / (2 pi sigma_x^2)
    kRSkipX,       // n_sigma_skip sigma_x
    kRSkipEta,     // n_sigma_skip sigma_eta
    kStrRec
};

//! HydroSourceStrings::getStringTransverseCoord
template <typename Real>
MD_FN Real stringTransverseCoord(Real xl, Real xr, Real etaFrac, Real tilt) {
    const Real x_c = Real(0.5) * (xl + xr);
    const Real S_transverse_position =
        x_c + tilt * (Real(0.5) - etaFrac) * (xl - xr);
    Real xT;
    Real sub_nucleon_distance = Real(2.0);
    xT = S_transverse_position;
    if (tilt > Real(0.0)) {
        if ((S_transverse_position - x_c) > sub_nucleon_distance) {
            xT = x_c + sub_nucleon_distance;
        } else if ((S_transverse_position - x_c) < -sub_nucleon_distance) {
            xT = x_c - sub_nucleon_distance;
        }
    }
    return xT;
}

//! HydroSourceStrings::get_hydro_energy_source
template <typename Real>
MD_FN void stringEnergySource(
    MD_GLOBAL const Real *S, Real tau, Real x, Real y, Real eta_s,
    MD_THREAD Real *j_mu) {
    for (int k = 0; k < 4; k++) j_mu[k] = Real(0);
    if (S[kStrNEnergy] == Real(0) && S[kStrNRemnant] == Real(0)) return;
    const int nStr = int(S[kStrNStrings]);
    MD_GLOBAL const Real *R0 = S + kStrHead;
    const Real small_eps = Real(1e-16);
    const Real preEqFlowFactor_ = S[kStrPreEqFlow];
    const Real dtau = S[kStrDtau];
    const Real exp_tau = Real(1.) / tau;
    for (int k = 0; k < nStr; k++) {
        MD_GLOBAL const Real *it = R0 + k * kStrRec;
        if (it[kRInEnergy] == Real(0)) continue;
        const Real sigma_x = it[kRSigmaX];
        const Real sigma_eta = it[kRSigmaEta];
        const Real prefactor_prep = it[kRPrefPrep];
        const Real prefactor_etas = it[kRPrefEtas];
        const Real skip_dis_x = it[kRSkipX];
        const Real skip_dis_eta = it[kRSkipEta];
        const Real tau_0 = it[kRTau0];
        const Real delta_tau = it[kRTauForm];
        if (eta_s < it[kREtaL] - skip_dis_eta
            || eta_s > it[kREtaR] + skip_dis_eta)
            continue;
        Real eta_frac =
            ((eta_s - it[kREtaL]) / dmax(small_eps, it[kREtaR] - it[kREtaL]));
        eta_frac = dmax(Real(0.), dmin(Real(1.), eta_frac));
        const Real x_perp =
            stringTransverseCoord(it[kRXl], it[kRXr], eta_frac, it[kRTilt]);
        Real x_dis = x - x_perp;
        if (fabs(x_dis) > skip_dis_x) continue;
        const Real y_perp =
            stringTransverseCoord(it[kRYl], it[kRYr], eta_frac, it[kRTilt]);
        Real y_dis = y - y_perp;
        if (fabs(y_dis) > skip_dis_x) continue;
        Real eta_s_shift = Real(0.0);
        Real tau_L = tau - dtau / Real(2.);
        if (tau_L > tau_0 + delta_tau) {
            eta_s_shift = acosh(
                (tau_L * tau_L + tau_0 * tau_0 - delta_tau * delta_tau)
                / dmax(small_eps, Real(2.) * tau_L * tau_0));
        }
        Real eta_s_L = dmin(it[kREtaR], it[kREta0] - eta_s_shift);
        Real eta_s_R = dmax(it[kREtaL], it[kREta0] + eta_s_shift);
        Real eta_s_next_shift = Real(0.0);
        Real tau_next = tau + dtau / Real(2.);
        if (tau_next > tau_0 + delta_tau) {
            eta_s_next_shift = acosh(
                (tau_next * tau_next + tau_0 * tau_0 - delta_tau * delta_tau)
                / dmax(small_eps, Real(2.) * tau_next * tau_0));
        }
        Real eta_s_L_next = dmax(it[kREtaL], it[kREta0] - eta_s_next_shift);
        Real eta_s_R_next = dmin(it[kREtaR], it[kREta0] + eta_s_next_shift);
        bool flag_left = true;
        if (eta_s_L_next > eta_s_L) flag_left = false;
        bool flag_right = true;
        if (eta_s_R_next < eta_s_R) flag_right = false;
        Real exp_eta_s = Real(0.);
        if (flag_left) {
            if (eta_s > eta_s_L_next - skip_dis_eta
                && eta_s < eta_s_L + skip_dis_eta) {
                exp_eta_s +=
                    Real(0.5)
                    * (-derf((eta_s_L_next - eta_s) / (sqrt(Real(2.)) * sigma_eta))
                       + derf((eta_s_L - eta_s) / (sqrt(Real(2.)) * sigma_eta)));
            }
        }
        if (flag_right) {
            if (eta_s > eta_s_R - skip_dis_eta
                && eta_s < eta_s_R_next + skip_dis_eta) {
                exp_eta_s +=
                    Real(0.5)
                    * (-derf((eta_s_R - eta_s) / (sqrt(Real(2.)) * sigma_eta))
                       + derf((eta_s_R_next - eta_s)
                              / (sqrt(Real(2.)) * sigma_eta)));
            }
        }
        Real exp_xperp =
            exp(-(x_dis * x_dis + y_dis * y_dis) / (Real(2.) * sigma_x * sigma_x));
        Real cosh_perp =
            (cosh(preEqFlowFactor_ * sqrt(x_dis * x_dis + y_dis * y_dis)
                  / sigma_x));
        Real sinh_perp =
            (sinh(preEqFlowFactor_ * sqrt(x_dis * x_dis + y_dis * y_dis)
                  / sigma_x));
        Real phi_perp = atan2(y_dis, x_dis);
        Real e_local = exp_tau * exp_xperp * exp_eta_s * it[kRNorm];
        Real Delta_eta = it[kREtaR] - it[kREtaL];
        Real denorm_safe =
            copysign(dmax(small_eps, fabs(Delta_eta)), Delta_eta);
        Real y_string =
            (it[kRYL]
             + (it[kRYR] - it[kRYL]) / denorm_safe * (eta_s - it[kREtaL]));
        Real cosh_long = cosh(y_string - eta_s);
        Real sinh_long = sinh(y_string - eta_s);
        Real local_eperp = prefactor_etas * prefactor_prep * e_local;
        j_mu[0] += local_eperp * cosh_long * cosh_perp;
        j_mu[1] += local_eperp * sinh_perp * cos(phi_perp);
        j_mu[2] += local_eperp * sinh_perp * sin(phi_perp);
        j_mu[3] += local_eperp * sinh_long * cosh_perp;
    }
    for (int k = 0; k < nStr; k++) {
        MD_GLOBAL const Real *it = R0 + k * kStrRec;
        if (it[kRInRemnant] == Real(0)) continue;
        const Real sigma_x = it[kRSigmaX];
        const Real sigma_eta = it[kRSigmaEta];
        const Real prefactors = it[kRPrefactors];
        const Real skip_dis_x = it[kRSkipX];
        const Real skip_dis_eta = it[kRSkipEta];
        bool flag_left = false;
        if (it[kRTauEndL] >= tau - dtau / Real(2.)
            && it[kRTauEndL] < tau + dtau / Real(2.) && it[kRRemL] > Real(0.)) {
            flag_left = true;
        }
        bool flag_right = false;
        if (it[kRTauEndR] >= tau - dtau / Real(2.)
            && it[kRTauEndR] < tau + dtau / Real(2.) && it[kRRemR] > Real(0.)) {
            flag_right = true;
        }
        Real exp_eta_s_left = Real(0.0);
        if (flag_left) {
            Real eta_dis_left = fabs(eta_s - it[kREtaL]);
            if (eta_dis_left < skip_dis_eta) {
                exp_eta_s_left =
                    (exp(-eta_dis_left * eta_dis_left
                         / (Real(2.) * sigma_eta * sigma_eta)));
            }
        }
        Real exp_eta_s_right = Real(0.0);
        if (flag_right) {
            Real eta_dis_right = fabs(eta_s - it[kREtaR]);
            if (eta_dis_right < skip_dis_eta) {
                exp_eta_s_right =
                    (exp(-eta_dis_right * eta_dis_right
                         / (Real(2.) * sigma_eta * sigma_eta)));
            }
        }
        for (int side = 0; side < 2; side++) {
            if (side == 0 ? !flag_left : !flag_right) continue;
            const Real xEnd = it[side == 0 ? kRXl : kRXr];
            const Real yEnd = it[side == 0 ? kRYl : kRYr];
            const Real expEta = (side == 0 ? exp_eta_s_left : exp_eta_s_right);
            const Real rem = it[side == 0 ? kRRemL : kRRemR];
            const Real remNorm = it[side == 0 ? kRERemL : kRERemR];
            const Real yEndRap = it[side == 0 ? kRYL : kRYR];
            Real x_dis = x - xEnd;
            if (fabs(x_dis) <= skip_dis_x) {
                Real y_dis = y - yEnd;
                if (fabs(y_dis) <= skip_dis_x) {
                    Real cosh_long = exp_tau * expEta * (rem) * (remNorm)
                                     * cosh(yEndRap - eta_s);
                    Real sinh_long = exp_tau * expEta * (rem) * (remNorm)
                                     * sinh(yEndRap - eta_s);
                    Real cosh_perp = Real(1.);
                    Real sinh_perp = Real(0.);
                    Real phi_perp = Real(0.);
                    Real exp_xperp = Real(0.);
                    if (cosh_long > Real(0)) {
                        exp_xperp =
                            (prefactors
                             * exp(-(x_dis * x_dis + y_dis * y_dis)
                                   / (Real(2.) * sigma_x * sigma_x)));
                        cosh_perp =
                            (cosh(preEqFlowFactor_
                                  * sqrt(x_dis * x_dis + y_dis * y_dis)
                                  / sigma_x));
                        sinh_perp =
                            (sinh(preEqFlowFactor_
                                  * sqrt(x_dis * x_dis + y_dis * y_dis)
                                  / sigma_x));
                        phi_perp = atan2(y_dis, x_dis);
                    }
                    j_mu[0] += exp_xperp * cosh_long * cosh_perp;
                    j_mu[1] += exp_xperp * sinh_perp * cos(phi_perp);
                    j_mu[2] += exp_xperp * sinh_perp * sin(phi_perp);
                    j_mu[3] += exp_xperp * sinh_long * cosh_perp;
                }
            }
        }
    }
    const Real prefactor_tau = Real(1.) / dtau;
    const Real unit_convert = Real(1.0) / S[kStrHbarc];
    const Real prefactors = prefactor_tau * unit_convert;
    j_mu[0] *= prefactors;
    j_mu[1] *= prefactors;
    j_mu[2] *= prefactors;
    j_mu[3] *= prefactors;
}

//! the flow quantities of get_hydro_rhob_source and get_hydro_rhoQ_source
template <typename Real>
struct StrFlow {
    Real y_perp_flow, y_long_flow, cos_phi_flow, sin_phi_flow;
};

template <typename Real>
MD_FN StrFlow<Real> stringFlow(MD_THREAD const Real *u_mu, Real eta_s) {
    const Real small_eps = Real(1e-16);
    StrFlow<Real> f;
    const Real gamma_perp_flow = sqrt(Real(1.) + u_mu[1] * u_mu[1] + u_mu[2] * u_mu[2]);
    f.y_perp_flow = acosh(gamma_perp_flow);
    f.y_long_flow = asinh(u_mu[3] / gamma_perp_flow) + eta_s;
    const Real u_perp_flow = sqrt(u_mu[1] * u_mu[1] + u_mu[2] * u_mu[2]);
    f.cos_phi_flow = (u_perp_flow > small_eps ? u_mu[1] / u_perp_flow : Real(1.));
    f.sin_phi_flow = (u_perp_flow > small_eps ? u_mu[2] / u_perp_flow : Real(0.));
    return f;
}

//! the p.u of a deposited parton
template <typename Real>
MD_FN Real stringPDotU(
    MD_THREAD const Real *u_mu, Real eta_s, Real y_dump, Real y_dump_perp,
    MD_THREAD const StrFlow<Real> &f) {
    return (u_mu[0] * cosh(y_dump - eta_s) * cosh(y_dump_perp)
            - u_mu[1] * sinh(y_dump_perp) * f.cos_phi_flow
            - u_mu[2] * sinh(y_dump_perp) * f.sin_phi_flow
            - u_mu[3] * sinh(y_dump - eta_s) * cosh(y_dump_perp));
}

//! HydroSourceStrings::get_hydro_rhob_source
template <typename Real>
MD_FN Real stringRhobSource(
    MD_GLOBAL const Real *S, Real tau, Real x, Real y, Real eta_s,
    MD_THREAD const Real *u_mu) {
    Real res = Real(0.);
    if (S[kStrNBaryon] == Real(0)) return (res);
    const StrFlow<Real> fl = stringFlow(u_mu, eta_s);
    const Real dtau = S[kStrDtau];
    const Real parton_quench_factor = S[kStrQuench];
    const Real exp_tau = Real(1.0) / tau;
    const int nStr = int(S[kStrNStrings]);
    MD_GLOBAL const Real *R0 = S + kStrHead;
    for (int k = 0; k < nStr; k++) {
        MD_GLOBAL const Real *it = R0 + k * kStrRec;
        if (it[kRInBaryon] == Real(0)) continue;
        const Real sigma_x = it[kRSigmaX];
        const Real sigma_eta = it[kRSigmaEta];
        const Real prefactor_prep = it[kRPrefPrepB];
        const Real prefactor_etas = it[kRPrefEtas];
        const Real skip_dis_x = it[kRSkipX];
        const Real skip_dis_eta = it[kRSkipEta];
        int flag_left = 0;
        if (it[kRTauBL] >= tau - dtau / Real(2.)
            && it[kRTauBL] < tau + dtau / Real(2.) && it[kRFracL] > Real(0.)) {
            flag_left = 1;
        }
        int flag_right = 0;
        if (it[kRTauBR] >= tau - dtau / Real(2.)
            && it[kRTauBR] < tau + dtau / Real(2.) && it[kRFracR] > Real(0.)) {
            flag_right = 1;
        }
        if (flag_left == 0 && flag_right == 0) continue;
        if (eta_s < it[kREtaL] - skip_dis_eta
            || eta_s > it[kREtaR] + skip_dis_eta)
            continue;
        for (int side = 0; side < 2; side++) {
            if ((side == 0 ? flag_left : flag_right) != 1) continue;
            Real eta_dis = eta_s - it[side == 0 ? kREtaBL : kREtaBR];
            if (fabs(eta_dis) <= skip_dis_eta) {
                Real exp_eta =
                    exp(-(eta_dis * eta_dis) / (Real(2.) * sigma_eta * sigma_eta));
                Real x_dis = x - it[side == 0 ? kRXl : kRXr];
                if (fabs(x_dis) <= skip_dis_x) {
                    Real y_dis = y - it[side == 0 ? kRYl : kRYr];
                    if (fabs(y_dis) <= skip_dis_x) {
                        Real exp_xperp =
                            exp(-(x_dis * x_dis + y_dis * y_dis)
                                / (Real(2.) * sigma_x * sigma_x));
                        Real fsmear =
                            exp_tau
                            * (exp_xperp * exp_eta
                               * (it[side == 0 ? kRFracL : kRFracR]));
                        if (fsmear > Real(0.)) {
                            Real rapidity_local = it[side == 0 ? kRYBL : kRYBR];
                            Real y_dump =
                                ((Real(1.) - parton_quench_factor) * rapidity_local
                                 + parton_quench_factor * fl.y_long_flow);
                            Real y_dump_perp =
                                parton_quench_factor * fl.y_perp_flow;
                            Real p_dot_u = Real(1.);
                            if (parton_quench_factor < Real(1.)) {
                                p_dot_u = stringPDotU(
                                    u_mu, eta_s, y_dump, y_dump_perp, fl);
                            }
                            res += prefactor_etas * prefactor_prep * p_dot_u
                                   * fsmear;
                        }
                    }
                }
            }
        }
    }
    const Real prefactor_tau = Real(1.) / dtau;
    res *= prefactor_tau;
    return (res);
}

//! HydroSourceStrings::get_hydro_rhoQ_source
template <typename Real>
MD_FN Real stringRhoQSource(
    MD_GLOBAL const Real *S, Real tau, Real x, Real y, Real eta_s,
    MD_THREAD const Real *u_mu) {
    Real res = Real(0.);
    if (S[kStrUseBQ] != Real(0)) {
        Real nb = stringRhobSource(S, tau, x, y, eta_s, u_mu);
        res = Real(0.4) * nb;
        return (res);
    }
    if (S[kStrNElectric] == Real(0)) return (res);
    const StrFlow<Real> fl = stringFlow(u_mu, eta_s);
    const Real dtau = S[kStrDtau];
    const Real parton_quench_factor = S[kStrQuench];
    const Real y_dump = fl.y_long_flow;
    const Real y_dump_perp = parton_quench_factor * fl.y_perp_flow;
    Real p_dot_u = Real(1.);
    if (parton_quench_factor < Real(1.)) {
        p_dot_u = stringPDotU(u_mu, eta_s, y_dump, y_dump_perp, fl);
    }
    const Real exp_tau = Real(1.0) / tau;
    const int nStr = int(S[kStrNStrings]);
    MD_GLOBAL const Real *R0 = S + kStrHead;
    for (int k = 0; k < nStr; k++) {
        MD_GLOBAL const Real *it = R0 + k * kStrRec;
        if (it[kRInElectric] == Real(0)) continue;
        const Real sigma_x = it[kRSigmaX];
        const Real sigma_eta = it[kRSigmaEta];
        const Real prefactor_prep = it[kRPrefPrepB];
        const Real prefactor_etas = it[kRPrefEtas];
        const Real skip_dis_x = it[kRSkipX];
        const Real skip_dis_eta = it[kRSkipEta];
        const bool flag_left =
            (it[kRTauQL] >= tau - dtau / Real(2.)
             && it[kRTauQL] < tau + dtau / Real(2.)
             && fabs(it[kRQL]) > Real(1e-8));
        const bool flag_right =
            (it[kRTauQR] >= tau - dtau / Real(2.)
             && it[kRTauQR] < tau + dtau / Real(2.)
             && fabs(it[kRQR]) > Real(1e-8));
        if (!flag_left && !flag_right) continue;
        if (eta_s < it[kREtaL] - skip_dis_eta
            || eta_s > it[kREtaR] + skip_dis_eta)
            continue;
        for (int side = 0; side < 2; side++) {
            if (side == 0 ? !flag_left : !flag_right) continue;
            const Real eta_end = it[side == 0 ? kREtaQL : kREtaQR];
            const Real x_end = it[side == 0 ? kRXl : kRXr];
            const Real y_end = it[side == 0 ? kRYl : kRYr];
            const Real Qe = it[side == 0 ? kRQL : kRQR];
            const Real eta_dis = eta_s - eta_end;
            if (fabs(eta_dis) > skip_dis_eta) continue;
            const Real x_dis = x - x_end;
            if (fabs(x_dis) > skip_dis_x) continue;
            const Real y_dis = y - y_end;
            if (fabs(y_dis) > skip_dis_x) continue;
            const Real exp_eta =
                exp(-(eta_dis * eta_dis) / (Real(2.) * sigma_eta * sigma_eta));
            const Real exp_xperp =
                exp(-(x_dis * x_dis + y_dis * y_dis) / (Real(2.) * sigma_x * sigma_x));
            res += (prefactor_etas * prefactor_prep * p_dot_u * exp_tau
                    * exp_xperp * exp_eta * Qe);
        }
    }
    const Real prefactor_tau = Real(1.) / dtau;
    res *= prefactor_tau;
    return (res);
}

//! DeviceHydro::computeSource for cell i with the strings source: src from
//! the u of the stage's input slice sIn (tau_n, then the predictor)
template <typename Real>
MD_FN void kernelStringSource(
    MD_CONST const StepParams<Real> &P, MD_GLOBAL const Real *S,
    MD_GLOBAL const Real *sIn, MD_GLOBAL Real *src, int i, MD_FLAGS flags) {
    const int nx = P.nx, ny = P.ny;
    const int n = nx * ny * P.neta;
    const int ix = i % nx;
    const int iy = (i / nx) % ny;
    const int ieta = i / (nx * ny);
    const Real tauRk = P.tau + Real(P.rkFlag) * P.dtau;
    // as in Advance::AdvanceIt and FirstRKStepT
    Real eta_s_local = -S[kStrEtaSize] / Real(2.) + Real(ieta) * S[kStrDeta];
    Real x_local = -S[kStrXSize] / Real(2.) + Real(ix) * S[kStrDx];
    Real y_local = -S[kStrYSize] / Real(2.) + Real(iy) * S[kStrDy];
    Real qi_source[7] = {Real(0.), Real(0.), Real(0.), Real(0.),
                         Real(0.), Real(0.), Real(0.)};
    Real j_mu[4];
    Real u_local[4];
    for (int ii = 0; ii < 4; ii++) u_local[ii] = sIn[(4 + ii) * n + i];
    stringEnergySource(S, tauRk, x_local, y_local, eta_s_local, j_mu);
    bool isNan = false;
    for (int ii = 0; ii < 4; ii++) {
        qi_source[ii] = tauRk * j_mu[ii];
        if (isnan(qi_source[ii])) isNan = true;
    }
    if (S[kStrRhob] != Real(0)) {
        qi_source[4] =
            (tauRk
             * stringRhobSource(S, tauRk, x_local, y_local, eta_s_local, u_local));
    }
    if (S[kStrQS] != Real(0)) {
        qi_source[5] =
            (tauRk
             * stringRhoQSource(S, tauRk, x_local, y_local, eta_s_local, u_local));
    }
    for (int a = 0; a < 7; a++) src[a * n + i] = qi_source[a];
    if (isNan) flagOr(flags, kFlagSourceNan);
}

}  // namespace music_device

#endif  // SRC_DEVICE_DEVICE_STRINGS_H_
