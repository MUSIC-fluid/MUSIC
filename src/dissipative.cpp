// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include "dissipative.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "cell.h"
#include "data.h"
#include "util.h"

using Util::hbarc;
using Util::map_1d_idx_to_2d;
using Util::small_eps;

Diss::Diss(const EOS &eosIn, const InitData &Data_in)
    : DATA(Data_in), eos(eosIn), minmod(Data_in), transport_coeffs_(Data_in) {}

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Dissipative parts */
/* Sangyong Nov 18 2014 */
/* change: alpha first which is the case
for everywhere else. also, this change is necessary
to use Wmunu[rk_flag][4][mu] as the dissipative baryon current*/
/* this is the only one that is being subtracted in the rhs */
void Diss::MakeWSource(
    const double tau, const int ix, const int iy, const int ieta, TJbVec &dwmn,
    Fields &arenaCurr, Fields &arenaPrev, const int fieldIdx) {
    /* calculate d_m (tau W^{m,alpha}) + (geom source terms) */
    const double delta[4] = {0.0, DATA.delta_x, DATA.delta_y, DATA.delta_eta};
    const double tau_fac[4] = {0.0, tau, tau, 1.0};

    dwmn = {0.};
    for (int alpha = 0; alpha < 5; alpha++) {
        /* partial_tau W^tau alpha */
        /* this is partial_tau evaluated at tau */
        /* this is the first step. so rk_flag = 0 */
        /* Sangyong Nov 18 2014 */
        /* change: alpha first which is the case
                   for everywhere else. also, this change is necessary to use
                   Wmunu[rk_flag][4][mu] as the dissipative baryon current */
        // dW/dtau
        // backward time derivative (first order is more stable)
        int idx_1d_alpha0 = map_2d_idx_to_1d(alpha, 0);
        double dWdtau =
            ((arenaCurr.Wmunu_[idx_1d_alpha0][fieldIdx]
              - arenaPrev.Wmunu_[idx_1d_alpha0][fieldIdx])
             / DATA.delta_tau);

        /* bulk pressure term */
        double dPidtau = 0.0;
        double Pi_alpha0 = 0.0;
        // the total bulk pressure Pi = Pi_bulk + Pi_chem enters T^{mu nu}
        // whenever either part is evolved
        if (alpha < 4 && (DATA.turn_on_bulk || DATA.turn_on_bulk_chem)) {
            double gfac = (alpha == 0 ? -1.0 : 0.0);
            double PiCurr =
                arenaCurr.piBulk_[fieldIdx] + arenaCurr.piBulkChem_[fieldIdx];
            double PiPrev =
                arenaPrev.piBulk_[fieldIdx] + arenaPrev.piBulkChem_[fieldIdx];

            Pi_alpha0 =
                PiCurr
                * (gfac
                   + arenaCurr.u_[alpha][fieldIdx] * arenaCurr.u_[0][fieldIdx]);
            dPidtau = (Pi_alpha0
                       - PiPrev
                             * (gfac
                                + arenaPrev.u_[alpha][fieldIdx]
                                      * arenaPrev.u_[0][fieldIdx]))
                      / DATA.delta_tau;
        }

        double sf = tau * dWdtau + arenaCurr.Wmunu_[idx_1d_alpha0][fieldIdx];
        double bf = tau * dPidtau + Pi_alpha0;
        dwmn[alpha] += sf + bf;
        if (std::isnan(sf + bf)) {
            music_message << "[Error]Diss::MakeWSource: ";
            music_message << "sf=" << sf << " bf=" << bf << " Wmunu ="
                          << arenaCurr.Wmunu_[idx_1d_alpha0][fieldIdx]
                          << " pi_b =" << arenaCurr.piBulk_[fieldIdx]
                          << " prev_pi_b=" << arenaPrev.piBulk_[fieldIdx]
                          << " pi_b_chem =" << arenaCurr.piBulkChem_[fieldIdx]
                          << " prev_pi_b_chem="
                          << arenaPrev.piBulkChem_[fieldIdx];
            music_message.flush("error");
            music_message << "dWdtau = " << dWdtau;
            music_message.flush("error");
            music_message << "dPidtau = " << dPidtau
                          << ", Pi_alpha0 = " << Pi_alpha0;
            music_message.flush("error");
        }
    }

    EnergyFlowVec W_eta_p = {0.};  // save tau*W^{\eta \nu} at eta + deta/2
    EnergyFlowVec W_eta_m = {0.};  // save tau*W^{\eta \nu} at eta - deta/2
    FieldNeighbourLoop1(
        arenaCurr, ix, iy, ieta, FNLLAMBDAS1 {
            for (int alpha = 0; alpha < 5; alpha++) {
                int idx_1d = map_2d_idx_to_1d(alpha, direction);
                double sg = arenaCurr.Wmunu_[idx_1d][Ic] * tau_fac[direction];
                double sgp1 =
                    arenaCurr.Wmunu_[idx_1d][Ip1] * tau_fac[direction];
                double sgm1 =
                    arenaCurr.Wmunu_[idx_1d][Im1] * tau_fac[direction];
                // dWdx += minmod.minmod_dx(sgp1, sg, sgm1)/delta[direction];
                //  use central difference to preserve conservation law exactly
                double W_m = (sg + sgm1) * 0.5;
                double W_p = (sg + sgp1) * 0.5;

                double Pi_m = 0;
                double Pi_p = 0;
                if ((DATA.turn_on_bulk || DATA.turn_on_bulk_chem)
                    && alpha < 4) {
                    double gfac1 = (alpha == (direction) ? 1.0 : 0.0);
                    double PiIp1 =
                        arenaCurr.piBulk_[Ip1] + arenaCurr.piBulkChem_[Ip1];
                    double PiIc =
                        arenaCurr.piBulk_[Ic] + arenaCurr.piBulkChem_[Ic];
                    double PiIm1 =
                        arenaCurr.piBulk_[Im1] + arenaCurr.piBulkChem_[Im1];

                    double bgp1 = PiIp1
                                  * (gfac1
                                     + arenaCurr.u_[alpha][Ip1]
                                           * arenaCurr.u_[direction][Ip1])
                                  * tau_fac[direction];
                    double bg = PiIc
                                * (gfac1
                                   + arenaCurr.u_[alpha][Ic]
                                         * arenaCurr.u_[direction][Ic])
                                * tau_fac[direction];
                    double bgm1 = PiIm1
                                  * (gfac1
                                     + arenaCurr.u_[alpha][Im1]
                                           * arenaCurr.u_[direction][Im1])
                                  * tau_fac[direction];
                    // dPidx += minmod.minmod_dx(bgp1, bg,
                    // bgm1)/delta[direction];
                    //  use central difference to preserve conservation law
                    //  exactly
                    Pi_m = (bg + bgm1) * 0.5;
                    Pi_p = (bg + bgp1) * 0.5;
                }
                if (direction == 3 && (alpha == 0 || alpha == 3)) {
                    W_eta_p[alpha] = W_p + Pi_p;
                    W_eta_m[alpha] = W_m + Pi_m;
                } else {
                    dwmn[alpha] += (W_p - W_m + Pi_p - Pi_m) / delta[direction];
                }
            }
        });

    // add longitudinal flux with the discretized geometric terms
    // careful about the boost-invariant case when deta could be arbitary
    double cosh_deta =
        cosh(delta[3] / 2.) / std::max(delta[3], Util::small_eps);
    double sinh_deta =
        sinh(delta[3] / 2.) / std::max(delta[3], Util::small_eps);
    sinh_deta = std::max(0.5, sinh_deta);
    if (DATA.boost_invariant) {
        // if the simulation is boost-invariant,
        // we directly use the limiting value at \Delta eta = 0
        // Longitudinal derivatives should be 0, we set cosh_deta = 0 here
        cosh_deta = 0.0;
        sinh_deta = 0.5;
    }
    dwmn[0] +=
        ((W_eta_p[0] - W_eta_m[0]) * cosh_deta
         + (W_eta_p[3] + W_eta_m[3]) * sinh_deta);
    dwmn[3] +=
        ((W_eta_p[3] - W_eta_m[3]) * cosh_deta
         + (W_eta_p[0] + W_eta_m[0]) * sinh_deta);

    // sources due to coordinate transform this is added to partial_m W^mn
    // dwmn[0] += grid_pt.Wmunu[9];
    // dwmn[0] += grid_pt.pi_b*(1.0 + grid_pt.u[3]*grid_pt.u[3]);
    // dwmn[3] += grid_pt.Wmunu[3];
    // dwmn[3] += grid_pt.pi_b*(grid_pt.u[0]*grid_pt.u[3]);
}

void Diss::computeWSourceDivergence(
    const double tau, const Fields &arenaCurr, const Fields &arenaPrev,
    std::vector<double> &wSource) {
    // the same operations as MakeWSource, in the same order
    const double delta[4] = {0.0, DATA.delta_x, DATA.delta_y, DATA.delta_eta};
    const double tau_fac[4] = {0.0, tau, tau, 1.0};
    const bool withBulk = (DATA.turn_on_bulk || DATA.turn_on_bulk_chem);
    const int nx = arenaCurr.nX();
    const int ny = arenaCurr.nY();
    const int neta = arenaCurr.nEta();
    const size_t nStorage = arenaCurr.nStorage();
    wSource.resize(kWSourceComponents * nStorage);
    // the geometric factors of MakeWSource
    wSourceCoshDeta_ =
        cosh(delta[3] / 2.) / std::max(delta[3], Util::small_eps);
    wSourceSinhDeta_ =
        sinh(delta[3] / 2.) / std::max(delta[3], Util::small_eps);
    wSourceSinhDeta_ = std::max(0.5, wSourceSinhDeta_);
    if (DATA.boost_invariant) {
        wSourceCoshDeta_ = 0.0;
        wSourceSinhDeta_ = 0.5;
    }
    // face components: tau W^{i alpha} for alpha < 5, then the bulk part
    // Pi (g^{i alpha} + u^alpha u^i) tau for alpha < 4
    constexpr int nFaceComp = 9;

    for (int direction = 1; direction <= 3; direction++) {
        const int dx = (direction == 1);
        const int dy = (direction == 2);
        const int deta = (direction == 3);
        const int nfx = nx + dx;
        const int nfy = ny + dy;
        const int nfeta = neta + deta;
        const size_t nFaces = static_cast<size_t>(nfx) * nfy * nfeta;
        wSourceFaces_.resize(nFaceComp * nFaces);

        // the average over the cells L = f - 1 and R = f of face f
#pragma omp parallel for collapse(3) schedule(dynamic, 64)
        for (int feta = 0; feta < nfeta; feta++)
            for (int fy = 0; fy < nfy; fy++)
                for (int fx = 0; fx < nfx; fx++) {
                    const int IL = arenaCurr.getFieldIdxHalo(
                        fx - dx, fy - dy, feta - deta);
                    const int IR = arenaCurr.getFieldIdxHalo(fx, fy, feta);
                    const size_t f =
                        fx + nfx * (fy + static_cast<size_t>(nfy) * feta);
                    for (int alpha = 0; alpha < 5; alpha++) {
                        const int idx_1d = map_2d_idx_to_1d(alpha, direction);
                        const double sL =
                            arenaCurr.Wmunu_[idx_1d][IL] * tau_fac[direction];
                        const double sR =
                            arenaCurr.Wmunu_[idx_1d][IR] * tau_fac[direction];
                        wSourceFaces_[alpha * nFaces + f] = (sR + sL) * 0.5;
                    }
                    if (!withBulk) continue;
                    const double PiL =
                        arenaCurr.piBulk_[IL] + arenaCurr.piBulkChem_[IL];
                    const double PiR =
                        arenaCurr.piBulk_[IR] + arenaCurr.piBulkChem_[IR];
                    for (int alpha = 0; alpha < 4; alpha++) {
                        const double gfac1 = (alpha == direction ? 1.0 : 0.0);
                        const double bL =
                            PiL
                            * (gfac1
                               + arenaCurr.u_[alpha][IL]
                                     * arenaCurr.u_[direction][IL])
                            * tau_fac[direction];
                        const double bR =
                            PiR
                            * (gfac1
                               + arenaCurr.u_[alpha][IR]
                                     * arenaCurr.u_[direction][IR])
                            * tau_fac[direction];
                        wSourceFaces_[(5 + alpha) * nFaces + f] =
                            (bR + bL) * 0.5;
                    }
                }

#pragma omp parallel for collapse(3) schedule(dynamic, 64)
        for (int ieta = 0; ieta < neta; ieta++)
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++) {
                    const int fieldIdx = arenaCurr.getFieldIdx(ix, iy, ieta);
                    const size_t Imh =
                        ix + nfx * (iy + static_cast<size_t>(nfy) * ieta);
                    const size_t Iph =
                        (ix + dx)
                        + nfx
                              * ((iy + dy)
                                 + static_cast<size_t>(nfy) * (ieta + deta));
                    double *out = wSource.data() + fieldIdx;
                    if (direction == 1) {
                        // partial_tau (tau W^{tau alpha}) and its bulk part
                        for (int alpha = 0; alpha < 5; alpha++) {
                            const int idx_1d_alpha0 =
                                map_2d_idx_to_1d(alpha, 0);
                            double dWdtau =
                                ((arenaCurr.Wmunu_[idx_1d_alpha0][fieldIdx]
                                  - arenaPrev.Wmunu_[idx_1d_alpha0][fieldIdx])
                                 / DATA.delta_tau);
                            double dPidtau = 0.0;
                            double Pi_alpha0 = 0.0;
                            if (alpha < 4 && withBulk) {
                                double gfac = (alpha == 0 ? -1.0 : 0.0);
                                double PiCurr =
                                    arenaCurr.piBulk_[fieldIdx]
                                    + arenaCurr.piBulkChem_[fieldIdx];
                                double PiPrev =
                                    arenaPrev.piBulk_[fieldIdx]
                                    + arenaPrev.piBulkChem_[fieldIdx];
                                Pi_alpha0 =
                                    PiCurr
                                    * (gfac
                                       + arenaCurr.u_[alpha][fieldIdx]
                                             * arenaCurr.u_[0][fieldIdx]);
                                dPidtau =
                                    (Pi_alpha0
                                     - PiPrev
                                           * (gfac
                                              + arenaPrev.u_[alpha][fieldIdx]
                                                    * arenaPrev
                                                          .u_[0][fieldIdx]))
                                    / DATA.delta_tau;
                            }
                            double sf =
                                tau * dWdtau
                                + arenaCurr.Wmunu_[idx_1d_alpha0][fieldIdx];
                            double bf = tau * dPidtau + Pi_alpha0;
                            double d = 0.;
                            d += sf + bf;
                            out[alpha * nStorage] = d;
                            if (std::isnan(sf + bf)) {
                                music_message
                                    << "[Error]Diss::computeWSourceDivergence: "
                                    << "sf=" << sf << " bf=" << bf;
                                music_message.flush("error");
                            }
                        }
                    }
                    for (int alpha = 0; alpha < 5; alpha++) {
                        const double W_p = wSourceFaces_[alpha * nFaces + Iph];
                        const double W_m = wSourceFaces_[alpha * nFaces + Imh];
                        double Pi_m = 0;
                        double Pi_p = 0;
                        if (withBulk && alpha < 4) {
                            Pi_p = wSourceFaces_[(5 + alpha) * nFaces + Iph];
                            Pi_m = wSourceFaces_[(5 + alpha) * nFaces + Imh];
                        }
                        if (direction == 3 && (alpha == 0 || alpha == 3)) {
                            const int k = (alpha == 0 ? 0 : 1);
                            out[(5 + k) * nStorage] = W_p + Pi_p;
                            out[(7 + k) * nStorage] = W_m + Pi_m;
                        } else {
                            out[alpha * nStorage] +=
                                (W_p - W_m + Pi_p - Pi_m) / delta[direction];
                        }
                    }
                }
    }
}

void Diss::computeWDivergences(
    const double tau, const Fields &arenaCurr, const Fields &arenaPrev,
    const bool withFlux, std::vector<double> &wFluxDiv,
    std::vector<double> &wSource) {
    // the operations of computeWFluxDivergence and computeWSourceDivergence,
    // in the same order for every cell
    static const int piIdxArr[kWComponents] = {4, 5,  6,  7,  8,
                                               9, 14, 11, 12, 13};
    const Fields &arena = arenaCurr;
    const int nComp = nWComponents();
    const int nx = arena.nX();
    const int ny = arena.nY();
    const int neta = arena.nEta();
    const size_t nStorage = arena.nStorage();
    // the flux and the source use different eta spacings
    const double deltaFlux[4] = {
        0.0, DATA.delta_x, DATA.delta_y, DATA.delta_eta * tau};
    const double deltaSrc[4] = {
        0.0, DATA.delta_x, DATA.delta_y, DATA.delta_eta};
    const double tau_fac[4] = {0.0, tau, tau, 1.0};
    const double delta_tau = DATA.delta_tau;
    const bool withBulk = (DATA.turn_on_bulk || DATA.turn_on_bulk_chem);
    if (withFlux) wFluxDiv.resize(static_cast<size_t>(kWComponents) * nStorage);
    wSource.resize(kWSourceComponents * nStorage);
    wSourceCoshDeta_ =
        cosh(deltaSrc[3] / 2.) / std::max(deltaSrc[3], Util::small_eps);
    wSourceSinhDeta_ =
        sinh(deltaSrc[3] / 2.) / std::max(deltaSrc[3], Util::small_eps);
    wSourceSinhDeta_ = std::max(0.5, wSourceSinhDeta_);
    if (DATA.boost_invariant) {
        wSourceCoshDeta_ = 0.0;
        wSourceSinhDeta_ = 0.5;
    }

    const FieldComponent *fluxComp[kWComponents];
    for (int idx = 0; idx < kWComponents; idx++) {
        const int i = piIdxArr[idx];
        fluxComp[idx] =
            (i == 9 ? &arena.piBulk_
                    : (i == 14 ? &arena.piBulkChem_ : &arena.Wmunu_[i]));
    }

    // the values of one face: the KT fluxes of the W components, then the
    // averages of tau W^{d alpha} (alpha < 5) and of the bulk part (alpha < 4)
    constexpr int nFaceVal = kWComponents + 9;
    auto faceValues = [&](const int direction, const int fx, const int fy,
                          const int feta, double *v) {
        const int dx = (direction == 1);
        const int dy = (direction == 2);
        const int deta = (direction == 3);
        const int IL = arena.getFieldIdxHalo(fx - dx, fy - dy, feta - deta);
        const int IR = arena.getFieldIdxHalo(fx, fy, feta);
        if (withFlux) {
            const int ILL = arena.getFieldIdxHalo(
                fx - 2 * dx, fy - 2 * dy, feta - 2 * deta);
            const int IRR =
                arena.getFieldIdxHalo(fx + dx, fy + dy, feta + deta);
            const double aL =
                std::abs(arena.u_[direction][IL]) / arena.u_[0][IL];
            const double aR =
                std::abs(arena.u_[direction][IR]) / arena.u_[0][IR];
            const double ax = std::max(aL, aR);
            for (int idx = 0; idx < nComp; idx++) {
                const FieldComponent &W = *fluxComp[idx];
                const double fLL = W[ILL] * arena.u_[direction][ILL];
                const double fL = W[IL] * arena.u_[direction][IL];
                const double fR = W[IR] * arena.u_[direction][IR];
                const double fRR = W[IRR] * arena.u_[direction][IRR];
                const double gLL = W[ILL] * arena.u_[0][ILL];
                const double gL = W[IL] * arena.u_[0][IL];
                const double gR = W[IR] * arena.u_[0][IR];
                const double gRR = W[IRR] * arena.u_[0][IRR];
                const double uWR = fR - 0.5 * minmod.minmod_dx(fRR, fR, fL);
                const double uWL = fL + 0.5 * minmod.minmod_dx(fR, fL, fLL);
                const double WR = gR - 0.5 * minmod.minmod_dx(gRR, gR, gL);
                const double WL = gL + 0.5 * minmod.minmod_dx(gR, gL, gLL);
                v[idx] = ((uWR + uWL) - ax * (WR - WL)) * 0.5;
            }
        }
        double *s = v + kWComponents;
        for (int alpha = 0; alpha < 5; alpha++) {
            const int idx_1d = map_2d_idx_to_1d(alpha, direction);
            const double sL = arena.Wmunu_[idx_1d][IL] * tau_fac[direction];
            const double sR = arena.Wmunu_[idx_1d][IR] * tau_fac[direction];
            s[alpha] = (sR + sL) * 0.5;
        }
        if (!withBulk) return;
        const double PiL = arena.piBulk_[IL] + arena.piBulkChem_[IL];
        const double PiR = arena.piBulk_[IR] + arena.piBulkChem_[IR];
        for (int alpha = 0; alpha < 4; alpha++) {
            const double gfac1 = (alpha == direction ? 1.0 : 0.0);
            const double bL =
                PiL * (gfac1 + arena.u_[alpha][IL] * arena.u_[direction][IL])
                * tau_fac[direction];
            const double bR =
                PiR * (gfac1 + arena.u_[alpha][IR] * arena.u_[direction][IR])
                * tau_fac[direction];
            s[5 + alpha] = (bR + bL) * 0.5;
        }
    };

    // the partial_tau part of the source, written before the x faces
    auto timePart = [&](const int fieldIdx) {
        double *out = wSource.data() + fieldIdx;
        for (int alpha = 0; alpha < 5; alpha++) {
            const int idx_1d_alpha0 = map_2d_idx_to_1d(alpha, 0);
            double dWdtau =
                ((arenaCurr.Wmunu_[idx_1d_alpha0][fieldIdx]
                  - arenaPrev.Wmunu_[idx_1d_alpha0][fieldIdx])
                 / DATA.delta_tau);
            double dPidtau = 0.0;
            double Pi_alpha0 = 0.0;
            if (alpha < 4 && withBulk) {
                double gfac = (alpha == 0 ? -1.0 : 0.0);
                double PiCurr = arenaCurr.piBulk_[fieldIdx]
                                + arenaCurr.piBulkChem_[fieldIdx];
                double PiPrev = arenaPrev.piBulk_[fieldIdx]
                                + arenaPrev.piBulkChem_[fieldIdx];
                Pi_alpha0 = PiCurr
                            * (gfac
                               + arenaCurr.u_[alpha][fieldIdx]
                                     * arenaCurr.u_[0][fieldIdx]);
                dPidtau = (Pi_alpha0
                           - PiPrev
                                 * (gfac
                                    + arenaPrev.u_[alpha][fieldIdx]
                                          * arenaPrev.u_[0][fieldIdx]))
                          / DATA.delta_tau;
            }
            double sf =
                tau * dWdtau + arenaCurr.Wmunu_[idx_1d_alpha0][fieldIdx];
            double bf = tau * dPidtau + Pi_alpha0;
            double d = 0.;
            d += sf + bf;
            out[alpha * nStorage] = d;
            if (std::isnan(sf + bf)) {
                music_message << "[Error]Diss::computeWDivergences: "
                              << "sf=" << sf << " bf=" << bf;
                music_message.flush("error");
            }
        }
    };

    // the update of a cell from its faces at -1/2 (vm) and +1/2 (vp)
    auto cellUpdate = [&](const int direction, const int fieldIdx,
                          const double *vm, const double *vp) {
        if (withFlux) {
            for (int idx = 0; idx < nComp; idx++) {
                const double HW = (vp[idx] - vm[idx]) / deltaFlux[direction];
                double &div =
                    wFluxDiv[static_cast<size_t>(idx) * nStorage + fieldIdx];
                if (direction == 1) div = 0.;
                div += (-HW) * delta_tau;
            }
        }
        const double *sm = vm + kWComponents;
        const double *sp = vp + kWComponents;
        double *out = wSource.data() + fieldIdx;
        for (int alpha = 0; alpha < 5; alpha++) {
            const double W_p = sp[alpha];
            const double W_m = sm[alpha];
            double Pi_m = 0;
            double Pi_p = 0;
            if (withBulk && alpha < 4) {
                Pi_p = sp[5 + alpha];
                Pi_m = sm[5 + alpha];
            }
            if (direction == 3 && (alpha == 0 || alpha == 3)) {
                const int k = (alpha == 0 ? 0 : 1);
                out[(5 + k) * nStorage] = W_p + Pi_p;
                out[(7 + k) * nStorage] = W_m + Pi_m;
            } else {
                out[alpha * nStorage] +=
                    (W_p - W_m + Pi_p - Pi_m) / deltaSrc[direction];
            }
        }
    };

    // x: each thread walks whole x lines and keeps the last face
#pragma omp parallel for collapse(2) schedule(dynamic, 4)
    for (int ieta = 0; ieta < neta; ieta++)
        for (int iy = 0; iy < ny; iy++) {
            double faceA[nFaceVal], faceB[nFaceVal];
            double *vm = faceA;
            double *vp = faceB;
            faceValues(1, 0, iy, ieta, vm);
            for (int ix = 0; ix < nx; ix++) {
                faceValues(1, ix + 1, iy, ieta, vp);
                const int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                timePart(fieldIdx);
                cellUpdate(1, fieldIdx, vm, vp);
                std::swap(vm, vp);
            }
        }

    // y and eta: a block of x cells walks along the direction and keeps the
    // last row of faces, so the reads stay contiguous in x
    constexpr int xBlock = 32;
    const int nxBlocks = (nx + xBlock - 1) / xBlock;
    for (int direction = 2; direction <= 3; direction++) {
        const int nOuter = (direction == 2 ? neta : ny);
        const int nAlong = (direction == 2 ? ny : neta);
#pragma omp parallel
        {
            std::vector<double> rows(2 * xBlock * nFaceVal);
#pragma omp for collapse(2) schedule(dynamic, 1)
            for (int outer = 0; outer < nOuter; outer++)
                for (int xb = 0; xb < nxBlocks; xb++) {
                    const int x0 = xb * xBlock;
                    const int x1 = std::min(nx, x0 + xBlock);
                    double *rowM = rows.data();
                    double *rowP = rows.data() + xBlock * nFaceVal;
                    auto face = [&](const int f, const int ix, double *v) {
                        if (direction == 2) {
                            faceValues(2, ix, f, outer, v);
                        } else {
                            faceValues(3, ix, outer, f, v);
                        }
                    };
                    for (int ix = x0; ix < x1; ix++) {
                        face(0, ix, rowM + (ix - x0) * nFaceVal);
                    }
                    for (int i = 0; i < nAlong; i++) {
                        for (int ix = x0; ix < x1; ix++) {
                            double *vp = rowP + (ix - x0) * nFaceVal;
                            face(i + 1, ix, vp);
                            const int fieldIdx =
                                (direction == 2
                                     ? arena.getFieldIdx(ix, i, outer)
                                     : arena.getFieldIdx(ix, outer, i));
                            cellUpdate(
                                direction, fieldIdx,
                                rowM + (ix - x0) * nFaceVal, vp);
                        }
                        std::swap(rowM, rowP);
                    }
                }
        }
    }
}

void Diss::MakeWSourceFromDivergence(
    const std::vector<double> &wSource, const int nStorage, const int fieldIdx,
    TJbVec &dwmn) const {
    dwmn = {0.};
    for (int alpha = 0; alpha < 5; alpha++) {
        dwmn[alpha] = wSource[static_cast<size_t>(alpha) * nStorage + fieldIdx];
    }
    const double W_eta_p0 =
        wSource[5 * static_cast<size_t>(nStorage) + fieldIdx];
    const double W_eta_p3 =
        wSource[6 * static_cast<size_t>(nStorage) + fieldIdx];
    const double W_eta_m0 =
        wSource[7 * static_cast<size_t>(nStorage) + fieldIdx];
    const double W_eta_m3 =
        wSource[8 * static_cast<size_t>(nStorage) + fieldIdx];
    // the geometric terms, as in MakeWSource
    const double cosh_deta = wSourceCoshDeta_;
    const double sinh_deta = wSourceSinhDeta_;
    dwmn[0] +=
        ((W_eta_p0 - W_eta_m0) * cosh_deta + (W_eta_p3 + W_eta_m3) * sinh_deta);
    dwmn[3] +=
        ((W_eta_p3 - W_eta_m3) * cosh_deta + (W_eta_p0 + W_eta_m0) * sinh_deta);
}

void Diss::computeInverseReynoldsNumbers(
    const double enthalpy, const Cell_small &grid_pt, double &R_shear,
    double &R_bulk) const {
    double pi_00 = grid_pt.Wmunu[0];
    double pi_01 = grid_pt.Wmunu[1];
    double pi_02 = grid_pt.Wmunu[2];
    double pi_03 = grid_pt.Wmunu[3];
    double pi_11 = grid_pt.Wmunu[4];
    double pi_12 = grid_pt.Wmunu[5];
    double pi_13 = grid_pt.Wmunu[6];
    double pi_22 = grid_pt.Wmunu[7];
    double pi_23 = grid_pt.Wmunu[8];
    double pi_33 = grid_pt.Wmunu[9];

    double pisize =
        (pi_00 * pi_00 + pi_11 * pi_11 + pi_22 * pi_22 + pi_33 * pi_33
         - 2. * (pi_01 * pi_01 + pi_02 * pi_02 + pi_03 * pi_03)
         + 2. * (pi_12 * pi_12 + pi_13 * pi_13 + pi_23 * pi_23));

    const double h = std::max(enthalpy, Util::tiny_normal);
    R_shear = sqrt(pisize) / h;
    R_bulk = std::abs(grid_pt.pi_b) / h;
}

void Diss::Make_uWSource(
    const double tau, const Cell_small &grid_pt, const double theta_local,
    const DumuVec &a_local, const VelocityShearVec &sigma_1d,
    const VorticityVec &omega_1d, const ThermalVec &thermalVec,
    std::array<double, 5> &sourceTerms) {
    auto sigma = Util::UnpackVecToMatrix(sigma_1d);
    auto Wmunu = Util::UnpackVecToMatrix(grid_pt.Wmunu);

    double epsilon = thermalVec[0];
    double T = thermalVec[6];
    double muB = thermalVec[7];

    double shear_to_s = transport_coeffs_.get_eta_over_s(T, muB);

    bool include_WWterm = false;
    bool include_Wsigma_term = false;
    if (DATA.include_second_order_terms) {
        include_WWterm = true;
        include_Wsigma_term = true;
    }

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //                Defining transport coefficients                     //
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    double pressure = thermalVec[2];
    double shear = 0.;
    if (DATA.muB_dependent_shear_to_s == 0) {
        double entropy = thermalVec[12];
        shear = shear_to_s * entropy;
    } else {
        shear = shear_to_s * (epsilon + pressure) / std::max(T, small_eps);
    }
    double tau_pi =
        (transport_coeffs_.get_shear_relax_time_factor() * shear
         / std::max(epsilon + pressure, small_eps));
    tau_pi = std::min(10., std::max(3. * DATA.delta_tau, tau_pi));

    // transport coefficient for nonlinear terms -- shear only terms
    // transport coefficients of a massless gas of single component particles
    // phi7 pi pi enters as transport_coefficient * pi pi / tau_pi, with
    // transport_coefficient = phi7 (4/5) tau_pi / eta = phi7 (4/5) b_pi /
    // (e + P). Use the enthalpy form: tau_pi is clamped (>= 3 dtau) while
    // eta -> 0 in near-vacuum cells, and tau_pi / eta then diverges; a little
    // pi flowing into such a cell runs away to the regulator cap.
    double transport_coefficient =
        (transport_coeffs_.get_phi7_coeff()
         * transport_coeffs_.get_shear_relax_time_factor()
         / std::max(epsilon + pressure, small_eps) * (4. / 5.));
    double transport_coefficient2 =
        (transport_coeffs_.get_delta_pipi_coeff() * tau_pi);
    double transport_coefficient3 =
        (transport_coeffs_.get_tau_pipi_coeff() * tau_pi);

    // transport coefficient for nonlinear terms
    // -- coupling to bulk viscous pressure
    // transport coefficients not yet known -- fixed to zero
    double transport_coefficient_b =
        (transport_coeffs_.get_lambda_pibulkPi_coeff() * tau_pi);
    double transport_coefficient2_b = 0.;

    double resummedCorrection = 1;
    if (DATA.FlagResumTransportCoeff) {
        double R_shear = 0.;
        double R_bulk = 0.;
        computeInverseReynoldsNumbers(
            epsilon + pressure, grid_pt, R_shear, R_bulk);
        resummedCorrection =
            (transport_coeffs_.getResummedCorrFactor(R_shear, R_bulk));
        double resummedCorrSq = resummedCorrection * resummedCorrection;
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
        map_1d_idx_to_2d(idx_1d, mu, nu);
        /* This source has many terms */
        /* everything in the 1/(tau_pi) piece is here */
        /* third step in the split-operator time evol
           use Wmunu[rk_flag] and u[rk_flag] with rk_flag = 0 */

        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        //           Wmunu + transport_coefficient2*Wmunu*theta              //
        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////

        // full term is
        double tempf =
            (-(1.0 + transport_coefficient2 * theta_local) * Wmunu[mu][nu]);

        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        //           Navier-Stokes Term -- -2.*shear*sigma^munu              //
        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////

        // full Navier-Stokes term is
        // sign changes according to metric sign convention
        double NS_term = -2. * shear * sigma[mu][nu];

        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        //              Add coupling to bulk viscous pressure                //
        //             transport_coefficient_b*Bulk*sigma^mu nu              //
        //              transport_coefficient2_b*Bulk*W^mu nu                //
        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        double Coupling_to_Bulk = 0.0;
        if (DATA.include_second_order_terms) {
            double Bulk_Sigma = grid_pt.pi_b * sigma[mu][nu];
            double Bulk_W = grid_pt.pi_b * Wmunu[mu][nu];

            // multiply term by its respective transport coefficient
            double Bulk_Sigma_term = Bulk_Sigma * transport_coefficient_b;
            double Bulk_W_term = Bulk_W * transport_coefficient2_b;

            // full term is
            // first term: sign changes according to metric sign convention
            Coupling_to_Bulk = -Bulk_Sigma_term + Bulk_W_term;
        }

        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        //                      Vorticity Term                               //
        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        double Vorticity_term = 0.0;
        if (DATA.include_vorticity_terms) {
            double transport_coefficient4 = 2. * tau_pi;
            if (DATA.FlagResumTransportCoeff) {
                transport_coefficient4 *= resummedCorrection;
            }
            auto omega = Util::UnpackVecToMatrix(omega_1d);
            double term1_Vorticity =
                (-Wmunu[mu][0] * omega[nu][0] - Wmunu[nu][0] * omega[mu][0]
                 + Wmunu[mu][1] * omega[nu][1] + Wmunu[nu][1] * omega[mu][1]
                 + Wmunu[mu][2] * omega[nu][2] + Wmunu[nu][2] * omega[mu][2]
                 + Wmunu[mu][3] * omega[nu][3] + Wmunu[nu][3] * omega[mu][3])
                / 2.;
            // multiply term by its respective transport coefficient
            Vorticity_term = transport_coefficient4 * term1_Vorticity;
        }

        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        //              Add nonlinear term in shear-stress tensor            //
        //  transport_coefficient3                                           //
        //              *Delta(mu nu)(alpha beta)*Wmu gamma sigma nu gamma   //
        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        double Wsigma_term = 0.0;
        if (include_Wsigma_term) {
            double Wsigma =
                (Wmunu[0][0] * sigma[0][0] + Wmunu[1][1] * sigma[1][1]
                 + Wmunu[2][2] * sigma[2][2] + Wmunu[3][3] * sigma[3][3]
                 - 2.
                       * (Wmunu[0][1] * sigma[0][1] + Wmunu[0][2] * sigma[0][2]
                          + Wmunu[0][3] * sigma[0][3])
                 + 2.
                       * (Wmunu[1][2] * sigma[1][2] + Wmunu[1][3] * sigma[1][3]
                          + Wmunu[2][3] * sigma[2][3]));
            double term1_Wsigma =
                (-Wmunu[mu][0] * sigma[nu][0] - Wmunu[nu][0] * sigma[mu][0]
                 + Wmunu[mu][1] * sigma[nu][1] + Wmunu[nu][1] * sigma[mu][1]
                 + Wmunu[mu][2] * sigma[nu][2] + Wmunu[nu][2] * sigma[mu][2]
                 + Wmunu[mu][3] * sigma[nu][3] + Wmunu[nu][3] * sigma[mu][3])
                / 2.;

            double term2_Wsigma =
                (-(1. / 3.)
                 * (DATA.gmunu[mu][nu] + grid_pt.u[mu] * grid_pt.u[nu])
                 * Wsigma);
            // multiply term by its respective transport coefficient
            term1_Wsigma = transport_coefficient3 * term1_Wsigma;
            term2_Wsigma = transport_coefficient3 * term2_Wsigma;

            // full term is
            Wsigma_term = -term1_Wsigma - term2_Wsigma;
        }

        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        //              Add nonlinear term in shear-stress tensor            //
        // transport_coefficient*Delta(mu nu)(alpha beta)*Wmu gamma Wnu gamma//
        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        double WW_term = 0.0;
        if (include_WWterm) {
            double Wsquare =
                (Wmunu[0][0] * Wmunu[0][0] + Wmunu[1][1] * Wmunu[1][1]
                 + Wmunu[2][2] * Wmunu[2][2] + Wmunu[3][3] * Wmunu[3][3]
                 - 2.
                       * (Wmunu[0][1] * Wmunu[0][1] + Wmunu[0][2] * Wmunu[0][2]
                          + Wmunu[0][3] * Wmunu[0][3])
                 + 2.
                       * (Wmunu[1][2] * Wmunu[1][2] + Wmunu[1][3] * Wmunu[1][3]
                          + Wmunu[2][3] * Wmunu[2][3]));
            double term1_WW =
                (-Wmunu[mu][0] * Wmunu[nu][0] + Wmunu[mu][1] * Wmunu[nu][1]
                 + Wmunu[mu][2] * Wmunu[nu][2] + Wmunu[mu][3] * Wmunu[nu][3]);
            double term2_WW =
                (-(1. / 3.)
                 * (DATA.gmunu[mu][nu] + grid_pt.u[mu] * grid_pt.u[nu])
                 * Wsquare);

            // multiply term by its respective transport coefficient
            term1_WW = term1_WW * transport_coefficient;
            term2_WW = term2_WW * transport_coefficient;

            // full term is
            // sign changes according to metric sign convention
            WW_term = -term1_WW - term2_WW;
        }

        // final answer is
        sourceTerms[idx_1d - 4] =
            ((NS_term + tempf + Vorticity_term + Wsigma_term + WW_term
              + Coupling_to_Bulk)
             / tau_pi);
    }
}

int Diss::nWComponents() const {
    int loopIdx = 5;
    if (DATA.turn_on_bulk) {
        loopIdx = 6;
    }
    if (DATA.turn_on_bulk_chem) {
        loopIdx = 7;
    }
    if (DATA.turn_on_diff) {
        loopIdx = 10;
    }
    return (loopIdx);
}

/* Kurganov-Tadmor for Wmunu */
/* implement
   partial_tau (utau Wmn) + (1/tau)partial_eta (ueta Wmn)
   + partial_x (ux Wmn) + partial_y (uy Wmn) + utau Wmn/tau = SW
   or the right hand side of,
   partial_tau (utau Wmn) =
                    - (1/tau)partial_eta (ueta Wmn)
                    - partial_x (ux Wmn) - partial_y (uy Wmn)
                    - utau Wmn/tau + SW*/

/* the local velocity is just u_x/u_tau, u_y/u_tau, u_eta/tau/u_tau */
/* KT flux is given by
   H_{j+1/2} = (fRph + fLph)/2 - ax(uRph - uLph)
   Here fRph = ux WmnRph and ax uRph = |ux/utau|_max utau Wmn */
void Diss::computeWFluxDivergence(
    const double tau, const Fields &arena, std::vector<double> &wFluxDiv) {
    // the W components, with 9 for piBulk and 14 for piBulkChem
    static const int piIdxArr[kWComponents] = {4, 5,  6,  7,  8,
                                               9, 14, 11, 12, 13};
    const int nComp = nWComponents();
    const int nx = arena.nX();
    const int ny = arena.nY();
    const int neta = arena.nEta();
    const int nStorage = arena.nStorage();
    const double delta[4] = {
        0.0, DATA.delta_x, DATA.delta_y, DATA.delta_eta * tau};
    const double delta_tau = DATA.delta_tau;
    wFluxDiv.resize(static_cast<size_t>(kWComponents) * nStorage);

    auto component = [&](const int idx_1d) -> const FieldComponent & {
        if (idx_1d == 9) return (arena.piBulk_);
        if (idx_1d == 14) return (arena.piBulkChem_);
        return (arena.Wmunu_[idx_1d]);
    };

    for (int direction = 1; direction <= 3; direction++) {
        const int dx = (direction == 1);
        const int dy = (direction == 2);
        const int deta = (direction == 3);
        const int nfx = nx + dx;
        const int nfy = ny + dy;
        const int nfeta = neta + deta;
        const int nFaces = nfx * nfy * nfeta;
        wFaceFluxes_.resize(static_cast<size_t>(nComp) * nFaces);

        // the flux through the face between the cells L = f - 1 and R = f
#pragma omp parallel for collapse(3) schedule(dynamic, 64)
        for (int feta = 0; feta < nfeta; feta++)
            for (int fy = 0; fy < nfy; fy++)
                for (int fx = 0; fx < nfx; fx++) {
                    const int ILL = arena.getFieldIdxHalo(
                        fx - 2 * dx, fy - 2 * dy, feta - 2 * deta);
                    const int IL =
                        arena.getFieldIdxHalo(fx - dx, fy - dy, feta - deta);
                    const int IR = arena.getFieldIdxHalo(fx, fy, feta);
                    const int IRR =
                        arena.getFieldIdxHalo(fx + dx, fy + dy, feta + deta);
                    const int faceIdx = fx + nfx * (fy + nfy * feta);

                    const double aL =
                        std::abs(arena.u_[direction][IL]) / arena.u_[0][IL];
                    const double aR =
                        std::abs(arena.u_[direction][IR]) / arena.u_[0][IR];
                    const double ax = std::max(aL, aR);
                    for (int idx = 0; idx < nComp; idx++) {
                        const FieldComponent &W = component(piIdxArr[idx]);
                        // u^i W and u^tau W at the four cells
                        const double fLL = W[ILL] * arena.u_[direction][ILL];
                        const double fL = W[IL] * arena.u_[direction][IL];
                        const double fR = W[IR] * arena.u_[direction][IR];
                        const double fRR = W[IRR] * arena.u_[direction][IRR];
                        const double gLL = W[ILL] * arena.u_[0][ILL];
                        const double gL = W[IL] * arena.u_[0][IL];
                        const double gR = W[IR] * arena.u_[0][IR];
                        const double gRR = W[IRR] * arena.u_[0][IRR];

                        /* MakeuWmnHalfs */
                        const double uWR =
                            fR - 0.5 * minmod.minmod_dx(fRR, fR, fL);
                        const double uWL =
                            fL + 0.5 * minmod.minmod_dx(fR, fL, fLL);
                        const double WR =
                            gR - 0.5 * minmod.minmod_dx(gRR, gR, gL);
                        const double WL =
                            gL + 0.5 * minmod.minmod_dx(gR, gL, gLL);
                        wFaceFluxes_
                            [static_cast<size_t>(idx) * nFaces + faceIdx] =
                                ((uWR + uWL) - ax * (WR - WL)) * 0.5;
                    }
                }

        // make partial_i (u^i Wmn), accumulated over the directions in order
#pragma omp parallel for collapse(3) schedule(dynamic, 64)
        for (int ieta = 0; ieta < neta; ieta++)
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++) {
                    const int fieldIdx = arena.getFieldIdx(ix, iy, ieta);
                    const int Imh = ix + nfx * (iy + nfy * ieta);
                    const int Iph =
                        (ix + dx) + nfx * ((iy + dy) + nfy * (ieta + deta));
                    for (int idx = 0; idx < nComp; idx++) {
                        const size_t off = static_cast<size_t>(idx) * nFaces;
                        const double HW =
                            (wFaceFluxes_[off + Iph] - wFaceFluxes_[off + Imh])
                            / delta[direction];
                        double &div = wFluxDiv
                            [static_cast<size_t>(idx) * nStorage + fieldIdx];
                        if (direction == 1) div = 0.;
                        div += (-HW) * delta_tau;
                    }
                }
    }
}

void Diss::Make_uWRHS(
    const double tau, Fields &arena, const int fieldIdx,
    const std::vector<double> &wFluxDiv, std::array<double, 10> &w_rhs,
    const double theta_local, const DumuVec &a_local) {
    auto grid_pt = arena.getCell(fieldIdx);
    auto Wmunu_local = Util::UnpackVecToMatrix(grid_pt.Wmunu);

    /* This is the second step in the operator splitting. it uses
       rk_flag+1 as initial condition */
    // the KT flux part, from computeWFluxDivergence
    const int nStorage = arena.nStorage();
    const int nComp = nWComponents();
    for (int idx = 0; idx < nComp; idx++) {
        w_rhs[idx] += wFluxDiv[static_cast<size_t>(idx) * nStorage + fieldIdx];
    }

    /* add a source term -u^tau Wmn/tau
       due to the coordinate change to tau-eta */
    /* this is from udW = d(uW) - Wdu = RHS */
    /* or d(uW) = udW + Wdu */
    /* this term is being added to the rhs so that -4/3 + 1 = -1/3 */
    /* other source terms due to the coordinate change to tau-eta */

    // align gmunu in data for faster access - also changed **gmunu in data to
    // gmunu[4][4] moved two sums into w_rhs at top and bottom into one sum in
    // the end do not symmetrize in the end, just go through all nu's vectorized
    // innermost loop more efficiently by iterating over 4 indices instead of 3
    // to avoid masking

    int mu = 0;
    int nu = 0;
    for (int idx_1d = 4; idx_1d < 9; idx_1d++) {
        map_1d_idx_to_2d(idx_1d, mu, nu);
        double tempf =
            (-(DATA.gmunu[3][mu]) * (Wmunu_local[0][nu])
             - (DATA.gmunu[3][nu]) * (Wmunu_local[0][mu])
             + (DATA.gmunu[0][mu]) * (Wmunu_local[3][nu])
             + (DATA.gmunu[0][nu]) * (Wmunu_local[3][mu])
             + (Wmunu_local[3][nu]) * (grid_pt.u[mu]) * (grid_pt.u[0])
             + (Wmunu_local[3][mu]) * (grid_pt.u[nu]) * (grid_pt.u[0])
             - (Wmunu_local[0][nu]) * (grid_pt.u[mu]) * (grid_pt.u[3])
             - (Wmunu_local[0][mu]) * (grid_pt.u[nu]) * (grid_pt.u[3]))
            * (grid_pt.u[3] / tau);

        for (int ic = 0; ic < 4; ic++) {
            const double ic_fac = (ic == 0 ? -1.0 : 1.0);
            tempf +=
                ((Wmunu_local[ic][nu]) * (grid_pt.u[mu]) * (a_local[ic])
                     * ic_fac
                 + (Wmunu_local[ic][mu]) * (grid_pt.u[nu]) * (a_local[ic])
                       * ic_fac);
        }

        w_rhs[idx_1d - 4] +=
            (tempf * (DATA.delta_tau)
             + (-(grid_pt.u[0] * Wmunu_local[mu][nu]) / tau
                + (theta_local * Wmunu_local[mu][nu]))
                   * (DATA.delta_tau));
    }
    if (DATA.turn_on_bulk) {
        w_rhs[5] -= (grid_pt.pi_b) * (grid_pt.u[0]) / tau * DATA.delta_tau;
        w_rhs[5] += (grid_pt.pi_b) * theta_local * DATA.delta_tau;
    }

    if (DATA.turn_on_bulk_chem) {
        w_rhs[6] -= (grid_pt.pi_b_chem) * (grid_pt.u[0]) / tau * DATA.delta_tau;
        w_rhs[6] += (grid_pt.pi_b_chem) * theta_local * DATA.delta_tau;
    }
}

double Diss::Make_uPiSource(
    const double tau, const Cell_small &grid_pt, const double theta_local,
    const VelocityShearVec &sigma_1d, const ThermalVec &thermalVec) {
    double tempf;
    double bulk;
    double Bulk_Relax_time = 0.2;  // fm
    double transport_coeff1, transport_coeff2;
    double transport_coeff1_s, transport_coeff2_s;
    double NS_term, BB_term;
    double Final_Answer;

    // switch to include non-linear coupling terms in the bulk pi evolution
    bool include_BBterm = false;
    bool include_coupling_to_shear = false;
    if (DATA.include_second_order_terms) {
        include_BBterm = true;
        include_coupling_to_shear = true;
    }

    // defining bulk viscosity coefficient
    // shear viscosity = constant * entropy density
    // s_den = eos.get_entropy(epsilon, rhob);
    // shear = (DATA.shear_to_s)*s_den;
    // shear viscosity = constant * (e + P)/T
    double epsilon = thermalVec[0];
    double temperature = thermalVec[6];
    double mu_B = thermalVec[7];

    // cs2 is the velocity of sound squared
    // double cs2 = eos.get_cs2(epsilon, thermalVec[1]);
    double cs2 = thermalVec[5];
    double pressure = thermalVec[2];

    // T dependent bulk viscosity
    bulk = transport_coeffs_.get_zeta_over_s(temperature, mu_B);
    // guard T -> 0 (e.g. cells at the energy floor), as for the shear
    bulk = bulk * (epsilon + pressure) / std::max(temperature, small_eps);

    // defining bulk relaxation time and additional transport coefficients
    // Bulk relaxation time from kinetic theory
    if (DATA.bulk_relaxation_type == 0) {
        double csfactor = std::max(1. / 3. - cs2, small_eps);
        Bulk_Relax_time =
            (transport_coeffs_.get_bulk_relax_time_factor() * bulk
             / (csfactor * csfactor) / std::max(epsilon + pressure, small_eps));
    } else if (DATA.bulk_relaxation_type == 1) {
        double csfactor = std::max(1. / 3. - cs2, small_eps);
        Bulk_Relax_time =
            (transport_coeffs_.get_bulk_relax_time_factor() * bulk / csfactor
             / std::max(epsilon + pressure, small_eps));
    } else if (DATA.bulk_relaxation_type == 2) {
        Bulk_Relax_time =
            (transport_coeffs_.get_bulk_relax_time_factor() * bulk
             / std::max(epsilon + pressure, small_eps));
    }

    // avoid overflow or underflow of the bulk relaxation time
    Bulk_Relax_time =
        (std::min(10., std::max(3. * DATA.delta_tau, Bulk_Relax_time)));

    // from kinetic theory, small mass limit
    transport_coeff1 =
        (transport_coeffs_.get_delta_bulkPibulkPi_coeff() * Bulk_Relax_time);
    transport_coeff2 =
        (transport_coeffs_.get_tau_bulkPibulkPi_coeff() * Bulk_Relax_time);

    // from kinetic theory
    transport_coeff1_s =
        (transport_coeffs_.get_lambda_bulkPipi_coeff() * (1. / 3. - cs2)
         * Bulk_Relax_time);
    if (DATA.bulk_relaxation_type == 2) {
        transport_coeff1_s =
            (transport_coeffs_.get_lambda_bulkPipi_coeff() * 1. / 3.
             * Bulk_Relax_time);
    }
    transport_coeff2_s = 0.;  // not known;  put 0

    if (DATA.FlagResumTransportCoeff) {
        double resummedCorrection = 1;
        double R_shear = 0.;
        double R_bulk = 0.;
        computeInverseReynoldsNumbers(
            epsilon + pressure, grid_pt, R_shear, R_bulk);
        resummedCorrection =
            (transport_coeffs_.getResummedCorrFactor(R_shear, R_bulk));
        bulk *= resummedCorrection;
        transport_coeff1 *= resummedCorrection;
        transport_coeff2 *= resummedCorrection;
        transport_coeff1_s *= resummedCorrection;
        transport_coeff2_s *= resummedCorrection;
    }

    // Computing Navier-Stokes term (-bulk viscosity * theta)
    NS_term = -bulk * theta_local;

    // Computing relaxation term and nonlinear term:
    // - Bulk - transport_coeff1*Bulk*theta
    tempf = (-(grid_pt.pi_b) - transport_coeff1 * theta_local * (grid_pt.pi_b));

    // Computing nonlinear term: + transport_coeff2*Bulk*Bulk
    if (include_BBterm == 1) {
        BB_term = transport_coeff2 * (grid_pt.pi_b) * (grid_pt.pi_b);
    } else {
        BB_term = 0.0;
    }

    // Computing terms that Couple with shear-stress tensor
    double Wsigma, WW, Shear_Sigma_term, Shear_Shear_term, Coupling_to_Shear;

    if (include_coupling_to_shear) {
        auto sigma = Util::UnpackVecToMatrix(sigma_1d);
        auto Wmunu = Util::UnpackVecToMatrix(grid_pt.Wmunu);

        Wsigma =
            (Wmunu[0][0] * sigma[0][0] + Wmunu[1][1] * sigma[1][1]
             + Wmunu[2][2] * sigma[2][2] + Wmunu[3][3] * sigma[3][3]
             - 2.
                   * (Wmunu[0][1] * sigma[0][1] + Wmunu[0][2] * sigma[0][2]
                      + Wmunu[0][3] * sigma[0][3])
             + 2.
                   * (Wmunu[1][2] * sigma[1][2] + Wmunu[1][3] * sigma[1][3]
                      + Wmunu[2][3] * sigma[2][3]));

        WW =
            (Wmunu[0][0] * Wmunu[0][0] + Wmunu[1][1] * Wmunu[1][1]
             + Wmunu[2][2] * Wmunu[2][2] + Wmunu[3][3] * Wmunu[3][3]
             - 2.
                   * (Wmunu[0][1] * Wmunu[0][1] + Wmunu[0][2] * Wmunu[0][2]
                      + Wmunu[0][3] * Wmunu[0][3])
             + 2.
                   * (Wmunu[1][2] * Wmunu[1][2] + Wmunu[1][3] * Wmunu[1][3]
                      + Wmunu[2][3] * Wmunu[2][3]));
        // multiply term by its respective transport coefficient
        Shear_Sigma_term = Wsigma * transport_coeff1_s;
        Shear_Shear_term = WW * transport_coeff2_s;

        // full term that couples to shear is
        Coupling_to_Shear = -Shear_Sigma_term + Shear_Shear_term;
    } else {
        Coupling_to_Shear = 0.0;
    }

    // Final Answer
    Final_Answer = NS_term + tempf + BB_term + Coupling_to_Shear;

    return Final_Answer / (Bulk_Relax_time);
} /* Make_uPiSource */

/* Sangyong Nov 18 2014 */
/* baryon current parts */
/* this contains the source terms
   that is, all the terms that are not part of the current */
/* for the q part, we don't do tau*u*q we just do u*q
   this part contains
    -(1/tau_rho)(q[a] + kappa g[a][b]Dtildemu[b]
                 + kappa u[a] u[b]g[b][c]Dtildemu[c])
    +Delta[a][tau] u[eta] q[eta]/tau
    -Delta[a][eta] u[eta] q[tau]/tau
    -u[a]u[b]g[b][e] Dq[e]
*/

double Diss::Make_uPiChemSource(
    const double tau, const Cell_small &grid_pt, const double theta_local,
    const VelocityShearVec & /*sigma_1d*/, const ThermalVec &thermalVec) {
    // --- read thermal variables  ---
    const double epsilon = thermalVec[0];
    const double pressure = thermalVec[2];
    const double cs2 = thermalVec[5];
    const double temperature = thermalVec[6];  // 1/fm in MUSIC units

    // Solve Y_q using I_QCD and Pi_chem ---
    // I_QCD = e - 3P
    const double Iqcd = (epsilon - 3.0 * pressure);
    const double Iqcd_safe =
        std::copysign(std::max(std::abs(Iqcd), small_eps), Iqcd);

    // Eq.(14) with Ieq = I_QCD and I0 = 0:
    // Pi_chem = (1/3)(1 - sqrt(Y_q))[Ieq(e) - I0(e)]
    // => sqrt(Y_q) = 1 - 3 Pi_chem / Ieq(e), clamp to [0, 1]
    const double y_q_inner = 1.0 - 3.0 * grid_pt.pi_b_chem / Iqcd_safe;
    const double y_q_sqrt = std::min(1.0, std::abs(y_q_inner));
    const double y_q_raw = y_q_inner * y_q_inner;

    // Eq.(18)/(20): tau_Pi_chem = 2 / [R (1 + sqrt(Y_q))(sqrt(2) + sqrt(Y_q))]
    //  Provide DATA.chem_rate_C constant C in R = C*T
    const double chem_rate_C = std::max(DATA.chem_rate_C, small_eps);
    double tau_Pi_chem = 2.0
                         / std::max(
                             chem_rate_C * temperature * (1.0 + y_q_sqrt)
                                 * (std::sqrt(2.0) + y_q_sqrt),
                             small_eps);

    tau_Pi_chem = std::min(10.0, std::max(3.0 * DATA.delta_tau, tau_Pi_chem));

    // Eq.(19) with I0 = 0 and I(e,Y_q) = sqrt(Y_q) * Ieq(e):
    // zeta_chem = tau_Pi_chem * (1/3) * (1 - sqrt(Y_q)) * (dIeq/de) * (e + P)
    // with Ieq = I_QCD = e - 3P => dIeq/de = 1 - 3 cs^2
    const double dI3de = (1.0 - 3.0 * cs2);
    const double zeta_raw = tau_Pi_chem * (1.0 / 3.0) * (1.0 - y_q_sqrt) * dI3de
                            * (epsilon + pressure);
    const double zeta = std::max(0.0, zeta_raw);

    // Israel–Stewart form
    // D Pi_chem = [ - zeta*theta - Pi_chem ] / tau_Pi_chem
    const double NS_term = -zeta * theta_local;
    const double relax_term = -(grid_pt.pi_b_chem);

    if (y_q_raw > 1.0 + 1e-6 && echo_level_ > 5) {
        music_message << "[chem-bulk] Y_q>1 at tau=" << tau
                      << " Y_q_raw=" << y_q_raw << " sqrt(Y_q)=" << y_q_sqrt
                      << " Pi_chem=" << grid_pt.pi_b_chem << " Iqcd=" << Iqcd
                      << " T=" << temperature * hbarc << " GeV";
        music_message.flush("warning");
    }

    if (zeta_raw < -1e-12 && echo_level_ > 5) {
        music_message << "[chem-bulk] zeta_chem < 0 at tau=" << tau
                      << " zeta_chem_raw=" << zeta_raw << " cs2=" << cs2
                      << " sqrt(Y_q)=" << y_q_sqrt;
        music_message.flush("warning");
    }

    return (NS_term + relax_term) / tau_Pi_chem;
}

//! D(muB/T) = u^mu partial_mu (muB/T) from the ideal equations of motion,
//! De = -(e + P) theta and Dn = -n theta:
//! D(muB/T) = -theta [(e + P) d(muB/T)/de|_n + n d(muB/T)/dn|_e]
double Diss::D_muB_over_T_ideal(
    const double epsilon, const double rhob, const double pressure,
    const double theta) const {
    auto alpha = [this](const double e, const double n) {
        return (
            eos.get_muB(e, n) / std::max(eos.get_temperature(e, n), small_eps));
    };
    if (epsilon < small_eps) return 0.;
    const double de = 1e-4 * epsilon;
    const double dalpha_de =
        (alpha(epsilon + de, rhob) - alpha(epsilon - de, rhob)) / (2. * de);
    double dalpha_dn = 0.;
    if (std::abs(rhob) > small_eps) {
        const double dn = 1e-4 * std::abs(rhob);
        dalpha_dn =
            (alpha(epsilon, rhob + dn) - alpha(epsilon, rhob - dn)) / (2. * dn);
    }
    const double Dalpha =
        -theta * ((epsilon + pressure) * dalpha_de + rhob * dalpha_dn);
    return std::isfinite(Dalpha) ? Dalpha : 0.;
}

void Diss::Make_uqSource(
    const double tau, const Cell_small &grid_pt, const double theta_local,
    const DumuVec &a_local, const VelocityShearVec &sigma_1d,
    const VorticityVec &omega_1d, const DmuMuBoverTVec &baryon_diffusion_vec,
    const ThermalVec &thermalVec, std::array<double, 3> &sourceTerms) {
    double epsilon = thermalVec[0];
    double rhob = thermalVec[1];
    double pressure = thermalVec[2];
    double T = thermalVec[6];

    double kappa_coefficient = DATA.kappa_coefficient;
    double tau_rho = kappa_coefficient / std::max(T, small_eps);
    tau_rho = std::min(10., std::max(3. * DATA.delta_tau, tau_rho));

    double mub = thermalVec[7];
    double alpha = mub / std::max(T, small_eps);
    double denorm_safe = std::copysign(
        std::max(std::abs(3. * T * tanh(alpha)), small_eps),
        3. * T * tanh(alpha));
    double kappa = kappa_coefficient
                   * (rhob / denorm_safe
                      - rhob * rhob / std::max(epsilon + pressure, small_eps));

    if (DATA.Initial_profile == 1) {
        // for 1+1D numerical test
        double denorm_safe =
            std::copysign(std::max(std::abs(mub), small_eps), mub);
        kappa = kappa_coefficient * (rhob / denorm_safe);
    }
    // The kinetic-theory kappa turns negative where e + P < 3 n T tanh(muB/T)
    // (cold, baryon-rich cells with less than a nucleon mass of energy per
    // baryon, e.g. between the lumps of a lumpy TATB event). kappa < 0 is
    // anti-diffusion, which is ill-posed and blows up.
    kappa = std::max(0., kappa);

    // copy the value of \tilde{q^\mu}
    double q[4];
    for (int i = 0; i < 4; i++) {
        q[i] = grid_pt.Wmunu[10 + i];
    }

    /* -(1/tau_rho)(q[a] + kappa g[a][b]Dtildemu[b]
     *              + kappa u[a] u[b]g[b][c]Dtildemu[c])
     * + theta q[a] - q[a] u^\tau/tau
     * + Delta[a][tau] u[eta] q[eta]/tau
     * - Delta[a][eta] u[eta] q[tau]/tau
     * - u[a] u[b]g[b][e] Dq[e] -> u[a] q[e] g[e][b] Du[b]
     */

    // first: (1/tau_rho) part
    // recall that dUsup[4][i] = partial_i (muB/T)
    // and dUsup[4][0] = -partial_tau (muB/T) = partial^tau (muB/T)
    // and a[4] = u^a partial_a (muB/T) = DmuB/T
    // -(1/tau_rho)(q[a] + kappa g[a][b]DmuB/T[b]
    // + kappa u[a] u[b]g[b][c]DmuB/T[c])
    // a = nu

    double transport_coeff = transport_coeffs_.get_delta_qq_coeff() * tau_rho;
    double transport_coeff_2 =
        transport_coeffs_.get_lambda_qq_coeff() * tau_rho;

    // add a new non-linear term (-q^\mu \sigma_\mu\nu)
    auto sigma = Util::UnpackVecToMatrix(sigma_1d);

    // D(muB/T) from the ideal equations of motion instead of a_local[4],
    // which uses a backward time difference of muB/T. The time and spatial
    // parts of u^a partial_a (muB/T) nearly cancel where muB/T is large and
    // the lagged difference cannot track the limited spatial derivative:
    // q then carries time-step-scale noise that grows as dtau -> 0 (and
    // spoils N_B conservation in lumpy events). The two forms differ by
    // terms of second order in gradients (partial_mu q^mu, pi^{mu nu}
    // sigma_{mu nu}, ...).
    const double DmuBoverT =
        D_muB_over_T_ideal(epsilon, rhob, pressure, theta_local);

    for (int nu = 1; nu < 4; nu++) {
        double NS =
            kappa * (baryon_diffusion_vec[nu] + grid_pt.u[nu] * DmuBoverT);

        // add a new non-linear term (- q \theta)
        double Nonlinear1 = -transport_coeff * q[nu] * theta_local;

        double temptemp = 0.0;
        for (int i = 0; i < 4; i++) {
            temptemp += q[i] * sigma[i][nu] * DATA.gmunu[i][i];
        }
        double Nonlinear2 = -transport_coeff_2 * temptemp;

        // add a new non-linear term (-q_\mu \omega^{\mu\nu})
        double Nonlinear3 = 0.0;
        if (DATA.include_vorticity_terms) {
            double transport_coeff_3 = 1.0 * tau_rho;
            auto omega = Util::UnpackVecToMatrix(omega_1d);
            double temp3 = 0.0;
            for (int i = 0; i < 4; i++) {
                temp3 += q[i] * omega[i][nu] * DATA.gmunu[i][i];
            }
            Nonlinear3 = -transport_coeff_3 * temp3;
        }

        sourceTerms[nu - 1] =
            ((-q[nu] - NS + Nonlinear1 + Nonlinear2 + Nonlinear3) / tau_rho);
        if (DATA.Initial_profile == 1) {
            // for 1+1D numerical test
            sourceTerms[nu - 1] = (-q[nu] - NS) / tau_rho;
        }

        // all other geometric terms....
        // + theta q[a] - q[a] u^\tau/tau
        sourceTerms[nu - 1] += (theta_local - grid_pt.u[0] / tau) * q[nu];
        // if (isnan(SW)) {
        //     cout << "theta term is nan! " << endl;
        // }

        // +Delta[a][tau] u[eta] q[eta]/tau
        double tempf =
            ((DATA.gmunu[nu][0] + grid_pt.u[nu] * grid_pt.u[0]) * grid_pt.u[3]
                 * q[3] / tau
             - (DATA.gmunu[nu][3] + grid_pt.u[nu] * grid_pt.u[3]) * grid_pt.u[3]
                   * q[0] / tau);
        sourceTerms[nu - 1] += tempf;
        // if (isnan(tempf)) {
        //     cout << "Delta^{a \tau} and Delta^{a \eta} terms are nan!" <<
        //     endl;
        // }

        // -u[a] u[b]g[b][e] Dq[e] -> u[a] (q[e] g[e][b] Du[b])
        tempf = 0.0;
        for (int i = 0; i < 4; i++) {
            tempf += q[i] * Util::gmn(i) * a_local[i];
        }
        sourceTerms[nu - 1] += (grid_pt.u[nu]) * tempf;
        // if (isnan(tempf)) {
        //     cout << "u^a q_b Du^b term is nan! " << endl;
        // }
    }
}

//! this function outputs the T and muB dependence of the baryon diffusion
//! coefficient, kappa
void Diss::output_kappa_T_and_muB_dependence() {
    music_message.info("output kappa_B(T, mu_B) ...");
    std::ofstream of("kappa_B_T_and_muB_dependence.dat");
    // write out the header of the file
    of << "# e (GeV/fm^3)  rhob (1/fm^3) T (GeV)  mu_B (GeV)  kappa (1/fm^2)"
       << std::endl;

    // define the grid
    double e_min = 1e-5;   // fm^-4
    double e_max = 100.0;  // fm^-4
    int ne = 1000;
    double de = (e_max - e_min) / (ne - 1.);
    double rhob_min = 0.0;         // fm^-3
    double rhob_max = sqrt(10.0);  // fm^-3
    int nrhob = 1000;
    double drhob = (rhob_max - rhob_min) / (nrhob - 1.);

    for (int i = 0; i < ne; i++) {
        double e_local = e_min + i * de;
        for (int j = 0; j < nrhob; j++) {
            double rhob_local = rhob_min + j * drhob;
            rhob_local *= rhob_local;
            double mu_B_local = eos.get_muB(e_local, rhob_local);
            if (mu_B_local * hbarc > 0.78)
                continue;  // discard points out of the table
            double p_local = eos.get_pressure(e_local, rhob_local);
            double T_local = eos.get_temperature(e_local, rhob_local);
            double alpha_local =
                mu_B_local / std::max(T_local, Util::tiny_normal);

            double denorm_safe = std::copysign(
                std::max(small_eps, std::abs(3. * T_local * tanh(alpha_local))),
                3. * T_local * tanh(alpha_local));
            double kappa_local =
                (DATA.kappa_coefficient
                 * (rhob_local / denorm_safe
                    - rhob_local * rhob_local
                          / std::max(e_local + p_local, small_eps)));
            // output
            of << std::scientific << std::setw(18) << std::setprecision(8)
               << e_local * hbarc << "   " << rhob_local << "   "
               << T_local * hbarc << "   " << mu_B_local * hbarc << "   "
               << kappa_local << std::endl;
        }
    }
    of.close();  // close the file
}

//! this function outputs the T and muB dependence of the baryon diffusion
//! coefficient, kappa_B, along constant s/n_B trajectories
void Diss::output_kappa_along_const_sovernB() {
    music_message.info(
        "output kappa_B(T, mu_B) along constant s/n_B trajectories...");

    double sovernB[] = {10.0, 20.0, 30.0, 51.0, 70.0, 94.0, 144.0, 420.0};
    int array_length = sizeof(sovernB) / sizeof(double);
    double s_0 = 0.00;     // 1/fm^3
    double s_max = 100.0;  // 1/fm^3
    double ds = 0.005;     // 1/fm^3
    int ns = static_cast<int>((s_max - s_0) / ds) + 1;
    for (int i = 0; i < array_length; i++) {
        std::ostringstream file_name;
        file_name << "kappa_B_sovernB_" << sovernB[i] << ".dat";
        std::ofstream of(file_name.str().c_str());
        // write out the header of the file
        of << "# e (GeV/fm^3)  rhob (1/fm^3) s (1/fm^3)  "
           << "T (GeV)  mu_B (GeV)  kappa (1/fm^2)" << std::endl;
        for (int j = 0; j < ns; j++) {
            double s_local = s_0 + j * ds;
            double nB_local = s_local / sovernB[i];
            double e_local = eos.get_s2e(s_local, nB_local);
            double s_check = eos.get_entropy(e_local, nB_local);
            double p_local = eos.get_pressure(e_local, nB_local);
            double temperature = eos.get_temperature(e_local, nB_local);
            double mu_B = eos.get_muB(e_local, nB_local);
            if (mu_B * hbarc > 0.78)
                continue;  // discard points out of the table
            double alpha_local =
                mu_B / std::max(temperature, Util::tiny_normal);

            double denorm_safe = std::copysign(
                std::max(
                    small_eps, std::abs(3. * temperature * tanh(alpha_local))),
                3. * temperature * tanh(alpha_local));
            double kappa_local =
                (DATA.kappa_coefficient
                 * (nB_local / denorm_safe
                    - nB_local * nB_local
                          / std::max(e_local + p_local, small_eps)));
            // output
            of << std::scientific << std::setw(18) << std::setprecision(8)
               << e_local * hbarc << "   " << nB_local << "   " << s_check
               << "   " << temperature * hbarc << "   " << mu_B * hbarc << "   "
               << kappa_local << std::endl;
        }
        of.close();  // close the file
    }
}

//! this function outputs the T and muB dependence of the specific shear
//! viscosity eta/s
void Diss::output_eta_over_s_T_and_muB_dependence() {
    music_message.info("output eta/s(T, mu_B) ...");
    std::ofstream of("eta_over_s_T_and_muB_dependence.dat");
    // write out the header of the file
    of << "# e (GeV/fm^3)  rhob (1/fm^3) T (GeV)  mu_B (GeV)  eta/s"
       << std::endl;

    // define the grid
    double e_min = 1e-5;   // fm^-4
    double e_max = 100.0;  // fm^-4
    int ne = 1000;
    double de = (e_max - e_min) / (ne - 1.);
    double rhob_min = 0.0;         // fm^-3
    double rhob_max = sqrt(10.0);  // fm^-3
    int nrhob = 1000;
    double drhob = (rhob_max - rhob_min) / (nrhob - 1.);

    for (int i = 0; i < ne; i++) {
        double e_local = e_min + i * de;
        for (int j = 0; j < nrhob; j++) {
            double rhob_local = rhob_min + j * drhob;
            rhob_local *= rhob_local;
            double mu_B_local = eos.get_muB(e_local, rhob_local);
            if (mu_B_local * hbarc > 0.89)
                continue;  // discard points out of the table
            double p_local = eos.get_pressure(e_local, rhob_local);
            double s_local = eos.get_entropy(e_local, rhob_local);
            double T_local = eos.get_temperature(e_local, rhob_local);

            double shear_to_s = DATA.shear_to_s;
            shear_to_s = transport_coeffs_.get_eta_over_s(T_local, mu_B_local);

            double eta_over_s = shear_to_s;
            if (DATA.muB_dependent_shear_to_s != 0) {
                eta_over_s =
                    shear_to_s * (e_local + p_local) / std::max(T_local * s_local, Util::tiny_normal);
            }

            // output
            of << std::scientific << std::setw(18) << std::setprecision(8)
               << e_local * hbarc << "   " << rhob_local << "   "
               << T_local * hbarc << "   " << mu_B_local * hbarc << "   "
               << eta_over_s << std::endl;
        }
    }
    of.close();  // close the file
}

//! this function outputs the T and muB dependence of the specific shear
//! viscosity eta/s along constant s/n_B trajectories
void Diss::output_eta_over_s_along_const_sovernB() {
    music_message.info(
        "output eta/s(T, mu_B) along constant s/n_B trajectories...");

    double sovernB[] = {10.0, 20.0, 30.0, 51.0, 70.0, 94.0, 144.0, 420.0};
    int array_length = sizeof(sovernB) / sizeof(double);
    double s_0 = 0.00;     // 1/fm^3
    double s_max = 100.0;  // 1/fm^3
    double ds = 0.005;     // 1/fm^3
    int ns = static_cast<int>((s_max - s_0) / ds) + 1;
    for (int i = 0; i < array_length; i++) {
        std::ostringstream file_name;
        file_name << "eta_over_s_sovernB_" << sovernB[i] << ".dat";
        std::ofstream of(file_name.str().c_str());
        // write out the header of the file
        of << "# e (GeV/fm^3)  rhob (1/fm^3) s (1/fm^3)  "
           << "T (GeV)  mu_B (GeV)  eta/s" << std::endl;
        for (int j = 0; j < ns; j++) {
            double s_local = s_0 + j * ds;
            double nB_local = s_local / sovernB[i];
            double e_local = eos.get_s2e(s_local, nB_local);
            double mu_B = eos.get_muB(e_local, nB_local);
            if (mu_B * hbarc > 0.89)
                continue;  // discard points out of the table

            double p_local = eos.get_pressure(e_local, nB_local);
            double s_check = eos.get_entropy(e_local, nB_local);
            double T_local = eos.get_temperature(e_local, nB_local);

            double shear_to_s = DATA.shear_to_s;
            shear_to_s = transport_coeffs_.get_eta_over_s(T_local, mu_B);

            double eta_over_s = shear_to_s;
            if (DATA.muB_dependent_shear_to_s != 0) {
                eta_over_s =
                    shear_to_s * (e_local + p_local) / std::max(T_local * s_local, Util::tiny_normal);
            }

            // output
            of << std::scientific << std::setw(18) << std::setprecision(8)
               << e_local * hbarc << "   " << nB_local << "   " << s_check
               << "   " << T_local * hbarc << "   " << mu_B * hbarc << "   "
               << eta_over_s << std::endl;
        }
        of.close();  // close the file
    }
}

//! this function outputs the T and muB dependence of the specific bulk
//! viscosity zeta/s
void Diss::output_zeta_over_s_T_and_muB_dependence() {
    music_message.info("output zeta/s(T, mu_B) ...");
    std::ofstream of("zeta_over_s_T_and_muB_dependence.dat");
    // write out the header of the file
    of << "# e (GeV/fm^3)  rhob (1/fm^3) T (GeV)  mu_B (GeV)  zeta/s"
       << std::endl;

    // define the grid
    double e_min = 1e-5;   // fm^-4
    double e_max = 100.0;  // fm^-4
    int ne = 1000;
    double de = (e_max - e_min) / (ne - 1.);
    double rhob_min = 0.0;         // fm^-3
    double rhob_max = sqrt(10.0);  // fm^-3
    int nrhob = 1000;
    double drhob = (rhob_max - rhob_min) / (nrhob - 1.);

    for (int i = 0; i < ne; i++) {
        double e_local = e_min + i * de;
        for (int j = 0; j < nrhob; j++) {
            double rhob_local = rhob_min + j * drhob;
            rhob_local *= rhob_local;
            double mu_B_local = eos.get_muB(e_local, rhob_local);
            if (mu_B_local * hbarc > 0.89)
                continue;  // discard points out of the table
            double p_local = eos.get_pressure(e_local, rhob_local);
            double s_local = eos.get_entropy(e_local, rhob_local);
            double T_local = eos.get_temperature(e_local, rhob_local);
            double mu_B = eos.get_muB(e_local, rhob_local);

            double bulk = transport_coeffs_.get_zeta_over_s(T_local, mu_B);
            double zeta_over_s =
                bulk * (e_local + p_local) / std::max(T_local * s_local, Util::tiny_normal);

            // output
            of << std::scientific << std::setw(18) << std::setprecision(8)
               << e_local * hbarc << "   " << rhob_local << "   "
               << T_local * hbarc << "   " << mu_B_local * hbarc << "   "
               << zeta_over_s << std::endl;
        }
    }
    of.close();  // close the file
}

//! this function outputs the T and muB dependence of the specific bulk
//! viscosity zeta/s along constant s/n_B trajectories
void Diss::output_zeta_over_s_along_const_sovernB() {
    music_message.info(
        "output zeta/s(T, mu_B) along constant s/n_B trajectories...");

    double sovernB[] = {10.0, 20.0, 30.0, 51.0, 70.0, 94.0, 144.0, 420.0};
    int array_length = sizeof(sovernB) / sizeof(double);
    double s_0 = 0.00;     // 1/fm^3
    double s_max = 100.0;  // 1/fm^3
    double ds = 0.005;     // 1/fm^3
    int ns = static_cast<int>((s_max - s_0) / ds) + 1;
    for (int i = 0; i < array_length; i++) {
        std::ostringstream file_name;
        file_name << "zeta_over_s_sovernB_" << sovernB[i] << ".dat";
        std::ofstream of(file_name.str().c_str());
        // write out the header of the file
        of << "# e (GeV/fm^3)  rhob (1/fm^3) s (1/fm^3)  "
           << "T (GeV)  mu_B (GeV)  eta/s" << std::endl;
        for (int j = 0; j < ns; j++) {
            double s_local = s_0 + j * ds;
            double nB_local = s_local / sovernB[i];
            double e_local = eos.get_s2e(s_local, nB_local);
            double mu_B = eos.get_muB(e_local, nB_local);
            if (mu_B * hbarc > 0.89)
                continue;  // discard points out of the table

            double T_local = eos.get_temperature(e_local, nB_local);
            double p_local = eos.get_pressure(e_local, nB_local);
            double s_check = eos.get_entropy(e_local, nB_local);

            double bulk = transport_coeffs_.get_zeta_over_s(T_local, mu_B);
            double zeta_over_s =
                bulk * (e_local + p_local) / std::max(T_local * s_local, Util::tiny_normal);

            // output
            of << std::scientific << std::setw(18) << std::setprecision(8)
               << e_local * hbarc << "   " << nB_local << "   " << s_check
               << "   " << T_local * hbarc << "   " << mu_B * hbarc << "   "
               << zeta_over_s << std::endl;
        }
        of.close();  // close the file
    }
}
