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
#include "bulk_pi_chem.h"

using Util::hbarc;
using Util::map_1d_idx_to_2d;
using Util::small_eps;

extern BulkPiChem* g_bulkPiChem; // defined in music.cpp

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
        if (alpha < 4 && DATA.turn_on_bulk == 1) {
            double gfac = (alpha == 0 ? -1.0 : 0.0);
            Pi_alpha0 =
                (arenaCurr.piBulk_[fieldIdx]
                 * (gfac
                    + arenaCurr.u_[alpha][fieldIdx]
                          * arenaCurr.u_[0][fieldIdx]));
            dPidtau =
                ((Pi_alpha0
                  - arenaPrev.piBulk_[fieldIdx]
                        * (gfac
                           + arenaPrev.u_[alpha][fieldIdx]
                                 * arenaPrev.u_[0][fieldIdx]))
                 / DATA.delta_tau);
        }
        double sf = tau * dWdtau + arenaCurr.Wmunu_[idx_1d_alpha0][fieldIdx];
        double bf = tau * dPidtau + Pi_alpha0;
        dwmn[alpha] += sf + bf;
        if (std::isnan(sf + bf)) {
            music_message << "[Error]Diss::MakeWSource: ";
            music_message << "sf=" << sf << " bf=" << bf << " Wmunu ="
                          << arenaCurr.Wmunu_[idx_1d_alpha0][fieldIdx]
                          << " pi_b =" << arenaCurr.piBulk_[fieldIdx]
                          << " prev_pi_b=" << arenaPrev.piBulk_[fieldIdx];
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
                if (DATA.turn_on_bulk == 1 && alpha < 4) {
                    double gfac1 = (alpha == (direction) ? 1.0 : 0.0);
                    double bgp1 =
                        (arenaCurr.piBulk_[Ip1]
                         * (gfac1
                            + arenaCurr.u_[alpha][Ip1]
                                  * arenaCurr.u_[direction][Ip1])
                         * tau_fac[direction]);
                    double bg =
                        (arenaCurr.piBulk_[Ic]
                         * (gfac1
                            + arenaCurr.u_[alpha][Ic]
                                  * arenaCurr.u_[direction][Ic])
                         * tau_fac[direction]);
                    double bgm1 =
                        (arenaCurr.piBulk_[Im1]
                         * (gfac1
                            + arenaCurr.u_[alpha][Im1]
                                  * arenaCurr.u_[direction][Im1])
                         * tau_fac[direction]);
                    // dPidx += minmod.minmod_dx(bgp1, bg,
                    // bgm1)/delta[direction];
                    //  use central difference to preserve conservation law
                    //  exactly
                    Pi_m = (bg + bgm1) * 0.5;
                    Pi_p = (bg + bgp1) * 0.5;
                }
                if (direction == 3 && (alpha == 0 || alpha == 3)) {
                    W_eta_p[alpha] = W_p + Pi_m;
                    W_eta_m[alpha] = W_m + Pi_p;
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

    R_shear = sqrt(pisize) / enthalpy;
    R_bulk = std::abs(grid_pt.pi_b) / enthalpy;
}

void Diss::Make_uWSource(
    const double tau, const Cell_small &grid_pt, const double theta_local,
    const DumuVec &a_local, const VelocityShearVec &sigma_1d,
    const VorticityVec &omega_1d, const std::vector<double> &thermalVec,
    std::array<double, 5> &sourceTerms) {
    auto sigma = Util::UnpackVecToMatrix(sigma_1d);
    auto Wmunu = Util::UnpackVecToMatrix(grid_pt.Wmunu);

    double epsilon = thermalVec[0];
    double T = thermalVec[6];
    double muB = thermalVec[7];

    double shear_to_s = transport_coeffs_.get_eta_over_s(T, muB);

    bool include_WWterm = false;
    bool include_Wsigma_term = false;
    if (DATA.include_second_order_terms == 1) {
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
    double transport_coefficient =
        (transport_coeffs_.get_phi7_coeff() * tau_pi / shear * (4. / 5.));
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

        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        //              Add coupling to bulk viscous pressure                //
        //             transport_coefficient_b*Bulk*sigma^mu nu              //
        //              transport_coefficient2_b*Bulk*W^mu nu                //
        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        double Coupling_to_Bulk = 0.0;
        if (DATA.include_second_order_terms == 1) {
            double Bulk_Sigma = grid_pt.pi_b * sigma[mu][nu];
            double Bulk_W = grid_pt.pi_b * Wmunu[mu][nu];

            // multiply term by its respective transport coefficient
            double Bulk_Sigma_term = Bulk_Sigma * transport_coefficient_b;
            double Bulk_W_term = Bulk_W * transport_coefficient2_b;

            // full term is
            // first term: sign changes according to metric sign convention
            Coupling_to_Bulk = -Bulk_Sigma_term + Bulk_W_term;
        }

        // final answer is
        sourceTerms[idx_1d - 4] =
            ((NS_term + tempf + Vorticity_term + Wsigma_term + WW_term
              + Coupling_to_Bulk)
             / tau_pi);
    }
}

void Diss::Make_uWRHS(
    const double tau, Fields &arena, const int fieldIdx, const int ix,
    const int iy, const int ieta, std::array<double, 9> &w_rhs,
    const double theta_local, const DumuVec &a_local) {
    auto grid_pt = arena.getCell(fieldIdx);
    auto Wmunu_local = Util::UnpackVecToMatrix(grid_pt.Wmunu);

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
    /* This is the second step in the operator splitting. it uses
       rk_flag+1 as initial condition */
    double delta[4] = {0.0, DATA.delta_x, DATA.delta_y, DATA.delta_eta * tau};
    const double delta_tau = DATA.delta_tau;

    int piIdxArr[9] = {4, 5, 6, 7, 8, 9, 11, 12, 13};
    int loopIdx = 5;
    if (DATA.turn_on_bulk) {
        loopIdx = 6;
    }
    if (DATA.turn_on_diff) {
        loopIdx = 9;
    }

    // pi^\mu\nu is symmetric
    FieldNeighbourLoop2(
        arena, ix, iy, ieta, FNLLAMBDAS2 {
            double g, gp1, gp2, gm1, gm2;
            double f, fp1, fp2, fm1, fm2;
            for (int idx = 0; idx < loopIdx; idx++) {
                int idx_1d = piIdxArr[idx];
                /* Get_uWmns */
                if (idx_1d == 9) {
                    g = arena.piBulk_[Ic];
                    gp2 = arena.piBulk_[Ip2];
                    gp1 = arena.piBulk_[Ip1];
                    gm1 = arena.piBulk_[Im1];
                    gm2 = arena.piBulk_[Im2];
                } else {
                    g = arena.Wmunu_[idx_1d][Ic];
                    gp2 = arena.Wmunu_[idx_1d][Ip2];
                    gp1 = arena.Wmunu_[idx_1d][Ip1];
                    gm1 = arena.Wmunu_[idx_1d][Im1];
                    gm2 = arena.Wmunu_[idx_1d][Im2];
                }

                f = g * arena.u_[direction][Ic];
                g *= arena.u_[0][Ic];

                fp2 = gp2 * arena.u_[direction][Ip2];
                gp2 *= arena.u_[0][Ip2];

                fp1 = gp1 * arena.u_[direction][Ip1];
                gp1 *= arena.u_[0][Ip1];

                fm1 = gm1 * arena.u_[direction][Im1];
                gm1 *= arena.u_[0][Im1];

                fm2 = gm2 * arena.u_[direction][Im2];
                gm2 *= arena.u_[0][Im2];

                /* MakeuWmnHalfs */
                /* uWmn */
                double uWphR = fp1 - 0.5 * minmod.minmod_dx(fp2, fp1, f);
                double temp = 0.5 * minmod.minmod_dx(fp1, f, fm1);
                double uWphL = f + temp;
                double uWmhR = f - temp;
                double uWmhL = fm1 + 0.5 * minmod.minmod_dx(f, fm1, fm2);

                /* just Wmn */
                double WphR = gp1 - 0.5 * minmod.minmod_dx(gp2, gp1, g);
                temp = 0.5 * minmod.minmod_dx(gp1, g, gm1);
                double WphL = g + temp;
                double WmhR = g - temp;
                double WmhL = gm1 + 0.5 * minmod.minmod_dx(g, gm1, gm2);

                double a = std::abs(arena.u_[direction][Ic]) / arena.u_[0][Ic];
                double am1 =
                    std::abs(arena.u_[direction][Im1]) / arena.u_[0][Im1];
                double ap1 =
                    std::abs(arena.u_[direction][Ip1]) / arena.u_[0][Ip1];

                double ax = std::max(a, ap1);
                double HWph = ((uWphR + uWphL) - ax * (WphR - WphL)) * 0.5;

                ax = std::max(a, am1);
                double HWmh = ((uWmhR + uWmhL) - ax * (WmhR - WmhL)) * 0.5;

                double HW = (HWph - HWmh) / delta[direction];
                /* make partial_i (u^i Wmn) */
                w_rhs[idx] += (-HW) * delta_tau;
            }
        });

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
    if (DATA.turn_on_bulk == 1) {
        w_rhs[5] -= (grid_pt.pi_b) * (grid_pt.u[0]) / tau * DATA.delta_tau;
        w_rhs[5] += (grid_pt.pi_b) * theta_local * DATA.delta_tau;
    }
}

double Diss::Make_uPiSource(
    const double tau, const Cell_small &grid_pt, const double theta_local,
    const VelocityShearVec &sigma_1d, const std::vector<double> &thermalVec) {
    // ---- chem-bulk: prefer global instance, fallback to local default ----
    // We do NOT delete or alter your existing structure; we only add a pointer
    // selection that prefers g_bulkPiChem (constructed in music.cpp).
    BulkPiChem* rt = g_bulkPiChem;
    static BulkPiChem* s_bulk_rt = nullptr;
    if (!rt) {
        // If the global was not created, construct a safe local default
        // (chem_bulk_on=0 means runtime calls reduce to baseline behavior).
        if (!s_bulk_rt) {
            ChemBulkConfig cfg;
            cfg.chem_bulk_on = 0;    // default off unless enabled via input
            cfg.C_tauPi      = 5.0;
            cfg.tauPi_min    = 1e-3;
            cfg.Yq_eps       = 1e-16;
            s_bulk_rt = new BulkPiChem(cfg, &eos);
        }
        rt = s_bulk_rt;
    }
    // ----------------------------------------------------------------------

    double tempf;
    double bulk;                      // this is zeta (not zeta/s)
    double Bulk_Relax_time = 0.2;     // fm
    double transport_coeff1, transport_coeff2;
    double transport_coeff1_s, transport_coeff2_s;
    double NS_term, BB_term;
    double Final_Answer;

    bool include_BBterm = false;
    bool include_coupling_to_shear = false;
    if (DATA.include_second_order_terms == 1) {
        include_BBterm = true;
        include_coupling_to_shear = true;
    }

    double epsilon    = thermalVec[0];
    double rhob       = thermalVec[1];
    double temperature= thermalVec[6];
    double mu_B       = thermalVec[7];

    double cs2      = thermalVec[5];
    double pressure = thermalVec[2];

    // ---------- baseline zeta (from table) ----------
    bulk = transport_coeffs_.get_zeta_over_s(temperature, mu_B);
    bulk = bulk * (epsilon + pressure) / std::max(temperature, small_eps); // zeta_eq(e)

    // ---------- baseline tau_Pi (kinetic variants) ----------
    if (DATA.bulk_relaxation_type == 0) {
        double csfactor = std::max(1./3. - cs2, small_eps);
        Bulk_Relax_time =
            (transport_coeffs_.get_bulk_relax_time_factor() * bulk
             / (csfactor * csfactor) / std::max(epsilon + pressure, small_eps));
    } else if (DATA.bulk_relaxation_type == 1) {
        double csfactor = std::max(1./3. - cs2, small_eps);
        Bulk_Relax_time =
            (transport_coeffs_.get_bulk_relax_time_factor() * bulk / csfactor
             / std::max(epsilon + pressure, small_eps));
    } else if (DATA.bulk_relaxation_type == 2) {
        Bulk_Relax_time =
            (transport_coeffs_.get_bulk_relax_time_factor() * bulk
             / std::max(epsilon + pressure, small_eps));
    }
    Bulk_Relax_time =
        std::min(10., std::max(3.*DATA.delta_tau, Bulk_Relax_time));

    // ---------- resummation on baseline ----------
    double transport_coeff1_resc, transport_coeff2_resc;
    transport_coeff1_resc = transport_coeffs_.get_delta_bulkPibulkPi_coeff();
    transport_coeff2_resc = transport_coeffs_.get_tau_bulkPibulkPi_coeff();
    transport_coeff1 = transport_coeff1_resc * Bulk_Relax_time;
    transport_coeff2 = transport_coeff2_resc * Bulk_Relax_time;

    transport_coeff1_s =
        (transport_coeffs_.get_lambda_bulkPipi_coeff()
         * (DATA.bulk_relaxation_type == 2 ? 1./3. : (1./3. - cs2))
         * Bulk_Relax_time);
    transport_coeff2_s = 0.;

    if (DATA.FlagResumTransportCoeff) {
        double R_shear = 0., R_bulk = 0.;
        computeInverseReynoldsNumbers(epsilon + pressure, grid_pt, R_shear, R_bulk);
        double res = transport_coeffs_.getResummedCorrFactor(R_shear, R_bulk);
        bulk *= res;
        transport_coeff1 *= res;
        transport_coeff2 *= res;
        transport_coeff1_s *= res;
        transport_coeff2_s *= res;
    }

    // ---------- RUNTIME CHEMISTRY OVERRIDE (keeps same IS form) ----------
    // Always evaluate; if chem is off in config, rt->tau/zeta reduce to baselines.
    {
        const double Pi_tot = grid_pt.pi_b;
        const double tau_run  = rt->tau_Pi_runtime(epsilon, rhob, Pi_tot);
        const double zeta_run = rt->zeta_runtime(epsilon, rhob, Pi_tot, bulk);
        Bulk_Relax_time = std::min(10., std::max(3.*DATA.delta_tau, tau_run));
        bulk = zeta_run;  // use zeta(e,Pi)
    }
    // ---------------------------------------------------------------------

    // ---------- NS + nonlinear structure (unchanged) ----------
    NS_term = -bulk * theta_local;

    tempf = (-(grid_pt.pi_b) - transport_coeff1 * theta_local * (grid_pt.pi_b));

    if (include_BBterm) {
        BB_term = transport_coeff2 * (grid_pt.pi_b) * (grid_pt.pi_b);
    } else {
        BB_term = 0.0;
    }

    double Coupling_to_Shear = 0.0;
    if (include_coupling_to_shear) {
        auto sigma = Util::UnpackVecToMatrix(sigma_1d);
        auto Wmunu = Util::UnpackVecToMatrix(grid_pt.Wmunu);

        double Wsigma =
            (Wmunu[0][0]*sigma[0][0] + Wmunu[1][1]*sigma[1][1]
             + Wmunu[2][2]*sigma[2][2] + Wmunu[3][3]*sigma[3][3]
             - 2.*(Wmunu[0][1]*sigma[0][1] + Wmunu[0][2]*sigma[0][2] + Wmunu[0][3]*sigma[0][3])
             + 2.*(Wmunu[1][2]*sigma[1][2] + Wmunu[1][3]*sigma[1][3] + Wmunu[2][3]*sigma[2][3]));

        double WW =
            (Wmunu[0][0]*Wmunu[0][0] + Wmunu[1][1]*Wmunu[1][1]
             + Wmunu[2][2]*Wmunu[2][2] + Wmunu[3][3]*Wmunu[3][3]
             - 2.*(Wmunu[0][1]*Wmunu[0][1] + Wmunu[0][2]*Wmunu[0][2] + Wmunu[0][3]*Wmunu[0][3])
             + 2.*(Wmunu[1][2]*Wmunu[1][2] + Wmunu[1][3]*Wmunu[1][3] + Wmunu[2][3]*Wmunu[2][3]));

        double Shear_Sigma_term = Wsigma * transport_coeff1_s;
        double Shear_Shear_term = WW     * transport_coeff2_s;
        Coupling_to_Shear = -Shear_Sigma_term + Shear_Shear_term;
    }

    Final_Answer = NS_term + tempf + BB_term + Coupling_to_Shear;
    return Final_Answer / Bulk_Relax_time;
}


/* Make_uPiSource */

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
void Diss::Make_uqSource(
    const double tau, const Cell_small &grid_pt, const double theta_local,
    const DumuVec &a_local, const VelocityShearVec &sigma_1d,
    const VorticityVec &omega_1d, const DmuMuBoverTVec &baryon_diffusion_vec,
    const std::vector<double> &thermalVec, std::array<double, 3> &sourceTerms) {
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

    for (int nu = 1; nu < 4; nu++) {
        double NS =
            kappa * (baryon_diffusion_vec[nu] + grid_pt.u[nu] * a_local[4]);

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
            double alpha_local = mu_B_local / T_local;

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
            double alpha_local = mu_B / temperature;

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
                    shear_to_s * (e_local + p_local) / (T_local * s_local);
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
                    shear_to_s * (e_local + p_local) / (T_local * s_local);
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
                bulk * (e_local + p_local) / (T_local * s_local);

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
                bulk * (e_local + p_local) / (T_local * s_local);

            // output
            of << std::scientific << std::setw(18) << std::setprecision(8)
               << e_local * hbarc << "   " << nB_local << "   " << s_check
               << "   " << T_local * hbarc << "   " << mu_B * hbarc << "   "
               << zeta_over_s << std::endl;
        }
        of.close();  // close the file
    }
}
