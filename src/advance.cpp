// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <cassert>
#include <cmath>
#include <limits>
#include <memory>

#include "advance.h"
#include "cell.h"
#include "data.h"
#include "eos.h"
#include "evolve.h"
#include "reconst.h"
#include "util.h"

using Eigen::EigenSolver;
using Eigen::Matrix4d;
using Util::hbarc;
using Util::map_1d_idx_to_2d;
using Util::map_2d_idx_to_1d;

Advance::Advance(
    const EOS &eosIn, const InitData &DATA_in,
    std::shared_ptr<HydroSourceBase> hydro_source_ptr_in)
    : DATA(DATA_in),
      eos(eosIn),
      diss_helper(eosIn, DATA_in),
      minmod(DATA_in),
      reconst_helper(eos, DATA_in.echo_level, DATA_in.beastMode),
      transport_coeffs_(DATA_in) {
    hydro_source_terms_ptr = hydro_source_ptr_in;
    flag_add_hydro_source = false;
    if (hydro_source_terms_ptr) {
        if (DATA.Initial_profile == 42) {
            if (hydro_source_terms_ptr->get_number_of_sources() > 0) {
                flag_add_hydro_source = true;
            }
        } else {
            flag_add_hydro_source = true;
        }
    }
}

//! this function evolves one Runge-Kutta step in tau
void Advance::AdvanceIt(
    const double tau, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
    Fields &arenaFieldsNext, const int rk_flag) {
    const int grid_neta = arenaFieldsCurr.nEta();
    const int grid_nx = arenaFieldsCurr.nX();
    const int grid_ny = arenaFieldsCurr.nY();

    arenaFieldsCurr.updateGhosts();
    computePressure(arenaFieldsCurr, pressureCurr_);
    computeFaceFluxes(
        tau + rk_flag * DATA.delta_tau, arenaFieldsCurr, pressureCurr_,
        faceFluxes_);
    if (DATA.viscosity_flag == 1) {
        diss_helper.computeWFluxDivergence(
            tau + rk_flag * DATA.delta_tau, arenaFieldsCurr, wFluxDiv_);
    }

#pragma omp parallel for collapse(3) schedule(guided)
    for (int ieta = 0; ieta < grid_neta; ieta++)
        for (int ix = 0; ix < grid_nx; ix++)
            for (int iy = 0; iy < grid_ny; iy++) {
                int fieldIdx = arenaFieldsCurr.getFieldIdx(ix, iy, ieta);

                double eta_s_local =
                    -DATA.eta_size / 2. + ieta * DATA.delta_eta;
                double x_local = -DATA.x_size / 2. + ix * DATA.delta_x;
                double y_local = -DATA.y_size / 2. + iy * DATA.delta_y;

                TJbVec qiFinal = {0.};
                FirstRKStepT(
                    tau, x_local, y_local, eta_s_local, ix, iy, ieta, rk_flag,
                    fieldIdx, arenaFieldsCurr, arenaFieldsNext, arenaFieldsPrev,
                    pressureCurr_, qiFinal);

                if (DATA.viscosity_flag == 1) {
                    // arenaFieldsCurr is at tau_rk: tau_n in the first RK
                    // stage, tau_{n+1} (the predictor) in the second. The
                    // derivatives need its tau (u^tau / tau in theta, the
                    // Christoffel terms, tau delta_eta)
                    const double tau_rk = tau + rk_flag * DATA.delta_tau;
                    U_derivative u_derivative_helper(DATA, eos);
                    u_derivative_helper.MakedU(
                        tau_rk, arenaFieldsPrev, arenaFieldsCurr, fieldIdx, ix,
                        iy, ieta);
                    double theta_local =
                        u_derivative_helper.calculate_expansion_rate(
                            tau_rk, arenaFieldsCurr, fieldIdx);
                    DumuVec a_local;
                    u_derivative_helper.calculate_Du_supmu(
                        tau_rk, arenaFieldsCurr, fieldIdx, a_local);

                    VelocityShearVec sigma_local;
                    u_derivative_helper.calculate_velocity_shear_tensor(
                        tau_rk, arenaFieldsCurr, fieldIdx, theta_local, a_local,
                        sigma_local);

                    VorticityVec omega_local;
                    u_derivative_helper
                        .calculate_kinetic_vorticity_with_spatial_projector(
                            tau_rk, arenaFieldsCurr, fieldIdx, a_local,
                            omega_local);

                    DmuMuBoverTVec baryon_diffusion_vector;
                    u_derivative_helper.get_DmuMuBoverTVec(
                        baryon_diffusion_vector);

                    // the spatial part of the W update depends only on
                    // arenaFieldsCurr, theta and Du: compute it once, also
                    // for the iterations of the stage-2 fixed point
                    std::array<double, 10> w_rhs = {0.};
                    diss_helper.Make_uWRHS(
                        tau_rk, arenaFieldsCurr, fieldIdx, wFluxDiv_, w_rhs,
                        theta_local, a_local);

                    FirstRKStepW(
                        tau, arenaFieldsPrev, arenaFieldsCurr, arenaFieldsNext,
                        rk_flag, theta_local, a_local, sigma_local, omega_local,
                        baryon_diffusion_vector, w_rhs, ieta, ix, iy, fieldIdx);

                    if (rk_flag == 1) {
                        stage2ViscousFixedPoint(
                            tau, arenaFieldsPrev, arenaFieldsCurr,
                            arenaFieldsNext, qiFinal, theta_local, a_local,
                            sigma_local, omega_local, baryon_diffusion_vector,
                            w_rhs, ieta, ix, iy, fieldIdx);
                    }
                }
            }
    arenaFieldsNext.updateGhosts();
}

// The viscous part V of T^{tau alpha} (and q^tau of J^tau) enters through
// the backward difference (tau_rk + dtau) V_rk - tau_rk V_prev in
// Diss::MakeWSource. The first RK stage subtracts tau_{n+1} V_n - tau_n
// V_{n-1}, the second tau_{n+2} V* - tau_{n+1} V_n with the predictor V*,
// so the RK average telescopes only if V* = V_{n+1}; the remainder
// 1/2 tau_{n+2} (V_{n+1} - V*) drifts E, P and N_B wherever V changes
// quickly (Pi right after string deposition: 0.5% of E). The term involves
// only this cell, so V* is replaced by V_{n+1} with a local fixed point:
// shift the stage-2 result by s, reconstruct, and update W again (the
// fluxes and gradients do not depend on it), until
// s = 1/2 tau_{n+2} (V_{n+1}(s) - V*). The average then telescopes to the
// 1/2 (tau_{n+2} V_{n+1} + tau_{n+1} V_n) that
// Cell_info::check_conservation_law adds.
// The updates start as full steps, which converge in one or two steps in
// almost all cells. Where the miss grows (stiff baryon diffusion, cells
// with Pi ~ -(e + P)) the step is halved, and a cell that does not
// converge keeps the visited shift with the smallest miss (at worst s = 0,
// the predictor). A fixed number of full steps diverged in such cells, and
// a fixed under-relaxation did not converge right after string deposition;
// both made viscous string events blow up.
void Advance::stage2ViscousFixedPoint(
    const double tau, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
    Fields &arenaFieldsNext, const TJbVec &qiFinal, const double theta_local,
    const DumuVec &a_local, const VelocityShearVec &sigma_local,
    const VorticityVec &omega_local,
    const DmuMuBoverTVec &baryon_diffusion_vector,
    const std::array<double, 10> &w_rhs, const int ieta, const int ix,
    const int iy, const int fieldIdx) {
    const double tau_next = tau + DATA.delta_tau;
    const double tau_np2 = tau + 2. * DATA.delta_tau;
    const auto guess = arenaFieldsCurr.getCellIdeal(fieldIdx);
    const double scaleT = std::abs(qiFinal[0]);
    const double scaleJ = std::abs(qiFinal[4]) + 1e-12 * scaleT;

    auto applyShift = [&](const TJbVec &shift) {
        TJbVec qi = qiFinal;
        for (int alpha = 0; alpha < 5; alpha++) qi[alpha] -= shift[alpha];
        reconstructAndStore(tau_next, qi, guess, arenaFieldsNext, fieldIdx);
        FirstRKStepW(
            tau, arenaFieldsPrev, arenaFieldsCurr, arenaFieldsNext, 1,
            theta_local, a_local, sigma_local, omega_local,
            baryon_diffusion_vector, w_rhs, ieta, ix, iy, fieldIdx);
    };

    TJbVec shift = {0.};  // the state in arenaFieldsNext has this shift
    TJbVec best = shift;
    double bestMiss = std::numeric_limits<double>::max();
    double prevMiss = std::numeric_limits<double>::max();
    double lastMiss = 0.;
    double relaxation = 1.;
    for (int iter = 0; iter <= stage2MaxIterations_; iter++) {
        TJbVec target = {0.};
        double miss = 0.;
        bool converged = true;
        for (int alpha = 0; alpha < 5; alpha++) {
            target[alpha] =
                0.5 * tau_np2
                * (viscousTimeComponent(arenaFieldsNext, fieldIdx, alpha)
                   - viscousTimeComponent(arenaFieldsCurr, fieldIdx, alpha));
            const double scale = (alpha < 4 ? scaleT : scaleJ);
            const double d = std::abs(target[alpha] - shift[alpha]);
            miss += d / scale;
            if (d > stage2Tolerance_ * scale) converged = false;
        }
        lastMiss = miss;
        if (converged && iter > 0) return;
        if (miss < bestMiss) {
            bestMiss = miss;
            best = shift;
        }
        if (iter == stage2MaxIterations_) break;
        if (iter > 0 && miss > prevMiss) relaxation *= 0.5;
        prevMiss = miss;
        for (int alpha = 0; alpha < 5; alpha++) {
            shift[alpha] += relaxation * (target[alpha] - shift[alpha]);
        }
        applyShift(shift);
    }
    if (lastMiss > bestMiss) applyShift(best);
}

double Advance::viscousTimeComponent(
    const Fields &arena, const int fieldIdx, const int alpha) const {
    if (alpha == 4) return arena.Wmunu_[10][fieldIdx];  // q^tau
    double v = arena.Wmunu_[alpha][fieldIdx];           // W^{alpha tau}
    if (DATA.turn_on_bulk == 1 || DATA.turn_on_bulk_chem == 1) {
        const double g = (alpha == 0 ? -1. : 0.);
        const double Pi = arena.piBulk_[fieldIdx] + arena.piBulkChem_[fieldIdx];
        v += Pi * (g + arena.u_[alpha][fieldIdx] * arena.u_[0][fieldIdx]);
    }
    return v;
}

void Advance::reconstructAndStore(
    const double tau_next, const TJbVec &qi, const ReconstCell &guess,
    Fields &arenaFieldsNext, const int fieldIdx) {
    auto grid_rk_t = reconst_helper.ReconstIt_shell(tau_next, qi, guess);
    arenaFieldsNext.e_[fieldIdx] = grid_rk_t.e;
    arenaFieldsNext.rhob_[fieldIdx] = grid_rk_t.rhob;
    arenaFieldsNext.rhoq_[fieldIdx] = grid_rk_t.rhoq;
    arenaFieldsNext.rhos_[fieldIdx] = grid_rk_t.rhos;
    for (int ii = 0; ii < 4; ii++) {
        arenaFieldsNext.u_[ii][fieldIdx] = grid_rk_t.u[ii];
    }
}

void Advance::computePressure(
    const Fields &arena, std::vector<double> &pressure) const {
    // every cell, ghost cells included, so that the stencils can read the
    // pressure past the grid edge
    const int ncells = arena.nStorage();
    pressure.resize(ncells);
#pragma omp parallel for
    for (int idx = 0; idx < ncells; idx++) {
        pressure[idx] = eos.get_pressure(
            arena.e_[idx], arena.rhob_[idx], arena.rhoq_[idx],
            arena.rhos_[idx]);
    }
}

/* %%%%%%%%%%%%%%%%%%%%%% First steps begins here %%%%%%%%%%%%%%%%%% */
void Advance::FirstRKStepT(
    const double tau, const double x_local, const double y_local,
    const double eta_s_local, const int ix, const int iy, const int ieta,
    const int rk_flag, const int fieldIdx, Fields &arenaFieldsCurr,
    Fields &arenaFieldsNext, Fields &arenaFieldsPrev,
    const std::vector<double> &pressureCurr, TJbVec &qiFinal) {
    // this advances the ideal part
    double tau_rk = tau + rk_flag * (DATA.delta_tau);

    auto cellPrev = arenaFieldsPrev.getCellIdeal(fieldIdx);
    auto cellCurr = arenaFieldsCurr.getCellIdeal(fieldIdx);

    // Solve partial_a T^{a mu} = -partial_a W^{a mu}
    // Update T^{mu nu}
    // MakeDelatQI gets
    //   qi = q0 if rk_flag = 0 or
    //   qi = q0 + k1 if rk_flag = 1
    // rhs[alpha] is what MakeDeltaQI outputs.
    // It is the spatial derivative part of partial_a T^{a mu}
    // (including geometric terms)
    TJbVec qi = {0};

    const double pressure = pressureCurr[fieldIdx];

    for (int alpha = 0; alpha < 7; alpha++) {
        qi[alpha] = get_TJb(cellCurr, alpha, 0, pressure) * tau_rk;
    }

    MakeDeltaQI(faceFluxes_, ix, iy, ieta, qi);

    TJbVec qi_source = {0.0};

    if (flag_add_hydro_source) {
        EnergyFlowVec j_mu = {0};
        FlowVec u_local;
        for (int ii = 0; ii < 4; ii++)
            u_local[ii] = arenaFieldsCurr.u_[ii][fieldIdx];

        hydro_source_terms_ptr->get_hydro_energy_source(
            tau_rk, x_local, y_local, eta_s_local, u_local, j_mu);
        for (int ii = 0; ii < 4; ii++) {
            qi_source[ii] = tau_rk * j_mu[ii];
            if (isnan(qi_source[ii])) {
                music_message << "qi_source is nan. i = " << ii;
                music_message.flush("error");
                exit(1);
            }
        }

        if (DATA.turn_on_rhob == 1) {
            qi_source[4] =
                (tau_rk
                 * hydro_source_terms_ptr->get_hydro_rhob_source(
                     tau_rk, x_local, y_local, eta_s_local, u_local));
        }

        if (DATA.turn_on_QS == 1) {
            qi_source[5] =
                (tau_rk
                 * hydro_source_terms_ptr->get_hydro_rhoQ_source(
                     tau_rk, x_local, y_local, eta_s_local, u_local));
        }
    }

    // now MakeWSource returns partial_a W^{a mu}
    // (including geometric terms)

    TJbVec dwmn = {0.0};

    diss_helper.MakeWSource(
        tau_rk, ix, iy, ieta, dwmn, arenaFieldsCurr, arenaFieldsPrev, fieldIdx);

    double pressurePrev = eos.get_pressure(
        cellPrev.e, cellPrev.rhob, cellPrev.rhoq, cellPrev.rhos);
    for (int alpha = 0; alpha < 7; alpha++) {
        /* dwmn is the only one with the minus sign */
        qi[alpha] -= dwmn[alpha] * (DATA.delta_tau);

        // add energy momentum and net baryon density source terms
        qi[alpha] += qi_source[alpha] * DATA.delta_tau;

        // set baryon density back to zero if viscous correction made it
        // non-zero remove/modify if rho_b!=0
        // - this is only to remove the viscous correction that
        // can make rho_b negative which we do not want.
        // if (DATA.turn_on_rhob == 0) {
        //    if (alpha == 4 && std::abs(qi[alpha]) > 1e-12)
        //        qi[alpha] = 0.;
        //}

        /* if rk_flag > 0, we now have q0 + k1 + k2.
         * So add q0 and multiply by 1/2 */
        qi[alpha] += rk_flag * get_TJb(cellPrev, alpha, 0, pressurePrev) * tau;
        qi[alpha] *= 1. / (1. + rk_flag);
    }

    qiFinal = qi;
    reconstructAndStore(
        tau + DATA.delta_tau, qi, cellCurr, arenaFieldsNext, fieldIdx);
}

void Advance::FirstRKStepW(
    const double tau, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
    Fields &arenaFieldsNext, const int rk_flag, const double theta_local,
    const DumuVec &a_local, const VelocityShearVec &sigma_local,
    const VorticityVec &omega_local,
    const DmuMuBoverTVec &baryon_diffusion_vector,
    const std::array<double, 10> &w_rhs, const int ieta, const int ix,
    const int iy, const int fieldIdx) {
    auto grid_c = arenaFieldsCurr.getCell(fieldIdx);
    auto grid_f = arenaFieldsNext.getCell(fieldIdx);

    const double tau_now = tau + rk_flag * DATA.delta_tau;

    // Solve partial_a (u^a W^{mu nu}) = 0
    // Update W^{mu nu}
    // mu = 4 is the baryon current qmu

    // calculate delta uWmunu
    // need to use u[0][mu], remember rk_flag = 0 here
    // with the KT flux
    // solve partial_tau (u^0 W^{kl}) = -partial_i (u^i W^{kl}
    /* Advance uWmunu */

    // w_rhs: spatial gradients for all viscous quantities (Diss::Make_uWRHS
    // at tau_now, computed by the caller)

    // the sources are evaluated at the stage's own state (the predictor in
    // the second RK stage), like theta, sigma and Du; the transport
    // coefficients of tau_n there made the scheme first order in dtau
    ThermalVec thermalVec;
    eos.getThermalVariables(
        grid_c.epsilon, grid_c.rhob, grid_c.rhoq, grid_c.rhos, thermalVec);

    double tempf;
    double u0Prev = arenaFieldsPrev.u_[0][fieldIdx];
    if (DATA.turn_on_shear == 1) {
        std::array<double, 5> sourceTerms = {0.};
        diss_helper.Make_uWSource(
            tau_now, grid_c, theta_local, a_local, sigma_local, omega_local,
            thermalVec, sourceTerms);
        for (int idx_1d = 4; idx_1d < 9; idx_1d++) {
            double WmunuPrev = arenaFieldsPrev.Wmunu_[idx_1d][fieldIdx];
            tempf =
                ((1. - rk_flag) * (grid_c.Wmunu[idx_1d] * grid_c.u[0])
                 + rk_flag * (WmunuPrev * u0Prev));
            tempf += sourceTerms[idx_1d - 4] * (DATA.delta_tau);
            tempf += w_rhs[idx_1d - 4];
            tempf += rk_flag * ((grid_c.Wmunu[idx_1d]) * (grid_c.u[0]));
            tempf *= 1. / (1. + rk_flag);
            grid_f.Wmunu[idx_1d] = tempf / (grid_f.u[0]);
        }
    } else {
        for (int idx_1d = 4; idx_1d < 9; idx_1d++) {
            grid_f.Wmunu[idx_1d] = 0.0;
        }
    }

    if (DATA.turn_on_bulk == 1) {
        double piBulkPrev = arenaFieldsPrev.piBulk_[fieldIdx];
        tempf =
            ((1. - rk_flag) * (grid_c.pi_b * grid_c.u[0])
             + rk_flag * piBulkPrev * u0Prev);
        double temps = diss_helper.Make_uPiSource(
            tau_now, grid_c, theta_local, sigma_local, thermalVec);
        tempf += temps * (DATA.delta_tau);
        tempf += w_rhs[5];
        tempf += rk_flag * ((grid_c.pi_b) * (grid_c.u[0]));
        tempf *= 1. / (1. + rk_flag);
        grid_f.pi_b = tempf / (grid_f.u[0]);
    } else {
        grid_f.pi_b = 0.0;
        grid_f.pi_b_chem = 0.0;
    }

    // ---------------------------------------------------------------------
    // NEW: chemical-equilibration bulk pressure evolution (pi_b_chem)
    // ---------------------------------------------------------------------
    if (DATA.turn_on_bulk_chem == 1) {
        double piBulkChemPrev_local = arenaFieldsPrev.piBulkChem_[fieldIdx];

        tempf =
            ((1. - rk_flag) * (grid_c.pi_b_chem * grid_c.u[0])
             + rk_flag * piBulkChemPrev_local * u0Prev);

        double temps_chem = diss_helper.Make_uPiChemSource(
            tau_now, grid_c, theta_local, sigma_local, thermalVec);

        tempf += temps_chem * (DATA.delta_tau);
        tempf += w_rhs[6];
        tempf += rk_flag * ((grid_c.pi_b_chem) * (grid_c.u[0]));
        tempf *= 1. / (1. + rk_flag);

        grid_f.pi_b_chem = tempf / (grid_f.u[0]);
    }
    // ---------------------------------------------------------------------

    // CShen: add source term for baryon diffusion
    if (DATA.turn_on_diff == 1) {
        std::array<double, 3> sourceTerms = {0.};
        diss_helper.Make_uqSource(
            tau_now, grid_c, theta_local, a_local, sigma_local, omega_local,
            baryon_diffusion_vector, thermalVec, sourceTerms);
        for (int idx_1d = 11; idx_1d < 14; idx_1d++) {
            double WmunuPrev = arenaFieldsPrev.Wmunu_[idx_1d][fieldIdx];
            tempf =
                ((1. - rk_flag) * (grid_c.Wmunu[idx_1d] * grid_c.u[0])
                 + rk_flag * WmunuPrev * u0Prev);
            tempf += sourceTerms[idx_1d - 11] * (DATA.delta_tau);
            tempf += w_rhs[idx_1d - 4];
            tempf += rk_flag * (grid_c.Wmunu[idx_1d] * grid_c.u[0]);
            tempf *= 1. / (1. + rk_flag);

            grid_f.Wmunu[idx_1d] = tempf / (grid_f.u[0]);
        }
    } else {
        for (int idx_1d = 10; idx_1d < 14; idx_1d++) {
            grid_f.Wmunu[idx_1d] = 0.0;
        }
    }

    // re-make Wmunu[3][3] so that Wmunu[mu][nu] is traceless
    double u0sq = grid_f.u[0] * grid_f.u[0];
    grid_f.Wmunu[9] =
        ((2.
              * (grid_f.u[1] * grid_f.u[2] * grid_f.Wmunu[5]
                 + grid_f.u[1] * grid_f.u[3] * grid_f.Wmunu[6]
                 + grid_f.u[2] * grid_f.u[3] * grid_f.Wmunu[8])
          - (u0sq - grid_f.u[1] * grid_f.u[1]) * grid_f.Wmunu[4]
          - (u0sq - grid_f.u[2] * grid_f.u[2]) * grid_f.Wmunu[7])
         / (u0sq - grid_f.u[3] * grid_f.u[3]));

    // make Wmunu[i][0] using the transversality
    for (int mu = 1; mu < 4; mu++) {
        tempf = 0.0;
        for (int nu = 1; nu < 4; nu++) {
            int idx_1d = map_2d_idx_to_1d(mu, nu);
            tempf += grid_f.Wmunu[idx_1d] * grid_f.u[nu];
        }
        grid_f.Wmunu[mu] = tempf / (grid_f.u[0]);
    }

    // make Wmunu[0][0]
    tempf = 0.0;
    for (int nu = 1; nu < 4; nu++) tempf += grid_f.Wmunu[nu] * grid_f.u[nu];
    grid_f.Wmunu[0] = tempf / (grid_f.u[0]);

    // make qmu[0] using transversality
    tempf = 0.0;
    for (int idx_1d = 11; idx_1d < 14; idx_1d++) {
        tempf += grid_f.Wmunu[idx_1d] * grid_f.u[idx_1d - 10];
    }
    grid_f.Wmunu[10] = DATA.turn_on_diff * tempf / (grid_f.u[0]);

    // If the energy density of the fluid element is smaller than 0.01GeV
    // reduce Wmunu using the QuestRevert algorithm
    if (DATA.Initial_profile != 0 && DATA.Initial_profile != 1) {
        QuestRevert(grid_f, ieta, ix, iy);
        if (DATA.turn_on_diff == 1) {
            QuestRevert_qmu(grid_f, ieta, ix, iy);
        }
        if (DATA.FlagResumTransportCoeff) {
            QuestRevertResummedTransCoeff(grid_f);
        }
        if (DATA.FlagCausalityMethod == 1) {
            solveEigenvaluesWmunu(grid_f);
            nCausalityConstraints(grid_f);
        }
    }
    for (int idx_1d = 0; idx_1d < 14; idx_1d++) {
        arenaFieldsNext.Wmunu_[idx_1d][fieldIdx] = grid_f.Wmunu[idx_1d];
    }
    arenaFieldsNext.piBulk_[fieldIdx] = grid_f.pi_b;
    arenaFieldsNext.piBulkChem_[fieldIdx] = grid_f.pi_b_chem;
}

void Advance::QuestRevertResummedTransCoeff(Cell_small &grid_pt) {
    double e_local = grid_pt.epsilon;
    double rhob = grid_pt.rhob;
    double p_local =
        eos.get_pressure(e_local, rhob, grid_pt.rhoq, grid_pt.rhos);

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
    double pi_local = grid_pt.pi_b;

    double Rshear = sqrt(pisize) / (e_local + p_local + 1e-16);
    double Rbulk = std::abs(pi_local) / (e_local + p_local + 1e-16);
    double r_combined = 1.5 * (Rshear + Rbulk) + 1e-16;
    double renorm = std::min(1., (1. - 1e-6) / r_combined);
    for (int mu = 0; mu < 10; mu++) {
        grid_pt.Wmunu[mu] = renorm * grid_pt.Wmunu[mu];
    }
    grid_pt.pi_b = renorm * grid_pt.pi_b;
}

//! this function reduce the size of shear stress tensor and bulk pressure
//! in the dilute region to stablize numerical simulations
void Advance::QuestRevert(
    Cell_small &grid_pt, const int ieta, const int ix, const int iy) {
    double eps_scale = 0.1;  // 1/fm^4
    double e_local = grid_pt.epsilon;
    double rhob = grid_pt.rhob;
    double rhoq = grid_pt.rhoq;
    double rhos = grid_pt.rhos;

    // regulation factor in the default MUSIC
    // double factor = 300.*tanh(grid_pt.epsilon/eps_scale);
    double xi = 0.05;
    double factor = 10. * DATA.quest_revert_strength
                    * (1. / (exp(-(e_local - eps_scale) / xi) + 1.)
                       - 1. / (exp(eps_scale / xi) + 1.));
    double factor_bulk = factor;

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

    double pi_local = grid_pt.pi_b;
    double bulksize = 3. * pi_local * pi_local;

    double p_local = eos.get_pressure(e_local, rhob, rhoq, rhos);
    double eq_size = e_local * e_local + 3. * p_local * p_local;

    // In default MUSIC
    double rho_shear = sqrt(pisize / eq_size) / factor;
    double rho_bulk = sqrt(bulksize / eq_size) / factor_bulk;

    // Reducing the shear stress tensor
    double rho_shear_max = 0.1;
    if (std::isnan(rho_shear)) {
        for (int mu = 0; mu < 10; mu++) {
            grid_pt.Wmunu[mu] = 0.0;
        }
    } else if (rho_shear > rho_shear_max) {
        if (e_local > eps_scale && DATA.echo_level > 5) {
            music_message << "ieta = " << ieta << ", ix = " << ix
                          << ", iy = " << iy
                          << ", energy density = " << e_local * hbarc
                          << " GeV/fm^3, shear |pi/(epsilon+3*P)| = "
                          << rho_shear;
            music_message.flush("warning");
        }
        for (int mu = 0; mu < 10; mu++) {
            grid_pt.Wmunu[mu] = (rho_shear_max / rho_shear) * grid_pt.Wmunu[mu];
        }
    }

    // Reducing bulk viscous pressure
    double rho_bulk_max = 0.1;
    if (rho_bulk > rho_bulk_max) {
        if (e_local > eps_scale && DATA.echo_level > 5) {
            music_message << "ieta = " << ieta << ", ix = " << ix
                          << ", iy = " << iy
                          << ", energy density = " << e_local * hbarc
                          << " GeV/fm^3, bulk |Pi/(epsilon+3*P)| = "
                          << rho_bulk;
            music_message.flush("warning");
        }
        grid_pt.pi_b = (rho_bulk_max / rho_bulk) * grid_pt.pi_b;
    }

    // Reducing bulk_chem
    double rho_bulkChem =
        (sqrt(3. * grid_pt.pi_b_chem * grid_pt.pi_b_chem / eq_size)
         / factor_bulk);
    double rho_bulkChem_max = 0.1;
    if (rho_bulkChem > rho_bulkChem_max) {
        if (e_local > eps_scale && DATA.echo_level > 5) {
            music_message << "ieta = " << ieta << ", ix = " << ix
                          << ", iy = " << iy
                          << ", energy density = " << e_local * hbarc
                          << " GeV/fm^3, bulk_chem |Pi/(epsilon+3*P)| = "
                          << rho_bulkChem;
            music_message.flush("warning");
        }
        grid_pt.pi_b_chem =
            (rho_bulkChem_max / rho_bulkChem) * grid_pt.pi_b_chem;
    }
}

//! this function reduce the size of net baryon diffusion current
//! in the dilute region to stablize numerical simulations
void Advance::QuestRevert_qmu(
    Cell_small &grid_pt, const int ieta, const int ix, const int iy) {
    double eps_scale = 0.1;  // in 1/fm^4

    double xi = 0.05;
    double factor = 10. * DATA.quest_revert_strength
                    * (1. / (exp(-(grid_pt.epsilon - eps_scale) / xi) + 1.)
                       - 1. / (exp(eps_scale / xi) + 1.));

    double q_mu_local[4];
    for (int i = 0; i < 4; i++) {
        // copy the value from the grid
        q_mu_local[i] = grid_pt.Wmunu[10 + i];
    }

    // calculate the size of q^\mu
    double q_size = 0.0;
    for (int i = 0; i < 4; i++) {
        double gfac = (i == 0 ? -1.0 : 1.0);
        q_size += gfac * q_mu_local[i] * q_mu_local[i];
    }

    // first check the positivity of q^mu q_mu
    // (in the conversion of gmn = diag(-+++))
    if (q_size < 0.0) {
        music_message << "Advance::QuestRevert_qmu: q^mu q_mu = " << q_size
                      << " < 0!";
        music_message.flush("warning");
        music_message << "Reset it to zero!!!!";
        music_message.flush("warning");
        for (int i = 0; i < 4; i++) {
            int idx_1d = map_2d_idx_to_1d(4, i);
            grid_pt.Wmunu[idx_1d] = 0.0;
        }
    }

    // reduce the size of q^mu according to rhoB
    double e_local = grid_pt.epsilon;
    double rhob_local = grid_pt.rhob;
    double rho_q = sqrt(q_size / (rhob_local * rhob_local)) / factor;
    double rho_q_max = 0.1;
    if (rho_q > rho_q_max) {
        if (e_local > eps_scale && DATA.echo_level > 5) {
            music_message << "ieta = " << ieta << ", ix = " << ix
                          << ", iy = " << iy
                          << ", energy density = " << e_local * hbarc
                          << "GeV/fm^3"
                          << ", rhob = " << rhob_local << "1/fm^3"
                          << "-- diffusion |q/rhob| = " << rho_q;
            music_message.flush("warning");
        }
        for (int i = 0; i < 4; i++) {
            grid_pt.Wmunu[10 + i] = (rho_q_max / rho_q) * q_mu_local[i];
        }
    }
}

//! KT flux through every face of arena in the three directions
//! (see FaceFluxes). Each face is computed once, from the two cells on
//! either side and their outer neighbours (clamped at the grid edges).
void Advance::computeFaceFluxes(
    const double tau, const Fields &arena, const std::vector<double> &pressure,
    FaceFluxes &fluxes) {
    const int nx = arena.nX();
    const int ny = arena.nY();
    const int neta = arena.nEta();
    fluxes.resize(nx, ny, neta);
    for (int direction = 1; direction <= 3; direction++) {
        const int nfx = nx + (direction == 1);
        const int nfy = ny + (direction == 2);
        const int nfeta = neta + (direction == 3);
#pragma omp parallel for collapse(3) schedule(guided)
        for (int feta = 0; feta < nfeta; feta++)
            for (int fy = 0; fy < nfy; fy++)
                for (int fx = 0; fx < nfx; fx++) {
                    TJbVec flux;
                    computeFaceFlux(
                        tau, arena, pressure, direction, fx, fy, feta, flux);
                    const int faceIdx =
                        fluxes.getFaceIdx(direction, fx, fy, feta);
                    for (int alpha = 0; alpha < 7; alpha++) {
                        fluxes.flux(direction, alpha, faceIdx) = flux[alpha];
                    }
                }
    }
}

//! KT flux through the face (fx, fy, feta) normal to direction: the face
//! between the cells with coordinate f - 1 (L) and f (R) along direction
void Advance::computeFaceFlux(
    const double tau, const Fields &arena, const std::vector<double> &pressure,
    const int direction, const int fx, const int fy, const int feta,
    TJbVec &flux) {
    const double tau_fac[4] = {0.0, tau, tau, 1.0};
    const int dx = (direction == 1);
    const int dy = (direction == 2);
    const int deta = (direction == 3);
    const int ILL =
        arena.getFieldIdxHalo(fx - 2 * dx, fy - 2 * dy, feta - 2 * deta);
    const int IL = arena.getFieldIdxHalo(fx - dx, fy - dy, feta - deta);
    const int IR = arena.getFieldIdxHalo(fx, fy, feta);
    const int IRR = arena.getFieldIdxHalo(fx + dx, fy + dy, feta + deta);
    const auto cLL = arena.getCellIdeal(ILL);
    const auto cL = arena.getCellIdeal(IL);
    const auto cR = arena.getCellIdeal(IR);
    const auto cRR = arena.getCellIdeal(IRR);

    // the half-way states on the left and right side of the face
    TJbVec qL = {0.};
    TJbVec qR = {0.};
    for (int alpha = 0; alpha < 7; alpha++) {
        const double gLL = tau * get_TJb(cLL, alpha, 0, pressure[ILL]);
        const double gL = tau * get_TJb(cL, alpha, 0, pressure[IL]);
        const double gR = tau * get_TJb(cR, alpha, 0, pressure[IR]);
        const double gRR = tau * get_TJb(cRR, alpha, 0, pressure[IRR]);
        const double fL = 0.5 * minmod.minmod_dx(gR, gL, gLL);
        const double fR = -0.5 * minmod.minmod_dx(gRR, gR, gL);
        qL[alpha] = gL + fL;
        qR[alpha] = gR + fR;
    }

    // reconstruct e, rhob, and u[4] for the half-way states. Each uses the
    // cell on its own side of the face as the guess (and fallback)
    auto gridL = reconst_helper.ReconstIt_shell(tau, qL, cL);
    auto gridR = reconst_helper.ReconstIt_shell(tau, qR, cR);

    double pressureL, pressureR;  // of the half-way states, from MaxSpeed
    const double aL = MaxSpeed(tau, direction, gridL, pressureL);
    const double aR = MaxSpeed(tau, direction, gridR, pressureR);
    const double a = std::max(aL, aR);
    for (int alpha = 0; alpha < 7; alpha++) {
        const double FL =
            (tau_fac[direction] * get_TJb(gridL, alpha, direction, pressureL));
        const double FR =
            (tau_fac[direction] * get_TJb(gridR, alpha, direction, pressureR));
        // KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2}))/2
        //                  - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2
        flux[alpha] = 0.5 * ((FL + FR) - a * (qR[alpha] - qL[alpha]));
    }
}

//! This function computes the rhs array. It computes the spatial
//! derivatives of T^\mu\nu using the KT algorithm, from the face fluxes
//! of computeFaceFluxes
void Advance::MakeDeltaQI(
    const FaceFluxes &fluxes, const int ix, const int iy, const int ieta,
    TJbVec &qi) {
    const double delta[4] = {0.0, DATA.delta_x, DATA.delta_y, DATA.delta_eta};

    TJbVec rhs = {0.};
    EnergyFlowVec T_eta_m = {0.};
    EnergyFlowVec T_eta_p = {0.};
    for (int direction = 1; direction <= 3; direction++) {
        const int dx = (direction == 1);
        const int dy = (direction == 2);
        const int deta = (direction == 3);
        const int Imh = fluxes.getFaceIdx(direction, ix, iy, ieta);
        const int Iph =
            fluxes.getFaceIdx(direction, ix + dx, iy + dy, ieta + deta);
        for (int alpha = 0; alpha < 7; alpha++) {
            const double Fiph = fluxes.flux(direction, alpha, Iph);
            const double Fimh = fluxes.flux(direction, alpha, Imh);
            if (direction == 3 && (alpha == 0 || alpha == 3)) {
                T_eta_m[alpha] = Fimh;
                T_eta_p[alpha] = Fiph;
            } else {
                double DFmmp = (Fimh - Fiph) / delta[direction];
                rhs[alpha] += DFmmp * (DATA.delta_tau);
            }
        }
    }

    // add longitudinal flux with discretized geometric terms
    double cosh_deta =
        cosh(delta[3] / 2.) / std::max(delta[3], Util::small_eps);
    double sinh_deta =
        sinh(delta[3] / 2.) / std::max(delta[3], Util::small_eps);
    sinh_deta = std::max(0.5, sinh_deta);
    if (DATA.boost_invariant) {
        // if the simulation is boost-invariant,
        // we directly use the limiting value at \Delta eta = 0
        // Longitudinal derivatives should be 0, we set cosh_eta = 0 here
        cosh_deta = 0.0;
        sinh_deta = 0.5;
    }
    rhs[0] +=
        (((T_eta_m[0] - T_eta_p[0]) * cosh_deta
          - (T_eta_m[3] + T_eta_p[3]) * sinh_deta)
         * DATA.delta_tau);
    rhs[3] +=
        (((T_eta_m[3] - T_eta_p[3]) * cosh_deta
          - (T_eta_m[0] + T_eta_p[0]) * sinh_deta)
         * DATA.delta_tau);

    // geometric terms
    // rhs[0] -= get_TJb(arena_current(ix, iy, ieta), 3, 3)*DATA.delta_tau;
    // rhs[3] -= get_TJb(arena_current(ix, iy, ieta), 3, 0)*DATA.delta_tau;

    for (int i = 0; i < 7; i++) {
        qi[i] += rhs[i];
    }
}

// determine the maximum signal propagation speed at the given direction
double Advance::MaxSpeed(
    const double tau, const int direc, const ReconstCell &grid_p,
    double &pressure) {
    double g[] = {1., 1., 1. / tau};

    double utau = grid_p.u[0];
    double utau2 = utau * utau;
    double ux = std::abs(grid_p.u[direc]);
    double ut2mux2 = utau2 - ux * ux;

    double eps = grid_p.e;
    double rhob = grid_p.rhob;
    double rhoq = grid_p.rhoq;
    double rhos = grid_p.rhos;

    // double vs2 = eos.get_cs2(eps, rhob);
    double dpde, dpdrhob, dpdrhoq, dpdrhos, vs2;
    eos.get_pressure_with_gradients_and_cs2(
        eps, rhob, rhoq, rhos, pressure, dpde, dpdrhob, dpdrhoq, dpdrhos, vs2);
    double num_temp_sqrt = (ut2mux2 - (ut2mux2 - 1.) * vs2) * vs2;
    double num;
    if (num_temp_sqrt >= 0) {
        num = utau * ux * (1. - vs2) + sqrt(num_temp_sqrt);
    } else {
        // double dpde = eos.get_dpde(eps, rhob);
        double h = pressure + eps;
        if (dpde < 0.001) {
            num =
                (sqrt(-(h * dpde * h * (dpde * (-1.0 + ut2mux2) - ut2mux2)))
                 - h * (-1.0 + dpde) * utau * ux);
        } else {
            fprintf(stderr, "WARNING: in MaxSpeed. \n");
            fprintf(
                stderr, "Expression under sqrt in num=%lf. \n", num_temp_sqrt);
            fprintf(stderr, "at value e=%lf. \n", eps);
            fprintf(stderr, "at value p=%lf. \n", pressure);
            fprintf(stderr, "at value h=%lf. \n", h);
            fprintf(stderr, "at value rhob=%lf. \n", rhob);
            fprintf(stderr, "at value rhoq=%lf. \n", rhoq);
            fprintf(stderr, "at value rhos=%lf. \n", rhos);
            fprintf(stderr, "at value utau=%lf. \n", utau);
            fprintf(stderr, "at value uk=%lf. \n", ux);
            fprintf(stderr, "at value vs^2=%lf. \n", vs2);
            fprintf(stderr, "at value dpde=%lf. \n", dpde);
            fprintf(stderr, "at value dpdrhob=%lf. \n", dpdrhob);
            fprintf(stderr, "at value dpdrhoq=%lf. \n", dpdrhoq);
            fprintf(stderr, "at value dpdrhos=%lf. \n", dpdrhos);
            fprintf(stderr, "MaxSpeed: exiting.\n");
            exit(1);
        }
    }
    double den = utau2 * (1. - vs2) + vs2;
    double f = num / std::max(den, Util::small_eps);
    // check for problems
    if (f < 0.0) {
        fprintf(stderr, "SpeedMax = %e\n is negative.\n", f);
        fprintf(stderr, "Can't happen.\n");
        exit(1);
    } else if (f < ux / utau) {
        if (num != 0.0) {
            if (fabs(f - ux / utau) < 0.0001) {
                f = ux / utau;
            } else {
                fprintf(stderr, "SpeedMax-v = %lf\n", f - ux / utau);
                fprintf(
                    stderr, "SpeedMax = %e\n is smaller than v = %e.\n", f,
                    ux / utau);
                fprintf(stderr, "Can't happen.\n");
                exit(1);
            }
        }
    } else if (f > 1.0) {
        fprintf(stderr, "SpeedMax = %e\n is bigger than 1.\n", f);
        fprintf(stderr, "Can't happen.\n");
        fprintf(stderr, "SpeedMax = num/den, num = %e, den = %e \n", num, den);
        fprintf(stderr, "cs2 = %e \n", vs2);
        fprintf(stderr, "at value e=%e, rhob = %e \n", eps, rhob);
        fprintf(stderr, "utau = %e, ux = %e \n", utau, ux);
        f = 1.;
        exit(1);
    }
    f *= g[direc - 1];
    return f;
}

double Advance::get_TJb(
    const ReconstCell &grid_p, const int mu, const int nu,
    const double pressure) {
    // assert(mu < 5); assert(mu > -1);
    // assert(nu < 4); assert(nu > -1);

    const double u_nu = grid_p.u[nu];

    if (mu == 4) {
        return (grid_p.rhob) * u_nu;
    } else if (mu == 5) {
        return (grid_p.rhoq) * u_nu;
    } else if (mu == 6) {
        return (grid_p.rhos) * u_nu;
    }

    double e = grid_p.e;
    double gfac = 0.0;
    double u_mu = 0.0;
    if (mu == nu) {
        u_mu = u_nu;
        gfac = 1.0;
        if (mu == 0) {
            gfac = -1.0;
        }
    } else {
        u_mu = grid_p.u[mu];
    }
    const double T_munu = (e + pressure) * u_mu * u_nu + pressure * gfac;
    return (T_munu);
}

void Advance::solveEigenvaluesWmunu(Cell_small &grid_pt) {
    Matrix4d A = Matrix4d::Zero();

    A(0, 0) = -grid_pt.Wmunu[0];
    A(0, 1) = grid_pt.Wmunu[1];
    A(0, 2) = grid_pt.Wmunu[2];
    A(0, 3) = grid_pt.Wmunu[3];

    A(1, 0) = -grid_pt.Wmunu[1];
    A(1, 1) = grid_pt.Wmunu[4];
    A(1, 2) = grid_pt.Wmunu[5];
    A(1, 3) = grid_pt.Wmunu[6];

    A(2, 0) = -grid_pt.Wmunu[2];
    A(2, 1) = grid_pt.Wmunu[5];
    A(2, 2) = grid_pt.Wmunu[7];
    A(2, 3) = grid_pt.Wmunu[8];

    A(3, 0) = -grid_pt.Wmunu[3];
    A(3, 1) = grid_pt.Wmunu[6];
    A(3, 2) = grid_pt.Wmunu[8];
    A(3, 3) = grid_pt.Wmunu[9];

    EigenSolver<Matrix4d> es;
    es.compute(A, false);

    double min = es.eigenvalues().real()[0];
    double max = min;

    for (int i = 1; i < 4; i++) {
        double temp = es.eigenvalues().real()[i];
        if (temp > max) {
            max = temp;
        } else if (temp < min) {
            min = temp;
        }
    }
    grid_pt.Lambdas[0] = min;
    grid_pt.Lambdas[1] = -min - max;
    grid_pt.Lambdas[2] = max;
}

void Advance::nCausalityConstraints(Cell_small &grid_pt) {
    // check causality first, if violated, calculate alpha
    // and store them in an array. (alpha = 1 otherwise),
    // pick the min alpha between 0 and 1
    double eps = grid_pt.epsilon;
    double rhob = grid_pt.rhob;
    double cs2 = eos.get_cs2(eps, rhob, grid_pt.rhoq, grid_pt.rhos);
    double P = eos.get_pressure(eps, rhob, grid_pt.rhoq, grid_pt.rhos);
    double transportPart_n13 =
        2. * 1. / transport_coeffs_.get_shear_relax_time_factor();
    double viscousPart1_n13 = transport_coeffs_.get_lambda_pibulkPi_coeff();
    double viscousPart2_n13 = -1. / 2. * transport_coeffs_.get_tau_pipi_coeff();
    double transportPart_n56 =
        cs2 + 4. / 3. * 1. / transport_coeffs_.get_shear_relax_time_factor()
        + 1. / transport_coeffs_.get_bulk_relax_time_factor() * (1. / 3. - cs2)
              * (1. / 3. - cs2);
    double viscousPart1_n56 =
        2. / 3. * transport_coeffs_.get_lambda_pibulkPi_coeff()
        + transport_coeffs_.get_delta_bulkPibulkPi_coeff() + cs2;
    double viscousPart2_n56 =
        transport_coeffs_.get_delta_pipi_coeff()
        + 1. / 3. * transport_coeffs_.get_tau_pipi_coeff()
        + transport_coeffs_.get_lambda_bulkPipi_coeff() * (1. / 3. - cs2) + cs2;

    double n1 = transportPart_n13 + viscousPart1_n13 * grid_pt.pi_b / (eps + P)
                + viscousPart2_n13 * std::abs(grid_pt.Lambdas[0]) / (eps + P);
    double n3 = transportPart_n13 + viscousPart1_n13 * grid_pt.pi_b / (eps + P)
                + viscousPart2_n13 * grid_pt.Lambdas[2] / (eps + P);
    double n5 = transportPart_n56 + viscousPart1_n56 * grid_pt.pi_b / (eps + P)
                + viscousPart2_n56 * grid_pt.Lambdas[0] / (eps + P);
    double n6 = (1. - transportPart_n56)
                + (1. - viscousPart1_n56) * grid_pt.pi_b / (eps + P)
                + (1. - viscousPart2_n56) * grid_pt.Lambdas[2] / (eps + P);

    double minAlp = 1;
    double alp = 1.;
    if (n1 < 0) {
        alp = -transportPart_n13 / (n1 - transportPart_n13);
        minAlp = std::max(0., std::min(minAlp, alp));
    }

    if (n3 < 0) {
        alp = -transportPart_n13 / (n3 - transportPart_n13);
        minAlp = std::max(0., std::min(minAlp, alp));
    }

    if (n5 < 0) {
        alp = -transportPart_n56 / (n5 - transportPart_n56);
        minAlp = std::max(0., std::min(minAlp, alp));
    }

    if (n6 < 0) {
        alp = -(1. - transportPart_n56) / (n6 - (1. - transportPart_n56));
        minAlp = std::max(0., std::min(minAlp, alp));
    }

    for (double &pi : grid_pt.Wmunu) {
        pi = pi * minAlp;
    }
    grid_pt.pi_b = grid_pt.pi_b * minAlp;
}
