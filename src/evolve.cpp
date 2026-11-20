// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

#include "cornelius.h"
#include "emoji.h"
#include "evolve.h"
#include "u_derivative.h"
#include "util.h"

#ifndef _OPENMP
#define omp_get_thread_num() 0
#endif

using Util::hbarc;

Evolve::Evolve(
    const EOS &eosIn, InitData &DATA_in,
    std::shared_ptr<HydroSourceBase> hydro_source_ptr_in)
    : eos(eosIn),
      DATA(DATA_in),
      grid_info(DATA_in, eosIn),
      advance(eosIn, DATA_in, hydro_source_ptr_in) {
    rk_order = DATA_in.rk_order;
    if (DATA.freezeOutMethod == 4) {
        initialize_freezeout_surface_info();
    }
    hydro_source_terms_ptr = hydro_source_ptr_in;
}

// master control function for hydrodynamic evolution
int Evolve::EvolveIt(
    Fields &arenaFieldsPrev, Fields &arenaFieldsCurr, Fields &arenaFieldsNext,
    HydroinfoMUSIC &hydro_info_ptr) {
    // first pass some control parameters
    int freezeout_flag = DATA.doFreezeOut;
    int freezeout_lowtemp_flag = DATA.doFreezeOut_lowtemp;

    if (DATA.store_hydro_info_in_memory == 1) {
        hydro_info_ptr.set_grid_infomatioin(DATA);
    }

    FO_nBvsEta_.resize(arenaFieldsCurr.nEta(), 0);

    // main loop starts ...
    int itmax = DATA.nt;
    double tau0 = DATA.tau0;
    double tauMax = DATA.tau0 + (DATA.nt + 0.5) * DATA.delta_tau;

    int iFreezeStart = 0;
    double source_tau_max = 0.0;
    if (hydro_source_terms_ptr) {
        source_tau_max = hydro_source_terms_ptr->get_source_tau_max();
        if (freezeout_lowtemp_flag == 1) {
            double freezeOutTauStart =
                (hydro_source_terms_ptr->get_source_tauStart_max());
            freezeOutTauStart =
                std::min(DATA.freezeOutTauStartMax, freezeOutTauStart);
            iFreezeStart =
                (static_cast<int>((freezeOutTauStart - tau0) / DATA.delta_tau)
                 + 2);
        }
    }

    // beastMode >= 2 changes delta_tau during the evolution, so the
    // freeze-out start is tracked in tau, not in time-step index
    const double tauFreezeStart = tau0 + iFreezeStart * DATA.delta_tau;
    music_message << "Freeze-out surface starts at " << tauFreezeStart
                  << " fm/c.";
    music_message.flush("info");
    bool freezeStarted = false;
    int itFreezeStart = -1;     // time step at which freeze-out started
    double tauLastFreeze = 0.;  // tau of the stored freeze-out slice

    Fields *fpPrev = &arenaFieldsPrev;
    Fields *fpCurr = &arenaFieldsCurr;
    Fields *fpNext = &arenaFieldsNext;

    Fields freezeoutFieldPrev(
        arenaFieldsCurr.nX(), arenaFieldsCurr.nY(), arenaFieldsCurr.nEta());
    Fields freezeoutFieldCurr(
        arenaFieldsCurr.nX(), arenaFieldsCurr.nY(), arenaFieldsCurr.nEta());

    int it = 0;
    double eps_max_cur = -1.;
    int it_source_max =
        static_cast<int>((source_tau_max - tau0) / DATA.delta_tau + 0.5) + 1;
    const double max_allowed_e_increase_factor = 5.;
    double tau = tau0;
    int NtauBlock = 192;
    while (tau < tauMax) {
        const std::vector<Fields *> allFields = {
            fpPrev, fpCurr, fpNext, &freezeoutFieldPrev, &freezeoutFieldCurr};
        if (DATA.beastMode > 1 && it > 0 && (it % NtauBlock == 0)) {
            bool doubleDtau = false;
            bool coarsenGrid = false;
            if (DATA.delta_x / (2. * DATA.delta_tau) > 5) {
                doubleDtau = true;
            } else if (DATA.beastMode == 3 && DATA.delta_tau < 0.05) {
                // maximum delta_tau < 0.1 fm/c
                doubleDtau = true;
                coarsenGrid = true;
            }
            if (doubleDtau) {
                DATA.delta_tau = 2. * DATA.delta_tau;
                NtauBlock =
                    std::max(192, static_cast<int>(3.84 / DATA.delta_tau));
                DATA.facTau = std::max(1, static_cast<int>(DATA.facTau / 2.));
                DATA.output_evolution_every_N_timesteps = (std::max(
                    1, static_cast<int>(
                           DATA.output_evolution_every_N_timesteps / 2)));
                if (coarsenGrid) {
                    regridTransverse(allFields, 0., true);
                }
                // the previous slices are one old time step behind; shift
                // them to tau - delta_tau_new so that backward time
                // derivatives (curr - prev)/delta_tau stay consistent
                extrapolatePrevForNewDtau(*fpCurr, *fpPrev, 2.);
                advance.invalidateViscousCorrection();
                if (freezeStarted) {
                    extrapolatePrevForNewDtau(
                        freezeoutFieldCurr, freezeoutFieldPrev, 2.);
                }
                music_message
                    << "BeastMode: changed delta_tau = " << DATA.delta_tau
                    << " fm/c";
                if (coarsenGrid) {
                    music_message << ", dx = " << DATA.delta_x
                                  << " fm, dy = " << DATA.delta_y
                                  << " fm, nx = " << DATA.nx
                                  << ", ny = " << DATA.ny
                                  << ", x_size = " << DATA.x_size
                                  << " fm, y_size = " << DATA.y_size << " fm";
                }
                music_message.flush("info");
            }
        }
        // enlarge the transverse grid whenever the fluid above the
        // freeze-out energy density approaches the grid edge
        if (DATA.beastMode == 3 && DATA.gridPadding > 0.
            && fluidNearGridEdge(*fpCurr)) {
            // add at least one cell per side; a padding smaller than dx
            // would leave the grid unchanged and repeat this every step
            const double padding = std::max(
                DATA.gridPadding, std::max(DATA.delta_x, DATA.delta_y));
            regridTransverse(allFields, padding, false);
            advance.invalidateViscousCorrection();
            music_message << "BeastMode: enlarged the grid by " << padding
                          << " fm on each side, nx = " << DATA.nx
                          << ", ny = " << DATA.ny
                          << ", x_size = " << DATA.x_size
                          << " fm, y_size = " << DATA.y_size << " fm";
            music_message.flush("info");
        }
        int facTau = DATA.facTau;
        int Nskip_timestep = DATA.output_evolution_every_N_timesteps;
        double dt = DATA.delta_tau;

        if (hydro_source_terms_ptr) {
            hydro_source_terms_ptr->prepare_list_for_current_tau_frame(tau);
        }
        // store initial conditions
        const bool isFreezeStartStep =
            (!freezeStarted && tau > tauFreezeStart - 0.5 * dt);
        if (isFreezeStartStep) {
            freezeStarted = true;
            itFreezeStart = it;
            tauLastFreeze = tau;
            store_previous_step_for_freezeout(*fpPrev, freezeoutFieldPrev);
            store_previous_step_for_freezeout(*fpCurr, freezeoutFieldCurr);
        }
        const int itSinceFreezeStart = freezeStarted ? it - itFreezeStart : -1;

        // if (DATA.Initial_profile == 13) {
        //     if (tau >= source_tau_max + dt && tau < source_tau_max + 2*dt) {
        //         grid_info.output_energy_density_and_rhob_disitrubtion(
        //                     *fpCurr,
        //                     "energy_density_and_rhob_from_source_terms.dat");
        //     }
        // }

        int evoFileStatus = 0;
        if (it % Nskip_timestep == 0) {
            if (DATA.outputEvolutionData == 1) {
                grid_info.OutputEvolutionDataXYEta(*fpCurr, tau);
            } else if (DATA.outputEvolutionData == 2) {
                evoFileStatus =
                    grid_info.OutputEvolutionDataXYEta_chun(*fpCurr, tau);
            } else if (DATA.outputEvolutionData == 3) {
                grid_info.OutputEvolutionDataXYEta_photon(*fpCurr, tau);
            } else if (DATA.outputEvolutionData == 4) {
                grid_info.OutputEvolutionDataXYEta_vorticity(
                    *fpCurr, *fpPrev, tau);
            } else if (DATA.outputEvolutionData == 7) {
                grid_info.OutputEvolutionDataXYEta_MLtraining(*fpCurr, tau);
            }
            if (DATA.store_hydro_info_in_memory == 1) {
                grid_info.OutputEvolutionDataXYEta_memory(
                    *fpCurr, tau, hydro_info_ptr);
            }

            // if (DATA.output_movie_flag == 1) {
            //     grid_info.output_evolution_for_movie(*fpCurr, tau);
            // }
            if (DATA.output_outofequilibriumsize == 1) {
                grid_info.OutputEvolution_Knudsen_Reynoldsnumbers(*fpCurr, tau);
            }
        }
        if (evoFileStatus == -1) {
            DATA.reRunHydro = true;
            return (-1);
        }

        if (itSinceFreezeStart == 0 || itSinceFreezeStart == 10
            || itSinceFreezeStart == 30 || itSinceFreezeStart == 50) {
            grid_info.output_momentum_anisotropy_vs_etas(
                itSinceFreezeStart, tau, *fpCurr);
        }

        if (DATA.beastMode == 0) {
            if (DATA.Initial_profile == 0) {
                if (fabs(tau - 1.0) < 1e-8 || fabs(tau - 1.2) < 1e-8
                    || fabs(tau - 1.5) < 1e-8 || fabs(tau - 2.0) < 1e-8
                    || fabs(tau - 3.0) < 1e-8) {
                    grid_info.Gubser_flow_check_file(*fpCurr, tau);
                }
            } else if (DATA.Initial_profile == 1) {
                if (fabs(tau - 1.0) < 1e-8 || fabs(tau - 2.0) < 1e-8
                    || fabs(tau - 5.0) < 1e-8 || fabs(tau - 10.0) < 1e-8
                    || fabs(tau - 20.0) < 1e-8) {
                    grid_info.output_1p1D_check_file(*fpCurr, tau);
                }
            }

            grid_info.output_momentum_anisotropy_vs_tau(
                tau, -0.5, 0.5, *fpCurr);

            if (DATA.Initial_profile == 13 || DATA.Initial_profile == 131) {
                grid_info.output_average_phase_diagram_trajectory(
                    tau, -0.5, 0.5, *fpCurr);
                grid_info.output_average_phase_diagram_trajectory(
                    tau, 0.5, 2.0, *fpCurr);
                grid_info.output_average_phase_diagram_trajectory(
                    tau, 2.0, 3.0, *fpCurr);
                grid_info.output_average_phase_diagram_trajectory(
                    tau, 3.0, 4.0, *fpCurr);
                grid_info.output_average_phase_diagram_trajectory(
                    tau, 4.0, 5.0, *fpCurr);
            }
            // check energy conservation
            if (!DATA.boost_invariant) {
                grid_info.check_conservation_law(*fpCurr, *fpPrev, tau);
                if (DATA.output_vorticity == 1) {
                    if (std::abs(tau - 1.0) < 1e-8 || std::abs(tau - 2.0) < 1e-8
                        || std::abs(tau - 5.0) < 1e-8
                        || std::abs(tau - 10.0) < 1e-8) {
                        grid_info.output_vorticity_distribution(
                            *fpCurr, *fpPrev, tau, -0.5, 0.5);
                    }
                    grid_info.compute_angular_momentum(
                        *fpCurr, *fpPrev, tau, -0.5, 0.5);
                    grid_info.output_vorticity_time_evolution(
                        *fpCurr, *fpPrev, tau, -0.5, 0.5);
                    grid_info.compute_angular_momentum(
                        *fpCurr, *fpPrev, tau, -1.0, 1.0);
                    grid_info.output_vorticity_time_evolution(
                        *fpCurr, *fpPrev, tau, -1.0, 1.0);
                    grid_info.compute_angular_momentum(
                        *fpCurr, *fpPrev, tau, -DATA.eta_size / 2.,
                        DATA.eta_size / 2.);
                    grid_info.output_vorticity_time_evolution(
                        *fpCurr, *fpPrev, tau, -DATA.eta_size / 2.,
                        DATA.eta_size / 2.);
                }
            }

            if (DATA.output_OAM_density_Evolution == 2
                || (DATA.output_OAM_density_Evolution == 1
                    && it == it_source_max)) {
                grid_info.compute_Lmunu(*fpCurr, *fpPrev, tau);
            }

            double emax_loc = 0.;
            double Tmax_curr = 0.;
            double nB_max_curr = 0.;
            grid_info.get_maximum_energy_density(
                *fpCurr, emax_loc, nB_max_curr, Tmax_curr);
            if (tau > source_tau_max + dt && itSinceFreezeStart > 0) {
                if (eps_max_cur < 0.) {
                    eps_max_cur = emax_loc;
                } else {
                    if (emax_loc
                        > max_allowed_e_increase_factor * eps_max_cur) {
                        music_message
                            << "The maximum energy density increased by "
                            << "more than facotor of "
                            << max_allowed_e_increase_factor << "! ";
                        music_message << "This should not happen!";
                        music_message.flush("error");
                        exit(1);
                    } else {
                        eps_max_cur = std::min(emax_loc, eps_max_cur);
                    }
                }
            }

            if (DATA.output_hydro_debug_info == 1) {
                grid_info.monitor_a_fluid_cell(
                    *fpCurr, *fpPrev, 100, 100, 0, tau);
            }
        }

        // determine freeze-out surface
        int frozen = 0;
        if (freezeout_flag == 1) {
            if (freezeout_lowtemp_flag == 1 && isFreezeStartStep) {
                frozen = FreezeOut_equal_tau_Surface(tau, *fpCurr);
            }
            // avoid freeze-out at the first time step
            // Freeze-out is done every facTau*dt in tau. The interval is
            // measured from the stored slice, so it stays correct when
            // beastMode changes delta_tau and facTau mid-interval.
            if (itSinceFreezeStart > 0
                && tau - tauLastFreeze > (facTau - 0.5) * dt) {
                freezeoutDTau_ = tau - tauLastFreeze;
                tauLastFreeze = tau;
                if (!DATA.boost_invariant) {
                    frozen = FindFreezeOutSurface_Cornelius(
                        tau, *fpPrev, *fpCurr, freezeoutFieldPrev,
                        freezeoutFieldCurr);
                } else {
                    frozen = FindFreezeOutSurface_boostinvariant_Cornelius(
                        tau, *fpCurr, freezeoutFieldCurr);
                }
                store_previous_step_for_freezeout(*fpPrev, freezeoutFieldPrev);
                store_previous_step_for_freezeout(*fpCurr, freezeoutFieldCurr);
            }
            if (DATA.reRunHydro) {
                return (-1);
            }
        }

        /* execute rk steps */
        // all the evolution are at here !!!
        AdvanceRK(tau, fpPrev, fpCurr, fpNext);

        music_message << emoji::clock() << " Done time step " << it << "/"
                      << itmax << " tau = " << tau << " fm/c";
        music_message.flush("info");
        if (frozen == 1 && tau > source_tau_max) {
            if (DATA.outputEvolutionData > 1 && DATA.outputEvolutionData < 5) {
                double emax_loc = 0.;
                double Tmax_curr = 0.;
                double nB_max_curr = 0.;
                grid_info.get_maximum_energy_density(
                    *fpCurr, emax_loc, nB_max_curr, Tmax_curr);
                if (emax_loc < DATA.output_evolution_e_cut) {
                    music_message << "All cells e < "
                                  << DATA.output_evolution_e_cut
                                  << " GeV/fm^3.";
                    music_message.flush("info");
                    break;
                }
            } else {
                music_message.info("All cells frozen out. Exiting.");
                break;
            }
        }
        tau += dt;
        it++;
    }
    if (!DATA.boost_invariant) {
        std::ofstream FOinfo;
        FOinfo.open("FO_nBvseta.dat");
        for (unsigned int ieta = 0; ieta < FO_nBvsEta_.size(); ieta++) {
            double eta = (DATA.delta_eta) * ieta - (DATA.eta_size) / 2.0;
            FOinfo << std::scientific << std::setprecision(5) << eta << "  "
                   << FO_nBvsEta_[ieta] / DATA.delta_eta << std::endl;
        }
        FOinfo.close();
    }
    if (tau < tauMax) {
        music_message.info("Finished.");
    } else {
        music_message.warning("Maximum allowed time reached.");
    }
    return 1;
}

void Evolve::store_previous_step_for_freezeout(
    Fields &arenaCurr, Fields &arenaFreeze) {
    arenaFreeze.e_ = arenaCurr.e_;
    arenaFreeze.rhob_ = arenaCurr.rhob_;
    arenaFreeze.piBulk_ = arenaCurr.piBulk_;
    arenaFreeze.piBulkChem_ = arenaCurr.piBulkChem_;
    arenaFreeze.u_ = arenaCurr.u_;
    arenaFreeze.Wmunu_ = arenaCurr.Wmunu_;
}

void Evolve::regridTransverse(
    const std::vector<Fields *> &fieldList, double gridPadding, bool coarsen) {
    // update the transverse grid information in DATA and remap all fields
    // coarsen = true : double dx and dy (coarse grid stays centred at 0)
    // gridPadding > 0: add gridPadding [fm] of padding cells on both sides
    if (fieldList.empty()) return;
    const int nxOld = fieldList[0]->nX();
    const int nyOld = fieldList[0]->nY();
    // a dimension with a single cell is not coarsened (its spacing is
    // only a normalisation, and doubling it would double the content)
    const bool coarsenX = coarsen && nxOld > 1;
    const bool coarsenY = coarsen && nyOld > 1;
    if (coarsenX) DATA.delta_x = DATA.delta_x * 2;
    if (coarsenY) DATA.delta_y = DATA.delta_y * 2;
    const int gridOffset_x = static_cast<int>(gridPadding / DATA.delta_x);
    const int gridOffset_y = static_cast<int>(gridPadding / DATA.delta_y);
    // odd n -> (n+1)/2 points on the even sub-lattice,
    // even n -> n/2 points at the centres of cell pairs
    const int nxBase = coarsenX ? (nxOld + 1) / 2 : nxOld;
    const int nyBase = coarsenY ? (nyOld + 1) / 2 : nyOld;
    DATA.nx = nxBase + 2 * gridOffset_x;
    DATA.ny = nyBase + 2 * gridOffset_y;
    DATA.x_size = DATA.delta_x * (DATA.nx - 1);
    DATA.y_size = DATA.delta_y * (DATA.ny - 1);
    for (auto &field : fieldList) {
        coarseGrainAndEnlargeGrid(*field, gridPadding, coarsen);
    }
}

void Evolve::coarseGrainAndEnlargeGrid(
    Fields &arenaCurr, double gridPadding, bool coarsen) {
    // this function remaps the transverse grid onto DATA.nx x DATA.ny.
    // coarsen = true: dx and dy are doubled (DATA.delta_x/y must already be
    // the new spacings).
    //   Old grid: n points at x_i = -(n-1)dx/2 + i dx.
    //   New grid (before padding): (n+1)/2 points with spacing 2dx,
    //   centred at x = 0.
    //     n odd : new point j sits on old point 2j; use the 1/4-1/2-1/4
    //             filter over old points 2j-1, 2j, 2j+1. Neighbours outside
    //             the grid are dropped (not clamped), so every old cell
    //             gets a total weight of 1/2 and the content is conserved
    //     n = 1 : the direction is not coarsened
    //     n even: new point j sits at the midpoint of old points 2j and 2j+1;
    //             average the two
    //   Both keep the fluid at the same physical position. The filter is
    //   applied to the conserved densities T^{tau mu} (ideal + W + Pi) and
    //   J^tau, so their integrals over the plane are conserved; e, n_B and
    //   u of the coarse cell are reconstructed from them, and pi^{mu nu}
    //   and q^mu are projected to be transverse to (and pi traceless with
    //   respect to) the reconstructed u.
    // gridPadding > 0: gridPadding [fm] of padding cells is added on both
    // sides. Padding cells copy the nearest edge cell of the original fluid
    // with e, rhob, and u^i damped by a Gaussian in the distance from that
    // edge. They are set to local equilibrium (no dissipative fields), so
    // that pi^{mu nu} stays transverse to u and traceless.
    Fields arenaCopy(arenaCurr.nX(), arenaCurr.nY(), arenaCurr.nEta());
    store_previous_step_for_freezeout(arenaCurr, arenaCopy);
    const int gridOffset_x = static_cast<int>(gridPadding / DATA.delta_x);
    const int gridOffset_y = static_cast<int>(gridPadding / DATA.delta_y);
    const double sigma = std::max(1e-6, 0.2 * gridPadding);
    arenaCurr.resizeFields(DATA.nx, DATA.ny, arenaCurr.nEta());

    struct Stencil {
        int idx[3];
        double w[3];
        int n;
        double damp;
        bool isPadding;
    };
    // build the mapping from new to old indices along one direction
    auto buildStencils = [&](int nOld, int nNew, int offset, double dx) {
        const bool coarsenDir = coarsen && nOld > 1;
        const int nBase = coarsenDir ? (nOld + 1) / 2 : nOld;
        std::vector<Stencil> stencils(nNew);
        for (int i = 0; i < nNew; i++) {
            Stencil &st = stencils[i];
            const int j = i - offset;
            const int jClamped = std::max(0, std::min(nBase - 1, j));
            st.isPadding = (jClamped != j);
            st.damp = 1.0;
            if (st.isPadding) {
                double dis = std::abs(j - jClamped) * dx / sigma;
                st.damp = exp(-dis * dis / 2.);
            }
            if (!coarsenDir) {
                st.n = 1;
                st.idx[0] = jClamped;
                st.w[0] = 1.0;
            } else if (nOld % 2 == 1) {
                st.n = 0;
                const double w3[3] = {0.25, 0.5, 0.25};
                for (int k = 0; k < 3; k++) {
                    const int iOld = 2 * jClamped - 1 + k;
                    if (iOld < 0 || iOld >= nOld) continue;
                    st.idx[st.n] = iOld;
                    st.w[st.n] = w3[k];
                    st.n++;
                }
            } else {
                st.n = 2;
                st.idx[0] = 2 * jClamped;
                st.idx[1] = 2 * jClamped + 1;
                st.w[0] = 0.5;
                st.w[1] = 0.5;
            }
        }
        return stencils;
    };
    const auto stencilX = buildStencils(
        arenaCopy.nX(), arenaCurr.nX(), gridOffset_x, DATA.delta_x);
    const auto stencilY = buildStencils(
        arenaCopy.nY(), arenaCurr.nY(), gridOffset_y, DATA.delta_y);

    // W^{mu nu} index for (mu, nu) in tau, x, y, eta; q^mu is W[10 + mu]
    const int widx[4][4] = {
        {0, 1, 2, 3}, {1, 4, 5, 6}, {2, 5, 7, 8}, {3, 6, 8, 9}};
    const double gDiag[4] = {-1., 1., 1., 1.};
    // projects pi^{mu nu} and q^mu of Win onto the frame of u:
    // pi -> Delta^{mu nu}_{alpha beta} pi^{alpha beta}, q -> Delta^mu_a q^a
    auto projectDissipative = [&](const FlowVec &uu, const double *Win,
                                  double *Wout) {
        double Delta[4][4];  // Delta^mu_alpha = delta^mu_alpha + u^mu u_alpha
        for (int m = 0; m < 4; m++) {
            for (int a = 0; a < 4; a++) {
                Delta[m][a] = (m == a ? 1. : 0.) + uu[m] * uu[a] * gDiag[a];
            }
        }
        double A[4][4];
        for (int m = 0; m < 4; m++) {
            for (int n = m; n < 4; n++) {
                double sum = 0.;
                for (int a = 0; a < 4; a++) {
                    for (int b = 0; b < 4; b++) {
                        sum += Delta[m][a] * Delta[n][b] * Win[widx[a][b]];
                    }
                }
                A[m][n] = sum;
                A[n][m] = sum;
            }
        }
        double trace = 0.;
        for (int a = 0; a < 4; a++) trace += gDiag[a] * A[a][a];
        for (int m = 0; m < 4; m++) {
            for (int n = m; n < 4; n++) {
                const double DeltaUp = (m == n ? gDiag[m] : 0.) + uu[m] * uu[n];
                Wout[widx[m][n]] = A[m][n] - trace / 3. * DeltaUp;
            }
        }
        for (int m = 0; m < 4; m++) {
            double sum = 0.;
            for (int a = 0; a < 4; a++) sum += Delta[m][a] * Win[10 + a];
            Wout[10 + m] = sum;
        }
    };
    Reconst reconstHelper(eos, DATA.echo_level, DATA.beastMode);

    for (int ieta = 0; ieta < arenaCurr.nEta(); ieta++) {
        for (int iy = 0; iy < arenaCurr.nY(); iy++) {
            const Stencil &sy = stencilY[iy];
            for (int ix = 0; ix < arenaCurr.nX(); ix++) {
                const Stencil &sx = stencilX[ix];
                // densities are filtered with the stencil weights (their sum
                // is below 1 at an edge of an odd grid); u is averaged
                double e = 0., rhob = 0., piB = 0., piBChem = 0., wSum = 0.;
                double u[4] = {0., 0., 0., 0.};
                double W[14] = {0.};
                double Ttot[4] = {0., 0., 0., 0.};  // total T^{tau mu}
                double Jtot = 0.;                   // total J^tau
                for (int a = 0; a < sx.n; a++) {
                    for (int b = 0; b < sy.n; b++) {
                        const double w = sx.w[a] * sy.w[b];
                        int idxOld =
                            arenaCopy.getFieldIdx(sx.idx[a], sy.idx[b], ieta);
                        {
                            const double eOld = arenaCopy.e_[idxOld];
                            const double nOld = arenaCopy.rhob_[idxOld];
                            const double POld = eos.get_pressure(eOld, nOld);
                            const double PiOld =
                                arenaCopy.piBulk_[idxOld]
                                + arenaCopy.piBulkChem_[idxOld];
                            const double u0Old = arenaCopy.u_[0][idxOld];
                            for (int m = 0; m < 4; m++) {
                                const double g = (m == 0 ? -1. : 0.);
                                const double uu =
                                    u0Old * arenaCopy.u_[m][idxOld];
                                Ttot[m] += w
                                           * ((eOld + POld) * uu + POld * g
                                              + arenaCopy.Wmunu_[m][idxOld]
                                              + PiOld * (g + uu));
                            }
                            Jtot +=
                                w
                                * (nOld * u0Old + arenaCopy.Wmunu_[10][idxOld]);
                        }
                        e += w * arenaCopy.e_[idxOld];
                        rhob += w * arenaCopy.rhob_[idxOld];
                        piB += w * arenaCopy.piBulk_[idxOld];
                        piBChem += w * arenaCopy.piBulkChem_[idxOld];
                        for (int nu = 1; nu < 4; nu++) {
                            u[nu] += w * arenaCopy.u_[nu][idxOld];
                        }
                        wSum += w;
                        for (int nu = 0; nu < 14; nu++) {
                            W[nu] += w * arenaCopy.Wmunu_[nu][idxOld];
                        }
                    }
                }
                const bool isPadding = sx.isPadding || sy.isPadding;
                const double damp = sx.damp * sy.damp;
                const double dissFactor = isPadding ? 0. : 1.;

                int idxNew = arenaCurr.getFieldIdx(ix, iy, ieta);
                arenaCurr.e_[idxNew] = e * damp;
                if (isPadding) {
                    arenaCurr.e_[idxNew] =
                        std::max(Util::small_eps, arenaCurr.e_[idxNew]);
                }
                arenaCurr.rhob_[idxNew] = rhob * damp;
                for (int nu = 1; nu < 4; nu++) {
                    arenaCurr.u_[nu][idxNew] = u[nu] / wSum * damp;
                }
                arenaCurr.u_[0][idxNew] = sqrt(
                    1. + arenaCurr.u_[1][idxNew] * arenaCurr.u_[1][idxNew]
                    + arenaCurr.u_[2][idxNew] * arenaCurr.u_[2][idxNew]
                    + arenaCurr.u_[3][idxNew] * arenaCurr.u_[3][idxNew]);
                for (int nu = 0; nu < 14; nu++) {
                    arenaCurr.Wmunu_[nu][idxNew] = W[nu] * dissFactor;
                }
                if (coarsen && !isPadding && sx.n * sy.n > 1) {
                    // conservative remap: subtract the dissipative part of
                    // the averaged T^{tau mu}, J^tau at the coarse-cell u
                    // and reconstruct e, n_B, u from the ideal remainder;
                    // iterate since the projected W depends on u
                    const double PiTot = piB + piBChem;
                    ReconstCell cell;
                    cell.e = arenaCurr.e_[idxNew];
                    cell.rhob = arenaCurr.rhob_[idxNew];
                    for (int m = 0; m < 4; m++) {
                        cell.u[m] = arenaCurr.u_[m][idxNew];
                    }
                    double Wp[14];
                    for (int iter = 0; iter < 20; iter++) {
                        projectDissipative(cell.u, W, Wp);
                        TJbVec q;
                        for (int m = 0; m < 4; m++) {
                            const double g = (m == 0 ? -1. : 0.);
                            q[m] = Ttot[m] - Wp[m]
                                   - PiTot * (g + cell.u[0] * cell.u[m]);
                        }
                        q[4] = Jtot - Wp[10];
                        const ReconstCell next =
                            reconstHelper.ReconstIt_shell(1., q, cell);
                        double du = 0.;
                        for (int m = 0; m < 4; m++) {
                            du = std::max(du, std::abs(next.u[m] - cell.u[m]));
                        }
                        cell = next;
                        if (du < 1e-12) break;
                    }
                    projectDissipative(cell.u, W, Wp);
                    arenaCurr.e_[idxNew] = cell.e;
                    arenaCurr.rhob_[idxNew] = cell.rhob;
                    for (int m = 0; m < 4; m++) {
                        arenaCurr.u_[m][idxNew] = cell.u[m];
                    }
                    for (int nu = 0; nu < 14; nu++) {
                        arenaCurr.Wmunu_[nu][idxNew] = Wp[nu];
                    }
                }
                arenaCurr.piBulk_[idxNew] = piB * dissFactor;
                arenaCurr.piBulkChem_[idxNew] = piBChem * dissFactor;
            }
        }
    }
}

bool Evolve::fluidNearGridEdge(Fields &arena) {
    // returns true if any cell within a band near the transverse grid edges
    // has energy density above the lowest freeze-out energy density
    if (epsFO_list.empty()) return false;
    const double eThreshold =
        (*std::min_element(epsFO_list.begin(), epsFO_list.end())) / hbarc;
    const int nx = arena.nX();
    const int ny = arena.nY();
    // the freeze-out finder fails within 2*fac cells of the edge;
    // warn earlier, at least half of the padding size away from the edge
    const int bandX = std::min(
        nx / 2,
        std::max(
            2 * DATA.fac_x + 2, static_cast<int>(std::ceil(
                                    0.5 * DATA.gridPadding / DATA.delta_x))));
    const int bandY = std::min(
        ny / 2,
        std::max(
            2 * DATA.fac_y + 2, static_cast<int>(std::ceil(
                                    0.5 * DATA.gridPadding / DATA.delta_y))));
    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        for (int iy = 0; iy < ny; iy++) {
            const bool yEdge = (iy < bandY || iy >= ny - bandY);
            for (int ix = 0; ix < nx; ix++) {
                if (!yEdge && ix >= bandX && ix < nx - bandX) {
                    ix = nx - bandX - 1;  // skip the interior
                    continue;
                }
                if (arena.e_[arena.getFieldIdx(ix, iy, ieta)] > eThreshold) {
                    return true;
                }
            }
        }
    }
    return false;
}

void Evolve::extrapolatePrevForNewDtau(
    Fields &arenaCurr, Fields &arenaPrev, double dtauRatio) {
    // arenaPrev is at tau - dtau_old. Linearly extrapolate it to
    // tau - dtauRatio*dtau_old, so that (curr - prev)/dtau_new gives the
    // same backward time derivative as before the change of delta_tau.
    const int nCells = static_cast<int>(arenaCurr.e_.size());
    auto extrap = [dtauRatio](double curr, double prev) {
        return curr + dtauRatio * (prev - curr);
    };
    for (int i = 0; i < nCells; i++) {
        arenaPrev.e_[i] =
            std::max(Util::small_eps, extrap(arenaCurr.e_[i], arenaPrev.e_[i]));
        arenaPrev.rhob_[i] = extrap(arenaCurr.rhob_[i], arenaPrev.rhob_[i]);
        arenaPrev.piBulk_[i] =
            extrap(arenaCurr.piBulk_[i], arenaPrev.piBulk_[i]);
        arenaPrev.piBulkChem_[i] =
            extrap(arenaCurr.piBulkChem_[i], arenaPrev.piBulkChem_[i]);
        for (int nu = 1; nu < 4; nu++) {
            arenaPrev.u_[nu][i] =
                extrap(arenaCurr.u_[nu][i], arenaPrev.u_[nu][i]);
        }
        arenaPrev.u_[0][i] = sqrt(
            1. + arenaPrev.u_[1][i] * arenaPrev.u_[1][i]
            + arenaPrev.u_[2][i] * arenaPrev.u_[2][i]
            + arenaPrev.u_[3][i] * arenaPrev.u_[3][i]);
        for (int nu = 0; nu < 14; nu++) {
            arenaPrev.Wmunu_[nu][i] =
                extrap(arenaCurr.Wmunu_[nu][i], arenaPrev.Wmunu_[nu][i]);
        }
    }
}

void Evolve::AdvanceRK(
    double tau, Fields *&fpPrev, Fields *&fpCurr, Fields *&fpNext) {
    // control function for Runge-Kutta evolution in tau
    // loop over Runge-Kutta steps
    for (int rk_flag = 0; rk_flag < rk_order; rk_flag++) {
        advance.AdvanceIt(tau, *fpPrev, *fpCurr, *fpNext, rk_flag);
        if (rk_flag == 0) {
            Fields *temp2 = fpPrev;
            fpPrev = fpCurr;
            fpCurr = fpNext;
            fpNext = temp2;
        } else {
            Fields *temp = fpCurr;
            fpCurr = fpNext;
            fpNext = temp;
        }
    } /* loop over rk_flag */
}

// Cornelius freeze out  (C. Shen, 11/2014)
int Evolve::FindFreezeOutSurface_Cornelius(
    double tau, Fields &arena_prev, Fields &arena_current,
    Fields &arena_freezeout_prev, Fields &arena_freezeout) {
    const int neta = arena_current.nEta();
    const int fac_eta = 1;
    int intersections = 0;
    bool unresolvedHotCells = false;
    for (unsigned int i_freezesurf = 0; i_freezesurf < epsFO_list.size();
         i_freezesurf++) {
        const double epsFO = epsFO_list[i_freezesurf] / hbarc;  // 1/fm^4

        bool hitGridBoundary = false;
        // crossings of this threshold; the retry check must not see the
        // crossings of the other thresholds
        int surfaceIntersections = 0;
#pragma omp parallel for reduction(+ : surfaceIntersections) \
    reduction(|| : hitGridBoundary)
        for (int ieta = 0; ieta < (neta - fac_eta); ieta += fac_eta) {
            int thread_id = omp_get_thread_num();
            int nCells = FindFreezeOutSurface_Cornelius_XY(
                tau, ieta, arena_prev, arena_current, arena_freezeout_prev,
                arena_freezeout, thread_id, epsFO);
            if (nCells < 0) {
                hitGridBoundary = true;
            } else {
                surfaceIntersections += nCells;
            }
        }
        intersections += surfaceIntersections;
        if (!hitGridBoundary && surfaceIntersections == 0
            && hotCellsRemain(arena_current, epsFO)) {
            if (hotCellsAtGridEdge(arena_current, epsFO)) {
                // the hot region has no boundary inside the grid
                music_message << "No freeze-out surface at tau = " << tau
                              << " fm/c, but cells at the grid edge are "
                              << "above e_FO = " << epsFO * hbarc
                              << " GeV/fm^3. The grid is too small!";
                music_message.flush("error");
                hitGridBoundary = true;
            } else {
                // interior hot cells missed by the sampling stride: a
                // larger grid does not help; keep evolving until they cool
                unresolvedHotCells = true;
                music_message << "No freeze-out surface at tau = " << tau
                              << " fm/c, but interior cells are above e_FO = "
                              << epsFO * hbarc << " GeV/fm^3 (missed by "
                              << "the freeze-out sampling stride). "
                              << "Continuing the evolution.";
                music_message.flush("warning");
            }
        }
        if (hitGridBoundary) {
            DATA.reRunHydro = true;
        }
        if (DATA.reRunHydro) {
            return (0);
        }
    }

    // 1 signals that the whole fireball is frozen out
    if (unresolvedHotCells) intersections++;
    return (intersections + 1);
}

bool Evolve::hotCellsRemain(const Fields &arena, const double epsFO) const {
    // true if any cell is above the freeze-out energy density
    for (const double e : arena.e_) {
        if (e > epsFO) return true;
    }
    return false;
}

bool Evolve::hotCellsAtGridEdge(const Fields &arena, const double epsFO) const {
    // true if any cell above the freeze-out energy density lies in the
    // transverse band where the Cornelius finders cannot place a surface
    // (cf. "Freeze-out cell at the boundary"). The finders start cubes at
    // multiples of fac and flag the cube at index 0, which spans [0, fac],
    // and every cube starting at index >= n - 2 fac; cells past the last
    // sampled corner are never sampled at all. The band is therefore
    // ix <= fac and ix >= the first multiple of fac that is >= n - 2 fac
    // (likewise for iy), i.e. every cell of a flagged cube, corners
    // included. The corners cannot be left out: the finders' intersection
    // test only compares each corner at the current time with the opposite
    // corner at the freeze-out time, so "no surface found" does not mean
    // that all corners at the current time are on the same side of e_FO.
    const int nx = arena.nX();
    const int ny = arena.nY();
    auto firstFlaggedCube = [](int n, int fac) {
        const int lowest = n - 2 * fac;
        if (lowest <= 0) return (0);
        return (((lowest + fac - 1) / fac) * fac);
    };
    const int lowX = std::min(nx, DATA.fac_x + 1);
    const int lowY = std::min(ny, DATA.fac_y + 1);
    const int highStartX = std::max(lowX, firstFlaggedCube(nx, DATA.fac_x));
    const int highStartY = std::max(lowY, firstFlaggedCube(ny, DATA.fac_y));
    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        for (int iy = 0; iy < ny; iy++) {
            const bool yEdge = (iy < lowY || iy >= highStartY);
            for (int ix = 0; ix < nx; ix++) {
                if (!yEdge && ix >= lowX && ix < highStartX) {
                    ix = highStartX - 1;  // skip the interior
                    continue;
                }
                if (arena.e_[arena.getFieldIdx(ix, iy, ieta)] > epsFO) {
                    return true;
                }
            }
        }
    }
    return false;
}

int Evolve::FindFreezeOutSurface_Cornelius_XY(
    double tau, int ieta, Fields &arena_prev, Fields &arena_current,
    Fields &arena_freezeout_prev, Fields &arena_freezeout, int thread_id,
    double epsFO) {
    const bool surface_in_binary = DATA.freeze_surface_in_binary;
    const int nx = arena_current.nX();
    const int ny = arena_current.nY();

    std::stringstream strs_name;
    strs_name << "surface_eps_" << std::setprecision(4) << epsFO * hbarc << "_"
              << thread_id << ".dat";
    std::ofstream s_file;
    std::ios_base::openmode modes;

    modes = std::ios::out | std::ios::app;
    if (surface_in_binary) {
        modes = modes | std::ios::binary;
    }

    s_file.open(strs_name.str().c_str(), modes);

    const int dim = 4;
    int intersections = 0;

    int fac_x = DATA.fac_x;
    int fac_y = DATA.fac_y;
    int fac_eta = 1;

    // tau distance between arena_freezeout and arena_current
    const double DTAU = freezeoutDTau_;
    const double DX = fac_x * DATA.delta_x;
    const double DY = fac_y * DATA.delta_y;
    const double DETA = fac_eta * DATA.delta_eta;

    U_derivative u_derivative_helper(DATA, eos);

    // initialize Cornelius
    double lattice_spacing[4] = {DTAU, DX, DY, DETA};
    std::shared_ptr<Cornelius> cornelius_ptr(new Cornelius());
    cornelius_ptr->init(dim, epsFO, lattice_spacing);

    // initialize the hyper-cube for Cornelius
    Cell_small ****fluid_cube = new Cell_small ***[2];
    Cell_aux ****fluid_aux_cube = new Cell_aux ***[2];
    double ****cube = new double ***[2];
    for (int i = 0; i < 2; i++) {
        fluid_cube[i] = new Cell_small **[2];
        fluid_aux_cube[i] = new Cell_aux **[2];
        cube[i] = new double **[2];
        for (int j = 0; j < 2; j++) {
            fluid_cube[i][j] = new Cell_small *[2];
            fluid_aux_cube[i][j] = new Cell_aux *[2];
            cube[i][j] = new double *[2];
            for (int k = 0; k < 2; k++) {
                fluid_cube[i][j][k] = new Cell_small[2];
                fluid_aux_cube[i][j][k] = new Cell_aux[2];
                cube[i][j][k] = new double[2];
                for (int l = 0; l < 2; l++) cube[i][j][k][l] = 0.0;
            }
        }
    }

    double x_fraction[2][4];
    double eta = (DATA.delta_eta) * ieta - (DATA.eta_size) / 2.0;
    // set when a freeze-out cell touches the grid edge; this function runs
    // inside an OpenMP loop, so it must not write the shared DATA itself
    bool hitGridBoundary = false;
    for (int ix = 0; ix < nx - fac_x && !hitGridBoundary; ix += fac_x) {
        double x = ix * (DATA.delta_x) - (DATA.x_size / 2.0);
        for (int iy = 0; iy < ny - fac_y; iy += fac_y) {
            double y = iy * (DATA.delta_y) - (DATA.y_size / 2.0);
            int Idx000 = arena_current.getFieldIdx(ix, iy, ieta);
            int Idx100 = arena_current.getFieldIdx(ix + fac_x, iy, ieta);
            int Idx110 =
                arena_current.getFieldIdx(ix + fac_x, iy + fac_y, ieta);
            int Idx010 = arena_current.getFieldIdx(ix, iy + fac_y, ieta);
            int Idx001 = arena_current.getFieldIdx(ix, iy, ieta + fac_eta);
            int Idx101 =
                arena_current.getFieldIdx(ix + fac_x, iy, ieta + fac_eta);
            int Idx111 = arena_current.getFieldIdx(
                ix + fac_x, iy + fac_y, ieta + fac_eta);
            int Idx011 =
                arena_current.getFieldIdx(ix, iy + fac_y, ieta + fac_eta);

            // judge intersection (from Bjoern)
            int intersect = 1;
            if ((arena_current.e_[Idx111] - epsFO)
                    * (arena_freezeout.e_[Idx000] - epsFO)
                > 0.)
                if ((arena_current.e_[Idx100] - epsFO)
                        * (arena_freezeout.e_[Idx011] - epsFO)
                    > 0.)
                    if ((arena_current.e_[Idx010] - epsFO)
                            * (arena_freezeout.e_[Idx101] - epsFO)
                        > 0.)
                        if ((arena_current.e_[Idx001] - epsFO)
                                * (arena_freezeout.e_[Idx110] - epsFO)
                            > 0.)
                            if ((arena_current.e_[Idx110] - epsFO)
                                    * (arena_freezeout.e_[Idx001] - epsFO)
                                > 0.)
                                if ((arena_current.e_[Idx101] - epsFO)
                                        * (arena_freezeout.e_[Idx010] - epsFO)
                                    > 0.)
                                    if ((arena_current.e_[Idx011] - epsFO)
                                            * (arena_freezeout.e_[Idx100]
                                               - epsFO)
                                        > 0.)
                                        if ((arena_current.e_[Idx000] - epsFO)
                                                * (arena_freezeout.e_[Idx111]
                                                   - epsFO)
                                            > 0.)
                                            intersect = 0;

            if (intersect == 0) continue;

            if (ix == 0 || ix >= nx - 2 * fac_x || iy == 0
                || iy >= ny - 2 * fac_y) {
                music_message << "Freeze-out cell at the boundary! "
                              << "The grid is too small! ix = " << ix
                              << " iy = " << iy;
                music_message.flush("error");
                hitGridBoundary = true;
                break;
            }

            // if intersect, prepare for the hyper-cube
            intersections++;

            cube[0][0][0][0] = arena_freezeout.e_[Idx000];
            cube[0][0][1][0] = arena_freezeout.e_[Idx010];
            cube[0][1][0][0] = arena_freezeout.e_[Idx100];
            cube[0][1][1][0] = arena_freezeout.e_[Idx110];
            cube[1][0][0][0] = arena_current.e_[Idx000];
            cube[1][0][1][0] = arena_current.e_[Idx010];
            cube[1][1][0][0] = arena_current.e_[Idx100];
            cube[1][1][1][0] = arena_current.e_[Idx110];
            cube[0][0][0][1] = arena_freezeout.e_[Idx001];
            cube[0][0][1][1] = arena_freezeout.e_[Idx011];
            cube[0][1][0][1] = arena_freezeout.e_[Idx101];
            cube[0][1][1][1] = arena_freezeout.e_[Idx111];
            cube[1][0][0][1] = arena_current.e_[Idx001];
            cube[1][0][1][1] = arena_current.e_[Idx011];
            cube[1][1][0][1] = arena_current.e_[Idx101];
            cube[1][1][1][1] = arena_current.e_[Idx111];

            // Now, the magic will happen in the Cornelius ...
            cornelius_ptr->find_surface_4d(cube);

            // get positions of the freeze-out surface
            // and interpolating results
            for (int isurf = 0; isurf < cornelius_ptr->get_Nelements();
                 isurf++) {
                // surface normal vector d^3 \sigma_\mu
                double FULLSU[4];
                for (int ii = 0; ii < 4; ii++)
                    FULLSU[ii] = cornelius_ptr->get_normal_elem(isurf, ii);

                // check the size of the surface normal vector
                if (std::abs(FULLSU[0]) > (DX * DY * DETA + 0.01)) {
                    music_message << "problem: volume in tau direction "
                                  << std::abs(FULLSU[0])
                                  << "  > DX*DY*DETA = " << DX * DY * DETA;
                    music_message.flush("warning");
                }
                if (std::abs(FULLSU[1]) > (DTAU * DY * DETA + 0.01)) {
                    music_message << "problem: volume in x direction "
                                  << std::abs(FULLSU[1])
                                  << "  > DTAU*DY*DETA = " << DTAU * DY * DETA;
                    music_message.flush("warning");
                }
                if (std::abs(FULLSU[2]) > (DX * DTAU * DETA + 0.01)) {
                    music_message << "problem: volume in y direction "
                                  << std::abs(FULLSU[2])
                                  << "  > DX*DTAU*DETA = " << DX * DTAU * DETA;
                    music_message.flush("warning");
                }
                if (std::abs(FULLSU[3]) > (DX * DY * DTAU + 0.01)) {
                    music_message << "problem: volume in eta direction "
                                  << std::abs(FULLSU[3])
                                  << "  > DX*DY*DTAU = " << DX * DY * DTAU;
                    music_message.flush("warning");
                }

                // position of the freeze-out fluid cell
                for (int ii = 0; ii < 4; ii++) {
                    x_fraction[1][ii] =
                        cornelius_ptr->get_centroid_elem(isurf, ii);
                    x_fraction[0][ii] = lattice_spacing[ii] - x_fraction[1][ii];
                }
                const double tau_center = tau - DTAU + x_fraction[1][0];
                const double x_center = x + x_fraction[1][1];
                const double y_center = y + x_fraction[1][2];
                const double eta_center = eta + x_fraction[1][3];

                // perform 4-d linear interpolation for all fluid quantities
                for (int ii = 0; ii < 2; ii++)
                    for (int jj = 0; jj < 2; jj++)
                        for (int kk = 0; kk < 2; kk++) {
                            int fieldIdx = arena_current.getFieldIdx(
                                ix + ii * fac_x, iy + jj * fac_y,
                                ieta + kk * fac_eta);
                            fluid_cube[0][ii][jj][kk] =
                                (arena_freezeout.getCell(fieldIdx));
                            fluid_cube[1][ii][jj][kk] =
                                arena_current.getCell(fieldIdx);

                            if (DATA.output_vorticity == 0) continue;

                            // compute the vorticity tensors
                            Cell_aux aux_tmp;
                            double eta_local = eta + kk * DETA;
                            u_derivative_helper.compute_vorticity_shell(
                                tau, arena_prev, arena_current,
                                ieta + kk * fac_eta, ix + ii * fac_x,
                                iy + jj * fac_y, fieldIdx, eta_local,
                                aux_tmp.omega_kSP, aux_tmp.omega_k,
                                aux_tmp.omega_th, aux_tmp.omega_T,
                                aux_tmp.sigma_th, aux_tmp.DbetaMu);
                            fluid_aux_cube[1][ii][jj][kk] = aux_tmp;
                            u_derivative_helper.compute_vorticity_shell(
                                tau - DTAU, arena_freezeout_prev,
                                arena_freezeout, ieta + kk * fac_eta,
                                ix + ii * fac_x, iy + jj * fac_y, fieldIdx,
                                eta_local, aux_tmp.omega_kSP, aux_tmp.omega_k,
                                aux_tmp.omega_th, aux_tmp.omega_T,
                                aux_tmp.sigma_th, aux_tmp.DbetaMu);
                            fluid_aux_cube[0][ii][jj][kk] = aux_tmp;
                        }
                auto fluid_center = four_dimension_linear_interpolation(
                    lattice_spacing, x_fraction, fluid_cube);
                Cell_aux fluid_aux_center;
                if (DATA.output_vorticity == 1) {
                    fluid_aux_center = four_dimension_linear_interpolation(
                        lattice_spacing, x_fraction, fluid_aux_cube);
                }

                // reconstruct q^\tau from the transverality criteria
                FlowVec u_flow = fluid_center.u;
                double q_mu[4] = {
                    fluid_center.Wmunu[10], fluid_center.Wmunu[11],
                    fluid_center.Wmunu[12], fluid_center.Wmunu[13]};
                double q_regulated[4] = {0.0, 0.0, 0.0, 0.0};
                regulate_qmu(u_flow, q_mu, q_regulated);
                fluid_center.Wmunu[10] = q_regulated[0];
                fluid_center.Wmunu[11] = q_regulated[1];
                fluid_center.Wmunu[12] = q_regulated[2];
                fluid_center.Wmunu[13] = q_regulated[3];

                // regulate Wmunu according to transversality and traceless
                double Wmunu_input[4][4];
                double Wmunu_regulated[4][4];
                Wmunu_input[0][0] = fluid_center.Wmunu[0];
                Wmunu_input[0][1] = Wmunu_input[1][0] = fluid_center.Wmunu[1];
                Wmunu_input[0][2] = Wmunu_input[2][0] = fluid_center.Wmunu[2];
                Wmunu_input[0][3] = Wmunu_input[3][0] = fluid_center.Wmunu[3];
                Wmunu_input[1][1] = fluid_center.Wmunu[4];
                Wmunu_input[1][2] = Wmunu_input[2][1] = fluid_center.Wmunu[5];
                Wmunu_input[1][3] = Wmunu_input[3][1] = fluid_center.Wmunu[6];
                Wmunu_input[2][2] = fluid_center.Wmunu[7];
                Wmunu_input[2][3] = Wmunu_input[3][2] = fluid_center.Wmunu[8];
                Wmunu_input[3][3] = fluid_center.Wmunu[9];
                regulate_Wmunu(u_flow, Wmunu_input, Wmunu_regulated);
                fluid_center.Wmunu[0] = Wmunu_regulated[0][0];
                fluid_center.Wmunu[1] = Wmunu_regulated[0][1];
                fluid_center.Wmunu[2] = Wmunu_regulated[0][2];
                fluid_center.Wmunu[3] = Wmunu_regulated[0][3];
                fluid_center.Wmunu[4] = Wmunu_regulated[1][1];
                fluid_center.Wmunu[5] = Wmunu_regulated[1][2];
                fluid_center.Wmunu[6] = Wmunu_regulated[1][3];
                fluid_center.Wmunu[7] = Wmunu_regulated[2][2];
                fluid_center.Wmunu[8] = Wmunu_regulated[2][3];
                fluid_center.Wmunu[9] = Wmunu_regulated[3][3];

                // 4-dimension interpolation done
                const double TFO =
                    eos.get_temperature(epsFO, fluid_center.rhob);
                if (TFO < 0) {
                    music_message << "TFO=" << TFO << "<0. ERROR. exiting.";
                    music_message.flush("error");
                    exit(1);
                }
                const double muB = eos.get_muB(epsFO, fluid_center.rhob);
                const double muS = eos.get_muS(epsFO, fluid_center.rhob);
                const double muQ = eos.get_muQ(epsFO, fluid_center.rhob);

                const double pressure =
                    eos.get_pressure(epsFO, fluid_center.rhob);
                const double eps_plus_p_over_T_FO = (epsFO + pressure) / TFO;

                // finally output results !!!!
                if (surface_in_binary) {
                    const int FOsize = 34 + DATA.output_vorticity * (24 + 14);
                    std::vector<float> array(FOsize);
                    array[0] = static_cast<float>(tau_center);
                    array[1] = static_cast<float>(x_center);
                    array[2] = static_cast<float>(y_center);
                    array[3] = static_cast<float>(eta_center);
                    for (int ii = 0; ii < 4; ii++)
                        array[4 + ii] = static_cast<float>(FULLSU[ii]);
                    for (int ii = 0; ii < 4; ii++)
                        array[8 + ii] = static_cast<float>(fluid_center.u[ii]);
                    array[12] = static_cast<float>(epsFO);
                    array[13] = static_cast<float>(TFO);
                    array[14] = static_cast<float>(muB);
                    array[15] = static_cast<float>(muS);
                    array[16] = static_cast<float>(muQ);
                    array[17] = static_cast<float>(eps_plus_p_over_T_FO);
                    for (int ii = 0; ii < 10; ii++) {
                        array[18 + ii] =
                            static_cast<float>(fluid_center.Wmunu[ii]);
                    }
                    array[28] = static_cast<float>(fluid_center.pi_b);
                    array[29] = static_cast<float>(fluid_center.rhob);
                    for (int ii = 0; ii < 4; ii++) {
                        array[30 + ii] =
                            static_cast<float>(fluid_center.Wmunu[10 + ii]);
                    }
                    if (DATA.output_vorticity == 1) {
                        for (int ii = 0; ii < 6; ii++) {
                            // no minus sign because its definition is
                            // opposite to the kinetic vorticity
                            array[34 + ii] = static_cast<float>(
                                fluid_aux_center.omega_kSP[ii] / TFO);
                            // the extra minus sign is from metric
                            // output quantities for g = (1, -1, -1, -1)
                            array[40 + ii] = static_cast<float>(
                                -fluid_aux_center.omega_k[ii] / TFO);
                            array[46 + ii] = static_cast<float>(
                                -fluid_aux_center.omega_th[ii]);
                            array[52 + ii] = static_cast<float>(
                                -fluid_aux_center.omega_T[ii] / TFO / TFO);
                        }
                        // the extra minus sign is from metric
                        // output quantities for g = (1, -1, -1, -1)
                        for (int ii = 0; ii < 10; ii++) {
                            array[58 + ii] = static_cast<float>(
                                -fluid_aux_center.sigma_th[ii]);
                        }
                        for (int ii = 0; ii < 4; ii++) {
                            array[68 + ii] = static_cast<float>(
                                -fluid_aux_center.DbetaMu[ii]);
                        }
                    }
                    for (int i = 0; i < FOsize; i++)
                        s_file.write((char *)&(array[i]), sizeof(float));
                } else {
                    s_file << std::scientific << std::setprecision(10)
                           << tau_center << " " << x_center << " " << y_center
                           << " " << eta_center << " " << FULLSU[0] << " "
                           << FULLSU[1] << " " << FULLSU[2] << " " << FULLSU[3]
                           << " " << fluid_center.u[0] << " "
                           << fluid_center.u[1] << " " << fluid_center.u[2]
                           << " " << fluid_center.u[3] << " " << epsFO << " "
                           << TFO << " " << muB << " " << muS << " " << muQ
                           << " " << eps_plus_p_over_T_FO << " ";
                    for (int ii = 0; ii < 10; ii++)
                        s_file << std::scientific << std::setprecision(10)
                               << fluid_center.Wmunu[ii] << " ";
                    if (DATA.turn_on_bulk) s_file << fluid_center.pi_b << " ";
                    if (DATA.turn_on_rhob) s_file << fluid_center.rhob << " ";
                    if (DATA.turn_on_diff)
                        for (int ii = 10; ii < 14; ii++)
                            s_file << std::scientific << std::setprecision(10)
                                   << fluid_center.Wmunu[ii] << " ";
                    s_file << std::endl;
                }
                double u_dot_dsigma =
                    tau_center
                    * (FULLSU[0] * fluid_center.u[0]
                       + FULLSU[1] * fluid_center.u[1]
                       + FULLSU[2] * fluid_center.u[2]
                       + FULLSU[3] * fluid_center.u[3] / tau_center);
                if (u_dot_dsigma > 0.) {
                    // only count the time-like contribution
                    FO_nBvsEta_[ieta] += u_dot_dsigma * fluid_center.rhob;
                }
            }
        }
    }
    s_file.close();

    // clean up
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                delete[] cube[i][j][k];
                delete[] fluid_cube[i][j][k];
                delete[] fluid_aux_cube[i][j][k];
            }
            delete[] cube[i][j];
            delete[] fluid_cube[i][j];
            delete[] fluid_aux_cube[i][j];
        }
        delete[] cube[i];
        delete[] fluid_cube[i];
        delete[] fluid_aux_cube[i];
    }
    delete[] cube;
    delete[] fluid_cube;
    delete[] fluid_aux_cube;
    if (hitGridBoundary) return (-1);
    return (intersections);
}

// Cornelius freeze out (C. Shen, 11/2014)
int Evolve::FreezeOut_equal_tau_Surface(double tau, Fields &arena_current) {
    // this function freeze-out fluid cells between epsFO and epsFO_low
    // on an equal time hyper-surface at the first time step
    // this function will be trigged if freezeout_lowtemp_flag == 1
    const int neta = arena_current.nEta();
    const int fac_eta = 1;

    for (unsigned int i_freezesurf = 0; i_freezesurf < epsFO_list.size();
         i_freezesurf++) {
        double epsFO = epsFO_list[i_freezesurf] / hbarc;
        if (!DATA.boost_invariant) {
#pragma omp parallel for
            for (int ieta = 0; ieta < neta - fac_eta; ieta += fac_eta) {
                int thread_id = omp_get_thread_num();
                FreezeOut_equal_tau_Surface_XY(
                    tau, ieta, arena_current, thread_id, epsFO);
            }
        } else {
            FreezeOut_equal_tau_Surface_XY(tau, 0, arena_current, 0, epsFO);
        }
    }
    return (0);
}

void Evolve::FreezeOut_equal_tau_Surface_XY(
    double tau, int ieta, Fields &arena_current, int thread_id, double epsFO) {
    const bool surface_in_binary = DATA.freeze_surface_in_binary;
    const double epsFO_low = DATA.epsFO_low / hbarc;  // 1/fm^4

    const int nx = arena_current.nX();
    const int ny = arena_current.nY();

    std::stringstream strs_name;
    if (!DATA.boost_invariant) {
        strs_name << "surface_eps_" << std::setprecision(4) << epsFO * hbarc
                  << "_" << thread_id << ".dat";
    } else {
        strs_name << "surface_eps_" << std::setprecision(4) << epsFO * hbarc
                  << ".dat";
    }
    std::ofstream s_file;
    std::ios_base::openmode modes;

    modes = std::ios::out | std::ios::app;
    if (surface_in_binary) {
        modes = modes | std::ios::binary;
    }

    s_file.open(strs_name.str().c_str(), modes);

    const int fac_x = DATA.fac_x;
    const int fac_y = DATA.fac_y;
    const int fac_eta = 1;

    const double DX = fac_x * DATA.delta_x;
    const double DY = fac_y * DATA.delta_y;
    const double DETA = fac_eta * DATA.delta_eta;

    double eta = (DATA.delta_eta) * ieta - (DATA.eta_size) / 2.0;
    for (int ix = 0; ix < nx - fac_x; ix += fac_x) {
        double x = ix * (DATA.delta_x) - (DATA.x_size / 2.0);
        for (int iy = 0; iy < ny - fac_y; iy += fac_y) {
            double y = iy * (DATA.delta_y) - (DATA.y_size / 2.0);

            int fieldIdx = arena_current.getFieldIdx(ix, iy, ieta);

            // judge intersection
            if (arena_current.e_[fieldIdx] > epsFO) continue;
            if (arena_current.e_[fieldIdx] < epsFO_low) continue;

            // surface normal vector d^3 \sigma_\mu
            const double FULLSU[] = {DX * DY * DETA, 0.0, 0.0, 0.0};

            // get positions of the freeze-out surface
            const double tau_center = tau;
            const double x_center = x;
            const double y_center = y;
            const double eta_center = eta;

            // flow velocity
            const double ux_center = arena_current.u_[1][fieldIdx];
            const double uy_center = arena_current.u_[2][fieldIdx];
            const double ueta_center =
                arena_current.u_[3][fieldIdx];  // u^eta/tau
            // reconstruct u^tau from u^i
            const double utau_center = sqrt(
                1. + ux_center * ux_center + uy_center * uy_center
                + ueta_center * ueta_center);

            // baryon density rho_b
            const double rhob_center = arena_current.rhob_[fieldIdx];

            // baryon diffusion current
            double qtau_center = arena_current.Wmunu_[10][fieldIdx];
            double qx_center = arena_current.Wmunu_[11][fieldIdx];
            double qy_center = arena_current.Wmunu_[12][fieldIdx];
            double qeta_center = arena_current.Wmunu_[13][fieldIdx];

            // reconstruct q^\tau from the transverality criteria
            FlowVec u_flow = {utau_center, ux_center, uy_center, ueta_center};
            double q_mu[4] = {qtau_center, qx_center, qy_center, qeta_center};
            double q_regulated[4] = {0.0, 0.0, 0.0, 0.0};
            regulate_qmu(u_flow, q_mu, q_regulated);
            qtau_center = q_regulated[0];
            qx_center = q_regulated[1];
            qy_center = q_regulated[2];
            qeta_center = q_regulated[3];

            // bulk viscous pressure pi_b
            const double pi_b_center = arena_current.piBulk_[fieldIdx];

            // shear viscous tensor
            double Wtautau_center = arena_current.Wmunu_[0][fieldIdx];
            double Wtaux_center = arena_current.Wmunu_[1][fieldIdx];
            double Wtauy_center = arena_current.Wmunu_[2][fieldIdx];
            double Wtaueta_center = arena_current.Wmunu_[3][fieldIdx];
            double Wxx_center = arena_current.Wmunu_[4][fieldIdx];
            double Wxy_center = arena_current.Wmunu_[5][fieldIdx];
            double Wxeta_center = arena_current.Wmunu_[6][fieldIdx];
            double Wyy_center = arena_current.Wmunu_[7][fieldIdx];
            double Wyeta_center = arena_current.Wmunu_[8][fieldIdx];
            double Wetaeta_center = arena_current.Wmunu_[9][fieldIdx];
            // regulate Wmunu according to transversality and traceless
            double Wmunu_input[4][4];
            double Wmunu_regulated[4][4];
            Wmunu_input[0][0] = Wtautau_center;
            Wmunu_input[0][1] = Wmunu_input[1][0] = Wtaux_center;
            Wmunu_input[0][2] = Wmunu_input[2][0] = Wtauy_center;
            Wmunu_input[0][3] = Wmunu_input[3][0] = Wtaueta_center;
            Wmunu_input[1][1] = Wxx_center;
            Wmunu_input[1][2] = Wmunu_input[2][1] = Wxy_center;
            Wmunu_input[1][3] = Wmunu_input[3][1] = Wxeta_center;
            Wmunu_input[2][2] = Wyy_center;
            Wmunu_input[2][3] = Wmunu_input[3][2] = Wyeta_center;
            Wmunu_input[3][3] = Wetaeta_center;
            regulate_Wmunu(u_flow, Wmunu_input, Wmunu_regulated);
            Wtautau_center = Wmunu_regulated[0][0];
            Wtaux_center = Wmunu_regulated[0][1];
            Wtauy_center = Wmunu_regulated[0][2];
            Wtaueta_center = Wmunu_regulated[0][3];
            Wxx_center = Wmunu_regulated[1][1];
            Wxy_center = Wmunu_regulated[1][2];
            Wxeta_center = Wmunu_regulated[1][3];
            Wyy_center = Wmunu_regulated[2][2];
            Wyeta_center = Wmunu_regulated[2][3];
            Wetaeta_center = Wmunu_regulated[3][3];

            Cell_aux fluid_aux_center;

            // get other thermodynamical quantities
            double e_local = arena_current.e_[fieldIdx];
            double T_local = eos.get_temperature(e_local, rhob_center);
            if (T_local < 0) {
                music_message << "Evolve::FreezeOut_equal_tau_Surface: "
                              << "T_local = " << T_local
                              << " <0. ERROR. exiting.";
                music_message.flush("error");
                exit(1);
            }
            double muB_local = eos.get_muB(e_local, rhob_center);
            double muS_local = eos.get_muS(e_local, rhob_center);
            double muQ_local = eos.get_muQ(e_local, rhob_center);

            double pressure = eos.get_pressure(e_local, rhob_center);
            double eps_plus_p_over_T = (e_local + pressure) / T_local;

            // finally output results
            if (surface_in_binary) {
                const int FOsize = 34 + DATA.output_vorticity * (24 + 14);
                std::vector<float> array(FOsize);
                array[0] = static_cast<float>(tau_center);
                array[1] = static_cast<float>(x_center);
                array[2] = static_cast<float>(y_center);
                array[3] = static_cast<float>(eta_center);
                array[4] = static_cast<float>(FULLSU[0]);
                array[5] = static_cast<float>(FULLSU[1]);
                array[6] = static_cast<float>(FULLSU[2]);
                array[7] = static_cast<float>(FULLSU[3]);
                array[8] = static_cast<float>(utau_center);
                array[9] = static_cast<float>(ux_center);
                array[10] = static_cast<float>(uy_center);
                array[11] = static_cast<float>(ueta_center);
                array[12] = static_cast<float>(e_local);
                array[13] = static_cast<float>(T_local);
                array[14] = static_cast<float>(muB_local);
                array[15] = static_cast<float>(muS_local);
                array[16] = static_cast<float>(muQ_local);
                array[17] = static_cast<float>(eps_plus_p_over_T);
                array[18] = static_cast<float>(Wtautau_center);
                array[19] = static_cast<float>(Wtaux_center);
                array[20] = static_cast<float>(Wtauy_center);
                array[21] = static_cast<float>(Wtaueta_center);
                array[22] = static_cast<float>(Wxx_center);
                array[23] = static_cast<float>(Wxy_center);
                array[24] = static_cast<float>(Wxeta_center);
                array[25] = static_cast<float>(Wyy_center);
                array[26] = static_cast<float>(Wyeta_center);
                array[27] = static_cast<float>(Wetaeta_center);
                array[28] = static_cast<float>(pi_b_center);
                array[29] = static_cast<float>(rhob_center);
                array[30] = static_cast<float>(qtau_center);
                array[31] = static_cast<float>(qx_center);
                array[32] = static_cast<float>(qy_center);
                array[33] = static_cast<float>(qeta_center);
                if (DATA.output_vorticity == 1) {
                    for (int ii = 0; ii < 6; ii++) {
                        // no minus sign because its definition is opposite to
                        // the kinetic vorticity
                        array[34 + ii] = static_cast<float>(
                            fluid_aux_center.omega_kSP[ii] / T_local);
                        // the extra minus sign is from metric
                        // output quantities for g = (1, -1, -1, -1)
                        array[40 + ii] = static_cast<float>(
                            -fluid_aux_center.omega_k[ii] / T_local);
                        array[46 + ii] =
                            static_cast<float>(-fluid_aux_center.omega_th[ii]);
                        array[52 + ii] = static_cast<float>(
                            -fluid_aux_center.omega_T[ii] / T_local / T_local);
                    }
                    // the extra minus sign is from metric
                    // output quantities for g = (1, -1, -1, -1)
                    for (int ii = 0; ii < 10; ii++)
                        array[58 + ii] =
                            static_cast<float>(-fluid_aux_center.sigma_th[ii]);
                    for (int ii = 0; ii < 4; ii++)
                        array[68 + ii] =
                            static_cast<float>(-fluid_aux_center.DbetaMu[ii]);
                }
                for (int i = 0; i < FOsize; i++) {
                    s_file.write((char *)&(array[i]), sizeof(float));
                }
            } else {
                s_file << std::scientific << std::setprecision(10) << tau_center
                       << " " << x_center << " " << y_center << " "
                       << eta_center << " " << FULLSU[0] << " " << FULLSU[1]
                       << " " << FULLSU[2] << " " << FULLSU[3] << " "
                       << utau_center << " " << ux_center << " " << uy_center
                       << " " << ueta_center << " " << e_local << " " << T_local
                       << " " << muB_local << " " << muS_local << " "
                       << muQ_local << " " << eps_plus_p_over_T << " "
                       << Wtautau_center << " " << Wtaux_center << " "
                       << Wtauy_center << " " << Wtaueta_center << " "
                       << Wxx_center << " " << Wxy_center << " " << Wxeta_center
                       << " " << Wyy_center << " " << Wyeta_center << " "
                       << Wetaeta_center << " ";
                if (DATA.turn_on_bulk) s_file << pi_b_center << " ";
                if (DATA.turn_on_rhob) s_file << rhob_center << " ";
                if (DATA.turn_on_diff)
                    s_file << qtau_center << " " << qx_center << " "
                           << qy_center << " " << qeta_center << " ";
                s_file << std::endl;
            }
            FO_nBvsEta_[ieta] += FULLSU[0] * rhob_center;
        }
    }
    s_file.close();
}

int Evolve::FindFreezeOutSurface_boostinvariant_Cornelius(
    double tau, Fields &arena_current, Fields &arena_freezeout) {
    const bool surface_in_binary = DATA.freeze_surface_in_binary;

    // find boost-invariant hyper-surfaces
    int n_freeze_surf = epsFO_list.size();
    int *all_frozen = new int[n_freeze_surf];
    for (int i_freezesurf = 0; i_freezesurf < n_freeze_surf; i_freezesurf++) {
        double epsFO = epsFO_list[i_freezesurf] / hbarc;

        std::stringstream strs_name;
        strs_name << "surface_eps_" << std::setprecision(4) << epsFO * hbarc
                  << ".dat";

        std::ofstream s_file;
        std::ios_base::openmode modes;

        modes = std::ios::out | std::ios::app;
        if (surface_in_binary) {
            modes = modes | std::ios::binary;
        }

        s_file.open(strs_name.str().c_str(), modes);

        const int nx = arena_current.nX();
        const int ny = arena_current.nY();
        double FULLSU[4];  // d^3 \sigma_\mu

        int intersect;
        int intersections = 0;

        int fac_x = DATA.fac_x;
        int fac_y = DATA.fac_y;

        const double DX = fac_x * DATA.delta_x;
        const double DY = fac_y * DATA.delta_y;
        const double DETA = 1.0;
        // tau distance between arena_freezeout and arena_current
        const double DTAU = freezeoutDTau_;

        double lattice_spacing[3] = {DTAU, DX, DY};
        double x_fraction[2][3];

        // initialize Cornelius
        const int dim = 3;
        std::shared_ptr<Cornelius> cornelius_ptr(new Cornelius());
        cornelius_ptr->init(dim, epsFO, lattice_spacing);

        // initialize the hyper-cube for Cornelius
        Cell_small ***fluid_cube = new Cell_small **[2];
        double ***cube = new double **[2];
        for (int i = 0; i < 2; i++) {
            fluid_cube[i] = new Cell_small *[2];
            cube[i] = new double *[2];
            for (int j = 0; j < 2; j++) {
                fluid_cube[i][j] = new Cell_small[2];
                cube[i][j] = new double[2];
                for (int k = 0; k < 2; k++) cube[i][j][k] = 0.0;
            }
        }

        bool hitGridBoundary = false;
        for (int ix = 0; ix < nx - fac_x && !hitGridBoundary; ix += fac_x) {
            double x = ix * (DATA.delta_x) - (DATA.x_size / 2.0);
            for (int iy = 0; iy < ny - fac_y; iy += fac_y) {
                double y = iy * (DATA.delta_y) - (DATA.y_size / 2.0);

                int Idx00 = arena_current.getFieldIdx(ix, iy, 0);
                int Idx10 = arena_current.getFieldIdx(ix + fac_x, iy, 0);
                int Idx11 =
                    arena_current.getFieldIdx(ix + fac_x, iy + fac_y, 0);
                int Idx01 = arena_current.getFieldIdx(ix, iy + fac_y, 0);

                // judge intersection (from Bjoern)
                intersect = 1;
                if ((arena_current.e_[Idx11] - epsFO)
                        * (arena_freezeout.e_[Idx00] - epsFO)
                    > 0.)
                    if ((arena_current.e_[Idx10] - epsFO)
                            * (arena_freezeout.e_[Idx01] - epsFO)
                        > 0.)
                        if ((arena_current.e_[Idx01] - epsFO)
                                * (arena_freezeout.e_[Idx10] - epsFO)
                            > 0.)
                            if ((arena_current.e_[Idx00] - epsFO)
                                    * (arena_freezeout.e_[Idx11] - epsFO)
                                > 0.)
                                intersect = 0;
                if (intersect == 0) continue;

                if (ix == 0 || ix >= nx - 2 * fac_x || iy == 0
                    || iy >= ny - 2 * fac_y) {
                    music_message << "Freeze-out cell at the boundary! "
                                  << "The grid is too small! ix = " << ix
                                  << " iy = " << iy;
                    music_message.flush("error");
                    hitGridBoundary = true;
                    break;
                }

                // if intersect, prepare for the hyper-cube
                intersections++;
                cube[0][0][0] = arena_freezeout.e_[Idx00];
                cube[0][0][1] = arena_freezeout.e_[Idx01];
                cube[0][1][0] = arena_freezeout.e_[Idx10];
                cube[0][1][1] = arena_freezeout.e_[Idx11];
                cube[1][0][0] = arena_current.e_[Idx00];
                cube[1][0][1] = arena_current.e_[Idx01];
                cube[1][1][0] = arena_current.e_[Idx10];
                cube[1][1][1] = arena_current.e_[Idx11];

                // Now, the magic will happen in the Cornelius ...
                cornelius_ptr->find_surface_3d(cube);

                // get positions of the freeze-out surface
                // and interpolating results
                for (int isurf = 0; isurf < cornelius_ptr->get_Nelements();
                     isurf++) {
                    // surface normal vector d^3 \sigma_\mu
                    for (int ii = 0; ii < dim; ii++)
                        FULLSU[ii] = cornelius_ptr->get_normal_elem(isurf, ii);

                    FULLSU[3] = 0.0;  // rapidity direction is set to 0

                    // check the size of the surface normal vector
                    if (fabs(FULLSU[0]) > (DX * DY * DETA + 0.01)) {
                        music_message << "problem: volume in tau direction "
                                      << fabs(FULLSU[0])
                                      << "  > DX*DY*DETA = " << DX * DY * DETA;
                        music_message.flush("warning");
                    }
                    if (fabs(FULLSU[1]) > (DTAU * DY * DETA + 0.01)) {
                        music_message
                            << "problem: volume in x direction "
                            << fabs(FULLSU[1])
                            << "  > DTAU*DY*DETA = " << DTAU * DY * DETA;
                        music_message.flush("warning");
                    }
                    if (fabs(FULLSU[2]) > (DX * DTAU * DETA + 0.01)) {
                        music_message
                            << "problem: volume in y direction "
                            << fabs(FULLSU[2])
                            << "  > DX*DTAU*DETA = " << DX * DTAU * DETA;
                        music_message.flush("warning");
                    }

                    // position of the freeze-out fluid cell
                    for (int ii = 0; ii < dim; ii++) {
                        x_fraction[1][ii] =
                            (cornelius_ptr->get_centroid_elem(isurf, ii));
                        x_fraction[0][ii] =
                            (lattice_spacing[ii] - x_fraction[1][ii]);
                    }
                    const double tau_center = tau - DTAU + x_fraction[1][0];
                    const double x_center = x + x_fraction[1][1];
                    const double y_center = y + x_fraction[1][2];
                    const double eta_center = 0.0;

                    // perform 3-d linear interpolation for all fluid quantities
                    fluid_cube[0][0][0] = arena_freezeout.getCell(ix, iy, 0);
                    fluid_cube[0][0][1] =
                        arena_freezeout.getCell(ix, iy + fac_y, 0);
                    fluid_cube[0][1][0] =
                        arena_freezeout.getCell(ix + fac_x, iy, 0);
                    fluid_cube[0][1][1] =
                        arena_freezeout.getCell(ix + fac_x, iy + fac_y, 0);
                    fluid_cube[1][0][0] = arena_current.getCell(ix, iy, 0);
                    fluid_cube[1][0][1] =
                        arena_current.getCell(ix, iy + fac_y, 0);
                    fluid_cube[1][1][0] =
                        arena_current.getCell(ix + fac_x, iy, 0);
                    fluid_cube[1][1][1] =
                        arena_current.getCell(ix + fac_x, iy + fac_y, 0);
                    auto fluid_center = three_dimension_linear_interpolation(
                        lattice_spacing, x_fraction, fluid_cube);

                    // reconstruct q^\tau from the transverality criteria
                    FlowVec u_flow = fluid_center.u;
                    double q_mu[4] = {
                        fluid_center.Wmunu[10], fluid_center.Wmunu[11],
                        fluid_center.Wmunu[12], fluid_center.Wmunu[13]};
                    double q_regulated[4] = {0.0, 0.0, 0.0, 0.0};
                    regulate_qmu(u_flow, q_mu, q_regulated);
                    fluid_center.Wmunu[10] = q_regulated[0];
                    fluid_center.Wmunu[11] = q_regulated[1];
                    fluid_center.Wmunu[12] = q_regulated[2];
                    fluid_center.Wmunu[13] = q_regulated[3];

                    // regulate Wmunu according to transversality and traceless
                    double Wmunu_input[4][4];
                    double Wmunu_regulated[4][4];
                    Wmunu_input[0][0] = fluid_center.Wmunu[0];
                    Wmunu_input[0][1] = Wmunu_input[1][0] =
                        fluid_center.Wmunu[1];
                    Wmunu_input[0][2] = Wmunu_input[2][0] =
                        fluid_center.Wmunu[2];
                    Wmunu_input[0][3] = Wmunu_input[3][0] =
                        fluid_center.Wmunu[3];
                    Wmunu_input[1][1] = fluid_center.Wmunu[4];
                    Wmunu_input[1][2] = Wmunu_input[2][1] =
                        fluid_center.Wmunu[5];
                    Wmunu_input[1][3] = Wmunu_input[3][1] =
                        fluid_center.Wmunu[6];
                    Wmunu_input[2][2] = fluid_center.Wmunu[7];
                    Wmunu_input[2][3] = Wmunu_input[3][2] =
                        fluid_center.Wmunu[8];
                    Wmunu_input[3][3] = fluid_center.Wmunu[9];
                    regulate_Wmunu(u_flow, Wmunu_input, Wmunu_regulated);
                    fluid_center.Wmunu[0] = Wmunu_regulated[0][0];
                    fluid_center.Wmunu[1] = Wmunu_regulated[0][1];
                    fluid_center.Wmunu[2] = Wmunu_regulated[0][2];
                    fluid_center.Wmunu[3] = Wmunu_regulated[0][3];
                    fluid_center.Wmunu[4] = Wmunu_regulated[1][1];
                    fluid_center.Wmunu[5] = Wmunu_regulated[1][2];
                    fluid_center.Wmunu[6] = Wmunu_regulated[1][3];
                    fluid_center.Wmunu[7] = Wmunu_regulated[2][2];
                    fluid_center.Wmunu[8] = Wmunu_regulated[2][3];
                    fluid_center.Wmunu[9] = Wmunu_regulated[3][3];

                    // 3-dimension interpolation done
                    double TFO = eos.get_temperature(epsFO, fluid_center.rhob);
                    double muB = eos.get_muB(epsFO, fluid_center.rhob);
                    double muS = eos.get_muS(epsFO, fluid_center.rhob);
                    double muQ = eos.get_muQ(epsFO, fluid_center.rhob);
                    if (TFO < 0) {
                        music_message << "TFO=" << TFO << "<0. ERROR. exiting.";
                        music_message.flush("error");
                        exit(1);
                    }

                    double pressure =
                        eos.get_pressure(epsFO, fluid_center.rhob);
                    double eps_plus_p_over_T_FO = (epsFO + pressure) / TFO;

                    // finally output results !!!!
                    if (surface_in_binary) {
                        float array[34];
                        array[0] = static_cast<float>(tau_center);
                        array[1] = static_cast<float>(x_center);
                        array[2] = static_cast<float>(y_center);
                        array[3] = static_cast<float>(eta_center);
                        for (int ii = 0; ii < 4; ii++)
                            array[4 + ii] = static_cast<float>(FULLSU[ii]);
                        for (int ii = 0; ii < 4; ii++)
                            array[8 + ii] =
                                static_cast<float>(fluid_center.u[ii]);
                        array[12] = static_cast<float>(epsFO);
                        array[13] = static_cast<float>(TFO);
                        array[14] = static_cast<float>(muB);
                        array[15] = static_cast<float>(muS);
                        array[16] = static_cast<float>(muQ);
                        array[17] = static_cast<float>(eps_plus_p_over_T_FO);
                        for (int ii = 0; ii < 10; ii++)
                            array[18 + ii] =
                                static_cast<float>(fluid_center.Wmunu[ii]);
                        array[28] = fluid_center.pi_b;
                        array[29] = fluid_center.rhob;
                        for (int ii = 0; ii < 4; ii++)
                            array[30 + ii] =
                                static_cast<float>(fluid_center.Wmunu[10 + ii]);
                        for (int i = 0; i < 34; i++)
                            s_file.write((char *)&(array[i]), sizeof(float));
                    } else {
                        s_file << std::scientific << std::setprecision(10)
                               << tau_center << " " << x_center << " "
                               << y_center << " " << eta_center << " "
                               << FULLSU[0] << " " << FULLSU[1] << " "
                               << FULLSU[2] << " " << FULLSU[3] << " "
                               << fluid_center.u[0] << " " << fluid_center.u[1]
                               << " " << fluid_center.u[2] << " "
                               << fluid_center.u[3] << " " << epsFO << " "
                               << TFO << " " << muB << " " << muS << " " << muQ
                               << " " << eps_plus_p_over_T_FO << " ";
                        for (int ii = 0; ii < 10; ii++)
                            s_file << std::scientific << std::setprecision(10)
                                   << fluid_center.Wmunu[ii] << " ";
                        if (DATA.turn_on_bulk)
                            s_file << fluid_center.pi_b << " ";
                        if (DATA.turn_on_rhob)
                            s_file << fluid_center.rhob << " ";
                        if (DATA.turn_on_diff)
                            for (int ii = 10; ii < 14; ii++)
                                s_file << std::scientific
                                       << std::setprecision(10)
                                       << fluid_center.Wmunu[ii] << " ";
                        s_file << std::endl;
                    }
                }
            }
        }

        s_file.close();

        // clean up
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                delete[] cube[i][j];
                delete[] fluid_cube[i][j];
            }
            delete[] cube[i];
            delete[] fluid_cube[i];
        }
        delete[] cube;
        delete[] fluid_cube;

        bool unresolvedHotCells = false;
        if (!hitGridBoundary && intersections == 0
            && hotCellsRemain(arena_current, epsFO)) {
            if (hotCellsAtGridEdge(arena_current, epsFO)) {
                // the hot region has no boundary inside the grid
                music_message << "No freeze-out surface at tau = " << tau
                              << " fm/c, but cells at the grid edge are "
                              << "above e_FO = " << epsFO * hbarc
                              << " GeV/fm^3. The grid is too small!";
                music_message.flush("error");
                hitGridBoundary = true;
            } else {
                // interior hot cells missed by the sampling stride; do not
                // declare this surface frozen, keep evolving
                unresolvedHotCells = true;
                music_message << "No freeze-out surface at tau = " << tau
                              << " fm/c, but interior cells are above e_FO = "
                              << epsFO * hbarc << " GeV/fm^3 (missed by "
                              << "the freeze-out sampling stride). "
                              << "Continuing the evolution.";
                music_message.flush("warning");
            }
        }
        if (hitGridBoundary) {
            // like the 3D finder: the caller stops and main re-runs on a
            // larger grid where the initial profile allows it
            DATA.reRunHydro = true;
            delete[] all_frozen;
            return (0);
        }

        // judge whether the entire fireball is freeze-out
        all_frozen[i_freezesurf] = 0;
        if (intersections == 0 && !unresolvedHotCells) {
            all_frozen[i_freezesurf] = 1;
        }
    }

    int all_frozen_flag = 1;
    for (int ii = 0; ii < n_freeze_surf; ii++) {
        all_frozen_flag *= all_frozen[ii];
    }
    if (all_frozen_flag == 1) {
        music_message.info("All cells frozen out. Exiting.");
    }

    delete[] all_frozen;
    return (all_frozen_flag);
}

void Evolve::regulate_qmu(
    const FlowVec u, const double q[], double q_regulated[]) const {
    double u_dot_q = -u[0] * q[0] + u[1] * q[1] + u[2] * q[2] + u[3] * q[3];
    for (int i = 0; i < 4; i++) {
        q_regulated[i] = q[i] + u[i] * u_dot_q;
    }
}

void Evolve::regulate_Wmunu(
    const FlowVec u, const double Wmunu[4][4],
    double Wmunu_regulated[4][4]) const {
    const double gmunu[4][4] = {
        {-1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    double u_dot_pi[4];
    double u_mu[4];
    for (int i = 0; i < 4; i++) {
        u_dot_pi[i] =
            (-u[0] * Wmunu[0][i] + u[1] * Wmunu[1][i] + u[2] * Wmunu[2][i]
             + u[3] * Wmunu[3][i]);
        u_mu[i] = gmunu[i][i] * u[i];
    }
    double tr_pi = -Wmunu[0][0] + Wmunu[1][1] + Wmunu[2][2] + Wmunu[3][3];
    double u_dot_pi_dot_u = 0.0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            u_dot_pi_dot_u += u_mu[i] * Wmunu[i][j] * u_mu[j];
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Wmunu_regulated[i][j] =
                (Wmunu[i][j] + u[i] * u_dot_pi[j] + u[j] * u_dot_pi[i]
                 + u[i] * u[j] * u_dot_pi_dot_u
                 - 1. / 3. * (gmunu[i][j] + u[i] * u[j])
                       * (tr_pi + u_dot_pi_dot_u));
        }
    }
}

void Evolve::initialize_freezeout_surface_info() {
    if (DATA.useEpsFO == 0) {
        const double e_freeze = eos.get_T2e(DATA.TFO, 0.0) * Util::hbarc;
        epsFO_list.push_back(e_freeze);
        music_message << "Freeze out at a constant temperature T = " << DATA.TFO
                      << " GeV, e_fo = " << e_freeze << " GeV/fm^3";
        music_message.flush("info");
        return;
    }
    const int freeze_eps_flag = DATA.freeze_eps_flag;
    if (freeze_eps_flag == 0) {
        // constant spacing the energy density
        int n_freeze_surf = DATA.N_freeze_out;
        double freeze_max_ed = DATA.eps_freeze_max;
        double freeze_min_ed = DATA.eps_freeze_min;
        double d_epsFO =
            ((freeze_max_ed - freeze_min_ed)
             / (n_freeze_surf - 1 + Util::small_eps));
        for (int isurf = 0; isurf < n_freeze_surf; isurf++) {
            double temp_epsFO = freeze_min_ed + isurf * d_epsFO;
            epsFO_list.push_back(temp_epsFO);
        }
    } else if (freeze_eps_flag == 1) {
        // read in from a file
        std::string eps_freeze_list_filename = DATA.freeze_list_filename;
        music_message << "read in freeze out surface information from "
                      << eps_freeze_list_filename;
        music_message.flush("info");
        std::ifstream freeze_list_file(eps_freeze_list_filename.c_str());
        if (!freeze_list_file) {
            music_message << "Evolve::initialize_freezeout_surface_info: "
                          << "can not open freeze-out list file: "
                          << eps_freeze_list_filename;
            music_message.flush("error");
            exit(1);
        }
        std::string dummy;
        double temp_epsFO, dummyd;
        getline(freeze_list_file, dummy);  // get rid of the comment
        // keep every complete line, including a last line without a
        // trailing newline (checking eof() after reading dropped it)
        while (freeze_list_file >> temp_epsFO >> dummyd >> dummyd >> dummyd
               >> dummyd >> dummyd >> dummyd) {
            epsFO_list.push_back(temp_epsFO);
        }
        freeze_list_file.close();
        if (epsFO_list.empty()) {
            music_message << "Evolve::initialize_freezeout_surface_info: "
                          << "no freeze-out energy density found in "
                          << eps_freeze_list_filename;
            music_message.flush("error");
            exit(1);
        }
    } else {
        music_message << "Evolve::initialize_freezeout_surface_info: "
                      << "unrecoginze freeze_eps_flag = " << freeze_eps_flag;
        music_message.flush("error");
        exit(1);
    }
    // the list read from file can have a different length than N_freeze_out
    music_message << "In total, " << epsFO_list.size()
                  << " freeze-out surface will be generated ...";
    music_message.flush("info");
    for (unsigned int isurf = 0; isurf < epsFO_list.size(); isurf++) {
        music_message << "eps_freezeout = " << epsFO_list[isurf]
                      << " GeV/fm^3.";
        music_message.flush("info");
    }
}

Cell_small Evolve::three_dimension_linear_interpolation(
    double *lattice_spacing, double fraction[2][3], Cell_small ***cube) {
    double denorm = 1.0;
    for (int i = 0; i < 3; i++) {
        denorm *= lattice_spacing[i];
    }

    Cell_small results;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                results = results
                          + (cube[i][j][k] * fraction[i][0] * fraction[j][1]
                             * fraction[k][2]);
            }
        }
    }
    results = results * (1. / denorm);
    return (results);
}

Cell_small Evolve::four_dimension_linear_interpolation(
    double *lattice_spacing, double fraction[2][4], Cell_small ****cube) {
    double denorm = 1.0;
    Cell_small results;
    for (int i = 0; i < 4; i++) {
        denorm *= lattice_spacing[i];
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                for (int l = 0; l < 2; l++) {
                    results =
                        results
                        + (cube[i][j][k][l] * fraction[i][0] * fraction[j][1]
                           * fraction[k][2] * fraction[l][3]);
                }
            }
        }
    }
    results = results * (1. / denorm);
    return (results);
}

Cell_aux Evolve::three_dimension_linear_interpolation(
    double *lattice_spacing, double fraction[2][3], Cell_aux ***cube) {
    double denorm = 1.0;
    for (int i = 0; i < 3; i++) {
        denorm *= lattice_spacing[i];
    }

    Cell_aux results;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                results = results
                          + (cube[i][j][k] * fraction[i][0] * fraction[j][1]
                             * fraction[k][2]);
            }
        }
    }
    results = results * (1. / denorm);
    return (results);
}

Cell_aux Evolve::four_dimension_linear_interpolation(
    double *lattice_spacing, double fraction[2][4], Cell_aux ****cube) {
    double denorm = 1.0;
    Cell_aux results;
    for (int i = 0; i < 4; i++) {
        denorm *= lattice_spacing[i];
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                for (int l = 0; l < 2; l++) {
                    results =
                        results
                        + (cube[i][j][k][l] * fraction[i][0] * fraction[j][1]
                           * fraction[k][2] * fraction[l][3]);
                }
            }
        }
    }
    results = results * (1. / denorm);
    return (results);
}
