// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <algorithm>
#include <memory>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "evolve.h"
#include "cornelius.h"
#include "u_derivative.h"
#include "emoji.h"
#include "util.h"

#ifndef _OPENMP
  #define omp_get_thread_num() 0
#endif

using Util::hbarc;

Evolve::Evolve(const EOS &eosIn, const InitData &DATA_in,
               std::shared_ptr<HydroSourceBase> hydro_source_ptr_in) :
    eos(eosIn), DATA(DATA_in),
    grid_info(DATA_in, eosIn), advance(eosIn, DATA_in, hydro_source_ptr_in) {

    rk_order  = DATA_in.rk_order;
    if (DATA.freezeOutMethod == 4) {
        initialize_freezeout_surface_info();
    }
    hydro_source_terms_ptr = hydro_source_ptr_in;
}

// master control function for hydrodynamic evolution
int Evolve::EvolveIt(SCGrid &arena_prev, SCGrid &arena_current,
                     SCGrid &arena_future, HydroinfoMUSIC &hydro_info_ptr) {
    // first pass some control parameters
    facTau                      = DATA.facTau;
    int Nskip_timestep          = DATA.output_evolution_every_N_timesteps;
    int freezeout_flag          = DATA.doFreezeOut;
    int freezeout_lowtemp_flag  = DATA.doFreezeOut_lowtemp;

    // Output information about the hydro parameters 
    // in the format of a C header file
    //if (DATA.output_hydro_params_header || DATA.outputEvolutionData == 1)
    //    grid_info.Output_hydro_information_header();

    if (DATA.store_hydro_info_in_memory == 1) {
        hydro_info_ptr.set_grid_infomatioin(DATA);
    }

    // main loop starts ...
    int    itmax = DATA.nt;
    double tau0  = DATA.tau0;
    double dt    = DATA.delta_tau;

    double tau;
    int it_start = 3;
    double source_tau_max = 0.0;
    if (!Util::weak_ptr_is_uninitialized(hydro_source_terms_ptr)) {
        source_tau_max = hydro_source_terms_ptr.lock()->get_source_tau_max();
    }

    const auto closer = [](SCGrid* g) { /*Don't delete memory we don't own*/ };
    GridPointer ap_prev   (&arena_prev, closer);
    GridPointer ap_current(&arena_current, closer);
    GridPointer ap_future (&arena_future, closer);

    SCGrid arena_freezeout_prev(arena_current.nX(),
                                arena_current.nY(),
                                arena_current.nEta());
    SCGrid arena_freezeout(arena_current.nX(),
                           arena_current.nY(),
                           arena_current.nEta());

    int it = 0;
    double eps_max_cur = -1.;
    const double max_allowed_e_increase_factor = 2.;
    for (it = 0; it <= itmax; it++) {
        tau = tau0 + dt*it;

        if (!Util::weak_ptr_is_uninitialized(hydro_source_terms_ptr)) {
            hydro_source_terms_ptr.lock()->prepare_list_for_current_tau_frame(tau);
        }
        // store initial conditions
        if (it == it_start) {
            store_previous_step_for_freezeout(*ap_prev, arena_freezeout_prev);
            store_previous_step_for_freezeout(*ap_current, arena_freezeout);
        }

        if (DATA.Initial_profile == 0) {
            if (   fabs(tau - 1.0) < 1e-8 || fabs(tau - 1.2) < 1e-8
                || fabs(tau - 1.5) < 1e-8 || fabs(tau - 2.0) < 1e-8
                || fabs(tau - 3.0) < 1e-8) {
                grid_info.Gubser_flow_check_file(*ap_current, tau);
            }
        } else if (DATA.Initial_profile == 1) {
            if (   fabs(tau -  1.0) < 1e-8 || fabs(tau -  2.0) < 1e-8
                || fabs(tau -  5.0) < 1e-8 || fabs(tau - 10.0) < 1e-8
                || fabs(tau - 20.0) < 1e-8) {
                grid_info.output_1p1D_check_file(*ap_current, tau);
            }
        }

        //if (DATA.Initial_profile == 13) {
        //    if (tau >= source_tau_max + dt && tau < source_tau_max + 2*dt) {
        //        grid_info.output_energy_density_and_rhob_disitrubtion(
        //                    *ap_current,
        //                    "energy_density_and_rhob_from_source_terms.dat");
        //    }
        //}

        if (it % Nskip_timestep == 0) {
            if (DATA.outputEvolutionData == 1) {
                grid_info.OutputEvolutionDataXYEta(*ap_current, tau);
            } else if (DATA.outputEvolutionData == 2) {
                grid_info.OutputEvolutionDataXYEta_chun(*ap_current, tau);
            } else if (DATA.outputEvolutionData == 3) {
                grid_info.OutputEvolutionDataXYEta_photon(*ap_current, tau);
            } else if (DATA.outputEvolutionData == 4) {
                grid_info.OutputEvolutionDataXYEta_vorticity(
                                            *ap_current, *ap_prev, tau);
            }

            if (DATA.output_movie_flag == 1) {
                grid_info.output_evolution_for_movie(*ap_current, tau);
            }

            if (DATA.store_hydro_info_in_memory == 1) {
                grid_info.OutputEvolutionDataXYEta_memory(*ap_current, tau,
                                                          hydro_info_ptr);
            }

            if (DATA.output_outofequilibriumsize == 1) {
                grid_info.OutputEvolution_Knudsen_Reynoldsnumbers(*ap_current,
                                                                  tau);
            }
        }

        if (it == it_start)
            grid_info.output_momentum_anisotropy_vs_etas(tau, *ap_current);
        grid_info.output_momentum_anisotropy_vs_tau(
                                            tau, -0.5, 0.5, *ap_current);
        if (DATA.Initial_profile == 13) {
            grid_info.output_average_phase_diagram_trajectory(
                                            tau, -0.5, 0.5, *ap_current);
            grid_info.output_average_phase_diagram_trajectory(
                                            tau, 0.5, 2.0, *ap_current);
            grid_info.output_average_phase_diagram_trajectory(
                                            tau, 2.0, 3.0, *ap_current);
            grid_info.output_average_phase_diagram_trajectory(
                                            tau, 3.0, 4.0, *ap_current);
            grid_info.output_average_phase_diagram_trajectory(
                                            tau, 4.0, 5.0, *ap_current);
        }


        // check energy conservation
        if (!DATA.boost_invariant) {
            grid_info.check_conservation_law(*ap_current, *ap_prev, tau);
            if (DATA.output_vorticity) {
                if (   fabs(tau -  1.0) < 1e-8 || fabs(tau -  2.0) < 1e-8
                    || fabs(tau -  5.0) < 1e-8 || fabs(tau - 10.0) < 1e-8) {
                    grid_info.output_vorticity_distribution(
                                    *ap_current, *ap_prev, tau, -0.5, 0.5);
                }
                grid_info.compute_angular_momentum(
                                    *ap_current, *ap_prev, tau, -0.5, 0.5);
                grid_info.output_vorticity_time_evolution(
                                    *ap_current, *ap_prev, tau, -0.5, 0.5);
                grid_info.compute_angular_momentum(
                                    *ap_current, *ap_prev, tau, -1.0, 1.0);
                grid_info.output_vorticity_time_evolution(
                                    *ap_current, *ap_prev, tau, -1.0, 1.0);
                grid_info.compute_angular_momentum(
                                    *ap_current, *ap_prev, tau,
                                    -DATA.eta_size/2., DATA.eta_size/2.);
                grid_info.output_vorticity_time_evolution(
                                    *ap_current, *ap_prev, tau,
                                    -DATA.eta_size/2., DATA.eta_size/2.);
            }
        }

        double emax_loc = 0.;
        double Tmax_curr = 0.;
        double nB_max_curr = 0.;
        grid_info.get_maximum_energy_density(*ap_current, emax_loc,
                                             nB_max_curr, Tmax_curr);
        if (tau > source_tau_max && it > 0) {
            if (eps_max_cur < 0.) {
                eps_max_cur = emax_loc;
            } else {
                if (emax_loc > max_allowed_e_increase_factor*eps_max_cur) {
                    music_message << "The maximum energy density increased by "
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
            grid_info.monitor_a_fluid_cell(*ap_current, *ap_prev,
                                           100, 100, 0, tau);
        }

        /* execute rk steps */
        // all the evolution are at here !!!
        AdvanceRK(tau, ap_prev, ap_current, ap_future);

        //determine freeze-out surface
        int frozen = 0;
        if (freezeout_flag == 1) {
            if (freezeout_lowtemp_flag == 1 && it == it_start) {
                frozen = FreezeOut_equal_tau_Surface(tau, *ap_current);
            }
            // avoid freeze-out at the first time step
            if ((it - it_start)%facTau == 0 && it > it_start) {
                if (!DATA.boost_invariant) {
                    frozen = FindFreezeOutSurface_Cornelius(
                                tau, *ap_prev, *ap_current,
                                arena_freezeout_prev, arena_freezeout);
                } else {
                    frozen = FindFreezeOutSurface_boostinvariant_Cornelius(
                                tau, *ap_current, arena_freezeout);
                }
                store_previous_step_for_freezeout(*ap_prev,
                                                  arena_freezeout_prev);
                store_previous_step_for_freezeout(*ap_current,
                                                  arena_freezeout);
            }
        }
        music_message << emoji::clock()
                      << " Done time step " << it << "/" << itmax
                      << " tau = " << tau << " fm/c";
        music_message.flush("info");
        if (frozen == 1 && tau > source_tau_max) {
            if (   DATA.outputEvolutionData == 2
                || DATA.outputEvolutionData == 3) {
                if (eps_max_cur < DATA.output_evolution_e_cut) {
                    music_message << "All cells e < "
                                  << DATA.output_evolution_e_cut
                                  << " GeV/fm^3.";
                    music_message.flush("info");
                    break;
                }
            } else if (DATA.outputEvolutionData == 4) {
                if (Tmax_curr < DATA.output_evolution_T_cut) {
                    music_message << "All cells T < "
                                  << DATA.output_evolution_T_cut << " GeV.";
                    music_message.flush("info");
                    break;
                }
            } else {
                music_message.info("All cells frozen out. Exiting.");
                break;
            }
        }
    }
    if (it < itmax) {
        music_message.info("Finished.");
    } else {
        music_message.warning("Maximum allowed time reached.");
    }
    return 1;
}

void Evolve::store_previous_step_for_freezeout(SCGrid &arena_current,
                                               SCGrid &arena_freezeout) {
    const int nx   = arena_current.nX();
    const int ny   = arena_current.nY();
    const int neta = arena_current.nEta();
    #pragma omp parallel for collapse(3)
    for (int ieta = 0; ieta < neta; ieta++)
    for (int ix = 0;   ix   < nx;   ix++)
    for (int iy = 0;   iy   < ny;   iy++) {
        arena_freezeout(ix, iy, ieta) = arena_current(ix, iy, ieta);
    }
}

void Evolve::AdvanceRK(double tau, GridPointer &arena_prev, GridPointer &arena_current, GridPointer &arena_future) {
    // control function for Runge-Kutta evolution in tau
    // loop over Runge-Kutta steps
    for (int rk_flag = 0; rk_flag < rk_order; rk_flag++) {
        advance.AdvanceIt(tau, *arena_prev, *arena_current, *arena_future,
                          rk_flag);
        if (rk_flag == 0) {
            auto temp     = std::move(arena_prev);
            arena_prev    = std::move(arena_current);
            arena_current = std::move(arena_future);
            arena_future  = std::move(temp);
        } else {
            std::swap(arena_current, arena_future);
        }
    }  /* loop over rk_flag */
}

// Cornelius freeze out  (C. Shen, 11/2014)
int Evolve::FindFreezeOutSurface_Cornelius(double tau,
        SCGrid &arena_prev, SCGrid &arena_current,
        SCGrid &arena_freezeout_prev, SCGrid &arena_freezeout) {
    const int neta = arena_current.nEta();
    const int fac_eta = 1;
    int intersections = 0;
    for (int i_freezesurf = 0; i_freezesurf < n_freeze_surf; i_freezesurf++) {
        const double epsFO = epsFO_list[i_freezesurf]/hbarc;   // 1/fm^4

        #pragma omp parallel for reduction(+:intersections)
        for (int ieta = 0; ieta < (neta-fac_eta); ieta += fac_eta) {
            int thread_id = omp_get_thread_num();
            intersections += FindFreezeOutSurface_Cornelius_XY(
                tau, ieta, arena_prev, arena_current,
                arena_freezeout_prev, arena_freezeout, thread_id, epsFO);
        }
    }

    return(intersections + 1);
}

int Evolve::FindFreezeOutSurface_Cornelius_XY(double tau, int ieta,
                                              SCGrid &arena_prev,
                                              SCGrid &arena_current,
                                              SCGrid &arena_freezeout_prev,
                                              SCGrid &arena_freezeout,
                                              int thread_id, double epsFO) {
    const bool surface_in_binary = DATA.freeze_surface_in_binary;
    const int nx = arena_current.nX();
    const int ny = arena_current.nY();

    std::stringstream strs_name;
    strs_name << "surface_eps_" << std::setprecision(4) << epsFO*hbarc
              << "_" << thread_id << ".dat";
    std::ofstream s_file;
    std::ios_base::openmode modes;

    if (surface_in_binary) {
        modes=std::ios::out | std::ios::binary;
    } else {
        modes=std::ios::out;
    }

    // Only append at the end of the file if it's not the first timestep
    // (that is, overwrite file at first timestep)
    if (tau != DATA.tau0+DATA.delta_tau) {
        modes = modes | std::ios::app;
    }

    s_file.open(strs_name.str().c_str(), modes);

    const int dim = 4;
    int intersections = 0;

    facTau      = DATA.facTau;   // step to skip in tau direction
    int fac_x   = DATA.fac_x;
    int fac_y   = DATA.fac_y;
    int fac_eta = 1;

    const double DTAU = facTau*DATA.delta_tau;
    const double DX   = fac_x*DATA.delta_x;
    const double DY   = fac_y*DATA.delta_y;
    const double DETA = fac_eta*DATA.delta_eta;

    U_derivative u_derivative_helper(DATA, eos);

    // initialize Cornelius
    double lattice_spacing[4] = {DTAU, DX, DY, DETA};
    std::shared_ptr<Cornelius> cornelius_ptr(new Cornelius());
    cornelius_ptr->init(dim, epsFO, lattice_spacing);

    // initialize the hyper-cube for Cornelius
    Cell_small ****fluid_cube = new Cell_small*** [2];
    Cell_aux ****fluid_aux_cube = new Cell_aux*** [2];
    double ****cube = new double*** [2];
    for (int i = 0; i < 2; i++) {
        fluid_cube[i] = new Cell_small** [2];
        fluid_aux_cube[i] = new Cell_aux** [2];
        cube[i] = new double** [2];
        for (int j = 0; j < 2; j++) {
            fluid_cube[i][j] = new Cell_small* [2];
            fluid_aux_cube[i][j] = new Cell_aux* [2];
            cube[i][j] = new double* [2];
            for (int k = 0; k < 2; k++) {
                fluid_cube[i][j][k] = new Cell_small[2];
                fluid_aux_cube[i][j][k] = new Cell_aux[2];
                cube[i][j][k] = new double[2];
                for (int l = 0; l < 2; l++)
                    cube[i][j][k][l] = 0.0;
            }
        }
    }

    double x_fraction[2][4];
    double eta = (DATA.delta_eta)*ieta - (DATA.eta_size)/2.0;
    for (int ix = 0; ix < nx - fac_x; ix += fac_x) {
        double x = ix*(DATA.delta_x) - (DATA.x_size/2.0);
        for (int iy = 0; iy < ny - fac_y; iy += fac_y) {
            double y = iy*(DATA.delta_y) - (DATA.y_size/2.0);

            // judge intersection (from Bjoern)
            int intersect = 1;
            if ((arena_current(ix+fac_x,iy+fac_y,ieta+fac_eta).epsilon-epsFO)
                *(arena_freezeout(ix,iy,ieta).epsilon-epsFO)>0.)
                if((arena_current(ix+fac_x,iy,ieta).epsilon-epsFO)
                    *(arena_freezeout(ix,iy+fac_y,ieta+fac_eta).epsilon-epsFO)>0.)
                    if((arena_current(ix,iy+fac_y,ieta).epsilon-epsFO)
                        *(arena_freezeout(ix+fac_x,iy,ieta+fac_eta).epsilon-epsFO)>0.)
                        if((arena_current(ix,iy,ieta+fac_eta).epsilon-epsFO)
                            *(arena_freezeout(ix+fac_x,iy+fac_y,ieta).epsilon-epsFO)>0.)
                            if((arena_current(ix+fac_x,iy+fac_y,ieta).epsilon-epsFO)
                                *(arena_freezeout(ix,iy,ieta+fac_eta).epsilon-epsFO)>0.)
                                if((arena_current(ix+fac_x,iy,ieta+fac_eta).epsilon-epsFO)
                                    *(arena_freezeout(ix,iy+fac_y,ieta).epsilon-epsFO)>0.)
                                    if((arena_current(ix,iy+fac_y,ieta+fac_eta).epsilon-epsFO)
                                        *(arena_freezeout(ix+fac_x,iy,ieta).epsilon-epsFO)>0.)
                                        if((arena_current(ix,iy,ieta).epsilon-epsFO)
                                            *(arena_freezeout(ix+fac_x,iy+fac_y,ieta+fac_eta).epsilon-epsFO)>0.)
                                                intersect=0;

            if (intersect==0) continue;
                
            if (ix == 0 || ix >= nx - 2*fac_x
                    || iy == 0 || iy >= ny - 2*fac_y) {
                music_message << "Freeze-out cell at the boundary! "
                              << "The grid is too small!";
                music_message.flush("error");
                exit(1);
            }

            if (ix == 0 || ix >= nx - 2*fac_x
                    || iy == 0 || iy >= ny - 2*fac_y) {
                music_message << "Freeze-out cell at the boundary! "
                              << "The grid is too small!";
                music_message.flush("error");
                exit(1);
            }

            // if intersect, prepare for the hyper-cube
            intersections++;

            cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).epsilon;
            cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).epsilon;
            cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).epsilon;
            cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).epsilon;
            cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).epsilon;
            cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).epsilon;
            cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).epsilon;
            cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).epsilon;
            cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).epsilon;
            cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).epsilon;
            cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).epsilon;
            cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).epsilon;
            cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).epsilon;
            cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).epsilon;
            cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).epsilon;
            cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).epsilon;

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
                if (std::abs(FULLSU[0]) > (DX*DY*DETA+0.01)) {
                    music_message << "problem: volume in tau direction "
                                  << std::abs(FULLSU[0]) << "  > DX*DY*DETA = "
                                  << DX*DY*DETA;
                    music_message.flush("warning");
                }
                if (std::abs(FULLSU[1]) > (DTAU*DY*DETA+0.01)) {
                    music_message << "problem: volume in x direction "
                                  << std::abs(FULLSU[1])
                                  << "  > DTAU*DY*DETA = " << DTAU*DY*DETA;
                    music_message.flush("warning");
                }
                if (std::abs(FULLSU[2]) > (DX*DTAU*DETA+0.01)) {
                    music_message << "problem: volume in y direction "
                                  << std::abs(FULLSU[2])
                                  << "  > DX*DTAU*DETA = " << DX*DTAU*DETA;
                    music_message.flush("warning");
                }
                if (std::abs(FULLSU[3]) > (DX*DY*DTAU+0.01)) {
                    music_message << "problem: volume in eta direction "
                                  << std::abs(FULLSU[3]) << "  > DX*DY*DTAU = "
                                  << DX*DY*DTAU;
                    music_message.flush("warning");
                }

                // position of the freeze-out fluid cell
                for (int ii = 0; ii < 4; ii++) {
                    x_fraction[1][ii] =
                        cornelius_ptr->get_centroid_elem(isurf, ii);
                    x_fraction[0][ii] =
                        lattice_spacing[ii] - x_fraction[1][ii];
                }
                const double tau_center = tau - DTAU + x_fraction[1][0];
                const double x_center = x + x_fraction[1][1];
                const double y_center = y + x_fraction[1][2];
                const double eta_center = eta + x_fraction[1][3];


                // perform 4-d linear interpolation for all fluid quantities
                for (int ii = 0; ii < 2; ii++)
                for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++) {
                    fluid_cube[0][ii][jj][kk] = arena_freezeout(
                            ix + ii*fac_x, iy + jj*fac_y, ieta + kk*fac_eta);
                    fluid_cube[1][ii][jj][kk] = arena_current(
                            ix + ii*fac_x, iy + jj*fac_y, ieta + kk*fac_eta);

                    if (DATA.output_vorticity == 0) continue;

                    // compute the vorticity tensors
                    Cell_aux aux_tmp;
                    double eta_local = eta + kk*DETA;
                    u_derivative_helper.compute_vorticity_shell(
                        tau, arena_prev, arena_current,
                        ieta + kk*fac_eta, ix + ii*fac_eta, iy + jj*fac_eta,
                        eta_local,
                        aux_tmp.omega_kSP, aux_tmp.omega_k,
                        aux_tmp.omega_th, aux_tmp.omega_T);
                    fluid_aux_cube[1][ii][jj][kk] = aux_tmp;
                    u_derivative_helper.compute_vorticity_shell(
                        tau - DTAU, arena_freezeout_prev, arena_freezeout,
                        ieta + kk*fac_eta, ix + ii*fac_eta, iy + jj*fac_eta,
                        eta_local,
                        aux_tmp.omega_kSP, aux_tmp.omega_k,
                        aux_tmp.omega_th, aux_tmp.omega_T);
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
                const double TFO = eos.get_temperature(epsFO,
                                                       fluid_center.rhob);
                if (TFO < 0) {
                    music_message << "TFO=" << TFO
                                  << "<0. ERROR. exiting.";
                    music_message.flush("error");
                    exit(1);
                }
                const double muB = eos.get_muB(epsFO, fluid_center.rhob);
                const double muS = eos.get_muS(epsFO, fluid_center.rhob);
                const double muC = eos.get_muC(epsFO, fluid_center.rhob);

                const double pressure = eos.get_pressure(epsFO, fluid_center.rhob);
                const double eps_plus_p_over_T_FO = (epsFO + pressure)/TFO;

                // finally output results !!!!
                if (surface_in_binary) {
                    const int FOsize = 34 + DATA.output_vorticity*24;
                    float array[FOsize];
                    array[0] = static_cast<float>(tau_center);
                    array[1] = static_cast<float>(x_center);
                    array[2] = static_cast<float>(y_center);
                    array[3] = static_cast<float>(eta_center);
                    for (int ii = 0; ii < 4; ii++)
                        array[4+ii] = static_cast<float>(FULLSU[ii]);
                    for (int ii = 0; ii < 4; ii++)
                        array[8+ii] = static_cast<float>(fluid_center.u[ii]);
                    array[12] = static_cast<float>(epsFO);
                    array[13] = static_cast<float>(TFO);
                    array[14] = static_cast<float>(muB);
                    array[15] = static_cast<float>(muS);
                    array[16] = static_cast<float>(muC);
                    array[17] = static_cast<float>(eps_plus_p_over_T_FO);
                    for (int ii = 0; ii < 10; ii++)
                        array[18+ii] = static_cast<float>(fluid_center.Wmunu[ii]);
                    array[28] = fluid_center.pi_b;
                    array[29] = fluid_center.rhob;
                    for (int ii = 0; ii < 4; ii++)
                        array[30+ii] = static_cast<float>(fluid_center.Wmunu[10+ii]);
                    if (DATA.output_vorticity == 1) {
                        for (int ii = 0; ii < 6; ii++) {
                            array[34+ii] = fluid_aux_center.omega_kSP[ii]/TFO;  // no minus sign because its definition is opposite to the kinetic vorticity
                            // the extra minus sign is from metric
                            // output quantities for g = (1, -1, -1, -1)
                            array[40+ii] = -fluid_aux_center.omega_k[ii]/TFO;
                            array[46+ii] = -fluid_aux_center.omega_th[ii];
                            array[52+ii] = (-fluid_aux_center.omega_T[ii]
                                            /TFO/TFO);
                        }
                    }
                    for (int i = 0; i < FOsize; i++)
                        s_file.write((char*) &(array[i]), sizeof(float));
                } else {
                    s_file << std::scientific << std::setprecision(10)
                           << tau_center << " " << x_center << " "
                           << y_center << " " << eta_center << " "
                           << FULLSU[0] << " " << FULLSU[1] << " "
                           << FULLSU[2] << " " << FULLSU[3] << " "
                           << fluid_center.u[0] << " " << fluid_center.u[1] << " "
                           << fluid_center.u[2] << " " << fluid_center.u[3] << " "
                           << epsFO << " " << TFO << " " << muB << " "
                           << muS << " " << muC << " "
                           << eps_plus_p_over_T_FO << " ";
                    for (int ii = 0; ii < 10; ii++)
                        s_file << std::scientific << std::setprecision(10)
                               << fluid_center.Wmunu[ii] << " ";
                    if (DATA.turn_on_bulk)
                        s_file << fluid_center.pi_b << " ";
                    if (DATA.turn_on_rhob)
                        s_file << fluid_center.rhob << " ";
                    if (DATA.turn_on_diff)
                        for (int ii = 10; ii < 14; ii++)
                            s_file << std::scientific << std::setprecision(10)
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
            for (int k = 0; k < 2; k++) {
                delete [] cube[i][j][k];
                delete [] fluid_cube[i][j][k];
                delete [] fluid_aux_cube[i][j][k];
            }
            delete [] cube[i][j];
            delete [] fluid_cube[i][j];
            delete [] fluid_aux_cube[i][j];
        }
        delete [] cube[i];
        delete [] fluid_cube[i];
        delete [] fluid_aux_cube[i];
    }
    delete [] cube;
    delete [] fluid_cube;
    delete [] fluid_aux_cube;
    return(intersections);
}


// Cornelius freeze out (C. Shen, 11/2014)
int Evolve::FreezeOut_equal_tau_Surface(double tau,
                                        SCGrid &arena_current) {
    // this function freeze-out fluid cells between epsFO and epsFO_low
    // on an equal time hyper-surface at the first time step
    // this function will be trigged if freezeout_lowtemp_flag == 1
    const int neta = arena_current.nEta();
    const int fac_eta = 1;

    for (int i_freezesurf = 0; i_freezesurf < n_freeze_surf; i_freezesurf++) {
        double epsFO = epsFO_list[i_freezesurf]/hbarc;
        if (!DATA.boost_invariant) {
            #pragma omp parallel for
            for (int ieta = 0; ieta < neta - fac_eta; ieta += fac_eta) {
                int thread_id = omp_get_thread_num();
                FreezeOut_equal_tau_Surface_XY(tau,  ieta, arena_current,
                                               thread_id, epsFO);
            }
        } else {
            FreezeOut_equal_tau_Surface_XY(tau, 0, arena_current, 0, epsFO);
        }
    }
    return(0);
}


void Evolve::FreezeOut_equal_tau_Surface_XY(double tau, int ieta,
                                            SCGrid &arena_current,
                                            int thread_id, double epsFO) {
    const bool surface_in_binary = DATA.freeze_surface_in_binary;
    double epsFO_low = 0.05/hbarc;        // 1/fm^4

    const int nx = arena_current.nX();
    const int ny = arena_current.nY();

    std::stringstream strs_name;
    if (!DATA.boost_invariant) {
        strs_name << "surface_eps_" << std::setprecision(4) << epsFO*hbarc
                  << "_" << thread_id << ".dat";
    } else {
        strs_name << "surface_eps_" << std::setprecision(4) << epsFO*hbarc
                  << ".dat";
    }
    std::ofstream s_file;
    std::ios_base::openmode modes;

    if (surface_in_binary) {
        modes=std::ios::out | std::ios::binary;
    } else {
        modes=std::ios::out;
    }

    // Only append at the end of the file if it's not the first timestep
    // (that is, overwrite file at first timestep)
    if (tau != DATA.tau0+DATA.delta_tau) {
            modes = modes | std::ios::app;
    }

    s_file.open(strs_name.str().c_str(), modes);

    const int fac_x   = DATA.fac_x;
    const int fac_y   = DATA.fac_y;
    const int fac_eta = 1;

    const double DX   = fac_x*DATA.delta_x;
    const double DY   = fac_y*DATA.delta_y;
    const double DETA = fac_eta*DATA.delta_eta;

    double eta = (DATA.delta_eta)*ieta - (DATA.eta_size)/2.0;
    for (int ix = 0; ix < nx - fac_x; ix += fac_x) {
        double x = ix*(DATA.delta_x) - (DATA.x_size/2.0); 
        for (int iy = 0; iy < ny - fac_y; iy += fac_y) {
            double y = iy*(DATA.delta_y) - (DATA.y_size/2.0);

            // judge intersection
            if (arena_current(ix,iy,ieta).epsilon > epsFO) continue;
            if (arena_current(ix,iy,ieta).epsilon < epsFO_low) continue;

            // surface normal vector d^3 \sigma_\mu
            const double FULLSU[] = {DX*DY*DETA, 0.0, 0.0, 0.0};

            // get positions of the freeze-out surface
            const double tau_center = tau;
            const double x_center   = x;
            const double y_center   = y;
            const double eta_center = eta;

            // flow velocity
            const double ux_center   = arena_current(ix, iy, ieta).u[1];
            const double uy_center   = arena_current(ix, iy, ieta).u[2];
            const double ueta_center = arena_current(ix, iy, ieta).u[3];  // u^eta/tau
            // reconstruct u^tau from u^i
            const double utau_center = sqrt(1. + ux_center*ux_center 
                                               + uy_center*uy_center 
                                               + ueta_center*ueta_center);

            // baryon density rho_b
            const double rhob_center = arena_current(ix, iy, ieta).rhob;

            // baryon diffusion current
            double qtau_center = arena_current(ix, iy, ieta).Wmunu[10];
            double qx_center   = arena_current(ix, iy, ieta).Wmunu[11];
            double qy_center   = arena_current(ix, iy, ieta).Wmunu[12];
            double qeta_center = arena_current(ix, iy, ieta).Wmunu[13];

            // reconstruct q^\tau from the transverality criteria
            FlowVec u_flow = {utau_center, ux_center, uy_center, ueta_center};
            double q_mu[4] = {qtau_center, qx_center, qy_center, qeta_center};
            double q_regulated[4] = {0.0, 0.0, 0.0, 0.0};
            regulate_qmu(u_flow, q_mu, q_regulated);
            qtau_center = q_regulated[0];
            qx_center   = q_regulated[1];
            qy_center   = q_regulated[2];
            qeta_center = q_regulated[3];

            // bulk viscous pressure pi_b
            const double pi_b_center = arena_current(ix,iy,ieta).pi_b;

            // shear viscous tensor
            double Wtautau_center = arena_current(ix, iy, ieta).Wmunu[0];
            double Wtaux_center   = arena_current(ix, iy, ieta).Wmunu[1];
            double Wtauy_center   = arena_current(ix, iy, ieta).Wmunu[2];
            double Wtaueta_center = arena_current(ix, iy, ieta).Wmunu[3];
            double Wxx_center     = arena_current(ix, iy, ieta).Wmunu[4];
            double Wxy_center     = arena_current(ix, iy, ieta).Wmunu[5];
            double Wxeta_center   = arena_current(ix, iy, ieta).Wmunu[6];
            double Wyy_center     = arena_current(ix, iy, ieta).Wmunu[7];
            double Wyeta_center   = arena_current(ix, iy, ieta).Wmunu[8];
            double Wetaeta_center = arena_current(ix, iy, ieta).Wmunu[9];
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
            Wtaux_center   = Wmunu_regulated[0][1];
            Wtauy_center   = Wmunu_regulated[0][2];
            Wtaueta_center = Wmunu_regulated[0][3];
            Wxx_center     = Wmunu_regulated[1][1];
            Wxy_center     = Wmunu_regulated[1][2];
            Wxeta_center   = Wmunu_regulated[1][3];
            Wyy_center     = Wmunu_regulated[2][2];
            Wyeta_center   = Wmunu_regulated[2][3];
            Wetaeta_center = Wmunu_regulated[3][3];

            Cell_aux fluid_aux_center;

            // get other thermodynamical quantities
            double e_local   = arena_current(ix, iy, ieta).epsilon;
            double T_local   = eos.get_temperature(e_local, rhob_center);
            if (T_local < 0) {
                music_message << "Evolve::FreezeOut_equal_tau_Surface: "
                              << "T_local = " << T_local
                              << " <0. ERROR. exiting.";
                music_message.flush("error");
                exit(1);
            }
            double muB_local = eos.get_muB(e_local, rhob_center);
            double muS_local = eos.get_muS(e_local, rhob_center);
            double muC_local = eos.get_muC(e_local, rhob_center);

            double pressure = eos.get_pressure(e_local, rhob_center);
            double eps_plus_p_over_T = (e_local + pressure)/T_local;

            // finally output results
            if (surface_in_binary) {
                const int FOsize = 34 + DATA.output_vorticity*24;
                float array[FOsize];
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
                array[16] = static_cast<float>(muC_local);
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
                        array[34+ii] = fluid_aux_center.omega_kSP[ii]/T_local;  // no minus sign because its definition is opposite to the kinetic vorticity
                        // the extra minus sign is from metric
                        // output quantities for g = (1, -1, -1, -1)
                        array[40+ii] = -fluid_aux_center.omega_k[ii]/T_local;
                        array[46+ii] = -fluid_aux_center.omega_th[ii];
                        array[52+ii] = (-fluid_aux_center.omega_T[ii]
                                        /T_local/T_local);
                    }
                }
                for (int i = 0; i < FOsize; i++) {
                    s_file.write((char*) &(array[i]), sizeof(float));
                }
            } else {
                s_file << std::scientific << std::setprecision(10) 
                       << tau_center     << " " << x_center          << " "
                       << y_center       << " " << eta_center        << " "
                       << FULLSU[0]      << " " << FULLSU[1]         << " "
                       << FULLSU[2]      << " " << FULLSU[3]         << " "
                       << utau_center    << " " << ux_center         << " "
                       << uy_center      << " " << ueta_center       << " "
                       << e_local        << " " << T_local           << " "
                       << muB_local      << " " << muS_local         << " "
                       << muC_local      << " " << eps_plus_p_over_T << " "
                       << Wtautau_center << " " << Wtaux_center      << " "
                       << Wtauy_center   << " " << Wtaueta_center    << " "
                       << Wxx_center     << " " << Wxy_center        << " "
                       << Wxeta_center   << " " << Wyy_center        << " "
                       << Wyeta_center   << " " << Wetaeta_center    << " " ;
                if (DATA.turn_on_bulk)
                    s_file << pi_b_center << " " ;
                if (DATA.turn_on_rhob)
                    s_file << rhob_center << " " ;
                if (DATA.turn_on_diff)
                    s_file << qtau_center << " " << qx_center << " "
                           << qy_center << " " << qeta_center << " " ;
                s_file << std::endl;
            }
        }
    }
    s_file.close();
}


int Evolve::FindFreezeOutSurface_boostinvariant_Cornelius(
                double tau, SCGrid &arena_current, SCGrid &arena_freezeout) {
    const bool surface_in_binary = DATA.freeze_surface_in_binary;
    // find boost-invariant hyper-surfaces
    int *all_frozen = new int[n_freeze_surf];
    for (int i_freezesurf = 0; i_freezesurf < n_freeze_surf; i_freezesurf++) {
        double epsFO = epsFO_list[i_freezesurf]/hbarc;

        std::stringstream strs_name;
        strs_name << "surface_eps_" << std::setprecision(4) << epsFO*hbarc
                  << ".dat";

        std::ofstream s_file;
        std::ios_base::openmode modes;

        if (surface_in_binary) {
            modes=std::ios::out | std::ios::binary;
        } else {
            modes=std::ios::out;
        }

        // Only append at the end of the file if it's not the first timestep
        // (that is, overwrite file at first timestep)
        if (tau != DATA.tau0+DATA.delta_tau) {
                modes = modes | std::ios::app;
        }

        s_file.open(strs_name.str().c_str(), modes);

        const int nx = arena_current.nX();
        const int ny = arena_current.nY();
        double FULLSU[4];  // d^3 \sigma_\mu

        int intersect;
        int intersections = 0;

        facTau    = DATA.facTau;   // step to skip in tau direction
        int fac_x = DATA.fac_x;
        int fac_y = DATA.fac_y;

        const double DX   = fac_x*DATA.delta_x;
        const double DY   = fac_y*DATA.delta_y;
        const double DETA = 1.0;
        const double DTAU = facTau*DATA.delta_tau;

        double lattice_spacing[3] = {DTAU, DX, DY};
        double x_fraction[2][3];

        // initialize Cornelius
        const int dim = 3;
        std::shared_ptr<Cornelius> cornelius_ptr(new Cornelius());
        cornelius_ptr->init(dim, epsFO, lattice_spacing);

        // initialize the hyper-cube for Cornelius
        Cell_small ***fluid_cube = new Cell_small ** [2];
        double ***cube = new double ** [2];
        for (int i = 0; i < 2; i++) {
            fluid_cube[i] = new Cell_small * [2];
            cube[i] = new double * [2];
            for (int j = 0; j < 2; j++) {
                fluid_cube[i][j] = new Cell_small[2];
                cube[i][j] = new double[2];
                for (int k = 0; k < 2; k++)
                    cube[i][j][k] = 0.0;
            }
        }

        for (int ix=0; ix < nx - fac_x; ix += fac_x) {
            double x = ix*(DATA.delta_x) - (DATA.x_size/2.0);
            for (int iy=0; iy < ny - fac_y; iy += fac_y) {
                double y = iy*(DATA.delta_y) - (DATA.y_size/2.0);

                // judge intersection (from Bjoern)
                intersect=1;
                if ((arena_current(ix+fac_x,iy+fac_y,0).epsilon-epsFO)
                    *(arena_freezeout(ix,iy,0).epsilon-epsFO) > 0.)
                    if ((arena_current(ix+fac_x,iy,0).epsilon-epsFO)
                        *(arena_freezeout(ix,iy+fac_y,0).epsilon-epsFO) > 0.)
                        if ((arena_current(ix,iy+fac_y,0).epsilon-epsFO)
                            *(arena_freezeout(ix+fac_x,iy,0).epsilon-epsFO) > 0.)
                            if ((arena_current(ix,iy,0).epsilon-epsFO)
                                *(arena_freezeout(ix+fac_x,iy+fac_y,0).epsilon-epsFO) > 0.)
                                    intersect = 0;
                if (intersect == 0) continue;

                if (ix == 0 || ix >= nx - 2*fac_x
                        || iy == 0 || iy >= ny - 2*fac_y) {
                    music_message << "Freeze-out cell at the boundary! "
                                  << "The grid is too small!";
                    music_message.flush("error");
                    exit(1);
                }

                // if intersect, prepare for the hyper-cube
                intersections++;
                cube[0][0][0] = arena_freezeout(ix      , iy      , 0).epsilon;
                cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).epsilon;
                cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).epsilon;
                cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).epsilon;
                cube[1][0][0] = arena_current  (ix      , iy      , 0).epsilon;
                cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).epsilon;
                cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).epsilon;
                cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).epsilon;

                // Now, the magic will happen in the Cornelius ...
                cornelius_ptr->find_surface_3d(cube);

                // get positions of the freeze-out surface 
                // and interpolating results
                for (int isurf = 0; isurf < cornelius_ptr->get_Nelements(); 
                     isurf++) {
                    // surface normal vector d^3 \sigma_\mu
                    for (int ii = 0; ii < dim; ii++)
                        FULLSU[ii] = cornelius_ptr->get_normal_elem(isurf, ii);

                    FULLSU[3] = 0.0; // rapidity direction is set to 0

                    // check the size of the surface normal vector
                    if (fabs(FULLSU[0]) > (DX*DY*DETA + 0.01)) {
                       music_message << "problem: volume in tau direction "
                                     << fabs(FULLSU[0]) << "  > DX*DY*DETA = "
                                     << DX*DY*DETA;
                        music_message.flush("warning");
                    }
                    if (fabs(FULLSU[1]) > (DTAU*DY*DETA + 0.01)) {
                        music_message << "problem: volume in x direction "
                                      << fabs(FULLSU[1])
                                      << "  > DTAU*DY*DETA = " << DTAU*DY*DETA;
                        music_message.flush("warning");
                    }
                    if (fabs(FULLSU[2]) > (DX*DTAU*DETA+0.01)) {
                        music_message << "problem: volume in y direction "
                                      << fabs(FULLSU[2])
                                      << "  > DX*DTAU*DETA = " << DX*DTAU*DETA;
                        music_message.flush("warning");
                    }

                    // position of the freeze-out fluid cell
                    for (int ii = 0; ii < dim; ii++) {
                        x_fraction[1][ii] = (
                            cornelius_ptr->get_centroid_elem(isurf, ii));
                        x_fraction[0][ii] = (
                            lattice_spacing[ii] - x_fraction[1][ii]);
                    }
                    const double tau_center = tau - DTAU + x_fraction[1][0];
                    const double x_center = x + x_fraction[1][1];
                    const double y_center = y + x_fraction[1][2];
                    const double eta_center = 0.0;

                    // perform 3-d linear interpolation for all fluid quantities
                    fluid_cube[0][0][0] = arena_freezeout(ix      , iy      , 0);
                    fluid_cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0);
                    fluid_cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0);
                    fluid_cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0);
                    fluid_cube[1][0][0] = arena_current  (ix      , iy      , 0);
                    fluid_cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0);
                    fluid_cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0);
                    fluid_cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0);
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

                    // 3-dimension interpolation done
                    double TFO = eos.get_temperature(epsFO, fluid_center.rhob);
                    double muB = eos.get_muB(epsFO, fluid_center.rhob);
                    double muS = eos.get_muS(epsFO, fluid_center.rhob);
                    double muC = eos.get_muC(epsFO, fluid_center.rhob);
                    if (TFO < 0) {
                        music_message << "TFO=" << TFO
                                      << "<0. ERROR. exiting.";
                        music_message.flush("error");
                        exit(1);
                    }

                    double pressure = eos.get_pressure(epsFO, fluid_center.rhob);
                    double eps_plus_p_over_T_FO = (epsFO + pressure)/TFO;

                    // finally output results !!!!
                    if (surface_in_binary) {
                        float array[34];
                        array[0] = static_cast<float>(tau_center);
                        array[1] = static_cast<float>(x_center);
                        array[2] = static_cast<float>(y_center);
                        array[3] = static_cast<float>(eta_center);
                        for (int ii = 0; ii < 4; ii++)
                            array[4+ii] = static_cast<float>(FULLSU[ii]);
                        for (int ii = 0; ii < 4; ii++)
                            array[8+ii] = static_cast<float>(fluid_center.u[ii]);
                        array[12] = static_cast<float>(epsFO);
                        array[13] = static_cast<float>(TFO);
                        array[14] = static_cast<float>(muB);
                        array[15] = static_cast<float>(muS);
                        array[16] = static_cast<float>(muC);
                        array[17] = static_cast<float>(eps_plus_p_over_T_FO);
                        for (int ii = 0; ii < 10; ii++)
                            array[18+ii] = static_cast<float>(fluid_center.Wmunu[ii]);
                        array[28] = fluid_center.pi_b;
                        array[29] = fluid_center.rhob;
                        for (int ii = 0; ii < 4; ii++)
                            array[30+ii] = static_cast<float>(fluid_center.Wmunu[10+ii]);
                        for (int i = 0; i < 34; i++)
                            s_file.write((char*) &(array[i]), sizeof(float));
                    } else {
                        s_file << std::scientific << std::setprecision(10)
                               << tau_center << " " << x_center << " "
                               << y_center << " " << eta_center << " "
                               << FULLSU[0] << " " << FULLSU[1] << " "
                               << FULLSU[2] << " " << FULLSU[3] << " "
                               << fluid_center.u[0] << " " << fluid_center.u[1] << " "
                               << fluid_center.u[2] << " " << fluid_center.u[3] << " "
                               << epsFO << " " << TFO << " " << muB << " "
                               << muS << " " << muC << " "
                               << eps_plus_p_over_T_FO << " ";
                        for (int ii = 0; ii < 10; ii++)
                            s_file << std::scientific << std::setprecision(10)
                                   << fluid_center.Wmunu[ii] << " ";
                        if (DATA.turn_on_bulk)
                            s_file << fluid_center.pi_b << " ";
                        if (DATA.turn_on_rhob)
                            s_file << fluid_center.rhob << " ";
                        if (DATA.turn_on_diff)
                            for (int ii = 10; ii < 14; ii++)
                                s_file << std::scientific << std::setprecision(10)
                                       << fluid_center.Wmunu[ii] << " ";
                    }
                }
            }
        }

        s_file.close();

        // clean up
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                delete [] cube[i][j];
                delete [] fluid_cube[i][j];
            }
            delete [] cube[i];
            delete [] fluid_cube[i];
        }
        delete [] cube;
        delete [] fluid_cube;

        // judge whether the entire fireball is freeze-out
        all_frozen[i_freezesurf] = 0;
        if (intersections == 0)
            all_frozen[i_freezesurf] = 1;
    }

    int all_frozen_flag = 1;
    for (int ii = 0; ii < n_freeze_surf; ii++) {
        all_frozen_flag *= all_frozen[ii];
    }
    if (all_frozen_flag == 1) {
        music_message.info("All cells frozen out. Exiting.");
    }

    delete [] all_frozen;
    return(all_frozen_flag);
}

void Evolve::regulate_qmu(const FlowVec u, const double q[],
                          double q_regulated[]) const {
    double u_dot_q = - u[0]*q[0] + u[1]*q[1] + u[2]*q[2] + u[3]*q[3];
    for (int i = 0; i < 4; i++) {
        q_regulated[i] = q[i] + u[i]*u_dot_q;
    }
}

void Evolve::regulate_Wmunu(const FlowVec u, const double Wmunu[4][4],
                            double Wmunu_regulated[4][4]) const {
    const double gmunu[4][4] = {{-1, 0, 0, 0},
                                { 0, 1, 0, 0},
                                { 0, 0, 1, 0},
                                { 0, 0, 0, 1}};
    double u_dot_pi[4];
    double u_mu[4];
    for (int i = 0; i < 4; i++) {
        u_dot_pi[i] = (- u[0]*Wmunu[0][i] + u[1]*Wmunu[1][i]
                       + u[2]*Wmunu[2][i] + u[3]*Wmunu[3][i]);
        u_mu[i] = gmunu[i][i]*u[i];
    }
    double tr_pi = - Wmunu[0][0] + Wmunu[1][1] + Wmunu[2][2] + Wmunu[3][3];
    double u_dot_pi_dot_u = 0.0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            u_dot_pi_dot_u += u_mu[i]*Wmunu[i][j]*u_mu[j];
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Wmunu_regulated[i][j] = (
                Wmunu[i][j] + u[i]*u_dot_pi[j] + u[j]*u_dot_pi[i]
                + u[i]*u[j]*u_dot_pi_dot_u
                - 1./3.*(gmunu[i][j] + u[i]*u[j])*(tr_pi + u_dot_pi_dot_u));
        }
    }
}

void Evolve::initialize_freezeout_surface_info() {
    if (DATA.useEpsFO == 0) {
        const double e_freeze = eos.get_T2e(DATA.TFO, 0.0)*Util::hbarc;
        n_freeze_surf = 1;
        for (int isurf = 0; isurf < n_freeze_surf; isurf++) {
            epsFO_list.push_back(e_freeze);
            music_message << "Freeze out at a constant temperature T = " 
                          << DATA.TFO << " GeV, e_fo = "
                          << e_freeze << " GeV/fm^3";
            music_message.flush("info");
        }
    }
    const int freeze_eps_flag = DATA.freeze_eps_flag;
    if (freeze_eps_flag == 0) {
        // constant spacing the energy density
        n_freeze_surf = DATA.N_freeze_out;
        double freeze_max_ed = DATA.eps_freeze_max;
        double freeze_min_ed = DATA.eps_freeze_min;
        double d_epsFO = ((freeze_max_ed - freeze_min_ed)
                          /(n_freeze_surf - 1 + Util::small_eps));
        for (int isurf = 0; isurf < n_freeze_surf; isurf++) {
            double temp_epsFO = freeze_min_ed + isurf*d_epsFO;
            epsFO_list.push_back(temp_epsFO);
        }
    } else if(freeze_eps_flag == 1) {
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
        int temp_n_surf = 0;
        std::string dummy;
        double temp_epsFO, dummyd;
        getline(freeze_list_file, dummy);  // get rid of the comment
        while(1) {
            freeze_list_file >> temp_epsFO >> dummyd >> dummyd 
                             >> dummyd >> dummyd >> dummyd >> dummyd;  
            if (!freeze_list_file.eof()) {
                epsFO_list.push_back(temp_epsFO);
                temp_n_surf++;
            } else {
                break;
            }
        }
        freeze_list_file.close();
        n_freeze_surf = temp_n_surf;
        music_message << "totally " << n_freeze_surf 
                      << " freeze-out surface will be generated ...";
        music_message.flush("info");
    } else {
        music_message << "Evolve::initialize_freezeout_surface_info: "
                      << "unrecoginze freeze_eps_flag = " << freeze_eps_flag;
        music_message.flush("error");
        exit(1);
    }
}


Cell_small Evolve::three_dimension_linear_interpolation(
        double* lattice_spacing, double fraction[2][3], Cell_small*** cube) {
    double denorm = 1.0;
    for (int i = 0; i < 3; i++) {
        denorm *= lattice_spacing[i];
    }

    Cell_small results;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                results = results + (cube[i][j][k]*fraction[i][0]
                                     *fraction[j][1]*fraction[k][2]);
            }
        }
    }
    results = results*(1./denorm);
    return(results);
}


Cell_small Evolve::four_dimension_linear_interpolation(
        double* lattice_spacing, double fraction[2][4], Cell_small**** cube) {
    double denorm = 1.0;
    Cell_small results;
    for (int i = 0; i < 4; i++) {
        denorm *= lattice_spacing[i];
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                for (int l = 0; l < 2; l++) {
                    results = results + (
                        cube[i][j][k][l]*fraction[i][0]*fraction[j][1]
                        *fraction[k][2]*fraction[l][3]);
                }
            }
        }
    }
    results = results*(1./denorm);
    return (results);
}


Cell_aux Evolve::three_dimension_linear_interpolation(
        double* lattice_spacing, double fraction[2][3], Cell_aux*** cube) {
    double denorm = 1.0;
    for (int i = 0; i < 3; i++) {
        denorm *= lattice_spacing[i];
    }

    Cell_aux results;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                results = results + (cube[i][j][k]*fraction[i][0]
                                     *fraction[j][1]*fraction[k][2]);
            }
        }
    }
    results = results*(1./denorm);
    return(results);
}


Cell_aux Evolve::four_dimension_linear_interpolation(
        double* lattice_spacing, double fraction[2][4], Cell_aux**** cube) {
    double denorm = 1.0;
    Cell_aux results;
    for (int i = 0; i < 4; i++) {
        denorm *= lattice_spacing[i];
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                for (int l = 0; l < 2; l++) {
                    results = results + (
                        cube[i][j][k][l]*fraction[i][0]*fraction[j][1]
                        *fraction[k][2]*fraction[l][3]);
                }
            }
        }
    }
    results = results*(1./denorm);
    return (results);
}
