// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <omp.h>
#include <algorithm>

#include "./evolve.h"
#include "./util.h"
#include "./data.h"
#include "./cell.h"
#include "./grid.h"
#include "./eos.h"
#include "./advance.h"
#include "./cornelius.h"

#ifndef _OPENMP
  #define omp_get_thread_num() 0
#endif

using namespace std;

Evolve::Evolve(const EOS &eosIn, const InitData &DATA_in,
               hydro_source &hydro_source_in) :
    eos(eosIn), DATA(DATA_in), hydro_source_terms(hydro_source_in),
    grid_info(DATA_in, eosIn), advance(eosIn, DATA_in, hydro_source_in),
    u_derivative(DATA_in, eosIn) {
    rk_order  = DATA_in.rk_order;
    if (DATA.freezeOutMethod == 4) {
        initialize_freezeout_surface_info();
    }
}

// master control function for hydrodynamic evolution
int Evolve::EvolveIt(SCGrid &arena_prev, SCGrid &arena_current,
                     SCGrid &arena_future) {
    // first pass some control parameters
    facTau                      = DATA.facTau;
    int Nskip_timestep          = DATA.output_evolution_every_N_timesteps;
    int freezeout_flag          = DATA.doFreezeOut;
    int freezeout_lowtemp_flag  = DATA.doFreezeOut_lowtemp;

    // Output information about the hydro parameters 
    // in the format of a C header file
    if (DATA.output_hydro_params_header || DATA.outputEvolutionData == 1)
        grid_info.Output_hydro_information_header();

    // main loop starts ...
    int    itmax = DATA.nt;
    double tau0  = DATA.tau0;
    double dt    = DATA.delta_tau;

    double tau;
    int it_start = 0;
    if (DATA.Initial_profile == 13) {
        double source_tau_max = hydro_source_terms.get_source_tau_max();
        it_start = static_cast<int>((source_tau_max - tau0)/dt);
        if (it_start < 0) it_start = 0;
    } else if (DATA.Initial_profile == 30) {
        double source_tau_max = hydro_source_terms.get_source_tau_max();
        it_start = static_cast<int>((source_tau_max - tau0)/dt);
        if (it_start < 0) it_start = 0;
    }


    const auto closer = [](SCGrid* g) { /*Don't delete memory we don't own*/ };
    GridPointer ap_prev   (&arena_prev, closer);
    GridPointer ap_current(&arena_current, closer);
    GridPointer ap_future (&arena_future, closer);

    SCGrid arena_freezeout(arena_current.nX(),
                           arena_current.nY(),
                           arena_current.nEta());

    for (int it = 0; it <= itmax; it++) {
        tau = tau0 + dt*it;

        if (DATA->Initial_profile == 30) {
            hydro_source_ptr->prepare_list_for_current_tau_frame(tau);
        }
        // store initial conditions
        if (it == it_start) {
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
        
        if (DATA.Initial_profile == 13) {
            if (it == it_start) {
                grid_info.output_energy_density_and_rhob_disitrubtion(
                            *ap_current,
                            "energy_density_and_rhob_from_source_terms.dat");
            }
        }

        if (it % Nskip_timestep == 0) {
            if (DATA.outputEvolutionData == 1) {
                grid_info.OutputEvolutionDataXYEta(*ap_current, tau);
            } else if (DATA.outputEvolutionData == 2) {
                grid_info.OutputEvolutionDataXYEta_chun(*ap_current, tau);
            }
            if (DATA.output_movie_flag == 1) {
                grid_info.output_evolution_for_movie(*ap_current, tau);
            }
        }

        if (DATA.Initial_profile == 13) {
            grid_info.output_average_phase_diagram_trajectory(
                                            tau, -0.5, 0.5, *ap_current);
        }

        // check energy conservation
        if (DATA.boost_invariant == 0)
            grid_info.check_conservation_law(*ap_current, *ap_prev, tau);
        grid_info.get_maximum_energy_density(*ap_current);

        if (DATA.output_hydro_debug_info == 1) {
            grid_info.monitor_fluid_cell(*ap_current, 100, 100, 0, tau);
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
                if (DATA.boost_invariant == 0) {
                    frozen = FindFreezeOutSurface_Cornelius(
                                tau, *ap_current, arena_freezeout);
                } else {
                    frozen = FindFreezeOutSurface_boostinvariant_Cornelius(
                                tau, *ap_current, arena_freezeout);
                }
                store_previous_step_for_freezeout(*ap_current,
                                                  arena_freezeout);
            }
        }
        music_message << "Done time step " << it << "/" << itmax
                      << " tau = " << tau << " fm/c";
        music_message.flush("info");
        if (frozen) break;
    }
    music_message.info("Finished.");
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
                                           SCGrid &arena_current,
                                           SCGrid &arena_freezeout) {
    const int neta = arena_current.nEta();
    const int fac_eta = 1;
    int flag_all_frozen = 0;
    for (int i_freezesurf = 0; i_freezesurf < n_freeze_surf; i_freezesurf++) {
        const double epsFO = epsFO_list[i_freezesurf]/hbarc;   // 1/fm^4

        int intersections = 0;
        #pragma omp parallel for reduction(+:intersections)
        for (int ieta = 0; ieta < (neta-fac_eta); ieta += fac_eta) {
            int thread_id = omp_get_thread_num();
            intersections += FindFreezeOutSurface_Cornelius_XY(
                tau, ieta, arena_current, arena_freezeout, thread_id, epsFO);
        }
        if (intersections == 0)
            flag_all_frozen = 1;
    }

    if (flag_all_frozen == 1) {
        std::cout << "All cells frozen out. Exiting." << std::endl;
    }
    return(flag_all_frozen);
}

int Evolve::FindFreezeOutSurface_Cornelius_XY(double tau, int ieta,
                                              SCGrid &arena_current,
                                              SCGrid &arena_freezeout,
                                              int thread_id, double epsFO) {
    const int nx = arena_current.nX();
    const int ny = arena_current.nY();

    stringstream strs_name;
    strs_name << "surface_eps_" << setprecision(4) << epsFO*hbarc
              << "_" << thread_id << ".dat";
    ofstream s_file;
    s_file.open(strs_name.str().c_str(), ios::out | ios::app);
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

    // initialize Cornelius
    double lattice_spacing[4] = {DTAU, DX, DY, DETA};
    Cornelius* cornelius_ptr = new Cornelius();
    cornelius_ptr->init(dim, epsFO, lattice_spacing);

    // initialize the hyper-cube for Cornelius
    double ****cube = new double*** [2];
    for (int i = 0; i < 2; i++) {
        cube[i] = new double** [2];
        for (int j = 0; j < 2; j++) {
            cube[i][j] = new double* [2];
            for (int k = 0; k < 2; k++) {
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

                // perform 4-d linear interpolation for all fluid
                // quantities

                // flow velocity u^x
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).u[1];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).u[1];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).u[1];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).u[1];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).u[1];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).u[1];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).u[1];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).u[1];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).u[1];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).u[1];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).u[1];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).u[1];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).u[1];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).u[1];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).u[1];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).u[1];
                const double ux_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);

                // flow velocity u^y
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).u[2];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).u[2];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).u[2];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).u[2];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).u[2];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).u[2];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).u[2];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).u[2];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).u[2];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).u[2];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).u[2];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).u[2];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).u[2];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).u[2];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).u[2];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).u[2];
                const double uy_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);

                // flow velocity u^eta
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).u[3];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).u[3];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).u[3];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).u[3];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).u[3];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).u[3];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).u[3];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).u[3];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).u[3];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).u[3];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).u[3];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).u[3];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).u[3];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).u[3];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).u[3];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).u[3];
                const double ueta_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);

                // reconstruct u^tau from u^i
                const double utau_center = sqrt(1. + ux_center*ux_center 
                                   + uy_center*uy_center 
                                   + ueta_center*ueta_center);

                // baryon density rho_b
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).rhob;
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).rhob;
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).rhob;
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).rhob;
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).rhob;
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).rhob;
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).rhob;
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).rhob;
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).rhob;
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).rhob;
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).rhob;
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).rhob;
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).rhob;
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).rhob;
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).rhob;
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).rhob;
                const double rhob_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);
          
                // baryon diffusion current q^tau
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[10];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[10];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[10];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[10];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[10];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[10];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[10];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[10];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[10];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[10];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[10];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[10];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[10];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[10];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[10];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[10];
                double qtau_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);
          
                // baryon diffusion current q^x
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[11];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[11];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[11];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[11];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[11];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[11];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[11];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[11];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[11];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[11];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[11];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[11];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[11];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[11];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[11];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[11];
                double qx_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);

                // baryon diffusion current q^y
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[12];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[12];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[12];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[12];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[12];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[12];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[12];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[12];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[12];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[12];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[12];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[12];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[12];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[12];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[12];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[12];
                double qy_center = 
                    Util::four_dimension_linear_interpolation(
                            lattice_spacing, x_fraction, cube);
          
                // baryon diffusion current q^eta
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[13];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[13];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[13];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[13];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[13];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[13];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[13];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[13];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[13];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[13];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[13];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[13];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[13];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[13];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[13];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[13];
                double qeta_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);

                // reconstruct q^\tau from the transverality criteria
                double u_flow[4] = {utau_center, ux_center, uy_center, ueta_center};
                double q_mu[4]   = {qtau_center, qx_center, qy_center, qeta_center};
                double q_regulated[4] = {0.0, 0.0, 0.0, 0.0};
                
                regulate_qmu(u_flow, q_mu, q_regulated);
                
                qtau_center = q_regulated[0];
                qx_center = q_regulated[1];
                qy_center = q_regulated[2];
                qeta_center = q_regulated[3];
    
                // bulk viscous pressure pi_b
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).pi_b;
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).pi_b;
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).pi_b;
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).pi_b;
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).pi_b;
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).pi_b;
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).pi_b;
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).pi_b;
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).pi_b;
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).pi_b;
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).pi_b;
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).pi_b;
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).pi_b;
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).pi_b;
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).pi_b;
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).pi_b;
                const double pi_b_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);

                // shear viscous tensor W^\tau\tau
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[0];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[0];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[0];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[0];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[0];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[0];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[0];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[0];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[0];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[0];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[0];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[0];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[0];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[0];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[0];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[0];
                double Wtautau_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);
      
                // shear viscous tensor W^{\tau x}
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[1];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[1];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[1];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[1];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[1];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[1];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[1];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[1];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[1];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[1];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[1];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[1];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[1];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[1];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[1];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[1];
                double Wtaux_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);

                // shear viscous tensor W^{\tau y}
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[2];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[2];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[2];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[2];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[2];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[2];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[2];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[2];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[2];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[2];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[2];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[2];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[2];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[2];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[2];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[2];
                double Wtauy_center = Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);
      
                // shear viscous tensor W^{\tau \eta}
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[3];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[3];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[3];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[3];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[3];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[3];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[3];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[3];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[3];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[3];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[3];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[3];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[3];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[3];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[3];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[3];
                double Wtaueta_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);
      
                // shear viscous tensor W^{xx}
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[4];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[4];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[4];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[4];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[4];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[4];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[4];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[4];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[4];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[4];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[4];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[4];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[4];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[4];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[4];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[4];
                double Wxx_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);

                // shear viscous tensor W^{xy}
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[5];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[5];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[5];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[5];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[5];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[5];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[5];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[5];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[5];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[5];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[5];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[5];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[5];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[5];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[5];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[5];
                double Wxy_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);

                // shear viscous tensor W^{x\eta}
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[6];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[6];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[6];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[6];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[6];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[6];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[6];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[6];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[6];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[6];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[6];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[6];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[6];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[6];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[6];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[6];
                double Wxeta_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);
      
                // shear viscous tensor W^{yy}
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[7];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[7];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[7];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[7];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[7];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[7];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[7];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[7];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[7];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[7];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[7];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[7];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[7];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[7];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[7];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[7];
                double Wyy_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);
      
                // shear viscous tensor W^{y\eta}
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[8];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[8];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[8];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[8];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[8];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[8];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[8];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[8];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[8];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[8];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[8];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[8];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[8];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[8];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[8];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[8];
                double Wyeta_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);
      
                // shear viscous tensor W^{\eta\eta}
                cube[0][0][0][0] = arena_freezeout(ix      , iy      , ieta        ).Wmunu[9];
                cube[0][0][1][0] = arena_freezeout(ix      , iy+fac_y, ieta        ).Wmunu[9];
                cube[0][1][0][0] = arena_freezeout(ix+fac_x, iy      , ieta        ).Wmunu[9];
                cube[0][1][1][0] = arena_freezeout(ix+fac_x, iy+fac_y, ieta        ).Wmunu[9];
                cube[1][0][0][0] = arena_current  (ix      , iy      , ieta        ).Wmunu[9];
                cube[1][0][1][0] = arena_current  (ix      , iy+fac_y, ieta        ).Wmunu[9];
                cube[1][1][0][0] = arena_current  (ix+fac_x, iy      , ieta        ).Wmunu[9];
                cube[1][1][1][0] = arena_current  (ix+fac_x, iy+fac_y, ieta        ).Wmunu[9];
                cube[0][0][0][1] = arena_freezeout(ix      , iy      , ieta+fac_eta).Wmunu[9];
                cube[0][0][1][1] = arena_freezeout(ix      , iy+fac_y, ieta+fac_eta).Wmunu[9];
                cube[0][1][0][1] = arena_freezeout(ix+fac_x, iy      , ieta+fac_eta).Wmunu[9];
                cube[0][1][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[9];
                cube[1][0][0][1] = arena_current  (ix      , iy      , ieta+fac_eta).Wmunu[9];
                cube[1][0][1][1] = arena_current  (ix      , iy+fac_y, ieta+fac_eta).Wmunu[9];
                cube[1][1][0][1] = arena_current  (ix+fac_x, iy      , ieta+fac_eta).Wmunu[9];
                cube[1][1][1][1] = arena_current  (ix+fac_x, iy+fac_y, ieta+fac_eta).Wmunu[9];
                double Wetaeta_center = 
                    Util::four_dimension_linear_interpolation(
                                lattice_spacing, x_fraction, cube);

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

                // 4-dimension interpolation done
                const double TFO = eos.get_temperature(epsFO, rhob_center);
                const double muB = eos.get_mu(epsFO, rhob_center);
                if (TFO < 0) {
                    music_message << "TFO=" << TFO
                                  << "<0. ERROR. exiting.";
                    music_message.flush("error");
                    exit(1);
                }

                const double pressure = eos.get_pressure(epsFO, rhob_center);
                const double eps_plus_p_over_T_FO = (epsFO + pressure)/TFO;

                // finally output results !!!!
                s_file << scientific << setprecision(10)
                       << tau_center << " " << x_center << " "
                       << y_center << " " << eta_center << " "
                       << FULLSU[0] << " " << FULLSU[1] << " "
                       << FULLSU[2] << " " << FULLSU[3] << " "
                       << utau_center << " " << ux_center << " "
                       << uy_center << " " << ueta_center << " "
                       << epsFO << " " << TFO << " " << muB << " "
                       << eps_plus_p_over_T_FO << " "
                       << Wtautau_center << " " << Wtaux_center << " "
                       << Wtauy_center << " " << Wtaueta_center << " "
                       << Wxx_center << " " << Wxy_center << " "
                       << Wxeta_center << " "
                       << Wyy_center << " " << Wyeta_center << " "
                       << Wetaeta_center << " ";
                if (DATA.turn_on_bulk)
                    s_file << pi_b_center << " ";
                if (DATA.turn_on_rhob)
                    s_file << rhob_center << " ";
                if (DATA.turn_on_diff)
                    s_file << qtau_center << " " << qx_center << " "
                           << qy_center << " " << qeta_center << " ";
                s_file << endl;
            }
        }
    }
    s_file.close();

    // clean up
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++)
                delete [] cube[i][j][k];
            delete [] cube[i][j];
        }
        delete [] cube[i];
    }
    delete [] cube;
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
        #pragma omp parallel for
        for (int ieta = 0; ieta < neta - fac_eta; ieta += fac_eta) {
            int thread_id = omp_get_thread_num();
            FreezeOut_equal_tau_Surface_XY(tau,  ieta, arena_current,
                                           thread_id, epsFO);
        }
    }
    return(0);
}


void Evolve::FreezeOut_equal_tau_Surface_XY(double tau, int ieta,
                                            SCGrid &arena_current,
                                            int thread_id, double epsFO) {
    double epsFO_low = 0.05/hbarc;        // 1/fm^4

    const int nx = arena_current.nX();
    const int ny = arena_current.nY();

    stringstream strs_name;
    strs_name << "surface_eps_" << setprecision(4) << epsFO*hbarc
              << "_" << thread_id << ".dat";
    ofstream s_file;
    s_file.open(strs_name.str().c_str(), ios::out | ios::app);

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
            double u_flow[]       = {utau_center, ux_center, uy_center, ueta_center};
            double q_mu[]         = {qtau_center, qx_center, qy_center, qeta_center};
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

            // get other thermodynamical quantities
            double e_local   = arena_current(ix, iy, ieta).epsilon;
            double T_local   = eos.get_temperature(e_local, rhob_center);
            double muB_local = eos.get_mu(e_local, rhob_center);
            if (T_local < 0) {
                music_message << "Evolve::FreezeOut_equal_tau_Surface: "
                              << "T_local = " << T_local
                              << " <0. ERROR. exiting.";
                music_message.flush("error");
                exit(1);
            }

            double pressure = eos.get_pressure(e_local, rhob_center);
            double eps_plus_p_over_T = (e_local + pressure)/T_local;

            // finally output results
            s_file << scientific << setprecision(10) 
                   << tau_center     << " " << x_center          << " " 
                   << y_center       << " " << eta_center        << " " 
                   << FULLSU[0]      << " " << FULLSU[1]         << " " 
                   << FULLSU[2]      << " " << FULLSU[3]         << " " 
                   << utau_center    << " " << ux_center         << " " 
                   << uy_center      << " " << ueta_center       << " " 
                   << e_local        << " " << T_local           << " "
                   << muB_local      << " " << eps_plus_p_over_T << " " 
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
            s_file << endl;
        }
    }
    s_file.close();
}


int Evolve::FindFreezeOutSurface_boostinvariant_Cornelius(
                double tau, SCGrid &arena_current, SCGrid &arena_freezeout) {
    // find boost-invariant hyper-surfaces
    int *all_frozen = new int[n_freeze_surf];
    for (int i_freezesurf = 0; i_freezesurf < n_freeze_surf; i_freezesurf++) {
        double epsFO = epsFO_list[i_freezesurf]/hbarc;

        stringstream strs_name;
        strs_name << "surface_eps_" << setprecision(4) << epsFO*hbarc
                  << ".dat";

        ofstream s_file;
        s_file.open(strs_name.str().c_str(), ios::out | ios::app);

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
        Cornelius* cornelius_ptr = new Cornelius();
        cornelius_ptr->init(dim, epsFO, lattice_spacing);

        // initialize the hyper-cube for Cornelius
        double ***cube = new double ** [2];
        for (int i = 0; i < 2; i++) {
            cube[i] = new double * [2];
            for (int j = 0; j < 2; j++) {
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

                    // flow velocity u^\tau
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).u[0];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).u[0];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).u[0];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).u[0];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).u[0];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).u[0];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).u[0];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).u[0];
                    double utau_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));

                    // flow velocity u^x
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).u[1];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).u[1];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).u[1];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).u[1];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).u[1];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).u[1];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).u[1];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).u[1];
                    double ux_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));

                    // flow velocity u^y
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).u[2];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).u[2];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).u[2];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).u[2];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).u[2];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).u[2];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).u[2];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).u[2];
                    double uy_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));

                    // flow velocity u^eta
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).u[3];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).u[3];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).u[3];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).u[3];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).u[3];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).u[3];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).u[3];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).u[3];
                    double ueta_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));

                    // baryon density rho_b
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).rhob;
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).rhob;
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).rhob;
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).rhob;
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).rhob;
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).rhob;
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).rhob;
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).rhob;
                    double rhob_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));

                    // bulk viscous pressure pi_b
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).pi_b;
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).pi_b;
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).pi_b;
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).pi_b;
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).pi_b;
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).pi_b;
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).pi_b;
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).pi_b;
                    double pi_b_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));
               
                    // baryon diffusion current q^\tau
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[10];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[10];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[10];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[10];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[10];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[10];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[10];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[10];
                    double qtau_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));
               
                    // baryon diffusion current q^x
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[11];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[11];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[11];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[11];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[11];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[11];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[11];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[11];
                    double qx_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));
               
                    // baryon diffusion current q^y
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[12];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[12];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[12];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[12];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[12];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[12];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[12];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[12];
                    double qy_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));
               
                    // baryon diffusion current q^eta
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[13];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[13];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[13];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[13];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[13];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[13];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[13];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[13];
                    double qeta_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));

                    // reconstruct q^\tau from the transverality criteria
                    double u_flow[4] = {utau_center, ux_center, uy_center, ueta_center};
                    double q_mu[4]   = {qtau_center, qx_center, qy_center, qeta_center};
                    double q_regulated[4] = {0.0, 0.0, 0.0, 0.0};
                    regulate_qmu(u_flow, q_mu, q_regulated);
                    qtau_center = q_regulated[0];
                    qx_center = q_regulated[1];
                    qy_center = q_regulated[2];
                    qeta_center = q_regulated[3];

                    // shear viscous tensor W^\tau\tau
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[0];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[0];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[0];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[0];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[0];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[0];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[0];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[0];
                    double Wtautau_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));
                  
                    // shear viscous tensor W^{\tau x}
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[1];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[1];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[1];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[1];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[1];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[1];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[1];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[1];
                    double Wtaux_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));

                    // shear viscous tensor W^{\tau y}
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[2];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[2];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[2];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[2];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[2];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[2];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[2];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[2];
                    double Wtauy_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));
                  
                    // shear viscous tensor W^{\tau \eta}
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[3];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[3];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[3];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[3];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[3];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[3];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[3];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[3];
                    double Wtaueta_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));
                  
                    // shear viscous tensor W^{xx}
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[4];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[4];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[4];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[4];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[4];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[4];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[4];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[4];
                    double Wxx_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));

                    // shear viscous tensor W^{xy}
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[5];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[5];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[5];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[5];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[5];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[5];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[5];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[5];
                    double Wxy_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));

                    // shear viscous tensor W^{x \eta}
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[6];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[6];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[6];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[6];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[6];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[6];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[6];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[6];
                    double Wxeta_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));
                  
                    // shear viscous tensor W^{yy}
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[7];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[7];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[7];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[7];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[7];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[7];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[7];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[7];
                    double Wyy_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));
                  
                    // shear viscous tensor W^{yeta}
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[8];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[8];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[8];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[8];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[8];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[8];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[8];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[8];
                    double Wyeta_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));
                  
                    // shear viscous tensor W^{\eta\eta}
                    cube[0][0][0] = arena_freezeout(ix      , iy      , 0).Wmunu[9];
                    cube[0][0][1] = arena_freezeout(ix      , iy+fac_y, 0).Wmunu[9];
                    cube[0][1][0] = arena_freezeout(ix+fac_x, iy      , 0).Wmunu[9];
                    cube[0][1][1] = arena_freezeout(ix+fac_x, iy+fac_y, 0).Wmunu[9];
                    cube[1][0][0] = arena_current  (ix      , iy      , 0).Wmunu[9];
                    cube[1][0][1] = arena_current  (ix      , iy+fac_y, 0).Wmunu[9];
                    cube[1][1][0] = arena_current  (ix+fac_x, iy      , 0).Wmunu[9];
                    cube[1][1][1] = arena_current  (ix+fac_x, iy+fac_y, 0).Wmunu[9];
                    double Wetaeta_center = (
                        Util::three_dimension_linear_interpolation(
                                        lattice_spacing, x_fraction, cube));

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

                    // 3-dimension interpolation done
                    double TFO = eos.get_temperature(epsFO, rhob_center);
                    double muB = eos.get_mu(epsFO, rhob_center);
                    if (TFO < 0) {
                        music_message << "TFO=" << TFO
                                      << "<0. ERROR. exiting.";
                        music_message.flush("error");
                        exit(1);
                    }

                    double pressure = eos.get_pressure(epsFO, rhob_center);
                    double eps_plus_p_over_T_FO = (epsFO + pressure)/TFO;

                    // finally output results !!!!
                    s_file << scientific << setprecision(10) 
                           << tau_center << " " << x_center << " " 
                           << y_center << " " << eta_center << " " 
                           << FULLSU[0] << " " << FULLSU[1] << " " 
                           << FULLSU[2] << " " << FULLSU[3] << " " 
                           << utau_center << " " << ux_center << " " 
                           << uy_center << " " << ueta_center << " " 
                           << epsFO << " " << TFO << " " << muB << " " 
                           << eps_plus_p_over_T_FO << " " 
                           << Wtautau_center << " " << Wtaux_center << " " 
                           << Wtauy_center << " " << Wtaueta_center << " " 
                           << Wxx_center << " " << Wxy_center << " " 
                           << Wxeta_center << " " 
                           << Wyy_center << " " << Wyeta_center << " " 
                           << Wetaeta_center << " " ;
                    if(DATA.turn_on_bulk)   // 27th column
                        s_file << pi_b_center << " " ;
                    if(DATA.turn_on_rhob)   // 28th column
                        s_file << rhob_center << " " ;
                    if(DATA.turn_on_diff)   // 29-32th column
                        s_file << qtau_center << " " << qx_center << " " 
                               << qy_center << " " << qeta_center << " " ;
                    s_file << endl;
                }
            }
        }

        s_file.close();
        // clean up
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++)
                delete [] cube[i][j];
            delete [] cube[i];
        }
        delete [] cube;
        delete cornelius_ptr;

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

void Evolve::regulate_qmu(const double u[], const double q[],
                          double q_regulated[]) const {
    double u_dot_q = - u[0]*q[0] + u[1]*q[1] + u[2]*q[2] + u[3]*q[3];
    for (int i = 0; i < 4; i++) {
        q_regulated[i] = q[i] + u[i]*u_dot_q;
    }
}

void Evolve::regulate_Wmunu(const double u[], const double Wmunu[4][4],
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
    int freeze_eps_flag = DATA.freeze_eps_flag;
    if (freeze_eps_flag == 0) {
        // constant spacing the energy density
        n_freeze_surf = DATA.N_freeze_out;
        double freeze_max_ed = DATA.eps_freeze_max;
        double freeze_min_ed = DATA.eps_freeze_min;
        double d_epsFO = ((freeze_max_ed - freeze_min_ed)
                          /(n_freeze_surf - 1 + 1e-15));
        for(int isurf = 0; isurf < n_freeze_surf; isurf++)
        {
            double temp_epsFO = freeze_min_ed + isurf*d_epsFO;
            epsFO_list.push_back(temp_epsFO);
        }
    } else if(freeze_eps_flag == 1) {
        // read in from a file
        string eps_freeze_list_filename = DATA.freeze_list_filename;
        music_message << "read in freeze out surface information from " 
                      << eps_freeze_list_filename;
        music_message.flush("info");
        ifstream freeze_list_file(eps_freeze_list_filename.c_str());
        if (!freeze_list_file) {
            music_message << "Evolve::initialize_freezeout_surface_info: "
                          << "can not open freeze-out list file: " 
                          << eps_freeze_list_filename;
            music_message.flush("error");
            exit(1);
        }
        int temp_n_surf = 0;
        string dummy;
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
