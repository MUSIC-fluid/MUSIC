// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <iomanip>
#include <omp.h>
#include "./util.h"
#include "./cell.h"
#include "./grid.h"
#include "./init.h"
#include "./eos.h"

using namespace std;

Init::Init(EOS *eosIn, InitData *DATA_in, hydro_source *hydro_source_in) {
    eos = eosIn;
    DATA_ptr = DATA_in;
    if (DATA_ptr->Initial_profile == 12 || DATA_ptr->Initial_profile == 13) {
        hydro_source_ptr = hydro_source_in;
    } else if (DATA_ptr->Initial_profile == 30) {
        hydro_source_ptr = hydro_source_in;
    }
}

void Init::InitArena(InitData *DATA, Grid &arena) {
    music_message.info("initArena");
    if (DATA->Initial_profile == 0) {
        music_message << "Using Initial_profile=" << DATA->Initial_profile;
        music_message << "nx=" << DATA->nx << ", ny=" << DATA->ny;
        music_message << "dx=" << DATA->delta_x << ", dy=" << DATA->delta_y;
        music_message.flush("info");
    } else if (DATA->Initial_profile == 1) {
        music_message << "Using Initial_profile=" << DATA->Initial_profile;
        DATA->nx = 2;
        DATA->ny = 2;
        DATA->neta = 695;
        DATA->delta_x = 0.1;
        DATA->delta_y = 0.1;
        DATA->delta_eta = 0.02;
        music_message << "nx=" << DATA->nx+1 << ", ny=" << DATA->ny+1;
        music_message << "dx=" << DATA->delta_x << ", dy=" << DATA->delta_y;
        music_message << "neta=" << DATA->neta << ", deta=" << DATA->delta_eta;
        music_message.flush("info");
    } else if (DATA->Initial_profile == 8) {
        music_message.info(DATA->initName);
        ifstream profile(DATA->initName.c_str());
        string dummy;
        int nx, ny, neta;
        double deta, dx, dy, dummy2;
        // read the first line with general info
        profile >> dummy >> dummy >> dummy2
                >> dummy >> neta >> dummy >> nx >> dummy >> ny
                >> dummy >> deta >> dummy >> dx >> dummy >> dy;
        profile.close();
        music_message << "Using Initial_profile=" << DATA->Initial_profile
                      << ". Overwriting lattice dimensions:";

        DATA->nx = nx;
        DATA->ny = ny;
        DATA->delta_x = dx;
        DATA->delta_y = dy;

        music_message << "neta=" << neta << ", nx=" << nx << ", ny=" << ny;
        music_message << "deta=" << DATA->delta_eta << ", dx=" << DATA->delta_x
                      << ", dy=" << DATA->delta_y;
        music_message.flush("info");
    } else if (DATA->Initial_profile == 9) {
        music_message.info(DATA->initName);
        ifstream profile(DATA->initName.c_str());
        string dummy;
        int nx, ny, neta;
        double deta, dx, dy, dummy2;
        // read the first line with general info
        profile >> dummy >> dummy >> dummy2
                >> dummy >> neta >> dummy >> nx >> dummy >> ny
                >> dummy >> deta >> dummy >> dx >> dummy >> dy;
        profile.close();
        music_message << "Using Initial_profile=" << DATA->Initial_profile
                      << ". Overwriting lattice dimensions:";

        DATA->nx = nx;
        DATA->ny = ny;
        DATA->neta = neta;
        DATA->delta_x = dx;
        DATA->delta_y = dy;
        DATA->delta_eta = 0.1;

        music_message << "neta=" << neta << ", nx=" << nx << ", ny=" << ny;
        music_message << "deta=" << DATA->delta_eta << ", dx=" << DATA->delta_x
                      << ", dy=" << DATA->delta_y;
        music_message.flush("info");
    } else if (DATA->Initial_profile == 12 || DATA->Initial_profile == 13) {
        DATA->tau0 = hydro_source_ptr->get_source_tau_min();
    } else if (DATA->Initial_profile == 101) {
        cout << "Using Initial_profile=" << DATA->Initial_profile << endl;
        cout << "nx=" << DATA->nx << ", ny=" << DATA->ny << endl;
        cout << "dx=" << DATA->delta_x << ", dy=" << DATA->delta_y << endl;
    }

    // initialize arena
    arena = Grid(DATA->nx, DATA->ny,DATA->neta);
    music_message.info("Cell allocated.");

    InitTJb(DATA, arena);

    if (DATA->output_initial_density_profiles == 1) {
        output_initial_density_profiles(DATA, arena);
    }
}/* InitArena */

//! This is a shell function to initial hydrodynamic fields
int Init::InitTJb(InitData *DATA, Grid &arena) {
    if (DATA->Initial_profile == 0) {
        // Gubser flow test
        music_message.info(" Perform Gubser flow test ... ");
        music_message.info(" ----- information on initial distribution -----");
        
        #pragma omp parallel for
        for (int ieta = 0; ieta < arena.nEta(); ieta++) {
            //cout << "[Info] Thread " << omp_get_thread_num()
            //     << " executes loop iteraction ieta = " << ieta << endl;
            initial_Gubser_XY(DATA, ieta, arena);
        }/* ieta */
    } else if (DATA->Initial_profile == 1) {
        // code test in 1+1 D vs Monnai's results
        music_message.info(" Perform 1+1D test vs Monnai's results... ");
        initial_1p1D_eta(DATA, arena);
    } else if (DATA->Initial_profile == 8) {
        // read in the profile from file
        // - IPGlasma initial conditions with initial flow
        music_message.info(" ----- information on initial distribution -----");
        music_message << "file name used: " << DATA->initName;
        music_message.flush("info");
  
        #pragma omp parallel for
        for (int ieta = 0; ieta < arena.nEta(); ieta++) {
            //cout << "[Info] Thread " << omp_get_thread_num()
            //     << " executes loop iteraction ieta = " << ieta << endl;
            initial_IPGlasma_XY(DATA, ieta, arena);
        } /* ieta */
    } else if (DATA->Initial_profile == 9) {
        // read in the profile from file
        // - IPGlasma initial conditions with initial flow
        // and initial shear viscous tensor
        music_message.info(" ----- information on initial distribution -----");
        music_message << "file name used: " << DATA->initName;
        music_message.flush("info");
  
        #pragma omp parallel for
        for (int ieta = 0; ieta < arena.nEta(); ieta++) {
            //cout << "[Info] Thread " << omp_get_thread_num()
            //     << " executes loop iteraction ieta = " << ieta << endl;
            initial_IPGlasma_XY_with_pi(DATA, ieta, arena);
        } /* ieta */
    } else if (DATA->Initial_profile == 11) {
        // read in the transverse profile from file with finite rho_B
        // the initial entropy and net baryon density profile are
        // constructed by nuclear thickness function TA and TB.
        // Along the longitudinal direction an asymmetric contribution from
        // target and projectile thickness function is allowed
        music_message.info(" ----- information on initial distribution -----");
        music_message << "file name used: " << DATA->initName_TA << " and "
                      << DATA->initName_TB;
        music_message.flush("info");

        #pragma omp parallel for
        for (int ieta = 0; ieta < arena.nEta(); ieta++) {
            //cout << "[Info] Thread " << omp_get_thread_num()
            //     << " executes loop iteraction ieta = " << ieta << endl;
            initial_MCGlb_with_rhob_XY(DATA, ieta, arena);
        } /* ix, iy, ieta */
    } else if (DATA->Initial_profile == 12 || DATA->Initial_profile == 13) {
        music_message.info("Initialize hydro with source terms");
        #pragma omp parallel for
        for (int ieta = 0; ieta < arena.nEta(); ieta++) {
            //cout << "[Info] Thread " << omp_get_thread_num()
            //     << " executes loop iteraction ieta = " << ieta << endl;
            initial_MCGlbLEXUS_with_rhob_XY(DATA, ieta, arena);
        } /* ix, iy, ieta */
    } else if (DATA->Initial_profile == 30) {
        #pragma omp parallel for
        for (int ieta = 0; ieta < arena.nEta(); ieta++) {
            //cout << "[Info] Thread " << omp_get_thread_num()
            //     << " executes loop iteraction ieta = " << ieta << endl;
            initial_AMPT_XY(DATA, ieta, arena);
        }
    } else if (DATA->Initial_profile == 101) {
        initial_UMN_with_rhob(DATA, arena);
    }
    music_message.info("initial distribution done.");
    return 1;
}  /* InitTJb*/

void Init::initial_Gubser_XY(InitData *DATA, int ieta, Grid &arena) {
    string input_filename;
    string input_filename_prev;
    if (DATA->turn_on_shear == 1) {
        input_filename = "tests/Gubser_flow/Initial_Profile.dat";
    } else {
        input_filename = "tests/Gubser_flow/y=0_tau=1.00_ideal.dat";
        input_filename_prev = "tests/Gubser_flow/y=0_tau=0.98_ideal.dat";
    }
    
    ifstream profile(input_filename.c_str());
    if (!profile.good()) {
        music_message << "Init::InitTJb: "
                      << "Can not open the initial file: " << input_filename;
        music_message.flush("error");
        exit(1);
    }
    ifstream profile_prev;
    if (DATA->turn_on_shear == 0) {
        profile_prev.open(input_filename_prev.c_str());
        if (!profile_prev.good()) {
            music_message << "Init::InitTJb: "
                          << "Can not open the initial file: "
                          << input_filename_prev;
            music_message.flush("error");
            exit(1);
        }
    }

    const int nx = arena.nX();
    const int ny = arena.nY();
    double** temp_profile_ed        = new double* [nx];
    double** temp_profile_ux        = new double* [nx];
    double** temp_profile_uy        = new double* [nx];
    double **temp_profile_ed_prev   = NULL;
    double **temp_profile_rhob      = NULL;
    double **temp_profile_rhob_prev = NULL;
    double **temp_profile_ux_prev   = NULL;
    double **temp_profile_uy_prev   = NULL;
    double **temp_profile_pixx      = NULL;
    double **temp_profile_piyy      = NULL;
    double **temp_profile_pixy      = NULL;
    double **temp_profile_pi00      = NULL;
    double **temp_profile_pi0x      = NULL;
    double **temp_profile_pi0y      = NULL;
    double **temp_profile_pi33      = NULL;
    if (DATA->turn_on_shear == 1) {
        temp_profile_pixx = new double* [nx];
        temp_profile_piyy = new double* [nx];
        temp_profile_pixy = new double* [nx];
        temp_profile_pi00 = new double* [nx];
        temp_profile_pi0x = new double* [nx];
        temp_profile_pi0y = new double* [nx];
        temp_profile_pi33 = new double* [nx];
    } else {
        temp_profile_ed_prev   = new double* [nx];
        temp_profile_rhob      = new double* [nx];
        temp_profile_rhob_prev = new double* [nx];
        temp_profile_ux_prev   = new double* [nx];
        temp_profile_uy_prev   = new double* [nx];
    }
    for (int i = 0; i < nx; i++) {
        temp_profile_ed[i] = new double[ny];
        temp_profile_ux[i] = new double[ny];
        temp_profile_uy[i] = new double[ny];
        if (DATA->turn_on_shear == 1) {
            temp_profile_pixx[i] = new double[ny];
            temp_profile_pixy[i] = new double[ny];
            temp_profile_piyy[i] = new double[ny];
            temp_profile_pi00[i] = new double[ny];
            temp_profile_pi0x[i] = new double[ny];
            temp_profile_pi0y[i] = new double[ny];
            temp_profile_pi33[i] = new double[ny];
        } else {
            temp_profile_ed_prev[i]   = new double[ny];
            temp_profile_rhob[i]      = new double[ny];
            temp_profile_rhob_prev[i] = new double[ny];
            temp_profile_ux_prev[i]   = new double[ny];
            temp_profile_uy_prev[i]   = new double[ny];
        }
    }

    double dummy;
    double u[4];
    int rk_order = DATA->rk_order;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            if (DATA->turn_on_shear == 1) {
                profile >> dummy >> dummy >> temp_profile_ed[ix][iy]
                        >> temp_profile_ux[ix][iy] >> temp_profile_uy[ix][iy];
                profile >> temp_profile_pixx[ix][iy]
                        >> temp_profile_piyy[ix][iy]
                        >> temp_profile_pixy[ix][iy]
                        >> temp_profile_pi00[ix][iy]
                        >> temp_profile_pi0x[ix][iy]
                        >> temp_profile_pi0y[ix][iy]
                        >> temp_profile_pi33[ix][iy];
            } else {
                profile >> dummy >> dummy >> temp_profile_ed[ix][iy]
                        >> temp_profile_rhob[ix][iy]
                        >> temp_profile_ux[ix][iy] >> temp_profile_uy[ix][iy];
                profile_prev >> dummy >> dummy >> temp_profile_ed_prev[ix][iy]
                             >> temp_profile_rhob_prev[ix][iy]
                             >> temp_profile_ux_prev[ix][iy]
                             >> temp_profile_uy_prev[ix][iy];
            }
        }
    }
    profile.close();
    if (DATA->turn_on_shear == 0) {
        profile_prev.close();
    }

    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy< ny; iy++) {
            double rhob = 0.0;
            if (DATA->turn_on_shear == 0) {
                if (DATA->turn_on_rhob == 1) {
                    rhob = temp_profile_rhob[ix][iy];
                }
            }

            double epsilon = temp_profile_ed[ix][iy];
            
            // set all values in the grid element:
            arena(ix,iy,ieta).epsilon      = epsilon;
            arena(ix,iy,ieta).epsilon_t    = epsilon;
            arena(ix,iy,ieta).prev_epsilon = epsilon;
            arena(ix,iy,ieta).rhob         = rhob;
            arena(ix,iy,ieta).rhob_t       = rhob;
            arena(ix,iy,ieta).prev_rhob    = rhob;
            
            /* for HIC */
            double utau_local = sqrt(1.
                          + temp_profile_ux[ix][iy]*temp_profile_ux[ix][iy]
                          + temp_profile_uy[ix][iy]*temp_profile_uy[ix][iy]);
            arena(ix,iy,ieta).u[0][0] = utau_local;
            arena(ix,iy,ieta).u[0][1] = temp_profile_ux[ix][iy];
            arena(ix,iy,ieta).u[0][2] = temp_profile_uy[ix][iy];
            arena(ix,iy,ieta).u[0][3] = 0.0;

            u[0] = utau_local;
            u[1] = temp_profile_ux[ix][iy];
            u[2] = temp_profile_uy[ix][iy];
            u[3] = 0.0;
            for (int rk_i = 0; rk_i < 1; rk_i++) {
                arena(ix,iy,ieta).prev_u[rk_i][0] = u[0];
                arena(ix,iy,ieta).prev_u[rk_i][1] = u[1];
                arena(ix,iy,ieta).prev_u[rk_i][2] = u[2];
                arena(ix,iy,ieta).prev_u[rk_i][3] = u[3];
                arena(ix,iy,ieta).prev_pi_b[rk_i] = 0.0;
            }

            if (DATA->turn_on_shear == 0) {
                double utau_prev = sqrt(1.
                    + temp_profile_ux_prev[ix][iy]*temp_profile_ux_prev[ix][iy]
                    + temp_profile_uy_prev[ix][iy]*temp_profile_uy_prev[ix][iy]
                );
                for (int rk_i = 0; rk_i < 1; rk_i++) {
                    arena(ix,iy,ieta).prev_u[rk_i][0] = utau_prev;
                    arena(ix,iy,ieta).prev_u[rk_i][1] =
                                                temp_profile_ux_prev[ix][iy];
                    arena(ix,iy,ieta).prev_u[rk_i][2] =
                                                temp_profile_uy_prev[ix][iy];
                    arena(ix,iy,ieta).prev_u[rk_i][3] = 0.0;
                }
            }
            arena(ix,iy,ieta).pi_b[0] = 0.0;

            if (DATA->turn_on_shear == 1) {
                arena(ix,iy,ieta).Wmunu[0][0] = temp_profile_pi00[ix][iy];
                arena(ix,iy,ieta).Wmunu[0][1] = temp_profile_pi0x[ix][iy];
                arena(ix,iy,ieta).Wmunu[0][2] = temp_profile_pi0y[ix][iy];
                arena(ix,iy,ieta).Wmunu[0][3] = 0.0;
                arena(ix,iy,ieta).Wmunu[0][4] = temp_profile_pixx[ix][iy];
                arena(ix,iy,ieta).Wmunu[0][5] = temp_profile_pixy[ix][iy];
                arena(ix,iy,ieta).Wmunu[0][6] = 0.0;
                arena(ix,iy,ieta).Wmunu[0][7] = temp_profile_piyy[ix][iy];
                arena(ix,iy,ieta).Wmunu[0][8] = 0.0;
                arena(ix,iy,ieta).Wmunu[0][9] = temp_profile_pi33[ix][iy];
                for (int mu = 10; mu < 14; mu++) {
                        arena(ix,iy,ieta).Wmunu[0][mu] = 0.0;
                }
            } else {
                for (int mu = 0; mu < 14; mu++) {
                        arena(ix,iy,ieta).Wmunu[0][mu] = 0.0;
                }
            }
            for (int rkstep = 0; rkstep < 1; rkstep++) {
                for (int ii = 0; ii < 14; ii++) {
                    arena(ix,iy,ieta).prevWmunu[rkstep][ii] = 
                                        arena(ix,iy,ieta).Wmunu[0][ii];
                }
            }
        }
    }
    // clean up
    for (int i = 0; i < nx; i++) {
        delete[] temp_profile_ed[i];
        delete[] temp_profile_ux[i];
        delete[] temp_profile_uy[i];
        if (DATA->turn_on_shear == 1) {
            delete[] temp_profile_pixx[i];
            delete[] temp_profile_piyy[i];
            delete[] temp_profile_pixy[i];
            delete[] temp_profile_pi00[i];
            delete[] temp_profile_pi0x[i];
            delete[] temp_profile_pi0y[i];
            delete[] temp_profile_pi33[i];
        } else {
            delete[] temp_profile_ed_prev[i];
            delete[] temp_profile_rhob[i];
            delete[] temp_profile_rhob_prev[i];
            delete[] temp_profile_ux_prev[i];
            delete[] temp_profile_uy_prev[i];
        }
    }
    delete[] temp_profile_ed;
    delete[] temp_profile_ux;
    delete[] temp_profile_uy;
    if (DATA->turn_on_shear == 1) {
        delete[] temp_profile_pixx;
        delete[] temp_profile_piyy;
        delete[] temp_profile_pixy;
        delete[] temp_profile_pi00;
        delete[] temp_profile_pi0x;
        delete[] temp_profile_pi0y;
        delete[] temp_profile_pi33;
    } else {
        delete[] temp_profile_ed_prev;
        delete[] temp_profile_rhob;
        delete[] temp_profile_rhob_prev;
        delete[] temp_profile_ux_prev;
        delete[] temp_profile_uy_prev;
    }
}

void Init::initial_1p1D_eta(InitData *DATA, Grid &arena) {
    string input_ed_filename;
    string input_rhob_filename;
    input_ed_filename = "tests/test_1+1D_with_Akihiko/e_baryon_init.dat";
    input_rhob_filename = "tests/test_1+1D_with_Akihiko/rhoB_baryon_init.dat";

    ifstream profile_ed(input_ed_filename.c_str());
    if (!profile_ed.good()) {
        music_message << "Init::InitTJb: "
                      << "Can not open the initial file: "
                      << input_ed_filename;
        music_message.flush("error");
        exit(1);
    }
    ifstream profile_rhob;
    profile_rhob.open(input_rhob_filename.c_str());
    if (!profile_rhob.good()) {
        music_message << "Init::InitTJb: "
                      << "Can not open the initial file: "
                      << input_rhob_filename;
        music_message.flush("error");
        exit(1);
    }

    const int neta = arena.nEta();
    double *temp_profile_ed = new double[neta];
    double *temp_profile_rhob = new double[neta];

    double dummy;
    double u[4];
    int rk_order = DATA->rk_order;
    for (int ieta = 0; ieta < neta; ieta++) {
        profile_ed >> dummy >> temp_profile_ed[ieta];
        profile_rhob >> dummy >> temp_profile_rhob[ieta];
    }
    profile_ed.close();
    profile_rhob.close();

    const int nx = arena.nX();
    const int ny = arena.nY();
    for (int ieta = 0; ieta < neta; ieta++) {
        double rhob = temp_profile_rhob[ieta];
        double epsilon = temp_profile_ed[ieta]/hbarc;   // fm^-4
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy< ny; iy++) {
                // set all values in the grid element:
                arena(ix,iy,ieta).epsilon      = epsilon;
                arena(ix,iy,ieta).epsilon_t    = epsilon;
                arena(ix,iy,ieta).prev_epsilon = epsilon;
                arena(ix,iy,ieta).rhob         = rhob;
                arena(ix,iy,ieta).rhob_t       = rhob;
                arena(ix,iy,ieta).prev_rhob    = rhob;
            
                /* for HIC */
                arena(ix,iy,ieta).u[0][0] = 1.0;
                arena(ix,iy,ieta).u[0][1] = 0.0;
                arena(ix,iy,ieta).u[0][2] = 0.0;
                arena(ix,iy,ieta).u[0][3] = 0.0;

                u[0] = 1.0;
                u[1] = 0.0;
                u[2] = 0.0;
                u[3] = 0.0;
                for (int rk_i = 0; rk_i < 1; rk_i++) {
                    arena(ix,iy,ieta).prev_u[rk_i][0] = u[0];
                    arena(ix,iy,ieta).prev_u[rk_i][1] = u[1];
                    arena(ix,iy,ieta).prev_u[rk_i][2] = u[2];
                    arena(ix,iy,ieta).prev_u[rk_i][3] = u[3];
                    arena(ix,iy,ieta).prev_pi_b[rk_i] = 0.0;
                }
                arena(ix,iy,ieta).pi_b[0] = 0.0;

                for (int mu = 0; mu < 14; mu++) {
                    arena(ix,iy,ieta).Wmunu[0][mu] = 0.0;
                }
            
                for (int rkstep = 0; rkstep < 1; rkstep++) {
                    for (int ii = 0; ii < 14; ii++) {
                        arena(ix,iy,ieta).prevWmunu[rkstep][ii] = 
                                        arena(ix,iy,ieta).Wmunu[0][ii];
                    }
                }
            }
        }
    }
    // clean up
    delete[] temp_profile_ed;
    delete[] temp_profile_rhob;
}

void Init::initial_IPGlasma_XY(InitData *DATA, int ieta, Grid &arena) {
    ifstream profile(DATA->initName.c_str());

    string dummy;
    // read the information line
    std::getline(profile, dummy);

    const int nx = arena.nX();
    const int ny = arena.nY();
    double** temp_profile_ed   = new double* [nx];
    double** temp_profile_utau = new double* [nx];
    double** temp_profile_ux   = new double* [nx];
    double** temp_profile_uy   = new double* [nx];
    for (int i = 0; i < nx; i++) {
        temp_profile_ed[i]   = new double[ny];
        temp_profile_utau[i] = new double[ny];
        temp_profile_ux[i]   = new double[ny];
        temp_profile_uy[i]   = new double[ny];
    }

    // read the one slice
    double density, dummy1, dummy2, dummy3;
    double ux, uy, utau;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            profile >> dummy1 >> dummy2 >> dummy3
                    >> density >> utau >> ux >> uy
                    >> dummy  >> dummy  >> dummy  >> dummy;
            temp_profile_ed[ix][iy] = density;
            temp_profile_ux[ix][iy] = ux;
            temp_profile_uy[ix][iy] = uy;
            temp_profile_utau[ix][iy] = sqrt(1. + ux*ux + uy*uy);
            if (ix == 0 && iy == 0) {
                DATA->x_size = -dummy2*2;
                DATA->y_size = -dummy3*2;
                if (omp_get_thread_num() == 0) {
                    music_message << "eta_size=" << DATA->eta_size
                                  << ", x_size=" << DATA->x_size
                                  << ", y_size=" << DATA->y_size;
                    music_message.flush("info");
                }
            }
        }
    }
    profile.close();

    double eta = (DATA->delta_eta)*ieta - (DATA->eta_size)/2.0;
    double eta_envelop_ed = eta_profile_normalisation(DATA, eta);
    int entropy_flag = DATA->initializeEntropy;
    double u[4];
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy< ny; iy++) {
            double rhob = 0.0;
            double epsilon = 0.0;
            if (entropy_flag == 0) {
                epsilon = (temp_profile_ed[ix][iy]*eta_envelop_ed
                           *DATA->sFactor/hbarc);  // 1/fm^4
            } else {
                double local_sd = (temp_profile_ed[ix][iy]*DATA->sFactor
                                   *eta_envelop_ed);
                epsilon = eos->get_s2e(local_sd, rhob);
            }
            if (epsilon < 0.00000000001)
                epsilon = 0.00000000001;

            // set all values in the grid element:
            arena(ix,iy,ieta).epsilon      = epsilon;
            arena(ix,iy,ieta).epsilon_t    = epsilon;
            arena(ix,iy,ieta).prev_epsilon = epsilon;
            arena(ix,iy,ieta).rhob         = rhob;
            arena(ix,iy,ieta).rhob_t       = rhob;
            arena(ix,iy,ieta).prev_rhob    = rhob;

            /* for HIC */
            arena(ix,iy,ieta).u[0][0] = temp_profile_utau[ix][iy];
            arena(ix,iy,ieta).u[0][1] = temp_profile_ux[ix][iy];
            arena(ix,iy,ieta).u[0][2] = temp_profile_uy[ix][iy];
            arena(ix,iy,ieta).u[0][3] = 0.0;

            u[0] = temp_profile_utau[ix][iy];
            u[1] = temp_profile_ux[ix][iy];
            u[2] = temp_profile_uy[ix][iy];
            u[3] = 0.0;

            for (int ii = 0; ii < 1; ii++) {
                arena(ix,iy,ieta).prev_u[ii][0] = u[0];
                arena(ix,iy,ieta).prev_u[ii][1] = u[1];
                arena(ix,iy,ieta).prev_u[ii][2] = u[2];
                arena(ix,iy,ieta).prev_u[ii][3] = u[3];
                arena(ix,iy,ieta).prev_pi_b[ii] = 0.0;
            }

            arena(ix,iy,ieta).pi_b[0] = 0.0;

            for (int ii = 0; ii < 14; ii++) {
                arena(ix,iy,ieta).prevWmunu[0][ii] = 0.0;
                arena(ix,iy,ieta).Wmunu[0][ii] = 0.0;
            }
        }
    }
    // clean up
    for (int i = 0; i < nx; i++) {
        delete[] temp_profile_ed[i];
        delete[] temp_profile_utau[i];
        delete[] temp_profile_ux[i];
        delete[] temp_profile_uy[i];
    }
    delete[] temp_profile_ed;
    delete[] temp_profile_utau;
    delete[] temp_profile_ux;
    delete[] temp_profile_uy;
}

void Init::initial_IPGlasma_XY_with_pi(InitData *DATA, int ieta,
                                       Grid &arena) {
    double tau0 = DATA->tau0;
    ifstream profile(DATA->initName.c_str());

    string dummy;
    // read the information line
    std::getline(profile, dummy);

    const int nx = arena.nX();
    const int ny = arena.nY();
    double** temp_profile_ed       = new double* [nx];
    double** temp_profile_utau     = new double* [nx];
    double** temp_profile_ux       = new double* [nx];
    double** temp_profile_uy       = new double* [nx];
    double** temp_profile_ueta     = new double* [nx];
    double** temp_profile_pitautau = new double* [nx];
    double** temp_profile_pitaux   = new double* [nx];
    double** temp_profile_pitauy   = new double* [nx];
    double** temp_profile_pitaueta = new double* [nx];
    double** temp_profile_pixx     = new double* [nx];
    double** temp_profile_pixy     = new double* [nx];
    double** temp_profile_pixeta   = new double* [nx];
    double** temp_profile_piyy     = new double* [nx];
    double** temp_profile_piyeta   = new double* [nx];
    double** temp_profile_pietaeta = new double* [nx];
    for (int i = 0; i < nx; i++) {
        temp_profile_ed[i]       = new double[ny];
        temp_profile_utau[i]     = new double[ny];
        temp_profile_ux[i]       = new double[ny];
        temp_profile_uy[i]       = new double[ny];
        temp_profile_ueta[i]     = new double[ny];
        temp_profile_pitautau[i] = new double[ny];
        temp_profile_pitaux[i]   = new double[ny];
        temp_profile_pitauy[i]   = new double[ny];
        temp_profile_pitaueta[i] = new double[ny];
        temp_profile_pixx[i]     = new double[ny];
        temp_profile_pixy[i]     = new double[ny];
        temp_profile_pixeta[i]   = new double[ny];
        temp_profile_piyy[i]     = new double[ny];
        temp_profile_piyeta[i]   = new double[ny];
        temp_profile_pietaeta[i] = new double[ny];
    }

    // read the one slice
    double density, dummy1, dummy2, dummy3;
    double ux, uy, utau, ueta;
    double pitautau, pitaux, pitauy, pitaueta;
    double pixx, pixy, pixeta, piyy, piyeta, pietaeta;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            profile >> dummy1 >> dummy2 >> dummy3
                    >> density >> utau >> ux >> uy >> ueta
                    >> pitautau >> pitaux >> pitauy >> pitaueta
                    >> pixx >> pixy >> pixeta >> piyy >> piyeta >> pietaeta;
            temp_profile_ed      [ix][iy] = density;
            temp_profile_ux      [ix][iy] = ux;
            temp_profile_uy      [ix][iy] = uy;
            temp_profile_ueta    [ix][iy] = ueta*tau0;
            temp_profile_utau    [ix][iy] = sqrt(1. + ux*ux + uy*uy + ueta*ueta);
            temp_profile_pitautau[ix][iy] = pitautau*DATA->sFactor;
            temp_profile_pitaux  [ix][iy] = pitaux*DATA->sFactor;
            temp_profile_pitauy  [ix][iy] = pitauy*DATA->sFactor;
            temp_profile_pitaueta[ix][iy] = pitaueta*tau0*DATA->sFactor;
            temp_profile_pixx    [ix][iy] = pixx*DATA->sFactor;
            temp_profile_pixy    [ix][iy] = pixy*DATA->sFactor;
            temp_profile_pixeta  [ix][iy] = pixeta*tau0*DATA->sFactor;
            temp_profile_piyy    [ix][iy] = piyy*DATA->sFactor;
            temp_profile_piyeta  [ix][iy] = piyeta*tau0*DATA->sFactor;
            temp_profile_pietaeta[ix][iy] = piyeta*tau0*tau0*DATA->sFactor;
            if (ix == 0 && iy == 0) {
                DATA->x_size = -dummy2*2;
                DATA->y_size = -dummy3*2;
                if (omp_get_thread_num() == 0) {
                    music_message << "eta_size=" << DATA->eta_size
                                  << ", x_size=" << DATA->x_size
                                  << ", y_size=" << DATA->y_size;
                    music_message.flush("info");
                }
            }
        }
    }
    profile.close();

    double eta = (DATA->delta_eta)*(ieta) - (DATA->eta_size)/2.0;
    double eta_envelop_ed = eta_profile_normalisation(DATA, eta);
    int entropy_flag = DATA->initializeEntropy;
    int rk_order = DATA->rk_order;
    double u[4];
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy< ny; iy++) {
            double rhob = 0.0;
            double epsilon = 0.0;
            if (entropy_flag == 0) {
                epsilon = (temp_profile_ed[ix][iy]*eta_envelop_ed
                           *DATA->sFactor/hbarc);  // 1/fm^4
            } else {
                double local_sd = (temp_profile_ed[ix][iy]*DATA->sFactor
                                   *eta_envelop_ed);
                epsilon = eos->get_s2e(local_sd, rhob);
            }
            if (epsilon < 0.00000000001)
                epsilon = 0.00000000001;

            double pressure = eos->get_pressure(epsilon, rhob);

            // set all values in the grid element:
            arena(ix,iy,ieta).epsilon      = epsilon;
            arena(ix,iy,ieta).epsilon_t    = epsilon;
            arena(ix,iy,ieta).prev_epsilon = epsilon;
            arena(ix,iy,ieta).rhob         = rhob;
            arena(ix,iy,ieta).rhob_t       = rhob;
            arena(ix,iy,ieta).prev_rhob    = rhob;


            /* for HIC */
            arena(ix,iy,ieta).u[0][0] = temp_profile_utau[ix][iy];
            arena(ix,iy,ieta).u[0][1] = temp_profile_ux[ix][iy];
            arena(ix,iy,ieta).u[0][2] = temp_profile_uy[ix][iy];
            arena(ix,iy,ieta).u[0][3] = temp_profile_ueta[ix][iy];

            u[0] = temp_profile_utau[ix][iy];
            u[1] = temp_profile_ux[ix][iy];
            u[2] = temp_profile_uy[ix][iy];
            u[3] = temp_profile_ueta[ix][iy];

            for (int ii = 0; ii < 1; ii++) {
                arena(ix,iy,ieta).prev_u[ii][0] = u[0];
                arena(ix,iy,ieta).prev_u[ii][1] = u[1];
                arena(ix,iy,ieta).prev_u[ii][2] = u[2];
                arena(ix,iy,ieta).prev_u[ii][3] = u[3];
                arena(ix,iy,ieta).prev_pi_b[ii] = epsilon/3. - pressure;
            }

            arena(ix,iy,ieta).pi_b[0] = epsilon/3. - pressure;

            arena(ix,iy,ieta).prevWmunu[0][0] = temp_profile_pitautau[ix][iy];
            arena(ix,iy,ieta).Wmunu    [0][0] = temp_profile_pitautau[ix][iy];
            arena(ix,iy,ieta).prevWmunu[0][1] = temp_profile_pitaux[ix][iy];
            arena(ix,iy,ieta).Wmunu    [0][1] = temp_profile_pitaux[ix][iy];
            arena(ix,iy,ieta).prevWmunu[0][2] = temp_profile_pitauy[ix][iy];
            arena(ix,iy,ieta).Wmunu    [0][2] = temp_profile_pitauy[ix][iy];
            arena(ix,iy,ieta).prevWmunu[0][3] = temp_profile_pitaueta[ix][iy];
            arena(ix,iy,ieta).Wmunu    [0][3] = temp_profile_pitaueta[ix][iy];
            arena(ix,iy,ieta).prevWmunu[0][4] = temp_profile_pixx[ix][iy];
            arena(ix,iy,ieta).Wmunu    [0][4] = temp_profile_pixx[ix][iy];
            arena(ix,iy,ieta).prevWmunu[0][5] = temp_profile_pixy[ix][iy];
            arena(ix,iy,ieta).Wmunu    [0][5] = temp_profile_pixy[ix][iy];
            arena(ix,iy,ieta).prevWmunu[0][6] = temp_profile_pixeta[ix][iy];
            arena(ix,iy,ieta).Wmunu    [0][6] = temp_profile_pixeta[ix][iy];
            arena(ix,iy,ieta).prevWmunu[0][7] = temp_profile_piyy[ix][iy];
            arena(ix,iy,ieta).Wmunu    [0][7] = temp_profile_piyy[ix][iy];
            arena(ix,iy,ieta).prevWmunu[0][8] = temp_profile_piyeta[ix][iy];
            arena(ix,iy,ieta).Wmunu    [0][8] = temp_profile_piyeta[ix][iy];
            arena(ix,iy,ieta).prevWmunu[0][9] = temp_profile_pietaeta[ix][iy];
            arena(ix,iy,ieta).Wmunu    [0][9] = temp_profile_pietaeta[ix][iy];
        }
    }
    // clean up
    for (int i = 0; i < nx; i++) {
        delete[] temp_profile_ed[i];
        delete[] temp_profile_utau[i];
        delete[] temp_profile_ux[i];
        delete[] temp_profile_uy[i];
        delete[] temp_profile_ueta[i];
        delete[] temp_profile_pitautau[i];
        delete[] temp_profile_pitaux[i];
        delete[] temp_profile_pitauy[i];
        delete[] temp_profile_pitaueta[i];
        delete[] temp_profile_pixx[i];
        delete[] temp_profile_pixy[i];
        delete[] temp_profile_pixeta[i];
        delete[] temp_profile_piyy[i];
        delete[] temp_profile_piyeta[i];
        delete[] temp_profile_pietaeta[i];
    }
    delete[] temp_profile_ed;
    delete[] temp_profile_utau;
    delete[] temp_profile_ux;
    delete[] temp_profile_uy;
    delete[] temp_profile_ueta;
    delete[] temp_profile_pitautau;
    delete[] temp_profile_pitaux;
    delete[] temp_profile_pitauy;
    delete[] temp_profile_pitaueta;
    delete[] temp_profile_pixx;
    delete[] temp_profile_pixy;
    delete[] temp_profile_pixeta;
    delete[] temp_profile_piyy;
    delete[] temp_profile_piyeta;
    delete[] temp_profile_pietaeta;
}

void Init::initial_MCGlb_with_rhob_XY(InitData *DATA, int ieta,
                                      Grid &arena) {
    // first load in the transverse profile
    ifstream profile_TA(DATA->initName_TA.c_str());
    ifstream profile_TB(DATA->initName_TB.c_str());
    ifstream profile_rhob_TA(DATA->initName_rhob_TA.c_str());
    ifstream profile_rhob_TB(DATA->initName_rhob_TB.c_str());

    const int nx = arena.nX();
    const int ny = arena.nY();
    double** temp_profile_TA = new double* [nx];
    double** temp_profile_TB = new double* [nx];
    double** temp_profile_rhob_TA = new double* [nx];
    double** temp_profile_rhob_TB = new double* [nx];
    for (int i = 0; i < nx; i++) {
        temp_profile_TA[i] = new double[ny];
        temp_profile_TB[i] = new double[ny];
        temp_profile_rhob_TA[i] = new double[ny];
        temp_profile_rhob_TB[i] = new double[ny];
    }
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            profile_TA >> temp_profile_TA[i][j];
            profile_TB >> temp_profile_TB[i][j];
            profile_rhob_TA >> temp_profile_rhob_TA[i][j];
            profile_rhob_TB >> temp_profile_rhob_TB[i][j];
        }
    }
    profile_TA.close();
    profile_TB.close();
    profile_rhob_TA.close();
    profile_rhob_TB.close();

    double eta = (DATA->delta_eta)*ieta - (DATA->eta_size)/2.0;
    double eta_envelop_left  = eta_profile_left_factor(DATA, eta);
    double eta_envelop_right = eta_profile_right_factor(DATA, eta);
    double eta_rhob_left     = eta_rhob_left_factor(DATA, eta);
    double eta_rhob_right    = eta_rhob_right_factor(DATA, eta);

    int entropy_flag = DATA->initializeEntropy;
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy< ny; iy++) {
            double rhob = 0.0;
            double epsilon = 0.0;
            if (DATA->turn_on_rhob == 1) {
                rhob = (
                    (temp_profile_rhob_TA[ix][iy]*eta_rhob_left
                     + temp_profile_rhob_TB[ix][iy]*eta_rhob_right));
            } else {
                rhob = 0.0;
            }
            if (entropy_flag == 0) {
                epsilon = (
                    (temp_profile_TA[ix][iy]*eta_envelop_left
                     + temp_profile_TB[ix][iy]*eta_envelop_right)
                    *DATA->sFactor/hbarc);   // 1/fm^4
            } else {
                double local_sd = (
                    (temp_profile_TA[ix][iy]*eta_envelop_left
                     + temp_profile_TB[ix][iy]*eta_envelop_right)
                    *DATA->sFactor);         // 1/fm^3
                epsilon = eos->get_s2e(local_sd, rhob);
            }
            if (epsilon < 0.00000000001)
                epsilon = 0.00000000001;

            // set all values in the grid element:
            arena(ix,iy,ieta).epsilon      = epsilon;
            arena(ix,iy,ieta).epsilon_t    = epsilon;
            arena(ix,iy,ieta).prev_epsilon = epsilon;
            arena(ix,iy,ieta).rhob         = rhob;
            arena(ix,iy,ieta).rhob_t       = rhob;
            arena(ix,iy,ieta).prev_rhob    = rhob;

            /* for HIC */
            arena(ix,iy,ieta).u[0][0] = 1.0;
            arena(ix,iy,ieta).u[0][3] = 0.0;
            arena(ix,iy,ieta).u[0][1] = 0.0;
            arena(ix,iy,ieta).u[0][2] = 0.0;
            
            for (int ii = 0; ii < 1; ii++) {
                arena(ix,iy,ieta).prev_u[ii][0] = 1.0;
                arena(ix,iy,ieta).prev_u[ii][3] = 0.0;
                arena(ix,iy,ieta).prev_u[ii][1] = 0.0;
                arena(ix,iy,ieta).prev_u[ii][2] = 0.0;
                arena(ix,iy,ieta).prev_pi_b[ii] = 0.0;
            }

            arena(ix,iy,ieta).pi_b[0] = 0.0;

            for (int ii = 0; ii < 14; ii++) {
                arena(ix,iy,ieta).prevWmunu[0][ii] = 0.0;
                arena(ix,iy,ieta).Wmunu[0][ii] = 0.0;
            }
        }
    }
    // clean up
    for (int i = 0; i < nx; i++) {
        delete[] temp_profile_TA[i];
        delete[] temp_profile_TB[i];
        delete[] temp_profile_rhob_TA[i];
        delete[] temp_profile_rhob_TB[i];
    }
    delete[] temp_profile_TA;
    delete[] temp_profile_TB;
    delete[] temp_profile_rhob_TA;
    delete[] temp_profile_rhob_TB;
}

void Init::initial_MCGlbLEXUS_with_rhob_XY(InitData *DATA, int ieta,
                                           Grid &arena) {
    const int nx = arena.nX();
    const int ny = arena.nY();
    double *u = new double[4];
    u[0] = 1.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = 0.0;
    double *j_mu = new double[4];
    for (int i = 0; i < 3; i++) {
        j_mu[i] = 0.0;
    }
    int entropy_flag = DATA->initializeEntropy;
    for (int ix = 0; ix < nx; ix++) {
        // double x_local = - DATA->x_size/2. + ix*DATA->delta_x;
        for (int iy = 0; iy < ny; iy++) {
            // double y_local = - DATA->y_size/2. + iy*DATA->delta_y;
            double rhob = 0.0;
            double epsilon = 0.0;
            if (DATA->turn_on_rhob == 1) {
                //rhob = hydro_source_ptr->get_hydro_rhob_source(
                //                            tau0, x_local, y_local, eta, u);
                rhob = 0.0;
            } else {
                rhob = 0.0;
            }

            //hydro_source_ptr->get_hydro_energy_source(
            //                        tau0, x_local, y_local, eta, u, j_mu);

            if (entropy_flag == 0) {
                // epsilon = j_mu[0];           // 1/fm^4
                epsilon = 0.0;           // 1/fm^4
            } else {
                //double local_sd = j_mu[0]*DATA->sFactor;         // 1/fm^3
                //epsilon = eos->get_s2e(local_sd, rhob);
                epsilon = 0.0;           // 1/fm^4
            }

            if (epsilon < 0.00000000001)
                epsilon = 0.00000000001;

            // set all values in the grid element:
            arena(ix,iy,ieta).epsilon      = epsilon;
            arena(ix,iy,ieta).epsilon_t    = epsilon;
            arena(ix,iy,ieta).prev_epsilon = epsilon;
            arena(ix,iy,ieta).rhob         = rhob;
            arena(ix,iy,ieta).rhob_t       = rhob;
            arena(ix,iy,ieta).prev_rhob    = rhob;

            /* for HIC */
            arena(ix,iy,ieta).u[0][0] = u[0];
            arena(ix,iy,ieta).u[0][1] = u[1];
            arena(ix,iy,ieta).u[0][2] = u[2];
            arena(ix,iy,ieta).u[0][3] = u[3];
            
            for (int ii = 0; ii < 1; ii++) {
                arena(ix,iy,ieta).prev_u[ii][0] = 1.0;
                arena(ix,iy,ieta).prev_u[ii][1] = 0.0;
                arena(ix,iy,ieta).prev_u[ii][2] = 0.0;
                arena(ix,iy,ieta).prev_u[ii][3] = 0.0;
                arena(ix,iy,ieta).prev_pi_b[ii] = 0.0;
            }

            arena(ix,iy,ieta).pi_b[0] = 0.0;

            for (int ii = 0; ii < 14; ii++) {
                arena(ix,iy,ieta).prevWmunu[0][ii] = 0.0;
                arena(ix,iy,ieta).Wmunu[0][ii] = 0.0;
            }
        }
    }
    delete[] u;
    delete[] j_mu;
}


void Init::initial_UMN_with_rhob(InitData *DATA, Grid &arena) {
    // first load in the transverse profile
    ifstream profile(DATA->initName.c_str());
    ifstream profile_TA(DATA->initName_TA.c_str());
    ifstream profile_TB(DATA->initName_TB.c_str());

    const int nx = arena.nX();
    const int ny = arena.nY();
    double** temp_profile_TA = new double* [nx];
    double** temp_profile_TB = new double* [nx];
    for (int i = 0; i < nx+1; i++) {
        temp_profile_TA[i] = new double[ny];
        temp_profile_TB[i] = new double[ny];
    }
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            profile_TA >> temp_profile_TA[i][j];
            profile_TB >> temp_profile_TB[i][j];
        }
    }
    profile_TA.close();
    profile_TB.close();

    double dummy;
    double ed_local, rhob_local;
    for (int ieta = 0; ieta < DATA->neta; ieta++) {
        double eta = (DATA->delta_eta)*ieta - (DATA->eta_size)/2.0;
        double eta_envelop_left = eta_profile_left_factor(DATA, eta);
        double eta_envelop_right = eta_profile_right_factor(DATA, eta);
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy< ny; iy++) {
                double rhob = 0.0;
                double epsilon = 0.0;
                profile >> dummy >> dummy >> dummy >> ed_local >> rhob_local;
                rhob = rhob_local;
                double sd_bg = (temp_profile_TA[ix][iy]*eta_envelop_left
                    + temp_profile_TB[ix][iy]*eta_envelop_right)*DATA->sFactor;
                double ed_bg = eos->get_s2e(sd_bg, 0.0);
                epsilon = ed_bg + ed_local/hbarc;    // 1/fm^4
                if (epsilon < 0.00000000001) {
                    epsilon = 0.00000000001;
                }

                // set all values in the grid element:
                arena(ix,iy,ieta).epsilon      = epsilon;
                arena(ix,iy,ieta).epsilon_t    = epsilon;
                arena(ix,iy,ieta).prev_epsilon = epsilon;
                arena(ix,iy,ieta).rhob         = rhob;
                arena(ix,iy,ieta).rhob_t       = rhob;
                arena(ix,iy,ieta).prev_rhob    = rhob;

                arena(ix,iy,ieta).u[0][0] = 1.0;
                arena(ix,iy,ieta).u[0][3] = 0.0;
                arena(ix,iy,ieta).u[0][1] = 0.0;
                arena(ix,iy,ieta).u[0][2] = 0.0;
            
                for (int ii = 0; ii < 1; ii++) {
                    arena(ix,iy,ieta).prev_u[ii][0] = 1.0;
                    arena(ix,iy,ieta).prev_u[ii][3] = 0.0;
                    arena(ix,iy,ieta).prev_u[ii][1] = 0.0;
                    arena(ix,iy,ieta).prev_u[ii][2] = 0.0;
                    arena(ix,iy,ieta).prev_pi_b[ii] = 0.0;
                }

                arena(ix,iy,ieta).pi_b[0] = 0.0;

                for (int ii = 0; ii < 14; ii++) {
                    arena(ix,iy,ieta).prevWmunu[0][ii] = 0.0;
                    arena(ix,iy,ieta).Wmunu[0][ii] = 0.0;
                }
            }
        }
    }
    profile.close();
    // clean up
    for (int i = 0; i < nx; i++) {
        delete[] temp_profile_TA[i];
        delete[] temp_profile_TB[i];
    }
    delete[] temp_profile_TA;
    delete[] temp_profile_TB;
}

void Init::initial_AMPT_XY(InitData *DATA, int ieta, Grid &arena) {
    double *u = new double[4];
    u[0] = 1.0;
    u[1] = 0.0;
    u[2] = 0.0;
    u[3] = 0.0;
    double *j_mu = new double[4];
    for (int i = 0; i < 3; i++) {
        j_mu[i] = 0.0;
    }

    double eta = (DATA->delta_eta)*ieta - (DATA->eta_size)/2.0;
    double tau0 = DATA->tau0;
    const int nx = arena.nX();
    const int ny = arena.nY();
    for (int ix = 0; ix < nx; ix++) {
        double x_local = - DATA->x_size/2. + ix*DATA->delta_x;
        for (int iy = 0; iy < ny; iy++) {
            double y_local = - DATA->y_size/2. + iy*DATA->delta_y;
            double rhob = 0.0;
            double epsilon = 0.0;
            if (DATA->turn_on_rhob == 1) {
                rhob = hydro_source_ptr->get_hydro_rhob_source_before_tau(
                                                tau0, x_local, y_local, eta);
            } else {
                rhob = 0.0;
            }

            hydro_source_ptr->get_hydro_energy_source_before_tau(
                                    tau0, x_local, y_local, eta, j_mu);

            epsilon = j_mu[0];           // 1/fm^4

            if (epsilon < 0.00000000001)
                epsilon = 0.00000000001;

            // set all values in the grid element:
            arena(ix,iy,ieta).epsilon      = epsilon;
            arena(ix,iy,ieta).epsilon_t    = epsilon;
            arena(ix,iy,ieta).prev_epsilon = epsilon;
            arena(ix,iy,ieta).rhob         = rhob;
            arena(ix,iy,ieta).rhob_t       = rhob;
            arena(ix,iy,ieta).prev_rhob    = rhob;

            /* for HIC */
            arena(ix,iy,ieta).u[0][0] = u[0];
            arena(ix,iy,ieta).u[0][3] = u[3];
            arena(ix,iy,ieta).u[0][1] = u[1];
            arena(ix,iy,ieta).u[0][2] = u[2];
            
            for (int ii = 0; ii < 1; ii++) {
                arena(ix,iy,ieta).prev_u[ii][0] = 1.0;
                arena(ix,iy,ieta).prev_u[ii][3] = 0.0;
                arena(ix,iy,ieta).prev_u[ii][1] = 0.0;
                arena(ix,iy,ieta).prev_u[ii][2] = 0.0;
                arena(ix,iy,ieta).prev_pi_b[ii] = 0.0;
            }

            arena(ix,iy,ieta).pi_b[0] = 0.0;

            for (int ii = 0; ii < 14; ii++) {
                arena(ix,iy,ieta).prevWmunu[0][ii] = 0.0;
                arena(ix,iy,ieta).Wmunu[0][ii] = 0.0;
            }
        }
    }
    delete[] u;
    delete[] j_mu;
}

double Init::eta_profile_normalisation(InitData *DATA, double eta) {
    // this function return the eta envelope profile for energy density
    double res;
    // Hirano's plateau + Gaussian fall-off
    if (DATA->initial_eta_profile == 1) {
        double exparg1, exparg;
        exparg1 = (fabs(eta) - DATA->eta_flat/2.0)/DATA->eta_fall_off;
        exparg = exparg1*exparg1/2.0;
        res = exp(-exparg*theta(exparg1));
    } else if (DATA->initial_eta_profile == 2) {
        // Woods-Saxon
        // The radius is set to be half of DATA->eta_flat
        // The diffusiveness is set to DATA->eta_fall_off
        double ws_R = DATA->eta_flat/2.0;
        double ws_a = DATA->eta_fall_off;
        res = (1.0 + exp(-ws_R/ws_a))/(1.0 + exp((abs(eta) - ws_R)/ws_a));
    } else {
        music_message.error("initial_eta_profile out of range.");
        exit(1);
    }
    return res;
}

double Init::eta_profile_left_factor(InitData *Data, double eta) {
    // this function return the eta envelope for projectile
    double res = eta_profile_normalisation(Data, eta);
    if (fabs(eta) < Data->beam_rapidity) {
        res = (1. - eta/Data->beam_rapidity)*res;
    } else {
        res = 0.0;
    }
    return(res);
}

double Init::eta_profile_right_factor(InitData *Data, double eta) {
    // this function return the eta envelope for target
    double res = eta_profile_normalisation(Data, eta);
    if (fabs(eta) < Data->beam_rapidity) {
        res = (1. + eta/Data->beam_rapidity)*res;
    } else {
        res = 0.0;
    }
    return(res);
}

double Init::eta_rhob_profile_normalisation(InitData *DATA, double eta) {
    // this function return the eta envelope profile for net baryon density
    double res;
    int profile_flag = DATA->initial_eta_rhob_profile;
    double eta_0 = DATA->eta_rhob_0;
    double tau0 = DATA->tau0;
    if (profile_flag == 1) {
        const double eta_width = DATA->eta_rhob_width;
        const double norm      = 1./(2.*sqrt(2*M_PI)*eta_width*tau0);
        const double exparg1   = (eta - eta_0)/eta_width;
        const double exparg2   = (eta + eta_0)/eta_width;
        res = norm*(exp(-exparg1*exparg1/2.0) + exp(-exparg2*exparg2/2.0));
    } else if (profile_flag == 2) {
        double eta_abs     = fabs(eta);
        double delta_eta_1 = DATA->eta_rhob_width_1;
        double delta_eta_2 = DATA->eta_rhob_width_2;
        double A           = DATA->eta_rhob_plateau_height;
        double exparg1     = (eta_abs - eta_0)/delta_eta_1;
        double exparg2     = (eta_abs - eta_0)/delta_eta_2;
        double theta;
        double norm = 1./(tau0*(sqrt(2.*M_PI)*delta_eta_1
                          + (1. - A)*sqrt(2.*M_PI)*delta_eta_2 + 2.*A*eta_0));
        if (eta_abs > eta_0)
            theta = 1.0;
        else
            theta = 0.0;
        res = norm*(theta*exp(-exparg1*exparg1/2.)
                    + (1. - theta)*(A + (1. - A)*exp(-exparg2*exparg2/2.)));
    } else {
        music_message << "initial_eta_rhob_profile = " << profile_flag
                      << " out of range.";
        music_message.flush("error");
        exit(1);
    }
    return res;
}

double Init::eta_rhob_left_factor(InitData *DATA, double eta) {
    double eta_0       = -fabs(DATA->eta_rhob_0);
    double tau0        = DATA->tau0;
    double delta_eta_1 = DATA->eta_rhob_width_1;
    double delta_eta_2 = DATA->eta_rhob_width_2;
    double norm        = 2./(sqrt(M_PI)*tau0*(delta_eta_1 + delta_eta_2));
    double exp_arg     = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_1;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_2;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}

double Init::eta_rhob_right_factor(InitData *DATA, double eta) {
    double eta_0       = fabs(DATA->eta_rhob_0);
    double tau0        = DATA->tau0;
    double delta_eta_1 = DATA->eta_rhob_width_1;
    double delta_eta_2 = DATA->eta_rhob_width_2;
    double norm        = 2./(sqrt(M_PI)*tau0*(delta_eta_1 + delta_eta_2));
    double exp_arg     = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_2;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_1;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}

void Init::output_initial_density_profiles(InitData *DATA, Grid &arena) {
    // this function outputs the 3d initial energy density profile
    // and net baryon density profile (if turn_on_rhob == 1)
    // for checking purpose
    music_message.info("output initial density profiles into a file... ");
    ofstream of("check_initial_density_profiles.dat");
    of << "# x(fm)  y(fm)  eta  ed(GeV/fm^3)";
    if (DATA->turn_on_rhob == 1)
        of << "  rhob(1/fm^3)";
    of << endl;
    for (int ieta = 0; ieta < arena.nEta(); ieta++) {
        double eta_local = (DATA->delta_eta)*ieta - (DATA->eta_size)/2.0;
        for(int ix = 0; ix < arena.nX(); ix++) {
            double x_local = -DATA->x_size/2. + ix*DATA->delta_x;
            for(int iy = 0; iy < arena.nY(); iy++) {
                double y_local = -DATA->y_size/2. + iy*DATA->delta_y;
                of << scientific << setw(18) << std::setprecision(8)
                   << x_local << "   " << y_local << "   "
                   << eta_local << "   " << arena(ix,iy,ieta).epsilon*hbarc;
                if (DATA->turn_on_rhob == 1) {
                    of << "   " << arena(ix,iy,ieta).rhob;
                }
                of << endl;
            }
        }
    }
    music_message.info("done!");
}
