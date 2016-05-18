// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <omp.h>
#include "./util.h"
#include "./grid.h"
#include "./init.h"
#include "./eos.h"

using namespace std;

Init::Init(EOS *eosIn) {
    eos = eosIn;
    util = new Util;
}

// destructor
Init::~Init() {
    delete util;
}

void Init::InitArena(InitData *DATA, Grid ****arena) {
    Grid *helperGrid;
    helperGrid = new Grid;
    cout << "initArena" << endl;
    if (DATA->Initial_profile == 8) {
        cout << DATA->initName <<endl;
        ifstream profile(DATA->initName.c_str());
        string dummy;
        int nx, ny, neta;
        double deta, dx, dy, dummy2;
        // read the first line with general info
        profile >> dummy >> dummy >> dummy2
                >> dummy >> neta >> dummy >> nx >> dummy >> ny
                >> dummy >> deta >> dummy >> dx >> dummy >> dy;
        profile.close();
        cout << "Using Initial_profile=" << DATA->Initial_profile
             << ". Overwriting lattice dimensions:" << endl;

        DATA->nx = nx-1;
        DATA->ny = ny-1;
        DATA->delta_x = dx;
        DATA->delta_y = dy;

        cout << "neta=" << neta << ", nx=" << nx << ", ny=" << ny << endl;
        cout << "deta=" << DATA->delta_eta << ", dx=" << DATA->delta_x
             << ", dy=" << DATA->delta_y << endl;
    }

    // initialize arena
    *arena = helperGrid->grid_c_malloc(DATA->neta, DATA->nx+1, DATA->ny+1);

    cout << "Grid allocated." << endl;
    InitTJb(DATA, arena);

    LinkNeighbors(DATA, arena);
    delete helperGrid;
}/* InitArena */


void Init::LinkNeighbors(InitData *DATA, Grid ****arena) {
    int nx = DATA->nx;
    int ny = DATA->ny;
    int neta = DATA->neta;

    /* allocate memory */
    for (int ieta = 0; ieta < neta; ieta++) {
        for (int ix = 0; ix <= nx; ix++) {
            for (int iy = 0; iy <= ny; iy++) {
                (*arena)[ieta][ix][iy].nbr_p_1 = new Grid *[4];
                (*arena)[ieta][ix][iy].nbr_m_1 = new Grid *[4];
                (*arena)[ieta][ix][iy].nbr_p_2 = new Grid *[4];
                (*arena)[ieta][ix][iy].nbr_m_2 = new Grid *[4];

                (*arena)[ieta][ix][iy].position[1] = ix;
                (*arena)[ieta][ix][iy].position[2] = iy;
                (*arena)[ieta][ix][iy].position[3] = ieta;
            }
        }
    }

    int ieta;
    #pragma omp parallel private(ieta)
    {
        #pragma omp for
        for (ieta = 0; ieta < neta; ieta++) {
            //printf("Thread %d executes loop iteraction %d\n",
            //       omp_get_thread_num(), ieta);
            for (int ix = 0; ix <= nx; ix++) {
                for (int iy = 0; iy <= ny; iy++) {
                    if (ix != nx)
                        (*arena)[ieta][ix][iy].nbr_p_1[1] =
                                                    &(*arena)[ieta][ix+1][iy];
                    else
                        (*arena)[ieta][ix][iy].nbr_p_1[1] = 
                                                    &(*arena)[ieta][nx][iy];
                    if (ix < nx - 1)
                        (*arena)[ieta][ix][iy].nbr_p_2[1] =
                                                    &(*arena)[ieta][ix+2][iy];
                    else
                        (*arena)[ieta][ix][iy].nbr_p_2[1] = 
                                                    &(*arena)[ieta][nx][iy];
                    if (ix != 0)
                        (*arena)[ieta][ix][iy].nbr_m_1[1] =
                                                    &(*arena)[ieta][ix-1][iy];
                    else
                        (*arena)[ieta][ix][iy].nbr_m_1[1] =
                                                    &(*arena)[ieta][0][iy];
                    if (ix > 1)
                        (*arena)[ieta][ix][iy].nbr_m_2[1] =
                                                    &(*arena)[ieta][ix-2][iy];
                    else
                        (*arena)[ieta][ix][iy].nbr_m_2[1] = 
                                                    &(*arena)[ieta][0][iy];
                    if (iy != ny)
                        (*arena)[ieta][ix][iy].nbr_p_1[2] =
                                                    &(*arena)[ieta][ix][iy+1];
                    else
                        (*arena)[ieta][ix][iy].nbr_p_1[2] = 
                                                    &(*arena)[ieta][ix][ny];
                    if (iy < ny - 1)
                        (*arena)[ieta][ix][iy].nbr_p_2[2] =
                                                    &(*arena)[ieta][ix][iy+2];
                    else
                        (*arena)[ieta][ix][iy].nbr_p_2[2] =
                                                    &(*arena)[ieta][ix][ny];
                    if (iy != 0)
                        (*arena)[ieta][ix][iy].nbr_m_1[2] =
                                                    &(*arena)[ieta][ix][iy-1];
                    else
                        (*arena)[ieta][ix][iy].nbr_m_1[2] =
                                                    &(*arena)[ieta][ix][0];
                    if (iy > 1)
                        (*arena)[ieta][ix][iy].nbr_m_2[2] =
                                                    &(*arena)[ieta][ix][iy-2];
                    else
                        (*arena)[ieta][ix][iy].nbr_m_2[2] =
                                                    &(*arena)[ieta][ix][0];

                    if (ieta != neta-1)
                        (*arena)[ieta][ix][iy].nbr_p_1[3] =
                                                    &(*arena)[ieta+1][ix][iy];
                    else
                        (*arena)[ieta][ix][iy].nbr_p_1[3] =
                                                    &(*arena)[neta-1][ix][iy];
                    if (ieta < neta-2)
                        (*arena)[ieta][ix][iy].nbr_p_2[3] =
                                                    &(*arena)[ieta+2][ix][iy];
                    else
                        (*arena)[ieta][ix][iy].nbr_p_2[3] =
                                                    &(*arena)[neta-1][ix][iy];
                    if (ieta != 0)
                        (*arena)[ieta][ix][iy].nbr_m_1[3] =
                                                    &(*arena)[ieta-1][ix][iy];
                    else
                        (*arena)[ieta][ix][iy].nbr_m_1[3] =
                                                    &(*arena)[0][ix][iy];
                    if (ieta > 1)
                        (*arena)[ieta][ix][iy].nbr_m_2[3] =
                                                    &(*arena)[ieta-2][ix][iy];
                    else
                        (*arena)[ieta][ix][iy].nbr_m_2[3] =
                                                    &(*arena)[0][ix][iy];
                }
            }
        }
    }
}  /* LinkNeighbors */

int Init::InitTJb(InitData *DATA, Grid ****arena) {
    double epsilon = 0.0;
    double rhob = 0.0;
    double u[4];
    int rk_order = DATA->rk_order;
    cout << "rk_order=" << rk_order << endl;
    if (DATA->Initial_profile == 0) {
    } else if (DATA->Initial_profile == 8) {
        // read in the profile from file
        // - IPGlasma initial conditions with initial flow
        cout << " ----- information on initial distribution -----" << endl;
        cout << "file name used: " << DATA->initName << endl;
  
        ifstream profile(DATA->initName.c_str());

        string dummy;
        int nx, ny, neta;
        double dx, dy, deta;
        // read the information line
        profile >> dummy >> dummy >> dummy >> dummy >> neta
                >> dummy >> nx >> dummy >> ny
                >> dummy >> deta >> dummy >> dx >> dummy >> dy;

        cout << "neta=" << DATA->neta << ", nx=" << nx << ", ny=" << ny
             << ", deta=" << DATA->delta_eta << ", dx=" << dx << ", dy=" << dy
             << endl;

        double density, dummy1, dummy2, dummy3;
        double ux, uy, utau;

        double** temp_profile_ed = new double* [nx+1];
        double** temp_profile_utau = new double* [nx+1];
        double** temp_profile_ux = new double* [nx+1];
        double** temp_profile_uy = new double* [nx+1];
        for (int i = 0; i < nx+1; i++) {
            temp_profile_ed[i] = new double[ny+1];
            temp_profile_utau[i] = new double[ny+1];
            temp_profile_ux[i] = new double[ny+1];
            temp_profile_uy[i] = new double[ny+1];
        }

        // read the one slice
        for (int ix = 0; ix <= DATA->nx; ix++) {
            for (int iy = 0; iy <= DATA->ny; iy++) {
                profile >> dummy1 >> dummy2 >> dummy3
                        >> density >> utau >> ux >> uy
                        >> dummy  >> dummy  >> dummy  >> dummy;
                temp_profile_ed[ix][iy] = density;
                temp_profile_utau[ix][iy] = utau;
                temp_profile_ux[ix][iy] = ux;
                temp_profile_uy[ix][iy] = uy;
                if (ix == 0 && iy == 0) {
                    DATA->x_size = -dummy2*2;
                    DATA->y_size = -dummy3*2;
                    cout << "eta_size=" << DATA->eta_size
                         << ", x_size=" << DATA->x_size
                         << ", y_size=" << DATA->y_size << endl;
                }
            }
        }
        profile.close();

        int entropy_flag = DATA->initializeEntropy;
        int ieta;
        #pragma omp parallel private(ieta)
        {
            #pragma omp for
            for (ieta = 0; ieta < DATA->neta; ieta++) {
                double eta = (DATA->delta_eta)*(ieta) - (DATA->eta_size)/2.0;
                double eta_envelop_ed = eta_profile_normalisation(DATA, eta);
                for (int ix = 0; ix <= DATA->nx; ix++) {
                    for (int iy = 0; iy<= DATA->ny; iy++) {
                        rhob = 0.0;
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

                        // initial pressure distribution
                        double p = eos->get_pressure(epsilon, rhob);
                        // set all values in the grid element:
                        (*arena)[ieta][ix][iy].epsilon = epsilon;
                        (*arena)[ieta][ix][iy].epsilon_t = epsilon;
                        (*arena)[ieta][ix][iy].prev_epsilon = epsilon;
                        (*arena)[ieta][ix][iy].rhob = rhob;
                        (*arena)[ieta][ix][iy].rhob_t = rhob;
                        (*arena)[ieta][ix][iy].prev_rhob = rhob;
                        (*arena)[ieta][ix][iy].p = p;
                        (*arena)[ieta][ix][iy].p_t = p;
                        (*arena)[ieta][ix][iy].trouble = 0;

                        (*arena)[ieta][ix][iy].T = eos->get_temperature(
                                                                epsilon, rhob);
                        (*arena)[ieta][ix][iy].mu = eos->get_mu(epsilon, rhob);
                    
                        (*arena)[ieta][ix][iy].TJb =
                                        util->cube_malloc(rk_order+1, 5, 4);
                        (*arena)[ieta][ix][iy].dUsup =
                                        util->cube_malloc(rk_order+1, 5, 4);
                        (*arena)[ieta][ix][iy].u = 
                                        util->mtx_malloc(rk_order+1, 4);
                        (*arena)[ieta][ix][iy].a = 
                                        util->mtx_malloc(rk_order+1, 5);
                        (*arena)[ieta][ix][iy].theta_u =
                                        util->vector_malloc(rk_order+1);
                        (*arena)[ieta][ix][iy].sigma =
                                        util->cube_malloc(rk_order+1, 4, 4);
                        (*arena)[ieta][ix][iy].pi_b =
                                        util->vector_malloc(rk_order+1);
                        (*arena)[ieta][ix][iy].prev_u = util->mtx_malloc(1, 4);
                        (*arena)[ieta][ix][iy].Wmunu =
                                        util->cube_malloc(rk_order+1, 5, 4);
                        (*arena)[ieta][ix][iy].prevWmunu =
                                        util->cube_malloc(rk_order, 5, 4);
                        (*arena)[ieta][ix][iy].Pimunu =
                                        util->cube_malloc(rk_order+1, 5, 4);
                        (*arena)[ieta][ix][iy].prevPimunu =
                                        util->cube_malloc(rk_order, 5, 4);
                        (*arena)[ieta][ix][iy].W_prev = util->mtx_malloc(5, 4);

                        /* for HIC */
                        (*arena)[ieta][ix][iy].u[0][0] =
                                                    temp_profile_utau[ix][iy];
                        (*arena)[ieta][ix][iy].u[0][1] =
                                                    temp_profile_ux[ix][iy];
                        (*arena)[ieta][ix][iy].u[0][2] =
                                                    temp_profile_uy[ix][iy];
                        (*arena)[ieta][ix][iy].u[0][3] = 0.0;

                        u[0] = temp_profile_utau[ix][iy];
                        u[1] = temp_profile_ux[ix][iy];
                        u[2] = temp_profile_uy[ix][iy];
                        u[3] = 0.0;

                        (*arena)[ieta][ix][iy].prev_u[0][0] = u[0];
                        (*arena)[ieta][ix][iy].prev_u[0][1] = u[1];
                        (*arena)[ieta][ix][iy].prev_u[0][2] = u[2];
                        (*arena)[ieta][ix][iy].prev_u[0][3] = u[3];

                        (*arena)[ieta][ix][iy].pi_b[0] = 0.0;

                        for (int mu = 0; mu < 4; mu++) {
                            /* baryon density */
                            (*arena)[ieta][ix][iy].TJb[0][4][mu] = rhob*u[mu];

                            // diffusion current
                            (*arena)[ieta][ix][iy].Wmunu[0][4][mu] = 0.0;
                            (*arena)[ieta][ix][iy].prevWmunu[0][4][mu] = 0.0;
                        
                            for (int nu = 0; nu < 4; nu++) {
                                (*arena)[ieta][ix][iy].TJb[0][nu][mu] = (
                                    (epsilon + p)*u[mu]*u[nu]
                                    + p*(DATA->gmunu)[mu][nu]);
                                (*arena)[ieta][ix][iy].Wmunu[0][nu][mu] = 0.0;
                                (*arena)[ieta][ix][iy].prevWmunu[0][nu][mu] = 0.0;

                                (*arena)[ieta][ix][iy].Pimunu[0][nu][mu] = 0.0;
                                (*arena)[ieta][ix][iy].prevPimunu[0][nu][mu] = 0.0;
                            }/* nu */
                        }/* mu */
                    }
                }
            }/* ix, iy, ieta */
        }
        // clean up
        for (int i = 0; i < nx+1; i++) {
            delete[] temp_profile_ed[i];
            delete[] temp_profile_utau[i];
            delete[] temp_profile_ux[i];
            delete[] temp_profile_uy[i];
        }
        delete[] temp_profile_ed;
        delete[] temp_profile_utau;
        delete[] temp_profile_ux;
        delete[] temp_profile_uy;
    } else if (DATA->Initial_profile == 11) {
        // read in the transverse profile from file with finite rho_B
        // the initial entropy and net baryon density profile are
        // constructed by nuclear thickness function TA and TB.
        // Along the longitudinal direction an asymmetric contribution from
        // target and projectile thickness function is allowed
        cout << " ----- information on initial distribution -----" << endl;
        cout << "file name used: " << DATA->initName_TA << " and "
             << DATA->initName_TB << endl;
        // first load in the transverse profile
        ifstream profile_TA(DATA->initName_TA.c_str());
        ifstream profile_TB(DATA->initName_TB.c_str());
        int nx = DATA->nx;
        int ny = DATA->ny;
        double** temp_profile_TA = new double* [nx+1];
        double** temp_profile_TB = new double* [nx+1];
        for (int i = 0; i < nx+1; i++) {
            temp_profile_TA[i] = new double[ny+1];
            temp_profile_TB[i] = new double[ny+1];
        }
        for (int i = 0; i < nx+1; i++) {
            for (int j = 0; j < ny+1; j++) {
                profile_TA >> temp_profile_TA[i][j];
                profile_TB >> temp_profile_TB[i][j];
            }
        }
        profile_TA.close();
        profile_TB.close();

        int entropy_flag = DATA->initializeEntropy;
        int ieta;
        #pragma omp parallel private(ieta)
        {
            #pragma omp for
            for (ieta = 0; ieta < DATA->neta; ieta++) {
                //printf("Thread %d executes loop iteraction %d\n",
                //       omp_get_thread_num(), ieta);
                double eta = (DATA->delta_eta)*ieta - (DATA->eta_size)/2.0;
                double eta_envelop_left = eta_profile_left_factor(DATA, eta);
                double eta_envelop_right = eta_profile_right_factor(DATA, eta);
                double eta_rhob_left = eta_rhob_left_factor(DATA, eta);
                double eta_rhob_right = eta_rhob_right_factor(DATA, eta);
                for (int ix = 0; ix < (DATA->nx+1); ix++) {
                    for (int iy = 0; iy< (DATA->ny+1); iy++) {
                        if (DATA->turn_on_rhob == 1) {
                            rhob = (
                                (temp_profile_TA[ix][iy]*eta_rhob_left
                                 + temp_profile_TB[ix][iy]*eta_rhob_right));
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

                        // initial pressure distribution
                        double p = eos->get_pressure(epsilon, rhob);

                        // set all values in the grid element:
                        (*arena)[ieta][ix][iy].epsilon = epsilon;
                        (*arena)[ieta][ix][iy].epsilon_t = epsilon;
                        (*arena)[ieta][ix][iy].prev_epsilon = epsilon;
                        (*arena)[ieta][ix][iy].rhob = rhob;
                        (*arena)[ieta][ix][iy].rhob_t = rhob;
                        (*arena)[ieta][ix][iy].prev_rhob = rhob;
                        (*arena)[ieta][ix][iy].p = p;
                        (*arena)[ieta][ix][iy].p_t = p;
                        (*arena)[ieta][ix][iy].trouble = 0;
                        (*arena)[ieta][ix][iy].T =
                                        eos->get_temperature(epsilon, rhob);
                        (*arena)[ieta][ix][iy].mu = eos->get_mu(epsilon, rhob);
                        (*arena)[ieta][ix][iy].TJb =
                                        util->cube_malloc(rk_order+1, 5, 4);
                        (*arena)[ieta][ix][iy].dUsup =
                                        util->cube_malloc(rk_order+1, 5, 4);
                        (*arena)[ieta][ix][iy].u =
                                        util->mtx_malloc(rk_order+1, 4);
                        (*arena)[ieta][ix][iy].a =
                                        util->mtx_malloc(rk_order+1, 5);
                        (*arena)[ieta][ix][iy].theta_u =
                                        util->vector_malloc(rk_order+1);
                        (*arena)[ieta][ix][iy].sigma =
                                        util->cube_malloc(rk_order+1, 4, 4);
                        (*arena)[ieta][ix][iy].pi_b =
                                        util->vector_malloc(rk_order+1);
                        (*arena)[ieta][ix][iy].prev_u = util->mtx_malloc(1, 4);
                        (*arena)[ieta][ix][iy].Wmunu =
                                        util->cube_malloc(rk_order+1, 5, 4);
                        (*arena)[ieta][ix][iy].prevWmunu =
                                        util->cube_malloc(rk_order, 5, 4);
                        (*arena)[ieta][ix][iy].Pimunu =
                                        util->cube_malloc(rk_order+1, 5, 4);
                        (*arena)[ieta][ix][iy].prevPimunu =
                                        util->cube_malloc(rk_order, 5, 4);
                        (*arena)[ieta][ix][iy].W_prev = util->mtx_malloc(5, 4);

                        /* for HIC */
                        u[0] = (*arena)[ieta][ix][iy].u[0][0] = 1.0;
                        u[3] = (*arena)[ieta][ix][iy].u[0][3] = 0.0;
                        u[1] = (*arena)[ieta][ix][iy].u[0][1] = 0.0;
                        u[2] = (*arena)[ieta][ix][iy].u[0][2] = 0.0;
                        (*arena)[ieta][ix][iy].prev_u[0][0] = 1.0;
                        (*arena)[ieta][ix][iy].prev_u[0][3] = 0.0;
                        (*arena)[ieta][ix][iy].prev_u[0][1] = 0.0;
                        (*arena)[ieta][ix][iy].prev_u[0][2] = 0.0;

                        (*arena)[ieta][ix][iy].pi_b[0] = 0.0;

                        for (int mu = 0; mu < 4; mu++) {
                            /* baryon density */
                            (*arena)[ieta][ix][iy].TJb[0][4][mu] = rhob*u[mu];

                            // diffusion current
                            (*arena)[ieta][ix][iy].Wmunu[0][4][mu] = 0.0;
                            (*arena)[ieta][ix][iy].prevWmunu[0][4][mu] = 0.0;
                            for (int nu = 0; nu < 4; nu++) {
                                (*arena)[ieta][ix][iy].TJb[0][nu][mu] = (
                                    (epsilon + p)*u[mu]*u[nu]
                                    + p*(DATA->gmunu)[mu][nu]);
                                (*arena)[ieta][ix][iy].Wmunu[0][nu][mu] = 0.0;
                                (*arena)[ieta][ix][iy].prevWmunu[0][nu][mu] = 0.0;

                                (*arena)[ieta][ix][iy].Pimunu[0][nu][mu] = 0.0;
                                (*arena)[ieta][ix][iy].prevPimunu[0][nu][mu] = 0.0;
                            }/* nu */
                        }/* mu */
                    }
                }
            } /* ix, iy, ieta */
        }
        // clean up
        for (int i = 0; i < nx+1; i++) {
            delete[] temp_profile_TA[i];
            delete[] temp_profile_TB[i];
        }
        delete[] temp_profile_TA;
        delete[] temp_profile_TB;
    }
    cout << "initial distribution done." << endl;
    return 1;
}  /* InitTJb*/

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
        fprintf(stderr, "initial_eta_profile out of range.\n");
        exit(0);
    }
    return res;
}

double Init::eta_profile_left_factor(InitData *Data, double eta) {
    // this function return the eta envelope for projectile
    double res = eta_profile_normalisation(Data, eta);
    res = (1. - eta/Data->beam_rapidity)*res;
    return(res);
}

double Init::eta_profile_right_factor(InitData *Data, double eta) {
    // this function return the eta envelope for target
    double res = eta_profile_normalisation(Data, eta);
    res = (1. + eta/Data->beam_rapidity)*res;
    return(res);
}

double Init::eta_rhob_profile_normalisation(InitData *DATA, double eta) {
    // this function return the eta envelope profile for net baryon density
    double res;
    int profile_flag = DATA->initial_eta_rhob_profile;
    double eta_0 = DATA->eta_rhob_0;
    double tau0 = DATA->tau0;
    if (profile_flag == 1) {
        double eta_width = DATA->eta_rhob_width;
        double norm = 1./(2.*sqrt(2*M_PI)*eta_width*tau0);
        double exparg1 = (eta - eta_0)/eta_width;
        double exparg2 = (eta + eta_0)/eta_width;
        res = norm*(exp(-exparg1*exparg1/2.0) + exp(-exparg2*exparg2/2.0));
    } else if (profile_flag == 2) {
        double eta_abs = fabs(eta);
        double delta_eta_1 = DATA->eta_rhob_width_1;
        double delta_eta_2 = DATA->eta_rhob_width_2;
        double A = DATA->eta_rhob_plateau_height;
        double exparg1 = (eta_abs - eta_0)/delta_eta_1;
        double exparg2 = (eta_abs - eta_0)/delta_eta_2;
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
        fprintf(stderr, "initial_eta_rhob_profile = %d out of range.\n",
                profile_flag);
        exit(0);
    }
    return res;
}

double Init::eta_rhob_left_factor(InitData *DATA, double eta) {
    double eta_0 = -fabs(DATA->eta_rhob_0);
    double tau0 = DATA->tau0;
    double delta_eta_1 = DATA->eta_rhob_width_1;
    double delta_eta_2 = DATA->eta_rhob_width_2;
    double norm = 2./(sqrt(M_PI)*tau0*(delta_eta_1 + delta_eta_2));
    double exp_arg = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_1;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_2;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}

double Init::eta_rhob_right_factor(InitData *DATA, double eta) {
    double eta_0 = fabs(DATA->eta_rhob_0);
    double tau0 = DATA->tau0;
    double delta_eta_1 = DATA->eta_rhob_width_1;
    double delta_eta_2 = DATA->eta_rhob_width_2;
    double norm = 2./(sqrt(M_PI)*tau0*(delta_eta_1 + delta_eta_2));
    double exp_arg = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_2;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_1;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}
