// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include "./util.h"
#include "./grid.h"
#include "./init.h"
#include "./eos.h"

using namespace std;

Init::Init(EOS *eosIn) {
    eos = eosIn;
    util = new Util;
    random = new Random;
}

// destructor
Init::~Init() {
    delete random;
    delete util;
}

void Init::InitArena(InitData *DATA, Grid ****arena, Grid ****Lneighbor,
                     Grid ****Rneighbor, int size, int rank) {
    Grid *helperGrid;
    helperGrid = new Grid;
    cout << "initArena" << endl;
    if (DATA->Initial_profile == 6 || DATA->Initial_profile == 7
        || DATA->Initial_profile == 8 || DATA->Initial_profile == 9) {
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
        if (DATA->Initial_profile == 6 || DATA->Initial_profile == 7) {
            DATA->neta = neta/size;
            DATA->delta_eta = deta;
        }

        cout << "neta=" << neta << ", nx=" << nx << ", ny=" << ny << endl;
        cout << "deta=" << DATA->delta_eta << ", dx=" << DATA->delta_x
             << ", dy=" << DATA->delta_y << endl;
    }

    // initialize arena
    *arena = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, DATA->neta);
    *Lneighbor = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, 2);
    *Rneighbor = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, 2);

    cout << "Grid allocated." << endl;
    InitTJb(DATA, arena, Lneighbor, Rneighbor, size, rank);

    cout << "rank " << rank << " ok." << endl;
    LinkNeighbors(DATA, arena, size, rank);
    delete helperGrid;
}/* InitArena */


void Init::LinkNeighbors(InitData *DATA, Grid ****arena, int size, int rank) {
    int nx = DATA->nx;
    int ny = DATA->ny;
    int neta = DATA->neta;

    /* allocate memory */
    for (int ix = 0; ix <= nx; ix++) {
        for (int iy = 0; iy <= ny; iy++) {
            for (int ieta = 0; ieta < neta; ieta++) {
                (*arena)[ix][iy][ieta].nbr_p_1 = new Grid *[4];
                (*arena)[ix][iy][ieta].nbr_m_1 = new Grid *[4];
                (*arena)[ix][iy][ieta].nbr_p_2 = new Grid *[4];
                (*arena)[ix][iy][ieta].nbr_m_2 = new Grid *[4];

                (*arena)[ix][iy][ieta].position[1] = ix;
                (*arena)[ix][iy][ieta].position[2] = iy;
                (*arena)[ix][iy][ieta].position[3] = ieta;
            }
        }
    }

    for (int ix = 0; ix <= nx; ix++) {
        for (int iy = 0; iy <= ny; iy++) {
            for (int ieta = 0; ieta < neta; ieta++) {
                if (ix != nx)
                    (*arena)[ix][iy][ieta].nbr_p_1[1] =
                                                &(*arena)[ix+1][iy][ieta];
                else
                    (*arena)[ix][iy][ieta].nbr_p_1[1] = NULL;
                if (ix < nx - 1)
                    (*arena)[ix][iy][ieta].nbr_p_2[1] =
                                                &(*arena)[ix+2][iy][ieta];
                else
                    (*arena)[ix][iy][ieta].nbr_p_2[1] = NULL;
                if (ix != 0)
                    (*arena)[ix][iy][ieta].nbr_m_1[1] =
                                                &(*arena)[ix-1][iy][ieta];
                else
                    (*arena)[ix][iy][ieta].nbr_m_1[1] = NULL;
                if (ix > 1)
                    (*arena)[ix][iy][ieta].nbr_m_2[1] =
                                                &(*arena)[ix-2][iy][ieta];
                else
                    (*arena)[ix][iy][ieta].nbr_m_2[1] = NULL;
                if (iy != ny)
                    (*arena)[ix][iy][ieta].nbr_p_1[2] =
                                                &(*arena)[ix][iy+1][ieta];
                else
                    (*arena)[ix][iy][ieta].nbr_p_1[2] = NULL;
                if (iy < ny - 1)
                    (*arena)[ix][iy][ieta].nbr_p_2[2] =
                                                &(*arena)[ix][iy+2][ieta];
                else
                    (*arena)[ix][iy][ieta].nbr_p_2[2] = NULL;
                if (iy != 0)
                    (*arena)[ix][iy][ieta].nbr_m_1[2] =
                                                &(*arena)[ix][iy-1][ieta];
                else
                    (*arena)[ix][iy][ieta].nbr_m_1[2] = NULL;
                if (iy > 1)
                    (*arena)[ix][iy][ieta].nbr_m_2[2] =
                                                &(*arena)[ix][iy-2][ieta];
                else
                    (*arena)[ix][iy][ieta].nbr_m_2[2] = NULL;

                // do not care which rank it is - that is dealt with in
                // evolve.cpp
                if (ieta != neta-1)
                    (*arena)[ix][iy][ieta].nbr_p_1[3] =
                                                &(*arena)[ix][iy][ieta+1];
                else
                    (*arena)[ix][iy][ieta].nbr_p_1[3] = NULL;
                if (ieta < neta-2)
                    (*arena)[ix][iy][ieta].nbr_p_2[3] =
                                                &(*arena)[ix][iy][ieta+2];
                else
                    (*arena)[ix][iy][ieta].nbr_p_2[3] = NULL;
                if (ieta != 0)
                    (*arena)[ix][iy][ieta].nbr_m_1[3] =
                                                &(*arena)[ix][iy][ieta-1];
                else
                    (*arena)[ix][iy][ieta].nbr_m_1[3] = NULL;
                if (ieta > 1)
                    (*arena)[ix][iy][ieta].nbr_m_2[3] =
                                                &(*arena)[ix][iy][ieta-2];
                else
                    (*arena)[ix][iy][ieta].nbr_m_2[3] = NULL;
            }
        }
    }
}/* LinkNeighbors */

int Init::InitTJb(InitData *DATA, Grid ****arena, int size, int rank) {
    double epsilon, p, u[4], x, y, eta, rho;
    double rhob = 0.0;
    double exparg1, exparg;
    int ix, iy, ieta, mu, nu;
    int initializeEntropy = DATA->initializeEntropy;

    int rk_order = DATA->rk_order;
    double eta_fall_off = DATA->eta_fall_off;
    double eta_flat = DATA->eta_flat;
    double eta0 = 0.0;  /* modify this to Hirano's transvere row-on-row shift */
    double epsilon0 = DATA->epsilon0;  /* this is in GeV/fm^3 */

    if (initializeEntropy == 0)
        epsilon0 /= hbarc;  /* everything is in fm now. eps is in fm^-4 */

    cout << "rk_order=" << rk_order << endl;
    if (DATA->Initial_profile == 0) {
    } else if (DATA->Initial_profile == 8) {
        // read in the profile from file
        // - IPGlasma initial conditions with initial flow
        size = DATA->size;

        cout << "size=" << size << endl;
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
        for (ix = 0; ix <= DATA->nx; ix++) {
            for (iy = 0; iy <= DATA->ny; iy++) {
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

        for (ix = 0; ix <= DATA->nx; ix++) {
            for (iy = 0; iy <= DATA->ny; iy++) {
                (*Lneighbor)[ix][iy][0].TJb =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][0].TJb =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][1].TJb =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][1].TJb =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][0].Wmunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][0].Wmunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][1].Wmunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][1].Wmunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][0].Pimunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][0].Pimunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][1].Pimunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][1].Pimunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
            }
        }
        int entropy_flag = DATA->initializeEntropy;
        for (int ieta = 0; ieta < DATA->neta; ieta++) {
            double eta = ((DATA->delta_eta)*(ieta + DATA->neta*rank)
                          - (DATA->eta_size)/2.0);
            double eta_envelop_ed = eta_profile_normalisation(DATA, eta);
            for (ix = 0; ix <= DATA->nx; ix++) {
                for (iy = 0; iy<= DATA->ny; iy++) {
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
                    p = eos->get_pressure(epsilon, rhob);
                    // set all values in the grid element:
                    (*arena)[ix][iy][ieta].epsilon = epsilon;
                    (*arena)[ix][iy][ieta].epsilon_t = epsilon;
                    (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
                    (*arena)[ix][iy][ieta].rhob = rhob;
                    (*arena)[ix][iy][ieta].rhob_t = rhob;
                    (*arena)[ix][iy][ieta].prev_rhob = rhob;
                    (*arena)[ix][iy][ieta].p = p;
                    (*arena)[ix][iy][ieta].p_t = p;
                    (*arena)[ix][iy][ieta].trouble = 0;

                    (*arena)[ix][iy][ieta].T = eos->get_temperature(epsilon,
                                                                    rhob);
                    (*arena)[ix][iy][ieta].mu = eos->get_mu(epsilon, rhob);
                
                    (*arena)[ix][iy][ieta].TJb =
                                            util->cube_malloc(rk_order+1, 5, 4);
                    (*arena)[ix][iy][ieta].dUsup =
                                            util->cube_malloc(rk_order+1, 5, 4);
                    (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 5);
                    (*arena)[ix][iy][ieta].theta_u =
                                                util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].sigma =
                                            util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].pi_b =
                                            util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
                    (*arena)[ix][iy][ieta].Wmunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                    (*arena)[ix][iy][ieta].prevWmunu =
                                            util->cube_malloc(rk_order, 5, 4);
                    (*arena)[ix][iy][ieta].Pimunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                    (*arena)[ix][iy][ieta].prevPimunu =
                                            util->cube_malloc(rk_order, 5, 4);
                    (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5, 4);

                    /* for HIC */
                    (*arena)[ix][iy][ieta].u[0][0] = temp_profile_utau[ix][iy];
                    (*arena)[ix][iy][ieta].u[0][1] = temp_profile_ux[ix][iy];
                    (*arena)[ix][iy][ieta].u[0][2] = temp_profile_uy[ix][iy];
                    (*arena)[ix][iy][ieta].u[0][3] = 0.0;

                    u[0] = temp_profile_utau[ix][iy];
                    u[1] = temp_profile_ux[ix][iy];
                    u[2] = temp_profile_uy[ix][iy];
                    u[3] = 0.0;

                    (*arena)[ix][iy][ieta].prev_u[0][0] = u[0];
                    (*arena)[ix][iy][ieta].prev_u[0][1] = u[1];
                    (*arena)[ix][iy][ieta].prev_u[0][2] = u[2];
                    (*arena)[ix][iy][ieta].prev_u[0][3] = u[3];

                    (*arena)[ix][iy][ieta].pi_b[0] = 0.0;

                    for (int mu = 0; mu < 4; mu++) {
                        /* baryon density */
                        (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];

                        // diffusion current
                        (*arena)[ix][iy][ieta].Wmunu[0][4][mu] = 0.0;
                        (*arena)[ix][iy][ieta].prevWmunu[0][4][mu] = 0.0;
                    
                        for (nu = 0; nu < 4; nu++) {
                            (*arena)[ix][iy][ieta].TJb[0][nu][mu] = (
                                (epsilon + p)*u[mu]*u[nu]
                                + p*(DATA->gmunu)[mu][nu]);
                            (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = 0.0;

                            (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = 0.0;
                        }/* nu */
                    }/* mu */
                }
            }
        }/* ix, iy, ieta */
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
        size = DATA->size;
        cout << "size=" << size << endl;
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

        for (int ix = 0; ix < (DATA->nx+1); ix++) {
            for (int iy = 0; iy < (DATA->ny+1); iy++) {
                (*Lneighbor)[ix][iy][0].TJb =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][0].TJb =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][1].TJb =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][1].TJb =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][0].Wmunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][0].Wmunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][1].Wmunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][1].Wmunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][0].Pimunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][0].Pimunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][1].Pimunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][1].Pimunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
            }
        }
        
        int entropy_flag = DATA->initializeEntropy;
        for (int ieta = 0; ieta < DATA->neta; ieta++) {
            double eta = ((DATA->delta_eta)*(ieta + DATA->neta*rank)
                          - (DATA->eta_size)/2.0);
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
                    p = eos->get_pressure(epsilon, rhob);

                    // set all values in the grid element:
                    (*arena)[ix][iy][ieta].epsilon = epsilon;
                    (*arena)[ix][iy][ieta].epsilon_t = epsilon;
                    (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
                    (*arena)[ix][iy][ieta].rhob = rhob;
                    (*arena)[ix][iy][ieta].rhob_t = rhob;
                    (*arena)[ix][iy][ieta].prev_rhob = rhob;
                    (*arena)[ix][iy][ieta].p = p;
                    (*arena)[ix][iy][ieta].p_t = p;
                    (*arena)[ix][iy][ieta].trouble = 0;
                    (*arena)[ix][iy][ieta].T =
                                            eos->get_temperature(epsilon, rhob);
                    (*arena)[ix][iy][ieta].mu = eos->get_mu(epsilon, rhob);
                    (*arena)[ix][iy][ieta].TJb =
                                            util->cube_malloc(rk_order+1, 5, 4);
                    (*arena)[ix][iy][ieta].dUsup =
                                            util->cube_malloc(rk_order+1, 5, 4);
                    (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 5);
                    (*arena)[ix][iy][ieta].theta_u =
                                            util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].sigma =
                                            util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].pi_b =
                                            util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
                    (*arena)[ix][iy][ieta].Wmunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                    (*arena)[ix][iy][ieta].prevWmunu =
                                            util->cube_malloc(rk_order, 5, 4);
                    (*arena)[ix][iy][ieta].Pimunu =
                                            util->cube_malloc(rk_order+1, 5, 4);
                    (*arena)[ix][iy][ieta].prevPimunu =
                                            util->cube_malloc(rk_order, 5, 4);
                    (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5, 4);

                    /* for HIC */
                    u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
                    u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
                    u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
                    u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;
                    (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
                    (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
                    (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
                    (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;

                    (*arena)[ix][iy][ieta].pi_b[0] = 0.0;

                    for (int mu = 0; mu < 4; mu++) {
                        /* baryon density */
                        (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];

                        // diffusion current
                        (*arena)[ix][iy][ieta].Wmunu[0][4][mu] = 0.0;
                        (*arena)[ix][iy][ieta].prevWmunu[0][4][mu] = 0.0;
                        for (int nu = 0; nu < 4; nu++) {
                            (*arena)[ix][iy][ieta].TJb[0][nu][mu] = (
                                (epsilon + p)*u[mu]*u[nu]
                                + p*(DATA->gmunu)[mu][nu]);
                            (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = 0.0;

                            (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = 0.0;
                        }/* nu */
                    }/* mu */
                }
            }
        } /* ix, iy, ieta */
        // clean up
        for (int i = 0; i < nx+1; i++) {
            delete[] temp_profile_TA[i];
            delete[] temp_profile_TB[i];
        }
        delete[] temp_profile_TA;
        delete[] temp_profile_TB;
    } else if (DATA->Initial_profile == 12) {
        // read in the 3d profile from file
        size = DATA->size;

        cout << "size=" << size << endl;
        cout << " ----- information on initial distribution -----" << endl;
        cout << "file name used: " << DATA->initName << " and "
             << DATA->initName_rhob << endl;
  
      // first load in the transverse profile
      ifstream profile_ed(DATA->initName.c_str());
      ifstream profile_rhob(DATA->initName_rhob.c_str());

      int input_grid_nx = DATA->input_grid_nx;
      int input_grid_ny = DATA->input_grid_ny;
      int input_grid_neta = DATA->input_grid_neta;
      double input_grid_dx = DATA->input_grid_dx;
      double input_grid_dy = DATA->input_grid_dy;
      double input_grid_deta = DATA->input_grid_deta;
      double input_eta_0 = - (input_grid_neta - 1)/2.*input_grid_deta;
      double input_x_0 = - (input_grid_nx - 1)/2.*input_grid_dx;
      double input_y_0 = - (input_grid_ny - 1)/2.*input_grid_dy;

      double*** temp_profile_ed = new double** [input_grid_neta];
      double*** temp_profile_rhob = new double** [input_grid_neta];
      for(int i = 0; i < input_grid_neta; i++)
      {
          temp_profile_ed[i] = new double* [input_grid_nx];
          temp_profile_rhob[i] = new double* [input_grid_nx];
          for(int j = 0; j < input_grid_nx; j++)
          {
              temp_profile_ed[i][j] = new double [input_grid_ny];
              temp_profile_rhob[i][j] = new double [input_grid_ny];
          }
      }
      double check_Npart = 0.0;
      for(int i = 0; i < input_grid_neta; i++)
      {
         for(int j = 0; j < input_grid_nx; j++)
         {
            for(int k = 0; k < input_grid_ny; k++)
            {
               profile_ed >> temp_profile_ed[i][j][k];
               profile_rhob >> temp_profile_rhob[i][j][k];
               if(rank == 0)
                   check_Npart += temp_profile_rhob[i][j][k];
            }
         }
      }
      profile_ed.close();
      profile_rhob.close();
      if(rank == 0)
      {
          check_Npart *= input_grid_dx*input_grid_dy*input_grid_deta;
          cout << "Net baryon number (= Npart) = " << check_Npart << endl;
      }
      
      for(int ix=0; ix< (DATA->nx+1); ix++)
      {
         for(int iy=0; iy< (DATA->ny+1); iy++)
         {
            (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
            (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
            (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
            (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
            (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
            (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
            (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
            (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
         }
      }
  
      int entropy_flag = DATA->initializeEntropy;
      for(int ieta = 0; ieta < DATA->neta; ieta++)
      {
         double eta = (
            (DATA->delta_eta)*(ieta + DATA->neta*rank) - (DATA->eta_size)/2.0);
         int eta_idx = (int)((eta - input_eta_0)/input_grid_deta);
         double eta_frac = (
            (eta - (input_eta_0 + eta_idx*input_grid_deta))/input_grid_deta);
         for(int ix = 0; ix < (DATA->nx+1); ix++)
         {
            double x_local = (DATA->delta_x)*ix - (DATA->x_size)/2.0;
            int x_idx = (int)((x_local - input_x_0)/input_grid_dx);
            double x_frac = (
               (x_local - (input_x_0 + x_idx*input_grid_dx))/input_grid_dx);
            for(int iy = 0; iy< (DATA->ny+1); iy++)
            {
               double y_local = (DATA->delta_y)*iy - (DATA->y_size)/2.0;
               int y_idx = (int)((y_local - input_y_0)/input_grid_dy);
               double y_frac = (
                  (y_local - (input_y_0 + y_idx*input_grid_dy))/input_grid_dy);
     
               double rhob_interp = 0.0;
               double ed_interp = 0.0;
     
               if(   eta_idx < 0 || eta_idx > input_grid_neta - 2
                  || x_idx < 0 || x_idx > input_grid_nx - 2 
                  || y_idx < 0 || y_idx > input_grid_ny - 2)
               {
                   rhob = 0.0;
                   epsilon = 0.0;
               }
               else
               {
                   rhob_interp = (
                        temp_profile_rhob[eta_idx][x_idx][y_idx]
                        *(1. - eta_frac)*(1. - x_frac)*(1. - y_frac)
                      + temp_profile_rhob[eta_idx+1][x_idx][y_idx]
                        *eta_frac*(1. - x_frac)*(1. - y_frac)
                      + temp_profile_rhob[eta_idx][x_idx+1][y_idx]
                        *(1. - eta_frac)*x_frac*(1. - y_frac)
                      + temp_profile_rhob[eta_idx+1][x_idx+1][y_idx]
                        *eta_frac*x_frac*(1. - y_frac)
                      + temp_profile_rhob[eta_idx][x_idx][y_idx+1]
                        *(1. - eta_frac)*(1. - x_frac)*y_frac
                      + temp_profile_rhob[eta_idx+1][x_idx][y_idx+1]
                        *eta_frac*(1. - x_frac)*y_frac
                      + temp_profile_rhob[eta_idx][x_idx+1][y_idx+1]
                        *(1. - eta_frac)*x_frac*y_frac
                      + temp_profile_rhob[eta_idx+1][x_idx+1][y_idx+1]
                        *eta_frac*x_frac*y_frac
                   );
                   ed_interp = (
                        temp_profile_ed[eta_idx][x_idx][y_idx]
                        *(1. - eta_frac)*(1. - x_frac)*(1. - y_frac)
                      + temp_profile_ed[eta_idx+1][x_idx][y_idx]
                        *eta_frac*(1. - x_frac)*(1. - y_frac)
                      + temp_profile_ed[eta_idx][x_idx+1][y_idx]
                        *(1. - eta_frac)*x_frac*(1. - y_frac)
                      + temp_profile_ed[eta_idx+1][x_idx+1][y_idx]
                        *eta_frac*x_frac*(1. - y_frac)
                      + temp_profile_ed[eta_idx][x_idx][y_idx+1]
                        *(1. - eta_frac)*(1. - x_frac)*y_frac
                      + temp_profile_ed[eta_idx+1][x_idx][y_idx+1]
                        *eta_frac*(1. - x_frac)*y_frac
                      + temp_profile_ed[eta_idx][x_idx+1][y_idx+1]
                        *(1. - eta_frac)*x_frac*y_frac
                      + temp_profile_ed[eta_idx+1][x_idx+1][y_idx+1]
                        *eta_frac*x_frac*y_frac
                   );
               }
               
     
               if(DATA->turn_on_rhob == 1)
                  rhob = rhob_interp;  // 1/fm^3
               else
                  rhob = 0.0;
     
               if(entropy_flag == 0)
                  epsilon= ed_interp*DATA->sFactor/hbarc;  // 1/fm^4
               else
               {
                  double local_sd = ed_interp*DATA->sFactor;
                  epsilon = eos->get_s2e(local_sd, rhob);
               }

               if (epsilon<0.00000000001)
                   epsilon = 0.00000000001;
        
               // initial pressure distribution
               p = eos->get_pressure(epsilon, rhob);
               
               // set all values in the grid element:
               (*arena)[ix][iy][ieta].epsilon = epsilon;
               (*arena)[ix][iy][ieta].epsilon_t = epsilon;
               (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
               (*arena)[ix][iy][ieta].rhob = rhob;
               (*arena)[ix][iy][ieta].rhob_t = rhob;
               (*arena)[ix][iy][ieta].prev_rhob = rhob;
               (*arena)[ix][iy][ieta].p = p;
               (*arena)[ix][iy][ieta].p_t = p;
               (*arena)[ix][iy][ieta].trouble = 0;
               
               (*arena)[ix][iy][ieta].T = eos->get_temperature(epsilon, rhob);
               (*arena)[ix][iy][ieta].mu = eos->get_mu(epsilon, rhob);
               
               (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
               (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 5,4);
               (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
               (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 5);
               (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
               (*arena)[ix][iy][ieta].sigma = util->cube_malloc(rk_order+1, 4, 4);
               (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
               //(*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
               //(*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
               (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
               //(*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
               (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
               (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
               //(*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
               (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
               (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
               //(*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);
               (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);

               /* for HIC */
               u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
               u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
               u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
               u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;
               
               (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
               (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
               (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
               (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
               //(*arena)[ix][iy][ieta].pprev_u[0][0] = 1.0;
               //(*arena)[ix][iy][ieta].pprev_u[0][3] = 0.0;
               //(*arena)[ix][iy][ieta].pprev_u[0][1] = 0.0;
               //(*arena)[ix][iy][ieta].pprev_u[0][2] = 0.0;
         
               (*arena)[ix][iy][ieta].pi_b[0] = 0.0;
               //(*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
               //(*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
               for(int mu=0; mu<4; mu++)
               {
                  /* baryon density */
                  (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
                  
                  // diffusion current
                  (*arena)[ix][iy][ieta].Wmunu[0][4][mu] = (double) 0.0;
                  (*arena)[ix][iy][ieta].prevWmunu[0][4][mu] = (double) 0.0;
                  
                  for(nu=0; nu<4; nu++)
                  {
                     (*arena)[ix][iy][ieta].TJb[0][nu][mu] = (epsilon + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];
                     (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = (double) 0.0;
                     (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (double) 0.0;
                     //(*arena)[ix][iy][ieta].pprevWmunu[0][nu][mu] = (double) 0.0;
                     
                     (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
                     (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
                     //(*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
                  }/* nu */
               }/* mu */
            }
         }
      }/* ix, iy, ieta */
        
      // clean up
      for(int i = 0; i < input_grid_neta; i++)
      {
          for(int j = 0; j < input_grid_nx; j++)
          {
              delete[] temp_profile_ed[i][j];
              delete[] temp_profile_rhob[i][j];
          }
          delete[] temp_profile_ed[i];
          delete[] temp_profile_rhob[i];
      }
      delete[] temp_profile_ed;
      delete[] temp_profile_rhob;
    } else if (DATA->Initial_profile==13) {
        //read in the transverse profile from file with initial flow velocity
        size = DATA->size;
       
        cout << "size=" << size << endl;
        cout << " ----- information on initial distribution -----" << endl;
        cout << "file name used: " << DATA->initName << " and "
             << DATA->initName_rhob << endl;
  
      // first load in the transverse profile
      ifstream profile_ed(DATA->initName.c_str());
      ifstream profile_rhob(DATA->initName_rhob.c_str());
      ifstream profile_ux(DATA->initName_ux.c_str());
      ifstream profile_uy(DATA->initName_uy.c_str());
      int nx = DATA->nx;
      int ny = DATA->ny;
      double** temp_profile_ed = new double* [nx+1];
      double** temp_profile_rhob = new double* [nx+1];
      double** temp_profile_ux = new double* [nx+1];
      double** temp_profile_uy = new double* [nx+1];
      for(int i = 0; i < nx+1; i++)
      {
          temp_profile_ed[i] = new double [ny+1];
          temp_profile_rhob[i] = new double [ny+1];
          temp_profile_ux[i] = new double [ny+1];
          temp_profile_uy[i] = new double [ny+1];
      }
      for(int i = 0; i < nx+1; i++)
      {
         for(int j = 0; j < ny+1; j++)
         {
            profile_ed >> temp_profile_ed[i][j];
            profile_rhob >> temp_profile_rhob[i][j];
            profile_ux >> temp_profile_ux[i][j];
            profile_uy >> temp_profile_uy[i][j];
         }
      }
      profile_ed.close();
      profile_rhob.close();

      for(int ix=0; ix< (DATA->nx+1); ix++)
      {
         for(int iy=0; iy< (DATA->ny+1); iy++)
         {
            (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
            (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
            (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
            (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
            (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
            (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
            (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
            (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
            (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
            (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
         }
      }
  
      int entropy_flag = DATA->initializeEntropy;
      for(int ieta = 0; ieta < DATA->neta; ieta++)
      {
         double eta = (DATA->delta_eta)*(ieta + DATA->neta*rank) - (DATA->eta_size)/2.0;
         double eta_envelop_ed = eta_profile_normalisation(DATA, eta);
         double eta_envelop_rhob = eta_rhob_profile_normalisation(DATA, eta);
         for(int ix = 0; ix < (DATA->nx+1); ix++)
         {
            for(int iy = 0; iy< (DATA->ny+1); iy++)
            {
               if(DATA->turn_on_rhob == 1)
                  rhob = temp_profile_rhob[ix][iy]*eta_envelop_rhob;  // 1/fm^3
               else
                  rhob = 0.0;
               if(entropy_flag == 0)
                  epsilon= temp_profile_ed[ix][iy]*eta_envelop_ed*DATA->sFactor/hbarc;  // 1/fm^4
               else
               {
                  double local_sd = temp_profile_ed[ix][iy]*DATA->sFactor*eta_envelop_ed;
                  epsilon = eos->get_s2e(local_sd, rhob);
               }
               if (epsilon<0.00000000001)
                  epsilon = 0.00000000001;
         
               // initial pressure distribution
               p = eos->get_pressure(epsilon, rhob);
            
               // set all values in the grid element:
               (*arena)[ix][iy][ieta].epsilon = epsilon;
               (*arena)[ix][iy][ieta].epsilon_t = epsilon;
               (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
               (*arena)[ix][iy][ieta].rhob = rhob;
               (*arena)[ix][iy][ieta].rhob_t = rhob;
               (*arena)[ix][iy][ieta].prev_rhob = rhob;
               (*arena)[ix][iy][ieta].p = p;
               (*arena)[ix][iy][ieta].p_t = p;
               (*arena)[ix][iy][ieta].trouble = 0;
             
               (*arena)[ix][iy][ieta].T = eos->get_temperature(epsilon, rhob);
               (*arena)[ix][iy][ieta].mu = eos->get_mu(epsilon, rhob);
            
               (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
               (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 5,4);
               (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
               (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 5);
               (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
               (*arena)[ix][iy][ieta].sigma = util->cube_malloc(rk_order+1, 4, 4);
               (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
               //(*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
               //(*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
               (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
               //(*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
               (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
               (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
               //(*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
               (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
               (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
               //(*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);
            
               (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
               /* for HIC */
               double temp_ux = temp_profile_ux[ix][iy];
               double temp_uy = temp_profile_uy[ix][iy];
               double temp_gamma = sqrt(1. + temp_ux*temp_ux + temp_uy*temp_uy);
               u[0] = (*arena)[ix][iy][ieta].u[0][0] = temp_gamma;
               u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
               u[1] = (*arena)[ix][iy][ieta].u[0][1] = temp_ux;
               u[2] = (*arena)[ix][iy][ieta].u[0][2] = temp_uy;
               
              (*arena)[ix][iy][ieta].prev_u[0][0] = temp_gamma;
              (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
              (*arena)[ix][iy][ieta].prev_u[0][1] = temp_ux;
              (*arena)[ix][iy][ieta].prev_u[0][2] = temp_uy;
               //(*arena)[ix][iy][ieta].pprev_u[0][0] = 1.0;
               //(*arena)[ix][iy][ieta].pprev_u[0][3] = 0.0;
               //(*arena)[ix][iy][ieta].pprev_u[0][1] = 0.0;
               //(*arena)[ix][iy][ieta].pprev_u[0][2] = 0.0;
           
               (*arena)[ix][iy][ieta].pi_b[0] = 0.0;
               //(*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
               //(*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
               
               for(int mu=0; mu<4; mu++)
               {
                  /* baryon density */
                  (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
                  
                  // diffusion current
                  (*arena)[ix][iy][ieta].Wmunu[0][4][mu] = (double) 0.0;
                  (*arena)[ix][iy][ieta].prevWmunu[0][4][mu] = (double) 0.0;
                  
                  for(nu=0; nu<4; nu++)
                  {
                     (*arena)[ix][iy][ieta].TJb[0][nu][mu] = (epsilon + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];
                     (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = (double) 0.0;
                     (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (double) 0.0;
                     //(*arena)[ix][iy][ieta].pprevWmunu[0][nu][mu] = (double) 0.0;
                     
                     (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
                     (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
                     //(*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
                        }/* nu */
                    }/* mu */
                }
            }
        }/* ix, iy, ieta */
        
        // clean up
        for (int i = 0; i < nx+1; i++) {
            delete[] temp_profile_ed[i];
            delete[] temp_profile_rhob[i];
            delete[] temp_profile_ux[i];
            delete[] temp_profile_uy[i];
        }
        delete[] temp_profile_ed;
        delete[] temp_profile_rhob;
        delete[] temp_profile_ux;
        delete[] temp_profile_uy;
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
