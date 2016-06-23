#include <omp.h>
#include "./util.h"
#include "./data.h"
#include "./grid.h"
#include "./minmod.h"
#include "./eos.h"
#include "./u_derivative.h"

using namespace std;

U_derivative::U_derivative(EOS *eosIn, InitData* DATA_in) {
   eos = eosIn;
   minmod = new Minmod(DATA_in);
}

// destructor
U_derivative::~U_derivative() {
   delete minmod;
}

int U_derivative::MakedU(double tau, InitData *DATA,
                         Grid ***arena, int rk_flag) {
    // ideal hydro: no need to evaluate any flow derivatives
    if (DATA->viscosity_flag == 0)
        return 1;  

    int neta = DATA->neta-1;

    int ieta;
    #pragma omp parallel private(ieta)
    {
        #pragma omp for
        for (ieta = 0; ieta <= neta; ieta++) {
            MakedUXY(tau, ieta, DATA, arena, rk_flag);
        }/* ieta */
        #pragma omp barrier

        #pragma omp for
        for (ieta = 0; ieta <= neta; ieta++) {
            Make_expansion_rate_XY(tau, ieta, DATA, arena, rk_flag);
        }/* ieta */

        // calculate the velocity shear tensor sigma^{\mu\nu}
        // immigrate the code from dissipative.cpp
        #pragma omp for
        for (ieta = 0; ieta <= neta; ieta++) {
            Make_sigma_XY(tau, ieta, DATA, arena, rk_flag);
        }/* ieta */
        #pragma omp barrier
    }

   return 1; /* successful */
}/* MakedU */

void U_derivative::MakedUXY(double tau, int ieta, InitData *DATA,
                            Grid ***arena, int rk_flag) {
    // This function is a shell function to calculate parital^\nu u^\mu
    int nx = DATA->nx;
    int ny = DATA->ny;
    for (int ix = 0; ix <= nx; ix++) {
        for (int iy = 0; iy <= ny; iy++) {
	       /* this calculates du/dx, du/dy, (du/deta)/tau */
           MakeDSpatial(tau, DATA, &(arena[ieta][ix][iy]), rk_flag);
           /* this calculates du/dtau */
           MakeDTau(tau, DATA, &(arena[ieta][ix][iy]), rk_flag); 
        }/* ieta */
    }/*iy */
}

void U_derivative::Make_expansion_rate_XY(double tau, int ieta, InitData *DATA,
                                          Grid ***arena, int rk_flag) {
    // this function calculates the local expansion rate in the transverse
    // plane
    int nx = DATA->nx;
    int ny = DATA->ny;
    for (int ix = 0; ix <= nx; ix++) {
        for (int iy = 0; iy <= ny; iy++) {
            double u_local[4];
            for(int i = 0; i < 4; i++)
                u_local[i] = arena[ieta][ix][iy].u[rk_flag][i];
            double dUsup_local[5][4];   // dUsup[m][n] = partial^n u^m
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 4; j++) {
                    dUsup_local[i][j] = (
                            arena[ieta][ix][iy].dUsup[0][i][j]);
                }
            }
            double partial_mu_u_supmu = 0.0;
            for (int mu=0; mu<4; mu++) {
                double gfac = (mu==0 ? -1.0 : 1.0);
                // for expansion rate: theta
                partial_mu_u_supmu += dUsup_local[mu][mu]*gfac;

                double u_supnu_partial_nu_u_supmu = 0.0;
	            for (int nu=0; nu<4; nu++) {
                    double tfac = (nu==0 ? -1.0 : 1.0);
                    u_supnu_partial_nu_u_supmu += (
                            tfac*u_local[nu]*dUsup_local[mu][nu]);
                }
                arena[ieta][ix][iy].a[0][mu] = (
                                            u_supnu_partial_nu_u_supmu);
	        }/* mu */

            int mu = 4;
            if (DATA->turn_on_diff == 1) {
                // derivatives of chemical potentials for diffusion
                double u_supnu_partial_nu_muoverT = 0.0;
	            for (int nu = 0; nu < 4; nu++) {
                    double tfac = (nu==0 ? -1.0 : 1.0);
                    u_supnu_partial_nu_muoverT += (
                                    tfac*u_local[nu]*dUsup_local[mu][nu]);
                }
                arena[ieta][ix][iy].a[0][mu] = u_supnu_partial_nu_muoverT;
            } else {
                arena[ieta][ix][iy].a[0][mu] = 0.0;
            }

            // expansion rate need to add the tau-eta coordinate source term
            arena[ieta][ix][iy].theta_u[0] = (
                                    partial_mu_u_supmu + u_local[0]/tau);
        }/*iy */
    }/* ix */
}

void U_derivative::Make_sigma_XY(double tau, int ieta, InitData *DATA,
                                 Grid ***arena, int rk_flag) {
    int nx = DATA->nx;
    int ny = DATA->ny;
    for (int ix=0; ix <= nx; ix++) {
        for (int iy=0; iy <= ny; iy++) {
            double dUsup_local[4][4];
            double u_local[4];
            double a_local[4];
            double theta_u_local = arena[ieta][ix][iy].theta_u[0];
            double sigma_local[4][4];
            for (int i = 0; i < 4; i++) {
                u_local[i] = arena[ieta][ix][iy].u[rk_flag][i];
                a_local[i] = arena[ieta][ix][iy].a[0][i];
                for (int j = 0; j < 4; j++) {
                    dUsup_local[i][j] = (
                                arena[ieta][ix][iy].dUsup[0][i][j]);
                }
            }
            for (int a = 1; a < 4; a++) {
                for (int b = a; b < 4; b++) {
                    sigma_local[a][b] = (
                        (dUsup_local[a][b] + dUsup_local[b][a])/2.
                        - (DATA->gmunu[a][b] + u_local[a]*u_local[b])
                          *theta_u_local/3.
                        + u_local[0]/tau*DATA->gmunu[a][3]*DATA->gmunu[b][3]
                        - u_local[3]/tau/2.
                          *(DATA->gmunu[a][3]*DATA->gmunu[b][0] 
                            + DATA->gmunu[b][3]*DATA->gmunu[a][0])
                        + u_local[3]*u_local[0]/tau/2.
                          *(DATA->gmunu[a][3]*u_local[b] 
                            + DATA->gmunu[b][3]*u_local[a])
                        - u_local[3]*u_local[3]/tau/2.
                          *(DATA->gmunu[a][0]*u_local[b] 
                            + DATA->gmunu[b][0]*u_local[a])
                        + (u_local[a]*a_local[b]
                           + u_local[b]*a_local[a])/2.);
                    sigma_local[b][a] = sigma_local[a][b];
                }
            }
            // make sigma[3][3] using traceless condition
            sigma_local[3][3] = (
                (  2.*(  u_local[1]*u_local[2]*sigma_local[1][2]
                       + u_local[1]*u_local[3]*sigma_local[1][3]
                       + u_local[2]*u_local[3]*sigma_local[2][3])
                 - (u_local[0]*u_local[0] - u_local[1]*u_local[1])
                   *sigma_local[1][1]
                 - (u_local[0]*u_local[0] - u_local[2]*u_local[2])
                   *sigma_local[2][2])
                /(u_local[0]*u_local[0] - u_local[3]*u_local[3]));
            // make sigma[0][i] using transversality
            for (int a = 1; a < 4; a++) {
                double temp = 0.0;
                for (int b = 1; b < 4; b++) {
                    temp += sigma_local[a][b]*u_local[b];
                }
                sigma_local[0][a] = temp/u_local[0];
                sigma_local[a][0] = sigma_local[0][a];
            }
            // make sigma[0][0]
            double temp = 0.0;
            for (int a = 1; a < 4; a++) {
                temp += sigma_local[0][a]*u_local[a];
            }
            sigma_local[0][0] = temp/u_local[0];

            // store the 10 essential components in memory
            arena[ieta][ix][iy].sigma[0][0] = sigma_local[0][0];
            arena[ieta][ix][iy].sigma[0][1] = sigma_local[0][1];
            arena[ieta][ix][iy].sigma[0][2] = sigma_local[0][2];
            arena[ieta][ix][iy].sigma[0][3] = sigma_local[0][3];
            arena[ieta][ix][iy].sigma[0][4] = sigma_local[1][1];
            arena[ieta][ix][iy].sigma[0][5] = sigma_local[1][2];
            arena[ieta][ix][iy].sigma[0][6] = sigma_local[1][3];
            arena[ieta][ix][iy].sigma[0][7] = sigma_local[2][2];
            arena[ieta][ix][iy].sigma[0][8] = sigma_local[2][3];
            arena[ieta][ix][iy].sigma[0][9] = sigma_local[3][3];
        }/*iy */
    }/* ix */
}


int U_derivative::MakeDSpatial(double tau, InitData *DATA, Grid *grid_pt,
                               int rk_flag) {
    double g, f, fp1, fm1, taufactor;
    double delta[4];
    //Sangyong Nov 18 2014: added these doubles
    double rhob, eps, muB, T;
 
    delta[1] = DATA->delta_x;
    delta[2] = DATA->delta_y;
    delta[3] = DATA->delta_eta;

    /* dUsup[m][n] = partial_n u_m */
    /* for u[i] */
    for (int m = 1; m <= 3; m++) {
        // partial_n u[m]
        for (int n = 1; n <= 3; n++) {
            taufactor = 1.0;
            if (n == 3)
                taufactor = tau;
            f = grid_pt->u[rk_flag][m];
            fp1 = grid_pt->nbr_p_1[n]->u[rk_flag][m];
            fm1 = grid_pt->nbr_m_1[n]->u[rk_flag][m];
            g = minmod->minmod_dx(fp1, f, fm1);
            g /= delta[n]*taufactor;
            grid_pt->dUsup[0][m][n] = g;
        }  // n = x, y, eta
    }  // m = x, y, eta
    /* for u[0], use u[0]u[0] = 1 + u[i]u[i] */
    /* u[0]_m = u[i]_m (u[i]/u[0]) */
    /* for u[0] */
    for (int n = 1; n <= 3; n++) {
        f = 0.0;
        for (int m = 1; m <= 3; m++) {
	        /* (partial_n u^m) u[m] */
	        f += (grid_pt->dUsup[0][m][n])*(grid_pt->u[rk_flag][m]);
        }
        f /= grid_pt->u[rk_flag][0];
        grid_pt->dUsup[0][0][n] = f;
    }
    // Sangyong Nov 18 2014
    // Here we make derivatives of muB/T
    // dUsup[rk_flag][4][n] = partial_n (muB/T)
    // partial_x (muB/T) and partial_y (muB/T) first
    int m = 4; // means (muB/T)
    for (int n = 1; n <= 3; n++) {
        taufactor = 1.0;
        if (n == 3)
            taufactor = tau;

        // f = grid_pt->rhob_t;
        if (rk_flag == 0) {
            rhob = grid_pt->rhob;
            eps = grid_pt->epsilon;
        } else {
            rhob = grid_pt->rhob_t;
            eps = grid_pt->epsilon_t;
        }
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        f = muB/T; 
    
        //fp1 = grid_pt->nbr_p_1[n]->rhob;
        if (rk_flag == 0) {
            rhob = grid_pt->nbr_p_1[n]->rhob;
            eps = grid_pt->nbr_p_1[n]->epsilon;
        } else {
            rhob = grid_pt->nbr_p_1[n]->rhob_t;
            eps = grid_pt->nbr_p_1[n]->epsilon_t;
        }
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        fp1 = muB/T; 
       
        // fm1 = grid_pt->nbr_m_1[n]->rhob;
        if (rk_flag == 0) {
            rhob = grid_pt->nbr_m_1[n]->rhob;
            eps = grid_pt->nbr_m_1[n]->epsilon;
        } else {
            rhob = grid_pt->nbr_m_1[n]->rhob_t;
            eps = grid_pt->nbr_m_1[n]->epsilon_t;
        }
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        fm1 = muB/T; 

        g = minmod->minmod_dx(fp1, f, fm1);
        g /= delta[n]*taufactor;
        grid_pt->dUsup[0][m][n] = g;
    }  // n = x, y, eta
    return 1;
}/* MakeDSpatial */

int U_derivative::MakeDTau(double tau, InitData *DATA, Grid *grid_pt,
                           int rk_flag) {
    int m;
    double f;
    // Sangyong Nov 18 2014: added these doubles 
    double tildemu, tildemu_prev, rhob, eps, muB, T;

    /* this makes dU[m][0] = partial^tau u^m */
    /* note the minus sign at the end because of g[0][0] = -1 */
    /* rk_flag is 0, 2, 4, ... */

    if (rk_flag == 0) {
        for (m=1; m<=3; m++) {
            /* first order is more stable */
            f = ((grid_pt->u[rk_flag][m] - grid_pt->prev_u[0][m])
                 /DATA->delta_tau);
            grid_pt->dUsup[0][m][0] = -f; /* g00 = -1 */
        }/* m */
    } else if (rk_flag > 0) {
        for (m=1; m<=3; m++) {
            /* first order */
            // this is from the prev full RK step 
            f = (grid_pt->u[rk_flag][m] - grid_pt->u[0][m])/(DATA->delta_tau);
            grid_pt->dUsup[0][m][0] = -f; /* g00 = -1 */
        }/* m */
    }

    /* I have now partial^tau u^i */
    /* I need to calculate (u^i partial^tau u^i) = u^0 partial^tau u^0 */
    /* u_0 d^0 u^0 + u_m d^0 u^m = 0 */
    /* -u^0 d^0 u^0 + u_m d^0 u^m = 0 */
    /* d^0 u^0 = u_m d^0 u^m/u^0 */

    f = 0.0;
    for (m=1; m<=3; m++) {
        /* (partial_0 u^m) u[m] */
        f += (grid_pt->dUsup[0][m][0])*(grid_pt->u[rk_flag][m]);
    }
    f /= grid_pt->u[rk_flag][0];
    grid_pt->dUsup[0][0][0] = f;

    // Sangyong Nov 18 2014
    // Here we make the time derivative of (muB/T)
    if (rk_flag == 0) {
        m = 4;  
        // first order is more stable 
        // backward derivative
        // current values
        // f = (grid_pt->rhob);
        rhob = grid_pt->rhob;
        eps = grid_pt->epsilon;
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        tildemu = muB/T;
        // f -= (grid_pt->rhob_prev);
        rhob = grid_pt->prev_rhob;
        eps = grid_pt->prev_epsilon;
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        tildemu_prev = muB/T;
        f = (tildemu - tildemu_prev)/(DATA->delta_tau);
        grid_pt->dUsup[0][m][0] = -f; /* g00 = -1 */
    } else if (rk_flag > 0) {
        m = 4;  
        // first order 
        // forward derivative
        // f = (grid_pt->rhob_t); // this is from the prev full RK step 
        // f -= (grid_pt->rhob_prev);
        // f /= (DATA->delta_tau);
        rhob = grid_pt->rhob_t;
        eps = grid_pt->epsilon_t;
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        tildemu = muB/T;
        // f -= (grid_pt->rhob_prev);
        rhob = grid_pt->rhob;
        eps = grid_pt->epsilon;
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        tildemu_prev = muB/T;
        f = (tildemu - tildemu_prev)/(DATA->delta_tau);
        grid_pt->dUsup[0][m][0] = -f; /* g00 = -1 */
    }
    // Ends Sangyong's addition Nov 18 2014
    return 1;
}/* MakeDTau */

