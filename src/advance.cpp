// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <omp.h>
#include "./util.h"
#include "./data.h"
#include "./grid.h"
#include "./reconst.h"
#include "./eos.h"
#include "./evolve.h"
#include "./advance.h"

using namespace std;

Advance::Advance(EOS *eosIn, InitData *DATA_in,
                 hydro_source *hydro_source_in) {
    DATA_ptr = DATA_in;
    eos = eosIn;
    grid = new Grid();
    util = new Util;
    reconst_ptr = new Reconst(eos, DATA_in->reconst_type);
    diss = new Diss(eosIn, DATA_in);
    minmod = new Minmod(DATA_in);
    u_derivative = new U_derivative(eosIn, DATA_in);
    if (DATA_in->Initial_profile == 12 || DATA_in->Initial_profile == 13
            || DATA_in->Initial_profile == 30) {
        flag_add_hydro_source = true;
        hydro_source_ptr = hydro_source_in;
    } else {
        flag_add_hydro_source = false;
        hydro_source_ptr = NULL;
    }

    grid_nx = DATA_in->nx;
    grid_ny = DATA_in->ny;
    grid_neta = DATA_in->neta;
    rk_order = DATA_in->rk_order;
}

// destructor
Advance::~Advance() {
    delete grid;
    delete util;
    delete diss;
    delete reconst_ptr;
    delete minmod;
    delete u_derivative;
}


// evolve Runge-Kutta step in tau
int Advance::AdvanceIt(double tau, InitData *DATA, Grid ***arena,
                       int rk_flag) {
    int ieta;
    #pragma omp parallel private(ieta)
    {
        #pragma omp for
        for (ieta=0; ieta < grid_neta; ieta++) {
	        AdvanceLocalT(tau, DATA, ieta, arena, rk_flag);
	    }/* ieta */
        #pragma omp barrier
    }
  
    if (DATA->viscosity_flag == 1) {
        #pragma omp parallel private(ieta)
        {
            #pragma omp for
	        for (ieta = 0; ieta < grid_neta; ieta++) {
		        AdvanceLocalW(tau, DATA, ieta, arena, rk_flag);
		    } /* ieta */
            #pragma omp barrier
        }
    }/* if viscosity flag is set */
    return 1;
}/* AdvanceIt */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%  Advance Local T %%%%%%%%%%%%%%%%%% */
int Advance::AdvanceLocalT(double tau, InitData *DATA, int ieta, Grid ***arena, 
                           int rk_flag) {
    // this function advances the ideal part in the transverse plane
    Grid grid_rk;
    double **qirk, *qi, *rhs, **w_rhs;

    qi = util->vector_malloc(5);
    rhs = util->vector_malloc(5);
    qirk = util->mtx_malloc(5, 4);
    w_rhs = util->mtx_malloc(5, 4);
    grid_rk.TJb = util->cube_malloc(rk_order, 5, 4);
    grid_rk.u = util->mtx_malloc(rk_order, 4);
    
    NbrQs NbrCells;
    InitNbrQs(&NbrCells);
    BdryCells HalfwayCells;
    InitTempGrids(&HalfwayCells, rk_order); 
 
    double eta_s_local = -DATA_ptr->eta_size/2. + ieta*DATA_ptr->delta_eta;
    for (int ix=0; ix <= grid_nx; ix++) {
        double x_local = -DATA_ptr->x_size/2. + ix*DATA_ptr->delta_x;
        for (int iy=0; iy <= grid_ny; iy++) {
            double y_local = -DATA_ptr->y_size/2. + iy*DATA_ptr->delta_y;
            FirstRKStepT(tau, x_local, y_local, eta_s_local,
                         DATA, &(arena[ieta][ix][iy]), rk_flag, qi, rhs, 
                         w_rhs, qirk, &grid_rk, &NbrCells, &HalfwayCells);
        }
    }

    util->cube_free(grid_rk.TJb, rk_order, 5, 4);
    util->mtx_free(grid_rk.u, rk_order, 4);
    util->mtx_free(qirk, 5, 4);
    util->mtx_free(w_rhs, 5, 4);
    util->vector_free(qi);
    util->vector_free(rhs);

    clean_Nbr_Qs(&NbrCells);
    clean_temp_grids(&HalfwayCells, rk_order);
    return 1; /* if successful */
}/* AdvanceLocalT */


/* %%%%%%%%%%%%%%%%% Advance Local W %%%%%%%%%% */
int Advance::AdvanceLocalW(double tau, InitData *DATA, int ieta, Grid ***arena,
                           int rk_flag) {
    Grid grid_rk;
    int flag = 0;
    double **qirk, *qi, *rhs, **w_rhs;

    qi = util->vector_malloc(5);
    rhs = util->vector_malloc(5);
    qirk = util->mtx_malloc(5, 4);
    w_rhs = util->mtx_malloc(5, 4);
    grid_rk.TJb = util->cube_malloc(rk_order, 5, 4);
    grid_rk.u = util->mtx_malloc(rk_order, 4);

    for (int ix=0; ix <= grid_nx; ix++) {
	    for (int iy=0; iy <= grid_ny; iy++) {
            flag = FirstRKStepW(tau, DATA, &(arena[ieta][ix][iy]), rk_flag, qi,
                                rhs, w_rhs, qirk, &grid_rk);
	    } /*iy */
	} /* ix */

    util->cube_free(grid_rk.TJb, rk_order, 5, 4);
    util->mtx_free(grid_rk.u, rk_order, 4);
    util->mtx_free(qirk, 5, 4);
    util->mtx_free(w_rhs, 5, 4);
    util->vector_free(qi);
    util->vector_free(rhs);
    return flag; 
}/* AdvanceLocalW */


/* %%%%%%%%%%%%%%%%%%%%%% First steps begins here %%%%%%%%%%%%%%%%%% */
int Advance::FirstRKStepT(double tau, double x_local, double y_local,
                          double eta_s_local,
                          InitData *DATA, Grid *grid_pt,
                          int rk_flag, double *qi, double *rhs, double **w_rhs,
                          double **qirk, Grid *grid_rk, NbrQs *NbrCells,
                          BdryCells *HalfwayCells) { 
    // this advances the ideal part
    double tau_now = tau;
    double tau_next = tau + (DATA->delta_tau);
    double tau_rk;
    if (rk_flag == 0) {
        tau_rk = tau_now;
    } else if (rk_flag == 1) {
        tau_rk = tau_next;
    } else {
        fprintf(stderr,"rk_flag = %d out of range.\n", rk_flag);
        exit(0);
    }

    // Solve partial_a T^{a mu} = -partial_a W^{a mu}
    // Update T^{mu nu}

    // MakeDelatQI gets
    //   qi = q0 if rk_flag = 0 or
    //   qi = q0 + k1 if rk_flag = 1
    // rhs[alpha] is what MakeDeltaQI outputs. 
    // It is the spatial derivative part of partial_a T^{a mu}
    // (including geometric terms)
    MakeDeltaQI(tau_rk, grid_pt, qi, rhs, DATA, rk_flag,
                NbrCells, HalfwayCells);

    double *j_mu = new double[4];
    for (int ii = 0; ii < 4; ii++) {
        j_mu[ii] = 0.0;
    }
    double rhob_source = 0.0;
    if (flag_add_hydro_source) {
        double *u_local = new double[4];
        for (int ii = 0; ii < 4; ii++) {
            u_local[ii] = grid_pt->u[rk_flag][ii];
        }
        hydro_source_ptr->get_hydro_energy_source(
                tau_rk, x_local, y_local, eta_s_local, u_local, j_mu);
        for (int ii = 0; ii < 4; ii++) {
            j_mu[ii] *= tau_rk;
        }
        if (DATA->turn_on_rhob == 1) {
            rhob_source = tau_rk*hydro_source_ptr->get_hydro_rhob_source(
                    tau_rk, x_local, y_local, eta_s_local, u_local);
        }
        delete[] u_local;
    }

    for (int alpha = 0; alpha < 5; alpha++) {
        qirk[alpha][0] = qi[alpha] + rhs[alpha];
        if (qirk[alpha][0] > LARGE) {
	        fprintf(stderr, "qirk[%d][0] = %e is a nan.\n", alpha,
                    qirk[alpha][0]);
	        fprintf(stderr, "qi[%d] = %e\n", alpha, qi[alpha]);
	        fprintf(stderr, "rhs[%d] = %e\n", alpha, rhs[alpha]);
        }
        // now MakeWSource returns partial_a W^{a mu}
        // (including geometric terms) 
        double dwmn = diss->MakeWSource(tau_rk, alpha, grid_pt, DATA, rk_flag);
        /* dwmn is the only one with the minus sign */
        qirk[alpha][0] -= dwmn*(DATA->delta_tau);

        if (flag_add_hydro_source) {
            // adding hydro_source terms
            if (alpha < 4) {
                qirk[alpha][0] += j_mu[alpha]*DATA->delta_tau;
            } else {
                qirk[alpha][0] += rhob_source*DATA->delta_tau;
            }
        }
     
        // set baryon density back to zero if viscous correction made it
        // non-zero remove/modify if rho_b!=0 
        // - this is only to remove the viscous correction that 
        // can make rho_b negative which we do not want.
        if (DATA->turn_on_rhob == 0) {
            if (alpha == 4 && fabs(qirk[alpha][0]) > 1e-12)
                qirk[alpha][0] = 0.;
        }

        /* if rk_flag > 0, we now have q0 + k1 + k2. 
         * So add q0 and multiply by 1/2 */
        if (rk_flag > 0) {
            qirk[alpha][0] += (grid_pt->TJb[0][alpha][0])*tau_now;
            qirk[alpha][0] *= 0.5;
        }
    }
    delete[] j_mu;

    int flag = 0;
    flag = reconst_ptr->ReconstIt_shell(grid_rk, 0, tau_next, qirk, grid_pt,
                                        DATA, rk_flag); 

    if (flag == 0) {
        reconst_ptr->ReconstError("grid_rk", 0, rk_flag+1, qi, qirk, grid_pt);
        return 0;
    } else {
        UpdateTJbRK(grid_rk, grid_pt, rk_flag); 
        /* TJb[rk_flag+1] is filled */
    }
    return 1;
}/* FirstRKStepT */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/*
   Done with T 
   Start W
*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int Advance::FirstRKStepW(double tau, InitData *DATA, Grid *grid_pt,
                          int rk_flag, double *qi, double *rhs, double **w_rhs,
                          double **qirk, Grid *grid_rk) { 
    double tau_now = tau;
    double tau_next = tau + (DATA->delta_tau);
  
    // Sangyong Nov 18 2014 implemented mu_max
    int mu_max;
    if (DATA->turn_on_rhob == 1)
        mu_max = 4;
    else 
        mu_max = 3;
 
    // Solve partial_a (u^a W^{mu nu}) = 0
    // Update W^{mu nu}
    // mu = 4 is the baryon current qmu

    // calculate delta uWmunu
    // need to use u[0][mu], remember rk_flag = 0 here
    // with the KT flux 
    // solve partial_tau (u^0 W^{kl}) = -partial_i (u^i W^{kl}
 
    /* Advance uWmunu */
    double tempf, temps;
    if (rk_flag == 0) {
        diss->Make_uWRHS(tau_now, grid_pt, w_rhs, DATA, rk_flag);
        for (int mu = 1; mu < 4; mu++) {
            for (int nu = mu; nu < 4; nu++) {
                int idx_1d = util->map_2d_idx_to_1d(mu, nu);
                tempf = ((grid_pt->Wmunu[rk_flag][idx_1d])
                         *(grid_pt->u[rk_flag][0]));
                temps = diss->Make_uWSource(tau_now, grid_pt, mu, nu, DATA,
                                            rk_flag); 
                tempf += temps*(DATA->delta_tau);
                tempf += w_rhs[mu][nu];
                grid_pt->Wmunu[rk_flag+1][idx_1d] = (
                                    tempf/(grid_pt->u[rk_flag+1][0]));
            }
        }
    } else if (rk_flag > 0) {
        diss->Make_uWRHS(tau_next, grid_pt, w_rhs, DATA, rk_flag);
        for (int mu = 1; mu < 4; mu++) {
            for (int nu = mu; nu < 4; nu++) {
                int idx_1d = util->map_2d_idx_to_1d(mu, nu);
                tempf = (grid_pt->Wmunu[0][idx_1d])*(grid_pt->u[0][0]);
                temps = diss->Make_uWSource(tau_next, grid_pt, mu, nu, DATA,
                                            rk_flag); 
                tempf += temps*(DATA->delta_tau);
                tempf += w_rhs[mu][nu];

                tempf += ((grid_pt->Wmunu[rk_flag][idx_1d])
                          *(grid_pt->u[rk_flag][0]));
                tempf *= 0.5;
       
                grid_pt->Wmunu[rk_flag+1][idx_1d] = (
                                            tempf/(grid_pt->u[rk_flag+1][0]));
            }
        }
    } /* rk_flag > 0 */

    if (DATA->turn_on_bulk == 1) {
        /* calculate delta u pi */
        double p_rhs;
        if (rk_flag == 0) {
            /* calculate delta u^0 pi */
            diss->Make_uPRHS(tau_now, grid_pt, &p_rhs, DATA, rk_flag);
   
            tempf = (grid_pt->pi_b[rk_flag])*(grid_pt->u[rk_flag][0]);
            temps = diss->Make_uPiSource(tau_now, grid_pt, DATA, rk_flag);
            tempf += temps*(DATA->delta_tau);
            tempf += p_rhs;
   
            grid_pt->pi_b[rk_flag+1] = tempf/(grid_pt->u[rk_flag+1][0]);
        } else if (rk_flag > 0) {
            /* calculate delta u^0 pi */
            diss->Make_uPRHS(tau_next, grid_pt, &p_rhs, DATA, rk_flag);
   
            tempf = (grid_pt->pi_b[0])*(grid_pt->u[0][0]);
            temps = diss->Make_uPiSource(tau_next, grid_pt, DATA, rk_flag);
            tempf += temps*(DATA->delta_tau);
            tempf += p_rhs;
  
            tempf += (grid_pt->pi_b[1])*(grid_pt->u[0][0]);
            tempf *= 0.5;

            grid_pt->pi_b[rk_flag+1] = tempf/(grid_pt->u[rk_flag+1][0]);
        }
    } else {
            grid_pt->pi_b[rk_flag+1] = 0.0;
    }

    // CShen: add source term for baryon diffusion
    if (DATA->turn_on_diff == 1) {
        if (rk_flag == 0) {
            diss->Make_uqRHS(tau_now, grid_pt, w_rhs, DATA, rk_flag);
            int mu = 4;
            for (int nu = 1; nu < 4; nu++) {
                int idx_1d = util->map_2d_idx_to_1d(mu, nu);
                tempf = ((grid_pt->Wmunu[rk_flag][idx_1d])
                         *(grid_pt->u[rk_flag][0]));
                temps = diss->Make_uqSource(tau_now, grid_pt, nu, DATA,
                                            rk_flag); 
                tempf += temps*(DATA->delta_tau);
                tempf += w_rhs[mu][nu];

                grid_pt->Wmunu[rk_flag+1][idx_1d] = (
                                            tempf/(grid_pt->u[rk_flag+1][0]));
            }
        } else if (rk_flag > 0) {
            diss->Make_uqRHS(tau_next, grid_pt, w_rhs, DATA, rk_flag);
            int mu = 4;
            for (int nu = 1; nu < 4; nu++) {
                int idx_1d = util->map_2d_idx_to_1d(mu, nu);
                tempf = (grid_pt->Wmunu[0][idx_1d])*(grid_pt->u[0][0]);
                temps = diss->Make_uqSource(tau_next, grid_pt, nu, DATA,
                                            rk_flag); 
                tempf += temps*(DATA->delta_tau);
                tempf += w_rhs[mu][nu];

                tempf += ((grid_pt->Wmunu[rk_flag][idx_1d])
                          *(grid_pt->u[rk_flag][0]));
                tempf *= 0.5;
       
                grid_pt->Wmunu[rk_flag+1][idx_1d] = (
                                        tempf/(grid_pt->u[rk_flag+1][0]));
            }
        } /* rk_flag > 0 */
    } else {
        for (int nu = 0; nu < 4; nu++) {
            int idx_1d = util->map_2d_idx_to_1d(4, nu);
            grid_pt->Wmunu[rk_flag+1][idx_1d] = 0.0;
        }
    }
   
    // re-make Wmunu[3][3] so that Wmunu[mu][nu] is traceless
    grid_pt->Wmunu[rk_flag+1][9] = (
            (2.*(grid_pt->u[rk_flag+1][1]*grid_pt->u[rk_flag+1][2]
                *grid_pt->Wmunu[rk_flag+1][5]
                + grid_pt->u[rk_flag+1][1]*grid_pt->u[rk_flag+1][3]
                  *grid_pt->Wmunu[rk_flag+1][6]
                + grid_pt->u[rk_flag+1][2]*grid_pt->u[rk_flag+1][3]
                  *grid_pt->Wmunu[rk_flag+1][8])
                - (grid_pt->u[rk_flag+1][0]*grid_pt->u[rk_flag+1][0] 
                   - grid_pt->u[rk_flag+1][1]*grid_pt->u[rk_flag+1][1])
                   *grid_pt->Wmunu[rk_flag+1][4] 
                - (grid_pt->u[rk_flag+1][0]*grid_pt->u[rk_flag+1][0] 
                   - grid_pt->u[rk_flag+1][2]*grid_pt->u[rk_flag+1][2])
                  *grid_pt->Wmunu[rk_flag+1][7])
            /(grid_pt->u[rk_flag+1][0]*grid_pt->u[rk_flag+1][0] 
              - grid_pt->u[rk_flag+1][3]*grid_pt->u[rk_flag+1][3]));

    // make Wmunu[i][0] using the transversality
    for (int mu = 1; mu < 4; mu++) {
        tempf = 0.0;
        for (int nu = 1; nu < 4; nu++) {
            int idx_1d = util->map_2d_idx_to_1d(mu, nu);
            tempf += (
                grid_pt->Wmunu[rk_flag+1][idx_1d]*grid_pt->u[rk_flag+1][nu]);
        }
        grid_pt->Wmunu[rk_flag+1][mu] = tempf/(grid_pt->u[rk_flag+1][0]);
    }

    // make Wmunu[0][0]
    tempf = 0.0;
    for (int nu=1; nu<4; nu++)
        tempf += grid_pt->Wmunu[rk_flag+1][nu]*grid_pt->u[rk_flag+1][nu]; 
    grid_pt->Wmunu[rk_flag+1][0] = tempf/(grid_pt->u[rk_flag+1][0]);
 
    if (DATA->turn_on_diff == 1) {
        // make qmu[0] using transversality
        for (int mu = 4; mu < mu_max + 1; mu++) {
            tempf = 0.0;
            for (int nu = 1; nu < 4; nu++) {
                int idx_1d = util->map_2d_idx_to_1d(mu, nu);
                tempf += (grid_pt->Wmunu[rk_flag+1][idx_1d]
                          *grid_pt->u[rk_flag+1][nu]);
            }
            grid_pt->Wmunu[rk_flag+1][10] = (
                                        tempf/(grid_pt->u[rk_flag+1][0]));
        }
    } else {
        grid_pt->Wmunu[rk_flag+1][10] = 0.0;
    }

    // If the energy density of the fluid element is smaller than 0.01GeV
    // reduce Wmunu using the QuestRevert algorithm
    int revert_flag = 0;
    int revert_q_flag = 0;
    if (DATA->Initial_profile != 0) {
        revert_flag = QuestRevert(tau, grid_pt, rk_flag, DATA);
        if (DATA->turn_on_diff == 1) {
            revert_q_flag = QuestRevert_qmu(tau, grid_pt, rk_flag, DATA);
        }
    }

    if (revert_flag == 1 || revert_q_flag == 1)
        return -1;
    else
        return 1;
}/* FirstRKStepW */

// update results after RK evolution to grid_pt
void Advance::UpdateTJbRK(Grid *grid_rk, Grid *grid_pt, int rk_flag) {
    int trk_flag = rk_flag+1;

    grid_pt->epsilon_t = grid_rk->epsilon;
    grid_pt->p_t = grid_rk->p;
    grid_pt->rhob_t = grid_rk->rhob;
    
    // reconstructed grid_rk uses rk_flag 0 only
    for (int mu=0; mu<4; mu++) {
        grid_pt->u[trk_flag][mu] = grid_rk->u[0][mu];
        for (int alpha=0; alpha<5; alpha++)  // alpha = 4 is the baryon current
            grid_pt->TJb[trk_flag][alpha][mu] = grid_rk->TJb[0][alpha][mu];
    }/* mu */
}/* UpdateTJbRK */

//! this function reduce the size of shear stress tensor and bulk pressure
//! in the dilute region to stablize numerical simulations
int Advance::QuestRevert(double tau, Grid *grid_pt, int rk_flag,
                         InitData *DATA) {
    int revert_flag = 0;
    const double energy_density_warning = 0.01;  // GeV/fm^3, T~100 MeV

    double eps_scale = 1.0;  // 1/fm^4
    double factor = 300.*tanh(grid_pt->epsilon/eps_scale);

    double pi_00 = grid_pt->Wmunu[rk_flag+1][0];
    double pi_01 = grid_pt->Wmunu[rk_flag+1][1];
    double pi_02 = grid_pt->Wmunu[rk_flag+1][2];
    double pi_03 = grid_pt->Wmunu[rk_flag+1][3];
    double pi_11 = grid_pt->Wmunu[rk_flag+1][4];
    double pi_12 = grid_pt->Wmunu[rk_flag+1][5];
    double pi_13 = grid_pt->Wmunu[rk_flag+1][6];
    double pi_22 = grid_pt->Wmunu[rk_flag+1][7];
    double pi_23 = grid_pt->Wmunu[rk_flag+1][8];
    double pi_33 = grid_pt->Wmunu[rk_flag+1][9];

    double pisize = (pi_00*pi_00 + pi_11*pi_11 + pi_22*pi_22 + pi_33*pi_33
                     - 2.*(pi_01*pi_01 + pi_02*pi_02 + pi_03*pi_03)
                     + 2.*(pi_12*pi_12 + pi_13*pi_13 + pi_23*pi_23));
  
    double pi_local = grid_pt->pi_b[rk_flag+1];
    double bulksize = 3.*pi_local*pi_local;

    double e_local = grid_pt->epsilon;
    double p_local = grid_pt->p;
    double eq_size = e_local*e_local + 3.*p_local*p_local;
       
    double rho_shear = sqrt(pisize/eq_size)/factor; 
    double rho_bulk  = sqrt(bulksize/eq_size)/factor;
 
    // Reducing the shear stress tensor 
    double rho_shear_max = 0.1;
    if (rho_shear > rho_shear_max) {
        if (e_local*hbarc > energy_density_warning) {
            printf("energy density = %lf -- |pi/(epsilon+3*P)| = %lf\n",
                   e_local*hbarc, rho_shear);
        }
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu; nu < 4; nu++) {
                int idx_1d = util->map_2d_idx_to_1d(mu, nu);
                grid_pt->Wmunu[rk_flag+1][idx_1d] = (
                    (rho_shear_max/rho_shear)
                    *grid_pt->Wmunu[rk_flag+1][idx_1d]);
            }
        }
        revert_flag = 1;
    }
   
    // Reducing bulk viscous pressure 
    double rho_bulk_max = 0.1;
    if (rho_bulk > rho_bulk_max) {
        if (e_local*hbarc > energy_density_warning) {
            printf("energy density = %lf --  |Pi/(epsilon+3*P)| = %lf\n",
                   e_local*hbarc, rho_bulk);
        }
        grid_pt->pi_b[rk_flag+1] = (
                (rho_bulk_max/rho_bulk)*grid_pt->pi_b[rk_flag+1]);
        revert_flag = 1;
    }

    return(revert_flag);
}/* QuestRevert */


//! this function reduce the size of net baryon diffusion current
//! in the dilute region to stablize numerical simulations
int Advance::QuestRevert_qmu(double tau, Grid *grid_pt, int rk_flag,
                             InitData *DATA) {
    int revert_flag = 0;
    const double energy_density_warning = 0.01;  // GeV/fm^3, T~100 MeV
    double eps_scale = 1.0;   // in 1/fm^4
    double factor = 300.*tanh(grid_pt->epsilon/eps_scale);

    double q_mu_local[4];
    for (int i = 0; i < 4; i++) {
        // copy the value from the grid
        int idx_1d = util->map_2d_idx_to_1d(4, i);
        q_mu_local[i] = grid_pt->Wmunu[rk_flag+1][idx_1d];
    }

    // calculate the size of q^\mu
    double q_size = 0.0;
    for (int i = 0; i < 4; i++) {
        double gfac = (i == 0 ? -1.0 : 1.0);
        q_size += gfac*q_mu_local[i]*q_mu_local[i];
    }

    // first check the positivity of q^mu q_mu 
    // (in the conversion of gmn = diag(-+++))
    if (q_size < 0.0) {
        cout << "Advance::QuestRevert_qmu: q^mu q_mu = " << q_size << " < 0!"
             << endl;
        cout << "Reset it to zero!!!!" << endl;
        for (int i = 0; i < 4; i++) {
            int idx_1d = util->map_2d_idx_to_1d(4, i);
            grid_pt->Wmunu[rk_flag+1][idx_1d] = 0.0;
        }
        revert_flag = 1;
    }

    // reduce the size of q^mu according to rhoB
    double e_local = grid_pt->epsilon;
    double rhob_local = grid_pt->rhob;
    double rho_q = sqrt(q_size/(rhob_local*rhob_local))/factor;
    double rho_q_max = 0.1;
    if (rho_q > rho_q_max) {
        if (e_local*hbarc > energy_density_warning) {
            printf("energy density = %lf, rhob = %lf -- |q/rhob| = %lf\n",
                   e_local*hbarc, rhob_local, rho_q);
        }
        for (int i = 0; i < 4; i++) {
            int idx_1d = util->map_2d_idx_to_1d(4, i);
            grid_pt->Wmunu[rk_flag+1][idx_1d] =
                                (rho_q_max/rho_q)*q_mu_local[i];
        }
        revert_flag = 1;
    }
    return(revert_flag);
}

void Advance::MakeDeltaQI(double tau, Grid *grid_pt, double *qi, double *rhs, 
			              InitData *DATA, int rk_flag, NbrQs *NbrCells,
                          BdryCells *HalfwayCells) {
    double delta[4];
    delta[1] = DATA->delta_x;
    delta[2] = DATA->delta_y;
    delta[3] = DATA->delta_eta;

    /* \partial_tau (tau Ttautau) + \partial_eta Tetatau 
            + \partial_x (tau Txtau) + \partial_y (tau Tytau) + Tetaeta = 0 */
    /* \partial_tau (tau Ttaueta) + \partial_eta Teteta 
            + \partial_x (tau Txeta) + \partial_y (tau Txeta) + Tetatau = 0 */
    /* \partial_tau (tau Txtau) + \partial_eta Tetax + \partial_x tau T_xx
            + \partial_y tau Tyx = 0 */

    // tau*Tmu0
    for (int alpha = 0; alpha < 5; alpha++) {
        qi[alpha] = grid_pt->TJb[rk_flag][alpha][0]*tau;
    }/* get qi first */

    double **DFmmp;
    DFmmp = util->mtx_malloc(5,4);
  
    /* implement Kurganov-Tadmor scheme */
    GetQIs(tau, grid_pt, qi, NbrCells, rk_flag, DATA);
 
    MakeQIHalfs(qi, NbrCells, HalfwayCells, grid_pt, DATA);
 
    ConstHalfwayCells(tau, HalfwayCells, qi, grid_pt, DATA, rk_flag);
 
    MakeKTCurrents(tau, DFmmp, grid_pt, HalfwayCells, rk_flag);
 
    for (int alpha = 0; alpha < 5; alpha++) {
        double sumf = 0.0; 
        for (int i = 1; i <= 3; i++) {
            sumf += DFmmp[alpha][i]/delta[i];
        } /* i */
        if (alpha == 0) {
            sumf -= grid_pt->TJb[rk_flag][3][3];
        } else if(alpha==3) {
            sumf -= grid_pt->TJb[rk_flag][3][0];
        }
        rhs[alpha] = sumf*(DATA->delta_tau);
    }/* alpha */

    util->mtx_free(DFmmp, 5, 4);
 
    return;
}/* MakeDeltaQI */


void Advance::GetQIs(double tau, Grid *grid_pt, double *qi, NbrQs *NbrCells,
                     int rk_flag, InitData *DATA) {
    double tempg = tau;
    for (int alpha = 0; alpha < 5; alpha++) {
        /* qs from the neighbors */
        /* implement outflow boundary condition - simply set the two outside
         * values the same as at the boundary */
        for (int i = 1; i <= 3; i++) {
            NbrCells->qip1[alpha][i] = (
                    tempg*grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0]);
            
            NbrCells->qip2[alpha][i] = (
                    tempg*grid_pt->nbr_p_2[i]->TJb[rk_flag][alpha][0]);
            
            NbrCells->qim1[alpha][i] = (
                    tempg*grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0]);
            
            NbrCells->qim2[alpha][i] = (
                    tempg*grid_pt->nbr_m_2[i]->TJb[rk_flag][alpha][0]);
        }/* i */    
    }/* alpha */
 
    return; 
}/* GetQIs */     


int Advance::MakeQIHalfs(double *qi, NbrQs *NbrCells, BdryCells *HalfwayCells,
			             Grid *grid_pt, InitData *DATA) {
    int alpha, direc;
    double fphL, fphR, fmhL, fmhR;
    double gphL, gphR, gmhL, gmhR;
    double tempf;

    for (alpha=0; alpha<5; alpha++) {
        for (direc=1; direc<=3; direc++) {
            gphL = qi[alpha];
            fphL = 0.5*minmod->minmod_dx(NbrCells->qip1[alpha][direc],
                                         qi[alpha],
                                         NbrCells->qim1[alpha][direc]);
      
            gphR = NbrCells->qip1[alpha][direc];
            fphR = -0.5*minmod->minmod_dx(NbrCells->qip2[alpha][direc],
                                          NbrCells->qip1[alpha][direc],
                                          qi[alpha]);
      
            gmhL = NbrCells->qim1[alpha][direc];
            fmhL = 0.5*minmod->minmod_dx(qi[alpha],
                                         NbrCells->qim1[alpha][direc],
                                         NbrCells->qim2[alpha][direc]);
      
            gmhR = qi[alpha];
            fmhR = -0.5*minmod->minmod_dx(NbrCells->qip1[alpha][direc],
                                          qi[alpha],
                                          NbrCells->qim1[alpha][direc]);

            tempf = HalfwayCells->qiphL[alpha][direc] = gphL + fphL;
            if (tempf > LARGE) {
                fprintf(stderr, "qiphL is not finite with g = %e, f = %e\n",
                        gphL, fphL);
                fprintf(stderr, "alpha = %d\n", alpha);
                fprintf(stderr, "direc = %d\n", direc);
                exit(0);
            }
            tempf = HalfwayCells->qiphR[alpha][direc] = gphR + fphR;
            if (tempf > LARGE) {
                fprintf(stderr, "qiphR is not finite with g = %e, f = %e\n",
                        gphR, fphR);
                fprintf(stderr, "alpha = %d\n", alpha);
                fprintf(stderr, "direc = %d\n", direc);
                exit(0);
            }
            tempf = HalfwayCells->qimhL[alpha][direc] = gmhL + fmhL;
            if (tempf > LARGE) {
                fprintf(stderr, "qimhL is not finite with g = %e, f = %e\n",
                        gmhL, fmhL);
                fprintf(stderr, "alpha = %d\n", alpha);
                fprintf(stderr, "direc = %d\n", direc);
                exit(0);
            }
            tempf = HalfwayCells->qimhR[alpha][direc] = gmhR + fmhR;
            if (tempf > LARGE) {
                fprintf(stderr, "qimhR is not finite with g = %e, f = %e\n",
                        gmhR, fmhR);
                fprintf(stderr, "alpha = %d\n", alpha);
                fprintf(stderr, "direc = %d\n", direc);
                exit(0);
            }
        }/* direc */
    }/* alpha */
    return 1; /* if successful */
}/* MakeQIHalfs */


int Advance::ConstHalfwayCells(
        double tau, BdryCells *HalfwayCells, double *qi, Grid *grid_pt,
        InitData *DATA, int rk_flag) {
    // this function reconstruct e, rhob, and u[4] for half way cells
    int flag = 0;
    for (int direc = 1; direc <= 3; direc++) {
        /* for each direction, reconstruct half-way cells */
        flag = reconst_ptr->ReconstIt_shell(
                    &(HalfwayCells->grid_p_h_L[direc]), direc, tau,
                    HalfwayCells->qiphL, grid_pt, DATA, rk_flag); 
        flag *= reconst_ptr->ReconstIt_shell(
                    &(HalfwayCells->grid_p_h_R[direc]), direc, tau,
                    HalfwayCells->qiphR, grid_pt, DATA, rk_flag); 
        flag *= reconst_ptr->ReconstIt_shell(
                    &(HalfwayCells->grid_m_h_L[direc]), direc, tau,
                    HalfwayCells->qimhL, grid_pt, DATA, rk_flag); 
        flag *= reconst_ptr->ReconstIt_shell(
                    &(HalfwayCells->grid_m_h_R[direc]), direc, tau,
                    HalfwayCells->qimhR, grid_pt, DATA, rk_flag); 
    } /* direc */

    return (flag);  /* upon successful execution */
}/* ConstHalfwayCells */


void Advance::MakeKTCurrents(double tau, double **DFmmp, Grid *grid_pt, 
			                 BdryCells *HalfwayCells, int rk_flag) {
    int i, alpha;
    double FiphL[5][4], FiphR[5][4], FimhL[5][4], FimhR[5][4];
    double Fiph[5][4], Fimh[5][4];
    double aiph[4], aimh[4], tau_fac[4], tempf;

    MakeMaxSpeedAs(tau, HalfwayCells, aiph, aimh, rk_flag);

    /* Current J^i is calculated from the halfway cells in the same
     * direction because what we need is the divergence d_i J^i */
    /* d_tau tau Jtau + d_x (tau Jx) + d_y (tau Jy) + d_eta Jeta */

    tau_fac[1] = tau;
    tau_fac[2] = tau;
    tau_fac[3] = 1.0;

    for (alpha=0; alpha<5; alpha++) {
        for (i=1; i<=3; i++) {
            // x_i current for TJb[0][alpha][0] 
            // from reconstructed halfway cells.
            // Reconst only uses TJb[0]
            FiphL[alpha][i] = 
                HalfwayCells->grid_p_h_L[i].TJb[0][alpha][i]*tau_fac[i];
            FiphR[alpha][i] = 
                HalfwayCells->grid_p_h_R[i].TJb[0][alpha][i]*tau_fac[i];
            FimhL[alpha][i] = 
                HalfwayCells->grid_m_h_L[i].TJb[0][alpha][i]*tau_fac[i]; 
            FimhR[alpha][i] = 
                HalfwayCells->grid_m_h_R[i].TJb[0][alpha][i]*tau_fac[i];
            // KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2})/2
            //                  - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2

            Fiph[alpha][i] = 0.5*(FiphL[alpha][i] + FiphR[alpha][i]);
            Fiph[alpha][i] -= 0.5*aiph[i]*(HalfwayCells->qiphR[alpha][i] 
                                           - HalfwayCells->qiphL[alpha][i]);
            Fimh[alpha][i] = 0.5*(FimhL[alpha][i] + FimhR[alpha][i]);
            Fimh[alpha][i] -= 0.5*aimh[i]*(HalfwayCells->qimhR[alpha][i] 
                                           - HalfwayCells->qimhL[alpha][i]);
    
            tempf= DFmmp[alpha][i] = (Fimh[alpha][i] - Fiph[alpha][i]);
            if (tempf > LARGE) {
                fprintf(stderr, "DFmmp[%d][%d] is not finite.\n", alpha, i);
                fprintf(stderr, "FimhL[%d][%d] is %e.\n",
                        alpha, i, FimhL[alpha][i]);
                fprintf(stderr, "FiphL[%d][%d] is %e.\n",
                        alpha, i, FiphL[alpha][i]);
                fprintf(stderr, "FimhR[%d][%d] is %e.\n",
                        alpha, i, FimhR[alpha][i]);
                fprintf(stderr, "FiphR[%d][%d] is %e.\n",
                        alpha, i, FiphR[alpha][i]);
                fprintf(stderr, "Fimh[%d][%d] is %e.\n",
                        alpha, i, Fimh[alpha][i]);
                fprintf(stderr, "Fiph[%d][%d] is %e.\n",
                        alpha, i, Fiph[alpha][i]);
                fprintf(stderr, "aimh[%d] = %e\n", i, aimh[i]);
                fprintf(stderr, "aiph[%d] = %e\n", i, aiph[i]);
                fprintf(stderr, "MakeKTCurrents: exiting.\n");
            }
        } /* i - loop over directions */
    } /* alpha */
} /* MakeKTCurrents */


/* Calculate the right-hand-side */
/*
 du/dt = (1/Delta x)(F_imh - F_iph) + D
 F_iph(a) = (1/2)(f_a(qiphR) + f_a(qiphL))- (1/2)aiph(qiphR - qiphL)
 F_imh(a) = (1/2)(f_a(qimhR) + f_a(qimhL))- (1/2)aimh(qimhR - qimhL)
 RK 2nd order Heun's rules:

 k1 = f(t, u);
 u1 = u + k1*h
 k2 = f(t+h, u1);
 u2 = u1 + k2*h = u + (k1+k2)*h;
 u_next = u + (1/2)(k1+k2)*h
        = (1/2)(u + u2);
*/
void Advance::MakeMaxSpeedAs(double tau, BdryCells *HalfwayCells,
                             double aiph[], double aimh[], int rk_flag) {
    /* Implement Kurganov-Tadmor */
    double aiphL[4], aiphR[4], aimhL[4], aimhR[4];
    for (int i = 1; i <= 3; i++) {
        aiphL[i] = MaxSpeed(tau, i, &(HalfwayCells->grid_p_h_L[i]), rk_flag);
        aiphR[i] = MaxSpeed(tau, i, &(HalfwayCells->grid_p_h_R[i]), rk_flag);
        aimhL[i] = MaxSpeed(tau, i, &(HalfwayCells->grid_m_h_L[i]), rk_flag);
        aimhR[i] = MaxSpeed(tau, i, &(HalfwayCells->grid_m_h_R[i]), rk_flag);
        
        aiph[i] = maxi(aiphL[i], aiphR[i]);
        aimh[i] = maxi(aimhL[i], aimhR[i]);
    }
    return;
}/* MakeMaxSpeedAs */

// determine the maximum signal propagation speed at the given direction
double Advance::MaxSpeed(double tau, int direc, Grid *grid_p, int rk_flag) {
    //grid_p = grid_p_h_L, grid_p_h_R, grid_m_h_L, grid_m_h_R
    //these are reconstructed by Reconst which only uses u[0] and TJb[0]
    double utau = (grid_p->u[0][0]);
    double utau2 = utau*utau;
    double ux = fabs((grid_p->u[0][direc]));
    double ux2 = ux*ux;
    double ut2mux2 = utau2 - ux2;
  
    double eps = grid_p->epsilon;
    double rhob = grid_p->rhob;
  
    double vs2 = eos->get_cs2(eps, rhob);

    double den = utau2*(1. - vs2) + vs2;
    double num_temp_sqrt = (ut2mux2 - (ut2mux2 - 1.)*vs2)*vs2;
    double num;
    if (num_temp_sqrt >= 0) {
        num = utau*ux*(1. - vs2) + sqrt(num_temp_sqrt);
    } else {
        double dpde = eos->p_e_func(eps, rhob);
        double p = grid_p->p;
        double h = p+eps;
        if (dpde < 0.001) {
            num = (sqrt(-(h*dpde*h*(dpde*(-1.0 + ut2mux2) - ut2mux2))) 
                   - h*(-1.0 + dpde)*utau*ux);
        } else {
            fprintf(stderr,"WARNING: in MaxSpeed. \n");
            fprintf(stderr, "Expression under sqrt in num=%lf. \n",
                    num_temp_sqrt);
            fprintf(stderr,"at value e=%lf. \n",eps);
            fprintf(stderr,"at value p=%lf. \n",p);
            fprintf(stderr,"at value h=%lf. \n",h);
            fprintf(stderr,"at value rhob=%lf. \n",rhob);
            fprintf(stderr,"at value utau=%lf. \n", utau);
            fprintf(stderr,"at value uk=%lf. \n", ux);
            fprintf(stderr,"at value vs^2=%lf. \n", vs2);
            fprintf(stderr,"at value dpde=%lf. \n", eos->p_e_func(eps, rhob));
            fprintf(stderr,"at value dpdrhob=%lf. \n",
                    eos->p_rho_func(eps, rhob));
            fprintf(stderr, "MaxSpeed: exiting.\n");
            exit(1);
        }
    }
    
    double f = num/(den + 1e-15);
    if (f < 0.0) {
        fprintf(stderr, "SpeedMax = %e\n is negative.\n", f);
        fprintf(stderr, "Can't happen.\n");
        exit(0);
    } else if (f <  ux/utau) {
        if (num != 0.0) {
            if (fabs(f-ux/utau)<0.0001) {
                f = ux/utau;
            } else {
	            fprintf(stderr, "SpeedMax-v = %lf\n", f-ux/utau);
	            fprintf(stderr, "SpeedMax = %e\n is smaller than v = %e.\n",
                        f, ux/utau);
	            fprintf(stderr, "Can't happen.\n");
	            exit(0);
            }
        }
    } else if (f >  1.0) {
        fprintf(stderr, "SpeedMax = %e\n is bigger than 1.\n", f);
        fprintf(stderr, "Can't happen.\n");
        fprintf(stderr, "SpeedMax = num/den, num = %e, den = %e \n", num, den);
        fprintf(stderr, "cs2 = %e \n", vs2);
        f =1.;
        exit(1);
    }
    if (direc == 3)
        f /= tau;

    return f;
}/* MaxSpeed */


void Advance::InitNbrQs(NbrQs *NbrCells) {
    (NbrCells->qip1) = util->mtx_malloc(5, 4);
    (NbrCells->qip2) = util->mtx_malloc(5, 4);
    (NbrCells->qim1) = util->mtx_malloc(5, 4);
    (NbrCells->qim2) = util->mtx_malloc(5, 4);
}/* InitNbrQs */

void Advance::clean_Nbr_Qs(NbrQs *NbrCells) {
    util->mtx_free((NbrCells->qip1), 5, 4);
    util->mtx_free((NbrCells->qip2), 5, 4);
    util->mtx_free((NbrCells->qim1), 5, 4);
    util->mtx_free((NbrCells->qim2), 5, 4);
}

void Advance::InitTempGrids(BdryCells *HalfwayCells, int rk_order) {
    int direc;

    (HalfwayCells->grid_p_h_L) = grid->grid_v_malloc(4);
    (HalfwayCells->grid_p_h_R) = grid->grid_v_malloc(4);
    (HalfwayCells->grid_m_h_L) = grid->grid_v_malloc(4);
    (HalfwayCells->grid_m_h_R) = grid->grid_v_malloc(4);

    HalfwayCells->qiphL = util->mtx_malloc(5, 4);
    HalfwayCells->qiphR = util->mtx_malloc(5, 4);
    HalfwayCells->qimhL = util->mtx_malloc(5, 4);
    HalfwayCells->qimhR = util->mtx_malloc(5, 4);
 
    for (direc=0; direc<4; direc++) {
        (HalfwayCells->grid_p_h_L)[direc].TJb =
                                            util->cube_malloc(rk_order, 5, 4);
        (HalfwayCells->grid_p_h_R)[direc].TJb =
                                            util->cube_malloc(rk_order, 5, 4);
        (HalfwayCells->grid_m_h_L)[direc].TJb =
                                            util->cube_malloc(rk_order, 5, 4);
        (HalfwayCells->grid_m_h_R)[direc].TJb =
                                            util->cube_malloc(rk_order, 5, 4);
    
        (HalfwayCells->grid_p_h_L)[direc].u = util->mtx_malloc(rk_order, 4);
        (HalfwayCells->grid_p_h_R)[direc].u = util->mtx_malloc(rk_order, 4);
        (HalfwayCells->grid_m_h_L)[direc].u = util->mtx_malloc(rk_order, 4);
        (HalfwayCells->grid_m_h_R)[direc].u = util->mtx_malloc(rk_order, 4);
    }
    return;
}/* InitTempGrids */

void Advance::clean_temp_grids(BdryCells *HalfwayCells, int rk_order) {
    for (int direc = 0; direc < 4; direc++) {
        util->cube_free((HalfwayCells->grid_p_h_L)[direc].TJb, rk_order, 5, 4);
        util->cube_free((HalfwayCells->grid_p_h_R)[direc].TJb, rk_order, 5, 4);
        util->cube_free((HalfwayCells->grid_m_h_L)[direc].TJb, rk_order, 5, 4);
        util->cube_free((HalfwayCells->grid_m_h_R)[direc].TJb, rk_order, 5, 4);
    
        util->mtx_free((HalfwayCells->grid_p_h_L)[direc].u, rk_order, 4);
        util->mtx_free((HalfwayCells->grid_p_h_R)[direc].u, rk_order, 4);
        util->mtx_free((HalfwayCells->grid_m_h_L)[direc].u, rk_order, 4);
        util->mtx_free((HalfwayCells->grid_m_h_R)[direc].u, rk_order, 4);
    }

    delete[] HalfwayCells->grid_p_h_L;
    delete[] HalfwayCells->grid_p_h_R;
    delete[] HalfwayCells->grid_m_h_L;
    delete[] HalfwayCells->grid_m_h_R;

    util->mtx_free(HalfwayCells->qiphL, 5, 4);
    util->mtx_free(HalfwayCells->qiphR, 5, 4);
    util->mtx_free(HalfwayCells->qimhL, 5, 4);
    util->mtx_free(HalfwayCells->qimhR, 5, 4);
    return;
}
