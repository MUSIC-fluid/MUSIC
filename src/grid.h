// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_GRID_H_
#define SRC_GRID_H_
#include <iostream>
#include <iomanip>
#include "./data.h"
#include "./eos.h"

class Grid {
 public:
    double epsilon;
    double p;
    double rhob;
    double ***TJb; 
    /* stress energy tensor plus baryon current  */
    /* TJb[flag][alpha][mu] */
    /* flag = 0 is the actual values. flag != 0 are the intermediate values
       for the Runge-Kutta step */
    /* alpha = 4 means Jb */
    // double correction;
        
    /* temporary values for the final RK update */
    double epsilon_t;
    double p_t;
    double rhob_t;
        
    // store the epsilon and rhob at previous time step
    double prev_epsilon;
    double prev_rhob;

    /* u[flag][mu]: flag=0 is the actual values. flag != are for RK steps */
    double **u;
        
    /* to include shear viscosity */
    /* we need to calculate partial_tau u[mu] */
    double **prev_u;
    /* u[mu] from the previous time step including the rk flag */
    //double **pprev_u; /* u[mu] from 2 time step ago including the rk flag */
    //double **dU; /* dU[m][n] = partial_m u_n at the current time */
    //double ***pimunu; /* Stress part of the TJb */
    
    Grid **nbr_p_1; 
    Grid **nbr_m_1; 
    Grid **nbr_p_2; 
    Grid **nbr_m_2; 
        
    /* This is u^mu partial_mu u^nu */
    double **a;
        
    /* this is the expansion rate partial_mu u^mu */
    double *theta_u;

    // the velocity shear tensor
    // only record 10 essential components in a 1-d array
    double **sigma;
        
    /* we need to calculate partial_tau u[mu] */
    /* dU[flag][m][n] = u^{m,n} = partial^n u^m with the rk flag */
    /* note that they are superscripted. So partial^t = -partial_t */
    double ***dUsup; 

    double ***Wmunu; /* shear part of the TJb with the rk_flag */
    double ***prevWmunu; 
    //double ***pprevWmunu; 
        
    double *pi_b; /* bulk pressure */
    double *prev_pi_b; /* bulk pressure */
    //double *pprev_pi_b; /* bulk pressure */

    // the following variables are for hyper-surface finder 
    // to determine freeze-out surface
    // they are only updated every freeze-out step not every evolution
    // time step
    // the one for the freeze-out surface finder for interpolation
    double epsilon_prev;
    double rhob_prev;
    double u_prev[4];
    double pi_b_prev;
    double **W_prev;
        
    Grid();//constructor
    ~Grid();//destructor
    Grid *grid_v_malloc(int );
    Grid **grid_m_malloc(int , int );
    Grid ***grid_c_malloc(int , int , int );
};

#endif  // SRC_GRID_H_
  
