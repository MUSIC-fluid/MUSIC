// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef _SRC_CELL_H_
#define _SRC_CELL_H_
#include <iostream>
#include <iomanip>
#include "./data.h"

class Cell {
 public:
    double epsilon;
    double rhob;
    /* stress energy tensor plus baryon current  */
    /* TJb[flag][alpha][mu] */
    /* flag = 0 is the actual values. flag != 0 are the intermediate values
       for the Runge-Kutta step */
    /* alpha = 4 means Jb */
    // double correction;
        
    /* temporary values for the final RK update */
    double epsilon_t;
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
    
    Cell **nbr_p_1; 
    Cell **nbr_m_1; 
    Cell **nbr_p_2; 
    Cell **nbr_m_2; 
        
    /* we need to calculate partial_tau u[mu] */
    /* dU[flag][m][n] = u^{m,n} = partial^n u^m with the rk flag */
    /* note that they are superscripted. So partial^t = -partial_t */
    double ***dUsup; 

    /* shear part of the TJb with the rk_flag */
    double **Wmunu;
    double **prevWmunu; 
        
    /* bulk pressure */
    double *pi_b;
    double *prev_pi_b;

    // the following variables are for hyper-surface finder 
    // to determine freeze-out surface
    // they are only updated every freeze-out step not every evolution
    // time step
    // the one for the freeze-out surface finder for interpolation
    double epsilon_prev;
    double rhob_prev;
    double u_prev[4];
    double pi_b_prev;
    double *W_prev;
        
    Cell();//constructor
    Cell ***grid_c_malloc(int , int , int );
};

#endif  // SRC_GRID_H_
  
