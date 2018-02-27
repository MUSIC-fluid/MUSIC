// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef _SRC_CELL_H_
#define _SRC_CELL_H_

#include <array>

class Cell {
 public:
    double epsilon = 0;
    double rhob    = 0;
    /* stress energy tensor plus baryon current  */
    /* TJb[flag][alpha][mu] */
    /* flag = 0 is the actual values. flag != 0 are the intermediate values
       for the Runge-Kutta step */
    /* alpha = 4 means Jb */
    // double correction;
        
    /* temporary values for the final RK update */
    double epsilon_t = 0;
    double rhob_t    = 0;
        
    // store the epsilon and rhob at previous time step
    double prev_epsilon = 0;
    double prev_rhob    = 0;

    /* u[flag][mu]: flag=0 is the actual values. flag != are for RK steps */
    std::array<double, 4> u_rk0;
    std::array<double, 4> u_rk1;
        
    /* to include shear viscosity */
    /* we need to calculate partial_tau u[mu] */
    std::array<double, 4> prev_u;
    /* u[mu] from the previous time step including the rk flag */
    //double **pprev_u; /* u[mu] from 2 time step ago including the rk flag */
    //double **dU; /* dU[m][n] = partial_m u_n at the current time */
    //double ***pimunu; /* Stress part of the TJb */
    
    /* we need to calculate partial_tau u[mu] */
    /* dU[flag][m][n] = u^{m,n} = partial^n u^m with the rk flag */
    /* note that they are superscripted. So partial^t = -partial_t */
    std::array<std::array<double, 5>, 4> dUsup;

    /* shear part of the TJb with the rk_flag */
    std::array<double, 14> Wmunu_rk0;
    std::array<double, 14> Wmunu_rk1;
    std::array<double, 14> prevWmunu; 
        
    /* bulk pressure */
    double pi_b_rk0  = 0.0;
    double pi_b_rk1  = 0.0;
    double prev_pi_b = 0.0;

    // the following variables are for hyper-surface finder 
    // to determine freeze-out surface
    // they are only updated every freeze-out step not every evolution
    // time step
    // the one for the freeze-out surface finder for interpolation
    double epsilon_prev;
    double rhob_prev;
    std::array<double, 4> u_prev;
    double pi_b_prev;
    std::array<double, 4> W_prev;
};

#endif  // SRC_GRID_H_
