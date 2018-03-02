// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_RECONST_H_
#define SRC_RECONST_H_

#include <array>
#include <iostream>
#include "util.h"
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "eos.h"
#include "pretty_ostream.h"
#include "data_struct.h"

class Reconst {
 private:
    const EOS &eos;
    const InitData &DATA;
    double eos_eps_max;
    pretty_ostream music_message;

    int max_iter;
    double rel_err, abs_err;

    double LARGE;

    int echo_level;
    double v_critical;

 public:
    Reconst(const EOS &eos, const InitData &DATA_in);

    ReconstCell ReconstIt_shell(double tau, TJbVec &q_vec, const Cell_small &grid_pt, int rk_flag);

    int ReconstIt(ReconstCell &grid_p, double tau, const TJbVec &q, const Cell &grid_pt,
                  int rk_flag);
    double GuessEps(double T00, double K00, double cs2);
    
    void revert_grid(ReconstCell &grid_current, const Cell &grid_prev, int rk_flag);
    void revert_grid(ReconstCell &grid_current, const Cell_small &grid_prev, int rk_flag);

    int ReconstIt_velocity_iteration(ReconstCell &grid_p, double tau,
                                     const TJbVec &q, const Cell &grid_pt, int rk_flag);
    double reconst_velocity_f(double v, double T00, double M, double J0);
    double reconst_u0_f(double u0, double T00, double K00, double M,
                        double J0);

    int ReconstIt_velocity_Newton(ReconstCell &grid_p, double tau,
                                  const TJbVec &q, const Cell &grid_pt, int rk_flag);
    int ReconstIt_velocity_Newton(ReconstCell &grid_p, double tau,
                                  const TJbVec &q, const Cell_small &grid_pt, int rk_flag);

#pragma omp declare simd
    void reconst_velocity_fdf(const double v, const double T00, const double M, const double J0, double &fv, double &dfdv) const{
        // this function returns f(v) = M/(M0 + P)
        const double epsilon = T00 - v*M;
        const double temp    = sqrt(1. - v*v);
        const double rho     = J0*temp;
   
        const double pressure = eos.get_pressure(epsilon, rho);
        const double temp1    = T00 + pressure;
        const double temp2    = v/temp;
        const double dPde     = eos.p_e_func(epsilon, rho);
        const double dPdrho   = eos.p_rho_func(epsilon, rho);

        fv   = v - M/temp1;
        dfdv = 1. - M/(temp1*temp1)*(M*dPde + J0*temp2*dPdrho);
    }

#pragma omp declare simd
    void reconst_u_fdf(const double u0, const double T00, const double K00, const double M, const double J0, double &fu, double &dfdu) const {
        const double v       = sqrt(1. - 1./(u0*u0));
        const double epsilon = T00 - v*M;
        const double rho     = J0/u0;

        const double pressure = eos.get_pressure(epsilon, rho);
        const double temp1 = T00 + pressure;
        const double temp0 = sqrt(temp1*temp1 - K00);

        const double dedu0   = - M/(u0*u0*u0*v);
        const double drhodu0 = - J0/(u0*u0);

        const double dPde     = eos.p_e_func(epsilon, rho);
        const double dPdrho   = eos.p_rho_func(epsilon, rho);

        const double denorm = temp0*temp0*temp0;
        
        fu = u0 - temp1/temp0;
        dfdu  = 1. + (dedu0*dPde + drhodu0*dPdrho)*K00/denorm;
    }

    double reconst_velocity_f_Newton(double v, double T00, double M,
                                     double J0);
    double reconst_u0_f_Newton(double u0, double T00, double K00,
                               double M, double J0);
    double reconst_velocity_df(double v, double T00, double M, double J0);
    double reconst_u0_df(double u0, double T00, double K00, double M,
                         double J0);

    void regulate_grid(ReconstCell &grid_cell, double elocal);
};

#endif  // SRC_RECONST_H_
