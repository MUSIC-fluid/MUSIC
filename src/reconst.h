// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_RECONST_H_
#define SRC_RECONST_H_

#include <iostream>
#include "./util.h"
#include "./data.h"
#include "./grid.h"
#include "./eos.h"
#include "./pretty_ostream.h"

class Reconst {
 private:
    EOS *eos;
    double eos_eps_max;
    Util *util;
    pretty_ostream music_message;
    InitData *DATA_ptr;

    int max_iter;
    double rel_err, abs_err;

    double LARGE;

    int echo_level;
    double v_critical;

 public:
    Reconst(EOS *eos, InitData *DATA_in);
    ~Reconst();
      
    int ReconstIt_shell(Grid *grid_p, double tau, double *uq,
                        Grid *grid_pt, int rk_flag);

    int ReconstIt(Grid *grid_p, double tau, double *uq, Grid *grid_pt,
                  int rk_flag);
    double GuessEps(double T00, double K00, double cs2);
    
    void revert_grid(Grid *grid_current, Grid *grid_prev, int rk_flag);

    int ReconstIt_velocity_iteration(Grid *grid_p, double tau,
                                     double *uq, Grid *grid_pt, int rk_flag);
    double reconst_velocity_f(double v, double T00, double M, double J0);
    double reconst_u0_f(double u0, double T00, double K00, double M,
                        double J0);

    int ReconstIt_velocity_Newton(Grid *grid_p, double tau,
                                  double *uq, Grid *grid_pt, int rk_flag);
    double reconst_velocity_f_Newton(double v, double T00, double M,
                                     double J0);
    double reconst_u0_f_Newton(double u0, double T00, double K00,
                               double M, double J0);
    double reconst_velocity_df(double v, double T00, double M, double J0);
    double reconst_u0_df(double u0, double T00, double K00, double M,
                         double J0);

    void regulate_grid(Grid *grid_cell, double elocal, int rk_flag);
};

#endif  // SRC_RECONST_H_
