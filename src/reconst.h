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
    int reconst_type;
    EOS *eos;
    double eos_eps_max;
    Util *util;
    pretty_ostream music_message;
    InitData *DATA_ptr;

    int max_iter;
    double rel_err, abs_err;

 public:
    Reconst(EOS *eos, int reconst_type_in);
    ~Reconst();
      
    void ReconstError(const char *str, int i, int rk_flag, double *qi,
                      double **qi2, Grid *grid_pt);

    int ReconstIt_shell(Grid *grid_p, int direc, double tau, double **uq,
                        Grid *grid_pt, InitData *DATA, int rk_flag);

    // reconst_type == 0
    int ReconstIt(Grid *grid_p, int i, double tau, double **uq, Grid *grid_pt,
                  InitData *DATA, int rk_flag);
    double GuessEps(double T00, double K00, double cs2);
    
    void revert_grid(Grid *grid_current, Grid *grid_prev, int rk_flag);

    // reconst_type == 1
    int ReconstIt_velocity_iteration(Grid *grid_p, int direc, double tau,
                                     double **uq, Grid *grid_pt,
                                     InitData *DATA, int rk_flag);
    double reconst_velocity_f(double v, double T00, double M, double J0);
    double reconst_u0_f(double u0, double T00, double K00, double M,
                        double J0);

    // reconst_type == 2
    int ReconstIt_velocity_Newton(Grid *grid_p, int direc, double tau,
                                  double **uq, Grid *grid_pt,
                                  InitData *DATA, int rk_flag);
    double reconst_velocity_f_Newton(double v, double T00, double M,
                                     double J0);
    double reconst_u0_f_Newton(double u0, double T00, double K00,
                               double M, double J0);
    double reconst_velocity_df(double v, double T00, double M, double J0);
    double reconst_u0_df(double u0, double T00, double K00, double M,
                         double J0);
};

#endif  // SRC_RECONST_H_
