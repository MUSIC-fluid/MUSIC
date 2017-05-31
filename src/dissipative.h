// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_DISSIPATIVE_H_
#define SRC_DISSIPATIVE_H_

#include <iostream>
#include "./util.h"
#include "./grid.h"
#include "./data.h"
#include "./minmod.h"

class Diss {
 private:
    EOS *eos;
    Minmod *minmod;
    Util *util;

 public:
    Diss(EOS *eosIn, InitData* DATA_in);
    ~Diss();
  
    double MakeWSource(double tau, int alpha, Grid *grid_pt, InitData *DATA,
                       int rk_flag);
    int Make_uWRHS(double tau, Grid *grid_pt, double **w_rhs, InitData *DATA,
                   int rk_flag, double theta_local, double *a_local);
    double Make_uWSource(double tau, Grid *grid_pt, int mu, int nu,
                         InitData *DATA, int rk_flag, double theta_local,
                         double *a_local, double *sigma_1d);
    
    int Make_uPRHS(double tau, Grid *grid_pt, double *p_rhs, InitData *DATA,
                   int rk_flag, double theta_local);

    double Make_uPiSource(double tau, Grid *grid_pt, InitData *DATA,
                          int rk_flag, double theta_local, double *sigma_1d);
    int Make_uqRHS(double tau, Grid *grid_pt, double **w_rhs, InitData *DATA,
                   int rk_flag);
    double Make_uqSource(double tau, Grid *grid_pt, int nu, InitData *DATA,
                         int rk_flag, double theta_local, double *a_local,
                         double *sigma_1d); 
    double get_temperature_dependent_eta_s(InitData *DATA, double T);
    double get_temperature_dependent_zeta_s(double temperature);

    void output_kappa_T_and_muB_dependence(InitData *DATA);
    void output_kappa_along_const_sovernB(InitData *DATA);
    void output_eta_over_s_T_and_muB_dependence(InitData *DATA);
    void output_eta_over_s_along_const_sovernB(InitData *DATA);
};

#endif  // SRC_DISSIPATIVE_H_
