// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_DISSIPATIVE_H_
#define SRC_DISSIPATIVE_H_

#include <array>
#include <iostream>
#include "util.h"
#include "cell.h"
#include "grid.h"
#include "data.h"
#include "minmod.h"

class Diss {
 private:
    EOS *eos;
    Minmod minmod;

 public:
    Diss(EOS *eosIn, InitData* DATA_in);
  
    double MakeWSource(double tau, int alpha, SCGrid &arena_current, SCGrid &arena_prev, int ix, int iy, int ieta,
                         InitData *DATA, int rk_flag);
    int Make_uWRHS(double tau, SCGrid &arena, int ix, int iy, int ieta,
                     std::array< std::array<double,4>, 5> &w_rhs, InitData *DATA,
                     double theta_local, DumuVec &a_local);
    double Make_uWSource(double tau, Cell_small *grid_pt, Cell_small *grid_pt_prev, int mu, int nu,
                         InitData *DATA, int rk_flag, double theta_local,
                         DumuVec &a_local, VelocityShearVec &sigma_1d);
    
    int Make_uPRHS(double tau, SCGrid &arena, int ix, int iy, int ieta,
                   double *p_rhs, InitData *DATA,
                   int rk_flag, double theta_local);

    double Make_uPiSource(double tau, Cell_small *grid_pt, Cell_small *grid_pt_prev, InitData *DATA,
                        int rk_flag, double theta_local, VelocityShearVec &sigma_1d);

    int Make_uqRHS(double tau, SCGrid &arena_current, int ix, int iy, int ieta,
                   std::array< std::array<double,4>, 5> &w_rhs, InitData *DATA, int rk_flag);
    double Make_uqSource(double tau, Cell_small *grid_pt, Cell_small *grid_pt_prev, int nu, InitData *DATA,
                         int rk_flag, double theta_local, DumuVec &a_local,
                         VelocityShearVec &sigma_1d); 
    double get_temperature_dependent_eta_s(InitData *DATA, double T);
    double get_temperature_dependent_zeta_s(double temperature);

    void output_kappa_T_and_muB_dependence(InitData *DATA);
    void output_kappa_along_const_sovernB(InitData *DATA);
    void output_eta_over_s_T_and_muB_dependence(InitData *DATA);
    void output_eta_over_s_along_const_sovernB(InitData *DATA);
};

#endif  // SRC_DISSIPATIVE_H_
