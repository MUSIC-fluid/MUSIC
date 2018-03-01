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
    const EOS &eos;
    const Minmod minmod;

 public:
    Diss(const EOS &eosIn, const InitData &DATA_in);
  
    double MakeWSource(double tau, int alpha, Grid &arena,
                       int ix, int iy, int ieta, const InitData &DATA, int rk_flag);
    int Make_uWRHS(double tau, Grid &arena, int ix, int iy, int ieta,
                   std::array< std::array<double,4>, 5> &w_rhs, const InitData &DATA,
                   int rk_flag, double theta_local, DumuVec &a_local);
    double Make_uWSource(double tau, Cell *grid_pt, int mu, int nu,
                         const InitData &DATA, int rk_flag, double theta_local,
                         DumuVec &a_local, VelocityShearVec &sigma_1d);
    
    int Make_uPRHS(double tau, Grid &arena, int ix, int iy, int ieta,
                   double *p_rhs, const InitData &DATA,
                   int rk_flag, double theta_local);

    double Make_uPiSource(double tau, Cell *grid_pt, const InitData &DATA,
                          int rk_flag, double theta_local, VelocityShearVec &sigma_1d);
    int Make_uqRHS(double tau, Grid &arena, int ix, int iy, int ieta,
                   std::array< std::array<double,4>, 5> &w_rhs, const InitData &DATA, int rk_flag);
    double Make_uqSource(double tau, Cell *grid_pt, int nu, const InitData &DATA,
                         int rk_flag, double theta_local, DumuVec &a_local,
                         VelocityShearVec &sigma_1d); 
    double get_temperature_dependent_eta_s(const InitData &DATA, double T);
    double get_temperature_dependent_zeta_s(double temperature);

    void output_kappa_T_and_muB_dependence(const InitData &DATA);
    void output_kappa_along_const_sovernB(const InitData &DATA);
    void output_eta_over_s_T_and_muB_dependence(const InitData &DATA);
    void output_eta_over_s_along_const_sovernB(const InitData &DATA);
};

#endif  // SRC_DISSIPATIVE_H_
