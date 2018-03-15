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
    const InitData &DATA;
    const EOS &eos;
    const Minmod minmod;
    int map_2d_idx_to_1d(int a, int b) {
        static const int index_map[5][4] = {{0,   1,  2,  3},
                                            {1,   4,  5,  6},
                                            {2,   5,  7,  8},
                                            {3,   6,  8,  9},
                                            {10, 11, 12, 13}};
        return index_map[a][b];
    }

 public:
    Diss(const EOS &eosIn, const InitData &DATA_in);
    double MakeWSource(double tau, int alpha,
                       SCGrid &arena_current, SCGrid &arena_prev,
                       int ix, int iy, int ieta);

    int Make_uWRHS(double tau, SCGrid &arena, int ix, int iy, int ieta,
                   std::array< std::array<double,4>, 5> &w_rhs,
                   double theta_local, DumuVec &a_local);
    double Make_uWSource(double tau, Cell_small *grid_pt, Cell_small *grid_pt_prev,
                         int mu, int nu, int rk_flag, double theta_local,
                         DumuVec &a_local, VelocityShearVec &sigma_1d);

    int Make_uWRHS(double tau, SCGrid &arena, int ix, int iy, int ieta,
                   int mu, int nu, double &w_rhs,
                   double theta_local, DumuVec &a_local);

    int Make_uPRHS(double tau, SCGrid &arena, int ix, int iy, int ieta,
                   double *p_rhs, double theta_local);
    double Make_uPiSource(double tau, Cell_small *grid_pt, Cell_small *grid_pt_prev,
                          int rk_flag, double theta_local, VelocityShearVec &sigma_1d);

    double Make_uqRHS(double tau, SCGrid &arena_current, int ix, int iy, int ieta,
                      int mu, int nu);
    double Make_uqSource(double tau, Cell_small *grid_pt, Cell_small *grid_pt_prev, int nu,
                         int rk_flag, double theta_local, DumuVec &a_local,
                         VelocityShearVec &sigma_1d,
                         DmuMuBoverTVec &baryon_diffusion_vec);

    double get_temperature_dependent_eta_s(double T);
    double get_temperature_dependent_zeta_s(double temperature);

    void output_kappa_T_and_muB_dependence();
    void output_kappa_along_const_sovernB();
    void output_eta_over_s_T_and_muB_dependence();
    void output_eta_over_s_along_const_sovernB();
};

#endif  // SRC_DISSIPATIVE_H_
