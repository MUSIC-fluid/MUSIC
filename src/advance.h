// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_ADVANCE_H_
#define SRC_ADVANCE_H_

#include <iostream>
#include "./data.h"
#include "./cell.h"
#include "./grid.h"
#include "./dissipative.h"
#include "./minmod.h"
#include "./u_derivative.h"
#include "./reconst.h"
#include "./hydro_source.h"
#include "./pretty_ostream.h"

class Advance {
 private:
    const InitData &DATA;
    const EOS &eos;
    Diss *diss;
    Reconst *reconst_ptr;
    Minmod minmod;
    hydro_source *hydro_source_ptr;
    pretty_ostream music_message;

    bool flag_add_hydro_source;

    int map_2d_idx_to_1d(int a, int b) {
        static const int index_map[5][4] = {{0,   1,  2,  3},
                                            {1,   4,  5,  6},
                                            {2,   5,  7,  8},
                                            {3,   6,  8,  9},
                                            {10, 11, 12, 13}};
        return index_map[a][b];
    }

 public:
    Advance(const EOS &eosIn, const InitData &DATA_in, hydro_source *hydro_source_in);
    ~Advance();

    void AdvanceIt(double tau_init,
                   SCGrid &arena_prev, SCGrid &arena_current, SCGrid &arena_future,
                   int rk_flag);

    void FirstRKStepT(const double tau, double x_local, double y_local,
                      double eta_s_local,  SCGrid &arena_current, SCGrid &arena_future, SCGrid &arena_prev, int ix, int iy, int ieta,
                      int rk_flag);

    void FirstRKStepW(double tau_it, SCGrid &arena_prev, SCGrid &arena_current, SCGrid &arena_future,
                      int rk_flag, double theta_local, DumuVec &a_local,
                      VelocityShearVec &sigma_local, DmuMuBoverTVec &baryon_diffusion_vector, int ieta, int ix, int iy);

    void UpdateTJbRK(const ReconstCell &grid_rk, Cell_small &grid_pt);
    void QuestRevert(double tau, Cell_small *grid_pt, int ieta, int ix, int iy);
    void QuestRevert_qmu(double tau, Cell_small *grid_pt,
                         int ieta, int ix, int iy);

    void MakeDeltaQI(double tau, SCGrid &arena_current,
                     int ix, int iy, int ieta, TJbVec &qi, int rk_flag);
    double MaxSpeed(double tau, int direc, const ReconstCell &grid_p);
    double get_TJb(const Cell &grid_p, const int rk_flag, const int mu, const int nu);
    double get_TJb(const ReconstCell &grid_p, const int rk_flag, const int mu, const int nu);
    double get_TJb(const Cell_small &grid_p, const int mu, const int nu);
};

#endif  // SRC_ADVANCE_H_
