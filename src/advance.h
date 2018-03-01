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

//! advance routines separate for
//! T^{0 nu} \del T^{i\nu} (T)
//! W
//! T^{0 nu} with W source (TS)
//! W with source (WS)
class Advance {
 private:
    InitData* DATA_ptr;
    Diss *diss;
    Reconst *reconst_ptr;
    EOS *eos;
    Minmod minmod;
    U_derivative *u_derivative_ptr;
    hydro_source *hydro_source_ptr;
    pretty_ostream music_message;

    bool flag_add_hydro_source;


 public:
    Advance(EOS *eosIn, InitData* DATA_in, hydro_source *hydro_source_in);
    ~Advance();

    int AdvanceIt(double tau_init, InitData *DATA, Grid &arena, int rk_flag);

    void update_small_cell_to_cell(Cell &c, const Cell_small &c_s, int rk_flag);
    void update_cell_to_small_cell(const Cell &c, Cell_small &c_s, int rk_flag);
    void FirstRKStepT(const double tau, double x_local, double y_local,
                     double eta_s_local, InitData *DATA, SCGrid &arena_current, SCGrid &arena_future, SCGrid &arena_prev, int ix, int iy, int ieta,
                     int rk_flag);

    void FirstRKStepW(double tau_it, InitData *DATA, Grid &arena,
                      SCGrid &arena_prev, SCGrid &arena_current, SCGrid &arena_future,
                     int rk_flag, double theta_local, DumuVec &a_local,
                     VelocityShearVec &sigma_local, int ieta, int ix, int iy);

    void UpdateTJbRK(const ReconstCell &grid_rk, Cell *grid_pt, int rk_flag);
    void UpdateTJbRK(const ReconstCell &grid_rk, Cell_small &grid_pt);
    int QuestRevert(double tau, Cell *grid_pt, int rk_flag, InitData *DATA,
                    int ieta, int ix, int iy);
    int QuestRevert_qmu(double tau, Cell *grid_pt, int rk_flag,
                        InitData *DATA, int ieta, int ix, int iy);

    void MakeDeltaQI(double tau, SCGrid &arena_current, int ix, int iy, int ieta, TJbVec &qi, int rk_flag);
    double MaxSpeed(double tau, int direc, const ReconstCell &grid_p);
    double get_TJb(const Cell &grid_p, const int rk_flag, const int mu, const int nu);
    double get_TJb(const ReconstCell &grid_p, const int rk_flag, const int mu, const int nu);
    double get_TJb(const Cell_small &grid_p, const int mu, const int nu);
};

#endif  // SRC_ADVANCE_H_
