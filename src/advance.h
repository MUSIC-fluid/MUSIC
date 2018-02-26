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
    Util *util;
    Diss *diss;
    Reconst *reconst_ptr;
    EOS *eos;
    Minmod *minmod;
    U_derivative *u_derivative_ptr;
    hydro_source *hydro_source_ptr;
    pretty_ostream music_message;

    int grid_nx, grid_ny, grid_neta;
    int rk_order;
    bool flag_add_hydro_source;


 public:
    Advance(EOS *eosIn, InitData* DATA_in, hydro_source *hydro_source_in);
    ~Advance();

    int AdvanceIt(double tau_init, InitData *DATA, Grid &arena, int rk_flag);

    int FirstRKStepT(double tau, double x_local, double y_local,
                     double eta_s_local, InitData *DATA, Cell *grid_pt,
                     int rk_flag);

    int FirstRKStepW(double tau_it, InitData *DATA, Cell *grid_pt,
                     int rk_flag, double theta_local, double* a_local,
                     double *sigma_local, int ieta, int ix, int iy);

    void UpdateTJbRK(Cell *grid_rk, Cell *grid_pt, int rk_flag);
    int QuestRevert(double tau, Cell *grid_pt, int rk_flag, InitData *DATA,
                    int ieta, int ix, int iy);
    int QuestRevert_qmu(double tau, Cell *grid_pt, int rk_flag,
                        InitData *DATA, int ieta, int ix, int iy);

    void MakeDeltaQI(double tau, Cell *grid_pt, double *qi, int rk_flag);
    double MaxSpeed(double tau, int direc, Cell *grid_p);
    double get_TJb(Cell *grid_p, int rk_flag, int mu, int nu);
};

#endif  // SRC_ADVANCE_H_
