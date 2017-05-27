// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_ADVANCE_H_
#define SRC_ADVANCE_H_

#include <iostream>
#include "./data.h"
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
    hydro_source *hydro_source_ptr;
    pretty_ostream music_message;

    int grid_nx, grid_ny, grid_neta;
    int rk_order;
    bool flag_add_hydro_source;


 public:
    Advance(EOS *eosIn, InitData* DATA_in, hydro_source *hydro_source_in);
    ~Advance();

    int AdvanceIt(double tau_init, InitData *DATA, Grid ***arena, int rk_flag);

    
    int AdvanceLocalT(double tau_init, InitData *DATA, int ieta, Grid ***arena,
                      int rk_flag);
    int AdvanceLocalW(double tau_init, InitData *DATA, int ieta, Grid ***arena,
                      int rk_flag);

    int FirstRKStepT(double tau, double x_local, double y_local,
                     double eta_s_local, InitData *DATA, Grid *grid_pt,
                     int rk_flag);

    int FirstRKStepW(double tau_it, InitData *DATA, Grid *grid_pt,
                     int rk_flag);

    void UpdateTJbRK(Grid *grid_rk, Grid *grid_pt, int rk_flag);
    int QuestRevert(double tau, Grid *grid_pt, int rk_flag, InitData *DATA);
    int QuestRevert_qmu(double tau, Grid *grid_pt, int rk_flag,
                        InitData *DATA);

    void MakeDeltaQI(double tau, Grid *grid_pt, double *qi, int rk_flag);
    double MaxSpeed(double tau, int direc, Grid *grid_p);
    double get_TJb(Grid *grid_p, int rk_flag, int mu, int nu);
};

#endif  // SRC_ADVANCE_H_
