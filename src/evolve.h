// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_EVOLVE_H_
#define SRC_EVOLVE_H_

#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "util.h"
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "grid_info.h"
#include "eos.h"
#include "advance.h"
#include "hydro_source.h"
#include "u_derivative.h"
#include "pretty_ostream.h"

// this is a control class for the hydrodynamic evolution
class Evolve {
 private:
    const EOS &eos;        // declare EOS object
    const InitData &DATA;

    Cell_info grid_info;
    Advance advance;
    U_derivative u_derivative;
    hydro_source *hydro_source_ptr = nullptr;
    pretty_ostream music_message;


    // simulation information
    int rk_order;

    double SUM;
    int warnings;
    int cells;
    int weirdCases;
    int facTau;

    // information about freeze-out surface
    // (only used when freezeout_method == 4)
    int n_freeze_surf;
    vector<double> epsFO_list;

 public:
    Evolve(EOS *eos, InitData *DATA_in, hydro_source *hydro_source_in);
    int EvolveIt(InitData *DATA, SCGrid &arena_prev, SCGrid &arena_current, SCGrid &arena_future);

    int AdvanceRK(double tau, InitData *DATA, SCGrid &arena_prev,
                  SCGrid &arena_current, SCGrid &arena_future);

    int FreezeOut_equal_tau_Surface(double tau, Grid &arena);
    void FreezeOut_equal_tau_Surface_XY(double tau,
                                        int ieta, Grid &arena,
                                        int thread_id, double epsFO);
    // void FindFreezeOutSurface(double tau, const InitData &DATA,
    //                          Grid &arena, int size, int rank);
    // void FindFreezeOutSurface2(double tau, const InitData &DATA,
    //                           Grid &arena, int size, int rank);
    // int FindFreezeOutSurface3(double tau, const InitData &DATA,
    //                          Grid &arena, int size, int rank);
    int FindFreezeOutSurface_Cornelius(double tau, Grid &arena);
    int FindFreezeOutSurface_Cornelius_XY(double tau,
                                          int ieta, Grid &arena,
                                          int thread_id, double epsFO);
    int FindFreezeOutSurface_boostinvariant_Cornelius(
                                    double tau, Grid &arena);

    void store_previous_step_for_freezeout(Grid &arena);

    void regulate_qmu(double* u, double* q, double* q_regulated);
    void regulate_Wmunu(double* u, double** Wmunu, double** Wmunu_regulated);

    void initialize_freezeout_surface_info();
};

#endif  // SRC_EVOLVE_H_

