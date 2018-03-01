// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_RECONST_H_
#define SRC_RECONST_H_

#include <array>
#include <iostream>
#include "util.h"
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "eos.h"
#include "pretty_ostream.h"
#include "data_struct.h"

class Reconst {
 private:
    const EOS &eos;
    const InitData &DATA;
    double eos_eps_max;
    pretty_ostream music_message;

    int max_iter;
    double rel_err, abs_err;

    double LARGE;

    int echo_level;
    double v_critical;

 public:
    Reconst(const EOS &eos, const InitData &DATA_in);

    ReconstCell ReconstIt_shell(double tau, TJbVec &q_vec, const Cell &grid_pt, int rk_flag);

    int ReconstIt(ReconstCell &grid_p, double tau, const TJbVec &q, const Cell &grid_pt,
                  int rk_flag);
    double GuessEps(double T00, double K00, double cs2);
    
    void revert_grid(ReconstCell &grid_current, const Cell &grid_prev, int rk_flag);

    int ReconstIt_velocity_iteration(ReconstCell &grid_p, double tau,
                                     const TJbVec &q, const Cell &grid_pt, int rk_flag);
    double reconst_velocity_f(double v, double T00, double M, double J0);
    double reconst_u0_f(double u0, double T00, double K00, double M,
                        double J0);

    int ReconstIt_velocity_Newton(ReconstCell &grid_p, double tau,
                                  const TJbVec &q, const Cell &grid_pt, int rk_flag);
    double reconst_velocity_f_Newton(double v, double T00, double M,
                                     double J0);
    double reconst_u0_f_Newton(double u0, double T00, double K00,
                               double M, double J0);
    double reconst_velocity_df(double v, double T00, double M, double J0);
    double reconst_u0_df(double u0, double T00, double K00, double M,
                         double J0);

    void regulate_grid(ReconstCell &grid_cell, double elocal);
};

#endif  // SRC_RECONST_H_
