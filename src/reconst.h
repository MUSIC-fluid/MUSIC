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

    const int max_iter;
    const double rel_err;
    const double abs_err;

    const double LARGE;

    int echo_level;
    const double v_critical;

 public:
    Reconst() = default;
    Reconst(const EOS &eos, const InitData &DATA_in);

    ReconstCell ReconstIt_shell(double tau, const TJbVec &tauq_vec,
                                const Cell_small &grid_pt);

    void revert_grid(ReconstCell &grid_current,
                     const Cell_small &grid_prev) const;

    int ReconstIt_velocity_Newton(ReconstCell &grid_p, double tau,
                                  const TJbVec &q, const Cell_small &grid_pt);
    
    void reconst_velocity_fdf(const double v, const double T00, const double M,
                              const double J0, double &fv, double &dfdv) const;

    void reconst_u0_fdf(const double u0, const double T00, const double K00,
                        const double M, const double J0,
                        double &fu0, double &dfdu0) const;

    int solve_velocity_Newton(const double v_guess, const double T00,
                              const double M, const double J0,
                              double &v_solution);

    int solve_u0_Newton(const double u0_guess, const double T00,
                        const double K00, const double M, const double J0,
                        double &u0_solution);

    void regulate_grid(ReconstCell &grid_cell, double elocal) const;
};

#endif  // SRC_RECONST_H_
