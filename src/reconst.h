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

    ReconstCell ReconstIt_shell(double tau, TJbVec &q_vec,
                                const Cell_small &grid_pt);

    void revert_grid(ReconstCell &grid_current,
                     const Cell_small &grid_prev) const;

    int ReconstIt_velocity_Newton(ReconstCell &grid_p, double tau,
                                  const TJbVec &q, const Cell_small &grid_pt);

#pragma omp declare simd
    void reconst_velocity_fdf(const double v, const double T00, const double M,
                              const double J0, double &fv, double &dfdv) const {
        // this function returns f(v) = M/(M0 + P)
        const double epsilon = T00 - v*M;
        const double temp    = sqrt(1. - v*v);
        const double rho     = J0*temp;

        const double pressure = eos.get_pressure(epsilon, rho);
        const double temp1    = T00 + pressure;
        const double temp2    = v/temp;
        const double dPde     = eos.p_e_func(epsilon, rho);
        const double dPdrho   = eos.p_rho_func(epsilon, rho);

        fv   = v - M/temp1;
        dfdv = 1. - M/(temp1*temp1)*(M*dPde + J0*temp2*dPdrho);
    }

    void regulate_grid(ReconstCell &grid_cell, double elocal) const;
};

#endif  // SRC_RECONST_H_
