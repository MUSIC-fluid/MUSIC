// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_RECONST_H_
#define SRC_RECONST_H_

#include <iostream>
#include "util.h"
#include "data.h"
#include "grid.h"
#include "eos.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

struct reconst_v_params {
    double T00;
    double K00; 
    double J0;
};

class Reconst;

struct CCallbackHolder {  
    Reconst* cls;
    void* params;
};

class Reconst {
 private:
    int reconst_type;
    EOS *eos;
    Util *util;

    int max_iter;
    double rel_err, abs_err;
    // initialize gsl root finding solver
    int gsl_rootfinding_max_iter;
    double gsl_rootfinding_abserr, gsl_rootfinding_relerr;
    const gsl_root_fsolver_type *gsl_solverType;
    gsl_root_fsolver *gsl_rootfinding_solver;
    gsl_function gslFunc;

 public:
    Reconst(EOS *eos, int reconst_type_in);
    ~Reconst();
      
    void ReconstError(const char *str, int i, int rk_flag, double *qi,
                      double **qi2, Grid *grid_pt);

    int ReconstIt_shell(Grid *grid_p, int direc, double tau, double **uq,
                        Grid *grid_pt, InitData *DATA, int rk_flag);

    // reconst_type == 0
    int ReconstIt(Grid *grid_p, int i, double tau, double **uq, Grid *grid_pt,
                  InitData *DATA, int rk_flag);
    double GuessEps(double T00, double K00, double cs2);

    // reconst_type == 1
    int ReconstIt_velocity_iteration(Grid *grid_p, int direc, double tau,
                                     double **uq, Grid *grid_pt,
                                     InitData *DATA, int rk_flag);
    double reconst_velocity_f(double v, double T00, double M, double J0);
    double reconst_u0_f(double u0, double T00, double K00, double M,
                        double J0);

    // reconst_type == 2
    int ReconstIt_velocity_Newton(Grid *grid_p, int direc, double tau,
                                  double **uq, Grid *grid_pt,
                                  InitData *DATA, int rk_flag);
    double reconst_velocity_f_Newton(double v, double T00, double M,
                                     double J0);
    double reconst_u0_f_Newton(double u0, double T00, double K00,
                               double M, double J0);
    double reconst_velocity_df(double v, double T00, double M, double J0);
    double reconst_u0_df(double u0, double T00, double K00, double M,
                         double J0);

    // reconst_type == 99
    // (does not support openmp because of the static functions)
    int ReconstIt_velocity_gsl(Grid *grid_p, int direc, double tau,
                               double **uq, Grid *grid_pt,
                               InitData *DATA, int rk_flag);
    double reconst_velocity_function(double v, void *params);
    static double CCallback_reconst_v(double x, void* params) {
        CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
        return h->cls->reconst_velocity_function(x, h->params);
    }
    double reconst_u0_function(double u0, void *params);
    static double CCallback_reconst_u0(double x, void* params) {
        CCallbackHolder* h = static_cast<CCallbackHolder*>(params);
        return h->cls->reconst_u0_function(x, h->params);
    }
};

#endif  // SRC_RECONST_H_
