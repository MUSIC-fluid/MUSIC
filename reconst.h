#ifndef RECONST_H
#define RECONST_H

#include <iostream>
#include "util.h"
#include "data.h"
#include "grid.h"
#include "eos.h"

class Reconst{
 private:
  EOS *eos;
  Grid *grid;
  Util *util;
 public:
  Reconst(EOS *eos, Grid *grid);//constructor
  ~Reconst();//destructor
  
  void ReconstError
    (char *str, int i, int rk_flag, double *qi, double **qi2, Grid *grid_pt);
  
  double GuessEps(double T00, double K00, double cs2);
  
  int ReconstIt
    (Grid *grid_p, int i, double tau, double **uq, Grid *grid_pt,
     double eps_init, double rhob_init, InitData *DATA, int rk_flag);
  double phi(double r);
};
#endif
