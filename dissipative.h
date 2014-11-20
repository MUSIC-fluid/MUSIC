#ifndef DISSIPATIVE_H
#define DISSIPATIVE_H

#include "util.h"
#include "grid.h"
#include "data.h"
#include "minmod.h"
#include <iostream>

class Diss{
 private:
  EOS *eos;
  Minmod *minmod;

 public:
  Diss(EOS *eosIn);
  ~Diss();
  
  double MakeWSource(double tau, int alpha, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		     Grid *Lneighbor2, Grid *Rneighbor2, InitData *DATA, int rk_flag, int size, int rank);
  
  int Make_uWRHS
    (double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
     Grid *Lneighbor2, Grid *Rneighbor2, double **w_rhs, InitData *DATA, int rk_flag, int size, int rank);
  
  void Get_uWmns(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			 Grid *Lneighbor2, Grid *Rneighbor2, 
		 int mu, int nu, int direc,
		 double *g, double *f, double *gp1, double *fp1,
		 double *gp2, double *fp2, double *gm1, double *fm1, 
		 double *gm2, double *fm2, InitData *DATA, int rk_flag, int size, int rank);
  
  double Make_uWSource
    (double tau, Grid *grid_pt, int mu, int nu, InitData *DATA, int rk_flag); 
  
  int Make_uPRHS
    (double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
     Grid *Lneighbor2, Grid *Rneighbor2, 
     double *p_rhs, InitData *DATA, int rk_flag, int size, int rank);
  void Get_uPis(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			 Grid *Lneighbor2, Grid *Rneighbor2, int direc,
		double *g, double *f, double *gp1, double *fp1,
		double *gp2, double *fp2, double *gm1, double *fm1, 
		double *gm2, double *fm2, InitData *DATA, int rk_flag, int size, int rank); 
  
  double Make_uPiSource
    (double tau, Grid *grid_pt, InitData *DATA, int rk_flag);

  // Sangyong Nov 18 2014
  int Make_uqRHS
    (double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
     Grid *Lneighbor2, Grid *Rneighbor2, double **w_rhs, InitData *DATA, int rk_flag, int size, int rank);
  
  double Make_uqSource
    (double tau, Grid *grid_pt, int nu, InitData *DATA, int rk_flag); 
};
#endif
