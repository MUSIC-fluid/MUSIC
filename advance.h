
#ifndef ADVANCE_H
#define ADVANCE_H

#include "data.h"
#include "grid.h"
#include "reconst.h"
#include "dissipative.h"
#include "minmod.h"
#include "u_derivative.h"
#include <iostream>
#include "random.h"

class Advance{
 private:
  
  Util *util;
  Reconst *reconst; // declare Reconst object
  Diss *diss; // dissipative object
  EOS *eos;
  Grid *grid;
  Minmod *minmod;
  U_derivative *u_derivative;
  
  typedef struct bdry_cells
  {
    Grid *grid_p_h_L;
    Grid *grid_p_h_R;
    Grid *grid_m_h_L;
    Grid *grid_m_h_R;
    
    double **qiphL;
    double **qiphR;
    double **qimhL;
    double **qimhR;
    
  } BdryCells;
  
  
  typedef struct nbrs
  {
    double **qip1;
    double **qip2;
    double **qim1;
    double **qim2;
  } NbrQs;

 public:
  Advance(EOS *eosIn, Grid *grid);
  ~Advance();

  int AdvanceIt(double tau_init, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int rk_flag, int size, int rank);
  //int AdvanceIt(double tau, InitData *DATA, Grid ***arena, int rk_flag, int size, int rank);


  // advance routines separate for 
  //T^{0 nu} \del T^{i\nu} (T)
  //W
  //T^{0 nu} with W source (TS)
  //W with source (WS)

  void MPISendReceiveW(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag);
  void MPISendReceiveT(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag);
  void MPISendReceiveXi(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag);
  void MPISendReceiveXi2(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag);

  int AdvanceLocalT(double tau_init, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, 
		    int rk_flag, int size, int rank);
  
  
  int AdvanceLocalW(double tau_init, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, 
		    int rk_flag, int size, int rank);
  
  int AdvanceLocalDeltaTAndW(double tau_init, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, 
		    int rk_flag, int size, int rank);
  
  int AdvanceLocalTS
    (double tau_init, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank);
  
  int AdvanceLocalWS
    (double tau_init, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank);
  
  int FirstRKStepTS
    (double tau_it, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, 
     double *qi, double *rhs, double **w_rhs, double **qirk, Grid *grid_rk, int size, int rank);
  
  int FirstRKStepWS
    (double tau_it, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, 
     double *qi, double *rhs, double **w_rhs, double **qirk, Grid *grid_rk, int size, int rank);
  
  int FirstRKStepT
    (double tau_it, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, 
     double *qi, double *rhs, double **w_rhs, double **qirk, Grid *grid_rk, int size, int rank);
  
  int FirstRKStepW
    (double tau_it, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, 
     double *qi, double *rhs, double **w_rhs, double **qirk, Grid *grid_rk, int size, int rank);
  
  int FirstRKStepDeltaTAndW(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank);

  void UpdateTJbRK
    (Grid *grid_rk, Grid *grid_pt, InitData *DATA, int rk_flag);
  
  int QuestRevert(double tau, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank);
  
  void RememberInits(Grid *grid_pt);
  void TestW(double tau, InitData *DATA, Grid *grid_pt, int rk_flag);
  void ProjectSpin2W(double tau, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank);
  void ProjectSpin2WS(double tau, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank);
  void MakeDeltaQI(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		   Grid *Lneighbor2, Grid *Rneighbor2, double *qi, double *rhs, 
		   InitData *DATA, int rk_flag, int size, int rank);
  void GetQIs(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2,
	      double *qi, NbrQs *NbrCells, int rk_flag, InitData *DATA, int size, int rank);
  int MakeQIHalfs(double *qi, NbrQs *NbrCells, BdryCells *HalfwayCells, 
		  Grid *grid_pt, InitData *DATA);
  int ConstHalfwayCells(double tau, BdryCells *HalfwayCells, double *qi, Grid *grid_pt,
			InitData *DATA, int rk_flag, int size, int rank);
  void MakeKTCurrents(double tau, double **DFmmp, Grid *grid_pt, 
		      BdryCells *HalfwayCells, int rk_flag);
  void MakeMaxSpeedAs(double tau, BdryCells *HalfwayCells, double aiph[], double aimh[], int rk_flag);
  double MaxSpeed (double tau, int direc, Grid *grid_p, int rk_flag);
  void InitNbrQs(NbrQs *NbrCells);
  void InitTempGrids(BdryCells *HalfwayCells, int rk_order);
  void AdvanceXi(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
  double getDXi(double tau, int alpha, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor,
		InitData *DATA, int rk_flag, int size, int rank);
  //These routines have been added for linearized thermal fluctuations in hydrodynamics. -CFY 11/4/2013

  double getDeltaUDWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  double getDMuDeltaWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int nu);
  double getDMuXiMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int nu);
  double getDAlphaWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, int alpha);
  double getDAlphaDeltaWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, int alpha);
  double getDAlphaXiMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, int alpha);

  double getUiDiDeltaWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  double getUiDiXiMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  double getDeltaMuNu(Grid *grid_pt, int rk_flag, int mu, int nu);
  double getDeltaDeltaMuNu(Grid *grid_pt, int rk_flag, int mu, int nu);
  double getDMuDeltaUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);

  double getDeltaW0Delta(Grid *grid_pt, int rk_flag, int mu, int nu);
  double getDeltaDeltaW0Delta(Grid *grid_pt, int rk_flag, int mu, int nu);
  double getDeltaDeltaWDelta(Grid *grid_pt, int rk_flag, int mu, int nu);
  double getDeltaW0DeltaDelta(Grid *grid_pt, int rk_flag, int mu, int nu);

  double getDMuUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  double getDMuUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);

  double getDMuDeltaUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  double getDMuDeltaUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);

  double getDeltaMuUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  double getDeltaMuDeltaUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);

  ////double getU0MuU0DDeltaUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  //
  //double getS(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  //double getDeltaS(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);

  //These definitions use previously calculated values of divergences to operate more efficiently:
  double getS(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, double dMuUMu);
  double getDeltaS(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, double dMuUMu, double dMuDeltaUMu);

  //double getDeltaS0Delta(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  //double getDeltaDeltaSDelta(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  //double getDeltaDeltaS0Delta(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  //double getDeltaS0DeltaDelta(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);

  //Also to operate more efficiently, these subroutines combine several terms for fewer matrix multiplications:
  double getDeltaDeltaWDeltaSDelta(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag,
				   int size, int rank, double eta, double deltaEta, double dMuUMu, double dMuDeltaUMu,
				   int ii, int ij);
  double getDeltaXiDelta(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag,
				   int size, int rank, double eta, double deltaEta, double dMuUMu, double dMuDeltaUMu,
				   int ii, int ij);
  //double getDeltaDeltaW0S0Delta(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag,
  //				int size, int rank, double eta, double dMuUMu, int mu, int nu);
  //double getDeltaW0S0DeltaDelta(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag,
  //				int size, int rank, double eta, double dMuUMu, int mu, int nu);

  void reconstructDeltaPAndDeltaU(double tau, InitData *DATA, Grid *grid_pt, int rk_flag);
  double getDUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
  double getDeltaUDUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
  double getDDeltaUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
  double getDWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
  double getDUWMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
  double getDeltaUDUWMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
  double getDDeltaUWMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
  double getDeltaWDUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
  double getXiDUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
  double getDeltaUWDU(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
  double getDUDeltaWMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
  double getProjectorSource(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, double *DUW, double *deltaUDUW, double *DDeltaUW, double *DUDeltaW);
  double getProjectorSourceXi(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, double *DUXi);
};
#endif

//double getDotProduct(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
//double getU0Dot(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
//double getU0jDjDeltaUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
//double getDMuU0Mu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
//double getDiU0i(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
//double getDiDeltaUi(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
//double getU0NuDNuU0Mu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
//double getU0iDiDU0(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
//double getU0iDiDeltaP(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
//double getU0jDjU0Mu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
//double getDeltaUNuDNuU0Mu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
//double getUjDjU0Mu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
//double getDeltaUjDjU0Mu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
//double getDeltaUMuDMuW0(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
//double getU0MuDMuW0(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank);
//double getDiDeltaP(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu);
//double getDeltaT(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
//double getDeltaUMuU0DDeltaUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
//double getU0MuDeltaUDDeltaUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
//void getInverseA(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, double **A);
//void getC(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, double *C);
//double getU0MuU0DU0Nu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
//double getDeltaDeltaMuU0Nu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
//double getDeltaUMuU0DU0Nu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
//double getU0MuDeltaUDU0Nu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
//double getU0MuDeltaUDUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu);
