#ifndef IDEAL_H
#define IDEAL_H

#include "util.h"
#include "data.h"
#include "grid.h"


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


void MakeDeltaQI
(double tau, Grid *grid_pt, double *qi, double *rhs, 
 InitData *DATA, int rk_flag);

void MakeMaxSpeedAs
(double tau, BdryCells *HalfwayCells, 
 double aiph[], double aimh[], int rk_flag);

double MaxSpeed (double tau, int direc, Grid *grid_p, int rk_flag);

void InitTempGrids(BdryCells *HalfwayCells, int rk_order);

void InitNbrQs(NbrQs *NbrCells);

int MakeQIHalfs
(double *qi, NbrQs *NbrCells, BdryCells *HalfwayCells, 
 Grid *grid_pt, InitData *DATA);

void GetQIs
(double tau, Grid *grid_pt, double *qi, NbrQs *NbrCells, 
 int rk_flag, InitData *DATA);

int ConstHalfwayCells
 (double tau, BdryCells *HalfwayCells, double *qi, Grid *grid_pt,
  InitData *DATA, int rk_flag);


void MakeKTCurrents
 (double tau, double **DFmmp, Grid *grid_pt, 
  BdryCells *HalfwayCells, int rk_flag);

void AdvanceQI
(double tau, double **qirk, double **Fiph, double **Fimh, double *qi,
 Grid *grid_pt, InitData *DATA, int rk_flag);


#endif
