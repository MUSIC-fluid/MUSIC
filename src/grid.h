#ifndef GRID_H
#define GRID_H
#include "data.h"
#include "eos.h"
#include <mpi.h>
#include <iostream>
#include <iomanip>

class Grid{

 public:

  double prev_T00;
  double prev_T33;
  double epsilon;
  double p;
  double rhob;
  double ***TJb; 
  /* stress energy tensor plus baryon current  */
  /* TJb[flag][alpha][mu] */
  /* flag = 0 is the actual values. flag != 0 are the intermediate values
     for the Runge-Kutta step */
  /* alpha = 4 means Jb */
  // double correction;
  
  /* temporary values for the final RK update */
  double epsilon_t;
  double p_t;
  double rhob_t;

  /* u[flag][mu]: flag=0 is the actual values. flag != are for RK steps */
  double **u;
  
  /* to include shear viscosity */
  /* we need to calculate partial_tau u[mu] */
  double **prev_u;  /* u[mu] from the previous time step including the rk flag */
//   double **pprev_u; /* u[mu] from 2 time step ago including the rk flag */
//  double **dU; /* dU[m][n] = partial_m u_n at the current time */
//  double ***pimunu; /* Stress part of the TJb */
    
  Grid **nbr_p_1; 
  Grid **nbr_m_1; 
  Grid **nbr_p_2; 
  Grid **nbr_m_2; 
  
  Grid *grid_v_malloc(int );
  Grid **grid_m_malloc(int , int );
  Grid ***grid_c_malloc(int , int , int );

  /* to include shear viscosity */
  
  /* This is u^mu partial_mu u^nu */
  double **a;
  
  /* this is the expansion rate partial_mu u^mu */
  double *theta_u;

  // the velocity shear tensor
  double ***sigma;
  
  /* we need to calculate partial_tau u[mu] */
  double ***dUsup; 
  /* dU[flag][m][n] = u^{m,n} = partial^n u^m with the rk flag */
  /* note that they are superscripted. So partial^t = -partial_t */
  double ***Wmunu; /* shear part of the TJb with the rk_flag */
  double ***prevWmunu; 
//   double ***pprevWmunu; 
  
  double ***Pimunu; /* bulk part of the TJb with the rk_flag */
  double ***prevPimunu; 
//   double ***pprevPimunu; 
  
  double *pi_b; /* bulk pressure */
//   double *prev_pi_b; /* bulk pressure */
//   double *pprev_pi_b; /* bulk pressure */

  double **W_prev; // the one for the freeze-out surface finder for interpolation

  int revert_flag;
  int trouble;
  double T;  //added
  double mu; //added
  //  double tauF; //added
  
  int position[4];

  Grid(){};//constructor
  ~Grid(){};//destructor

  double epsilon_prev; //added (previous 10 epsilons)
  double rhob_prev;
  double pi_b_prev;
  double u_prev[4]; //added (previous 10 epsilons)

  void ComputeV2(InitData *DATA, Grid ***arena, double tau); //added
  void ComputeEccentricity(InitData *DATA, Grid ***arena, double tau); //added
  void print_rhob_evolution(InitData *DATA, Grid ***arena, double tau, EOS* eos);
  void print_qmu_evolution(InitData *DATA, Grid ***arena, double tau, EOS* eos);
  void ComputeAnisotropy(InitData *DATA, Grid ***arena, double tau); //added
  void PrintGrid(Grid *grid_p, int rk_order);
  //void InitArena(InitData *DATA, GridElem ****arena);
  void LinkNeighbors(InitData *DATA, Grid ****arena);
  void InitTJb(InitData *DATA, Grid ****arena);
  void PrintAxy(InitData *DATA, Grid ***arena, double tau);
  void PrintAxy2(InitData *DATA, Grid ***arena, double tau);
  void PrintdEdEta(InitData *DATA, Grid ***arena);
  void OutputEvolutionDataXYZ(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
  void OutputEvolutionDataXYEta(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
  void OutputPlotDataXYZ(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
  void OutputEvolutionOSCAR(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
  void OutputXY(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
  void PrintArena(Grid ***arena, InitData *DATA, double tau);
  void PrintEtaEpsilon(Grid ***arena, InitData *DATA, double tau, int size, int rank);
  void PrintxEpsilon(Grid ***arena, InitData *DATA, double tau, int size, int rank);
  void ComputeEnergyConservation(InitData *DATA, Grid ***arena, double tau);
  void getAverageTandPlasmaEvolution(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
  void Output_hydro_information_header(InitData *DATA, EOS *eos);
  void Tmax_profile(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank);
};
#endif
  
