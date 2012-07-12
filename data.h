#ifndef HYDRO_DATA_H
#define HYDRO_DATA_H

#define SMALL (1.0e-8)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

struct ReturnValue 
{
  double x;
  double y;
  int collided;
  int acceptances;
};

typedef struct init_data
{
  double **gmunu; /* metric */
  string Target;
  string Projectile;
  string initName;


  int nx;
  int ny;
  int neta;
  int nt;
  
  double x_size; /* in fermi -x_size/2 < x < x_size/2 */
  double y_size; /* in fermi, ditto */
  double eta_size; /* ditto */
  double tau_size; /* tau_0 < tau < tau0+tau_size  */
  double tau0;
  
  double delta_x;
  double delta_y;
  double delta_eta;
  double delta_tau;
  
  double epsilon0;
  double rhoB0;
  double eta_fall_off;
  double eta_flat;
  double R_A;
  double a_A;
  double a_short;
  double a_long;
  int rk_order;
  int taufac;

  double minmod_theta;
  
  double SigmaNN;
  double b;
 




  int LexusImax;
  int facTau; // maximally 10

  double deltaY;
  double ymax;
  
  int NumberOfParticlesToInclude; // # of resonances to include. maximum=319 (all up to 2 GeV) 
  int whichEOS;         // 1: EOS-Q as in AZHYDRO, 2: Huovinen/Petrczky lattice EOS, 3:Huovinen/Petrczky lattice EOS with PCE
  int Initial_profile;  // 1: for a realistic profile in AA collision
  int initializeEntropy;// 1: make the entropy density proportional to 3 of BC and WN instead of the energy density, 0: use energy denisty
  double epsilonFreeze; // freeze-out energy density in GeV/fm^3
  int particleSpectrumNumber; // number of particle for which the spectrum is to be computed. 0: all particles
  int freezeOutMethod;  // 1: Hirano's simple method, 2: 4d-triangulation
  int mode;             // 1: do everything; 
                        // 2: do hydro evolution only; 
                        // 3: do calculation of thermal spectra only;
                        // 4: do resonance decays only
  double hard;          // fraction of binary collisons scaling in initial distribution
  int rotateBy45degrees;// if 1: rotate the initial configuration in the x-y plane by 45 degrees (to test dependence of v_4 on lattice)
  int outputEvolutionData; // decide whether to output the evolution data (1) or not (0)
  int includeJet;
  int includeTrigger;
  int include_deltaf;

  double tau_pi;
  double tau_b_pi;
  double shear_to_s;
  double bulk_to_s;
  int viscosity_flag;
  int verbose_flag;
  int turn_on_rhob;
  int alpha_max;
  int T_dependent_shear_to_s;
  
  int turn_on_shear;
  int turn_on_bulk;
  
  int zero_or_stop;
  
  double eps_limit;
  
  double local_y_max;

  // eccentricities and axes with respect to which we have to compute v_2 or v_3
  double ecc2;
  double ecc3;
  double ecc3r3;
  double Psi2;
  double Psi3;
  double Psi3r3;

  double bmin; // minimal and maximal impact parameter when using Initial_Distribution 3, for sampling
  double bmax;
  int doFreezeOut;  
  int seed;
  double sigma0;
  double TFO; // freeze-out temperature. Used if useEpsFO=0
  int useEpsFO; // if 1, use energy density value to define freeze out condition, if 0 use temperature in TFO
  int size;
  int rank;
  double avgT;
  int nSteps;
  double avgT2;
  int nSteps2;
  double plasmaEvolutionTime;
  double plasmaEvolutionTime2;
  double sFactor;
} InitData;


#endif