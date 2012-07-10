#ifndef HYDRO_DATA_H
#define HYDRO_DATA_H

#define SMALL (1.0e-8)

typedef struct init_data
{
  double **gmunu; /* metric */
  
  int nx;
  int ny;
  int neta;
  int nt;
  
  double x_size; /* in fermi -x_size/2 < x < x_size/2 */
  double y_size; /* in fermi, ditto */
  double eta_size; /* ditto */
  double tau_size; /* tau_0 < tau < tau0+tau_size  */
  double tau0;
  double epsilonFreeze; // freeze-out energy density in GeV
  
  double delta_x;
  double delta_y;
  double delta_eta;
  double delta_tau;
  
  double epsilon0;
  double eta_fall_off;
  double eta_flat;
  double R_A;
  double a_A;
  double a_short;
  double a_long;
  int rk_order;
  
  double minmod_theta;
  
  double SigmaNN;
  double b;
  char* Target;
  char* Projectile;
  int LexusImax;
  
  int Initial_profile;
} InitData;


#endif
