// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_DATA_H_
#define SRC_DATA_H_

#define SMALL (1.0e-16)

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

struct ReturnValue  {
    double x;
    double y;
    int collided;
    int acceptances;
};

typedef struct init_data {
    double **gmunu; /* metric */
    string Target;
    string Projectile;
    string initName;
    string initName_rhob;
    string initName_ux;
    string initName_uy;
    string initName_TA;
    string initName_TB;

    int nx;
    int ny;
    int neta;
    int nt;

    int input_grid_nx;
    int input_grid_ny;
    int input_grid_neta;
    double input_grid_dx;
    double input_grid_dy;
    double input_grid_deta;

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
    double beam_rapidity;
    double R_A;
    double a_A;
    double a_short;
    double a_long;
    int rk_order;
    int reconst_type;
    int taufac;

    double minmod_theta;

    double SigmaNN;
    double b;

    // envelope function parameter for rhoB in eta_s direction
    double eta_rhob_0;               // peak position
    double eta_rhob_width;           // Gaussian width for profile == 1
    double eta_rhob_plateau_height;  // central plateau height profile == 2
    double eta_rhob_width_1;         // outside tail Gaussian width profile == 2
    double eta_rhob_width_2;         // inside Gaussian width profile == 2

    int LexusImax;

    // skip step for freeze out surface finder
    int facTau;
    int fac_x;
    int fac_y;

    double deltaY;
    double ymax;

    // for calculation of spectra
    double max_pseudorapidity;
    int pseudo_steps;
    int phi_steps;
    double min_pt;
    double max_pt;
    int pt_steps;

    int pseudofreeze;

    // # of resonances to include. maximum=319 (all up to 2 GeV)
    int NumberOfParticlesToInclude;
    int whichEOS;
    int Initial_profile;
    int initializeEntropy;
    double epsilonFreeze;       // freeze-out energy density in GeV/fm^3

    int N_freeze_out;
    double eps_freeze_min;
    double eps_freeze_max;
    int freeze_eps_flag;
    string freeze_list_filename;

    // number of particle for which the spectrum is to be computed.
    // 0: all particles
    int particleSpectrumNumber;
    int freezeOutMethod;    // 1: Hirano's simple method, 2: 4d-triangulation
    int mode;               // 1: do everything;
                            // 2: do hydro evolution only;
                            // 3: do calculation of thermal spectra only;
                            // 4: do resonance decays only
    // fraction of binary collisons scaling in initial distribution
    double hard;
    // if 1: rotate the initial configuration in the x-y plane by 45 degrees
    // (to test dependence of v_4 on lattice)
    int rotateBy45degrees;

    // decide whether to output the evolution data (1) or not (0)
    int outputEvolutionData;
    // decide whether to output "evolution_xyeta.dat" and
    // "evolution_Wmunu_over_shear_xyeta.dat" in binary format (1)
    // or in text format (0)
    int outputBinaryEvolution;
    int includeJet;
    int includeTrigger;
    int include_deltaf;
    int include_deltaf_qmu;
    int include_deltaf_bulk;
    int deltaf_14moments;

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
    int turn_on_diff;

    int zero_or_stop;

    double eps_limit;

    double local_y_max;

    // eccentricities and axes with respect to which we have to compute
    // v_2 or v_3
    double ecc2;
    double ecc3;
    double ecc3r3;
    double Psi2;
    double Psi3;
    double Psi3r3;

    // minimal and maximal impact parameter when using Initial_Distribution 3,
    // for sampling
    double bmin;
    double bmax;
    int doFreezeOut;
    int seed;

    int Nbin_to_file;

    double sigma0;
    double TFO;      // freeze-out temperature. Used if useEpsFO=0
    // if 1, use energy density value to define freeze out condition
    // if 0 use temperature in TFO
    int useEpsFO;

    int size;
    int rank;
    double avgT;
    int nSteps;
    double avgT2;
    int nSteps2;
    double plasmaEvolutionTime;
    double plasmaEvolutionTime2;
    double sFactor;
    // For use with freeze_out_method=3.
    // Maximum size of freeze out hypersurface in eta,
    // implemented in freeze out finder.
    double max_delta_eta;
    // For use with freeze_out_method=3.
    // Maximum size of freeze out hypersurface in eta,
    // implemented in Cooper-Frye.
    double max_delta_eta2;

    int boost_invariant;    // set to 1 for rapidity-indendent solution.
    int boostInvariant;     // alias for boost_invariant
    int check_FO3_at_boundary_eta, check_FO3_at_boundary_xy;
    bool output_hydro_debug_info;
    int output_evolution_every_N_timesteps;
    int output_evolution_every_N_x;
    int output_evolution_every_N_y;
    int output_evolution_every_N_eta;
    bool output_hydro_params_header;
    int initial_eta_profile;
    int initial_eta_rhob_profile;

    // QuestRevert
    double QuestRevert_rho_shear_max, QuestRevert_rho_bulk_max;
    double QuestRevert_factor, QuestRevert_epsilon_min;
    double QuestRevert_prefactor, QuestRevert_eps_factor;
    double QuestRevert_rho_q_max;

    // parameters for mode 14
    double dNdy_y_min;     // rapidity range for dN/dy as a function of y
    double dNdy_y_max;
    double dNdy_eta_min;   // pseudo-rapidity range for dN/dy as a function of y
    double dNdy_eta_max;
    int dNdy_nrap;
    // the integrated rapidity range for dN/dypTdpT
    double dNdyptdpt_y_min;
    double dNdyptdpt_y_max;
    double dNdyptdpt_eta_min;
    double dNdyptdpt_eta_max;

    int echo_level;
} InitData;

#endif  // SRC_DATA_H_
