// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_DATA_H_
#define SRC_DATA_H_

#define SMALL (1.0e-16)

#include <array>
#include <stdio.h>
#include <stdlib.h>
#include <string>

//! This is a data structure contains all the parameters for simulation
typedef struct init_data {

    std::array<std::array<double, 4>, 4> gmunu = 
      {{{-1,0,0,0},
        { 0,1,0,0},
        { 0,0,1,0},
        { 0,0,0,1}}};

    int echo_level;
    int mode;               //!< 1: do everything;
    //!< 2: do hydro evolution only;
    //!< 3: do calculation of thermal spectra only;
    //!< 4: do resonance decays only
    std::string initName;
    std::string initName_rhob;
    std::string initName_ux;
    std::string initName_uy;
    std::string initName_TA;
    std::string initName_TB;
    std::string initName_rhob_TA;
    std::string initName_rhob_TB;
    std::string initName_AMPT;

    //! random seed
    int seed;
    double ecm;
    double beam_rapidity;

    int initial_eta_profile;
    int initial_eta_rhob_profile;
    // envelope function parameter for energy density in eta_s direction
    double eta_fall_off;
    double eta_flat;
    // envelope function parameter for rhoB in eta_s direction
    double eta_rhob_0;               //!< peak position
    double eta_rhob_width;           //!< Gaussian width for profile == 1
    double eta_rhob_plateau_height;  //!< central plateau height profile == 2
    double eta_rhob_width_1;         //!< outside tail Gaussian width profile == 2
    double eta_rhob_width_2;         //!< inside Gaussian width profile == 2

    int Initial_profile;    //! type of initial condition
    int initializeEntropy;  //! flag to initial entropy or energy density

    int string_dump_mode;
    double string_quench_factor;
    double parton_quench_factor;

    int nx;
    int ny;
    int neta;
    int nt;

    double x_size;      //!< in fermi -x_size/2 < x < x_size/2
    double y_size;      //!< in fermi, -y_size/2 < y < y_size/2
    double eta_size;    //!< -eta_size/2 < eta < eta_size/2
    double tau_size;    //!< tau_0 < tau < tau0+tau_size
    double tau0;

    double delta_x;
    double delta_y;
    double delta_eta;
    double delta_tau;

    int rk_order;
    double minmod_theta;

    double sFactor;     //!< overall normalization on energy density profile
    int whichEOS;       //!< type of EoS
    //! flag for boost invariant simulations
    bool boost_invariant;

    //! flag to output initial density profile
    int output_initial_density_profiles;

    //! skip step for freeze out surface finder
    int facTau;
    //! skip step for freeze out surface finder
    int fac_x;
    //! skip step for freeze out surface finder
    int fac_y;
    //! skip step for freeze out surface finder
    int fac_eta;

    int alpha_max;          //!< dimension of TJb
    int turn_on_rhob;       //!< flag to include baryon current

    int viscosity_flag;     //!< flag to include viscosity in the simulation
    int turn_on_shear;      //!< flag to include shear viscosity
    int turn_on_bulk;       //!< flag to include bulk viscosity
    int turn_on_diff;       //!< flag to include net baryon diffusion
    double shear_to_s;      //!< value of specific shear viscosity

    //! flag to include temperature dependent eta/s(T)
    int T_dependent_shear_to_s;

    //! flag to control the temperature dependence of eta/s(T) if "T_dependent_shear_to_s==2"
    double eta_over_s_min;
    double eta_over_s_slope;
    double eta_over_s_curv;

    //! flag to control the temperature dependence of eta/s(T) if "T_dependent_shear_to_s==3"
    double eta_over_s_T_kink_in_GeV;
    double eta_over_s_low_T_slope_in_GeV;
    double eta_over_s_high_T_slope_in_GeV;
    double eta_over_s_at_kink;

    //! flag to include temperature dependent zeta/s(T)
    int T_dependent_bulk_to_s;

    //! flag to control the temperature dependence of zeta/s(T) if "T_dependent_bulk_to_s==2"
    double bulk_viscosity_normalisation;
    double bulk_viscosity_width_in_GeV;
    double bulk_viscosity_peak_in_GeV;

    //! flag to control the temperature dependence of zeta/s(T) if "T_dependent_bulk_to_s==3"
    double zeta_over_s_max;
    double zeta_over_s_width_in_GeV;
    double zeta_over_s_T_peak_in_GeV;
    double zeta_over_s_lambda_asymm;

    //! multiplicative factors for the relaxation times
    double shear_relax_time_factor;
    double bulk_relax_time_factor;

    //! flag to include second order non-linear coupling terms
    int include_second_order_terms;

    //! coefficient related to the net baryon diff.
    double kappa_coefficient;

    //! decide whether to output the evolution data (1) or not (0)
    int outputEvolutionData;

    //! flag to store hydro evolution in memory for jetscape
    int store_hydro_info_in_memory;

    //! decide whether to output files for movie
    int output_movie_flag;

    //! decide whether to output files for R_pi and R_Pi
    int output_outofequilibriumsize;

    //! decide whether to output "evolution_xyeta.dat" and
    //! "evolution_Wmunu_over_shear_xyeta.dat" in binary format (1)
    //! or in text format (0)
    int outputBinaryEvolution;

    int output_hydro_debug_info;
    int output_evolution_every_N_timesteps;
    int output_evolution_every_N_x;
    int output_evolution_every_N_y;
    int output_evolution_every_N_eta;
    bool output_hydro_params_header;
    double output_evolution_T_cut;

    int doFreezeOut;            //!< flag to output freeze-out surface

    //! flag to include low temperature cell at the initial time
    int doFreezeOut_lowtemp;

    int freezeOutMethod;        //!< freeze-out method

    double TFO;      //!< freeze-out temperature. Used if useEpsFO=0

    int useEpsFO;    //!< if 1, use energy density value
                     //!<       to define freeze out condition
                     //!< if 0 use temperature in TFO

    double epsilonFreeze;       //!< freeze-out energy density in GeV/fm^3
    int N_freeze_out;
    double eps_freeze_min;
    double eps_freeze_max;
    int freeze_eps_flag;
    std::string freeze_list_filename;
    bool freeze_surface_in_binary;

    // for calculation of spectra
    int pseudofreeze;    //! flag to compute spectra in pseudorapdity
    double max_pseudorapidity;
    int pseudo_steps;
    int phi_steps;
    double min_pt;
    double max_pt;
    int pt_steps;

    //! # of resonances to include. maximum=319 (all up to 2 GeV)
    int NumberOfParticlesToInclude;

    //! number of particle for which the spectrum is to be computed.
    //! 0: all particles
    int particleSpectrumNumber;

    int include_deltaf;        //!< flag to include shear delta f
    int include_deltaf_qmu;    //!< flag to include diffusion delta f
    int include_deltaf_bulk;   //!< flag to include bulk delta f
    int deltaf_14moments;      //!< use delta f from 14 moment approxmation

    // parameters for mode 13 and mode 14
    //! rapidity range for dN/dy as a function of y
    double dNdy_y_min;
    //! rapidity range for dN/dy as a function of y
    double dNdy_y_max;
    //! pseudo-rapidity range for dN/dy as a function of y
    double dNdy_eta_min;
    //! pseudo-rapidity range for dN/dy as a function of y
    double dNdy_eta_max;
    int dNdy_nrap;

    // the integrated rapidity range for dN/dypTdpT
    double dNdyptdpt_y_min;
    double dNdyptdpt_y_max;
    double dNdyptdpt_eta_min;
    double dNdyptdpt_eta_max;

} InitData;

#endif  // SRC_DATA_H_
