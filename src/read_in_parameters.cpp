
#include <iostream>
#include <cstring>
#include "read_in_parameters.h"
#include "util.h"

using namespace std;
using Util::getParameter;

namespace ReadInParameters {

pretty_ostream music_message;

InitData read_in_parameters(std::string input_file) {
    InitData parameter_list;

    // this function reads in parameters
    string tempinput;

    // echo_level controls the mount of
    // warning message output during the evolution
    parameter_list.echo_level = getParameter(input_file, "echo_level", 9);

    parameter_list.beastMode = getParameter(input_file, "beastMode", 0);

    // mode:
    // 1: Does everything. Evolution. Computation of thermal spectra.
    //    Resonance decays. Observables.
    // 2: Evolution only.
    // 3: Compute all thermal spectra only.
    // 4: Resonance decays only.
    // 13: Compute observables from previously-computed thermal spectra
    // 14: Compute observables from post-decay spectra
    parameter_list.mode = getParameter(input_file, "mode", 2);

    // Initial_profile:
    parameter_list.Initial_profile = (
            getParameter(input_file, "Initial_profile", 1));

    // boost-invariant
    int temp_boost_invariant = getParameter(input_file, "boost_invariant", 1);
    if (temp_boost_invariant == 0) {
        parameter_list.boost_invariant = false;
    } else {
        parameter_list.boost_invariant = true;
    }

    parameter_list.output_initial_density_profiles = (
            getParameter(input_file, "output_initial_density_profiles", 0));

    //////////////////////////////////////////////////////////////////////////
    // parameters for 3D-Glauber initial condition
    //////////////////////////////////////////////////////////////////////////
    parameter_list.string_dump_mode = (
            getParameter(input_file, "string_dump_mode", 1));

    parameter_list.stringSourceSigmaX = (
            getParameter(input_file, "string_source_sigma_x", 0.5));

    parameter_list.stringSourceSigmaEta = (
            getParameter(input_file, "string_source_sigma_eta", 0.5));

    parameter_list.stringTransverseShiftFrac = (
            getParameter(input_file, "stringTransverseShiftFrac", 0.0));

    parameter_list.preEqFlowFactor = (
            getParameter(input_file, "stringPreEqFlowFactor", 0.0));

    parameter_list.string_quench_factor = (
            getParameter(input_file, "string_quench_factor", 0.0));

    parameter_list.parton_quench_factor = (
            getParameter(input_file, "parton_quench_factor", 1.0));

    parameter_list.gridPadding = getParameter(input_file, "gridPadding", 3.0);

    //////////////////////////////////////////////////////////////////////////
    // parameters for parameteric initial conditions
    //////////////////////////////////////////////////////////////////////////
    // Initial_Distribution_input_filename
    string tempinitName = "initial/initial_ed.dat";
    tempinput = Util::StringFind4(input_file,
                                  "Initial_Distribution_input_filename");
    if (tempinput != "empty")
        tempinitName.assign(tempinput);
    parameter_list.initName.assign(tempinitName);

    // Initial_Distribution_Filename for TA
    string tempinitName_TA = "initial/initial_TA.dat";
    tempinput = Util::StringFind4(input_file,
                                  "Initial_TA_Distribution_Filename");
    if (tempinput != "empty")
        tempinitName_TA.assign(tempinput);
    parameter_list.initName_TA.assign(tempinitName_TA);
    // Initial_Distribution_Filename for TB
    string tempinitName_TB = "initial/initial_TB.dat";
    tempinput = Util::StringFind4(input_file,
                                  "Initial_TB_Distribution_Filename");
    if (tempinput != "empty")
        tempinitName_TB.assign(tempinput);
    parameter_list.initName_TB.assign(tempinitName_TB);
    // Initial_Distribution_Filename for participant list
    string tempinitName_part = "initial/participantList.dat";
    tempinput = Util::StringFind4(input_file,
                                  "Initial_participantList_Filename");
    if (tempinput != "empty")
        tempinitName_part.assign(tempinput);
    parameter_list.initName_participants.assign(tempinitName_part);

    parameter_list.nucleonWidth = (
        getParameter(input_file, "nucleon_width", 0.5));
    // compute beam rapidity according to the collision energy
    parameter_list.ecm = getParameter(input_file, "ecm", 200.);
    double y_beam = acosh(parameter_list.ecm/(2.*Util::m_N));
    parameter_list.beam_rapidity = y_beam;

    // Eta_plateau_size:
    // width of the flat region symmetrical around eta=0
    parameter_list.eta_flat = getParameter(input_file, "Eta_plateau_size", 20.);
    // Eta_fall_off:
    // width of half-Gaussian on each side of a central pleateau in eta
    parameter_list.eta_fall_off = getParameter(input_file, "Eta_fall_off", 0.4);
    // yL_frac: the fraction of Y_CM in the longitudinal velocity
    parameter_list.yL_frac = getParameter(input_file, "yL_frac", 0.0);

    // eta envelope function parameter for rhob
    parameter_list.initial_eta_rhob_profile = (
            getParameter(input_file, "initial_eta_rhob_profile", 1));
    parameter_list.eta_rhob_0 = getParameter(input_file, "eta_rhob_0", 3.0);
    parameter_list.eta_rhob_width = (
            getParameter(input_file, "eta_rhob_width", 1.0));
    parameter_list.eta_rhob_plateau_height = (
            getParameter(input_file, "eta_rhob_plateau_height", 0.5));
    parameter_list.eta_rhob_width_1 = (
            getParameter(input_file, "eta_rhob_width_1", 1.0));
    parameter_list.eta_rhob_width_2 = (
            getParameter(input_file, "eta_rhob_width_2", 1.0));
    parameter_list.eta_rhob_asym = (
            getParameter(input_file, "eta_rhob_asym", 1.0));

    //initialize_with_entropy: 0: with energy density, 1: with entropy density
    parameter_list.initializeEntropy = (
            getParameter(input_file, "initialize_with_entropy", 0));

    //////////////////////////////////////////////////////////////////////////
    // Parameters for IP-Glasma initial conditions
    //////////////////////////////////////////////////////////////////////////
    // s_factor: for use with IP-Glasma initial conditions
    parameter_list.sFactor = getParameter(input_file, "s_factor", 1.0);

    // preEqVisFactor: for use with IP-Glasma initial conditions
    parameter_list.preEqVisFactor = (
            getParameter(input_file, "preEqVisFactor", 1.0));


    //////////////////////////////////////////////////////////////////////////
    // Parameters for AMPT initial conditions
    //////////////////////////////////////////////////////////////////////////
    // Initial_Distribution_AMPT_filename for AMPT
    string tempinitName_AMPT = "initial/initial_AMPT.dat";
    tempinput = Util::StringFind4(input_file,
                                  "Initial_Distribution_AMPT_filename");
    if (tempinput != "empty")
        tempinitName_AMPT.assign(tempinput);
    parameter_list.initName_AMPT.assign(tempinitName_AMPT);

    //////////////////////////////////////////////////////////////////////////
    // Hydro Parameters
    //////////////////////////////////////////////////////////////////////////
    // Grid_size_in_*
    // number of cells in x,y direction
    parameter_list.nx = getParameter(input_file, "Grid_size_in_x", 10);
    parameter_list.ny = getParameter(input_file, "Grid_size_in_y", 10);
    parameter_list.neta = getParameter(input_file, "Grid_size_in_eta", 1);

    // grid_size_in_fm:
    // total length of box in x,y direction in fm (minus delta_*)
    parameter_list.x_size = getParameter(input_file, "X_grid_size_in_fm", 25.);
    parameter_list.y_size = getParameter(input_file, "Y_grid_size_in_fm", 25.);
    parameter_list.eta_size = getParameter(input_file, "Eta_grid_size", 8.);

    parameter_list.delta_x =
            parameter_list.x_size/static_cast<double>(parameter_list.nx - 1);
    parameter_list.delta_y =
            parameter_list.y_size/static_cast<double>(parameter_list.ny - 1);
    if (parameter_list.neta > 1) {
        parameter_list.delta_eta = (parameter_list.eta_size
                                /static_cast<double>(parameter_list.neta - 1));
    } else {
        parameter_list.delta_eta = 0.1;
    }
    music_message << " DeltaX = " << parameter_list.delta_x << " fm";
    music_message.flush("info");
    music_message << " DeltaY = " << parameter_list.delta_y << " fm";
    music_message.flush("info");
    music_message << " DeltaETA = " << parameter_list.delta_eta;
    music_message.flush("info");

    // Initial_time_tau_0:  in fm
    parameter_list.tau0 = getParameter(input_file, "Initial_time_tau_0", 1.0);

    // Delta_Tau: time step to use in [fm].
    parameter_list.delta_tau = getParameter(input_file, "Delta_Tau", 0.02);
    parameter_list.dtaudxRatio = getParameter(input_file, "dtaudxRatio", 0.1);
    music_message << " dtaudxRatio = " << parameter_list.dtaudxRatio;
    music_message.flush("info");

    music_message << " DeltaTau = " << parameter_list.delta_tau << " fm";
    music_message.flush("info");

    // Total_evolution_time_tau
    // total evolution time in [fm]. in case of freeze_out_method = 2,3,4
    // evolution will halt earlier if all cells are frozen out.
    parameter_list.tau_size = (
        getParameter(input_file, "Total_evolution_time_tau", 50.));
    parameter_list.nt = static_cast<int>(
            parameter_list.tau_size/(parameter_list.delta_tau) + 0.5);
    music_message << "read_in_parameters: Time step size = "
                  << parameter_list.delta_tau;
    music_message.flush("info");
    music_message << "read_in_parameters: Number of time steps required = "
                  << parameter_list.nt;
    music_message.flush("info");

    // switch for baryon current propagation
    parameter_list.turn_on_rhob = (
        getParameter(input_file, "Include_Rhob_Yes_1_No_0", 0));
    if (parameter_list.turn_on_rhob == 1) {
       parameter_list.alpha_max = 5;
    } else {
       parameter_list.alpha_max = 4;
    }

    // Runge_Kutta_order:  must be 1 or 2
    parameter_list.rk_order = getParameter(input_file, "Runge_Kutta_order", 2);
    // Minmod_Theta: theta parameter in the min-mod like limiter
    parameter_list.minmod_theta = getParameter(input_file, "Minmod_Theta", 1.8);

    // EOS_to_use:
    // 0: ideal gas
    // 1: EOS-Q from azhydro
    // 2: lattice EOS from Huovinen and Petreczky
    // 3: lattice EOS from Huovinen and Petreczky
    //    with partial chemical equilibrium (PCE) at 150 MeV
    //    (see https://wiki.bnl.gov/TECHQM/index.php/QCD_Equation_of_State)
    // 4: PCE EOS with chemical freeze out at 155 MeV
    // 5: PCE EOS at 160 MeV
    // 6: PCE EOS at 165 MeV
    // 7: lattice EOS with CE match with UrQMD
    // 9: hotQCD (UrQMD); 91: hotQCD (SMASH)
    // 10: finite muB EOS from A. Monnai (up to mu_B^4)
    // 11: finite muB EOS from Pasi
    // 12: NEoS-B
    // 13: NEoS-BS
    // 14: NEoS-BQS
    // 17: BEST-EOS
    // 18: UH-EOS
    // 20: NEoS-4D
    parameter_list.whichEOS = getParameter(input_file, "EOS_to_use", 9);

    // Viscosity_Flag_Yes_1_No_0:   set to 0 for ideal hydro
    parameter_list.viscosity_flag = (
        getParameter(input_file, "Viscosity_Flag_Yes_1_No_0", 1));
    // Include_Shear_Visc_Yes_1_No_0
    parameter_list.turn_on_shear = (
        getParameter(input_file, "Include_Shear_Visc_Yes_1_No_0", 0));
    // Include_Bulk_Visc_Yes_1_No_0
    parameter_list.turn_on_bulk = (
        getParameter(input_file, "Include_Bulk_Visc_Yes_1_No_0", 0));
    // Include secord order terms
    parameter_list.include_second_order_terms = (
        getParameter(input_file, "Include_second_order_terms", 0));
    int tempturn_on_vorticity_terms = (
        getParameter(input_file, "Include_vorticity_terms", 0));
    if (tempturn_on_vorticity_terms == 0) {
        parameter_list.include_vorticity_terms = false;
    } else {
        parameter_list.include_vorticity_terms = true;
    }
    parameter_list.turn_on_diff = (
        getParameter(input_file, "turn_on_baryon_diffusion", 0));

    // the strength for the viscous regulation
    parameter_list.quest_revert_strength = (
        getParameter(input_file, "quest_revert_strength", 10.));

    //Shear_to_S_ratio:  constant eta/s
    parameter_list.shear_to_s = (
        getParameter(input_file, "Shear_to_S_ratio", 0.08));
    // Relaxation time factors
    parameter_list.shear_relax_time_factor = (
        getParameter(input_file, "shear_relax_time_factor", 5.0));

    parameter_list.T_dependent_shear_to_s = (
        getParameter(input_file, "T_dependent_Shear_to_S_ratio", 0));
    // If "T_dependent_Shear_to_S_ratio==2",
    // (eta/s)(T) = eta_over_s_min + eta_over_s_slope*(T − Tc)*(T/Tc)^{eta_over_s_curv}
    // with T_c=0.154 GeV
    parameter_list.shear_2_min = (
        getParameter(input_file, "shear_viscosity_2_min", 0.08));
    parameter_list.shear_2_slope = (
        getParameter(input_file, "shear_viscosity_2_slope", 1.0));
    parameter_list.shear_2_curv = (
        getParameter(input_file, "shear_viscosity_2_curv", 0.0));

    // If "T_dependent_Shear_to_S_ratio==3",
    parameter_list.shear_3_T_kink_in_GeV = (
        getParameter(input_file,
                     "shear_viscosity_3_eta_over_s_T_kink_in_GeV", 0.16));
    parameter_list.shear_3_low_T_slope_in_GeV = (
        getParameter(input_file,
                     "shear_viscosity_3_eta_over_s_low_T_slope_in_GeV", 0.0));
    parameter_list.shear_3_high_T_slope_in_GeV = (
        getParameter(input_file,
                     "shear_viscosity_3_eta_over_s_high_T_slope_in_GeV", 0.0));
    parameter_list.shear_3_at_kink = (
        getParameter(input_file, "shear_viscosity_3_eta_over_s_at_kink", 0.08));

    parameter_list.muB_dependent_shear_to_s = (
        getParameter(input_file, "muB_dependent_Shear_to_S_ratio", 0));
    parameter_list.shear_muBf0p4 = (
        getParameter(input_file, "shear_muBf0p4", 1.));
    parameter_list.shear_muBf0p2 = (
        getParameter(input_file, "shear_muBf0p2",
                     (1. + parameter_list.shear_muBf0p4)/2.));
    parameter_list.shear_muBDep_alpha = (
            getParameter(input_file, "shear_muBDep_alpha", 1.));
    parameter_list.shear_muBDep_slope = (
        getParameter(input_file, "shear_muBDep_slope", 1.0));
    parameter_list.shear_muBDep_scale = (
        getParameter(input_file, "shear_muBDep_scale", 0.6));

    // type of bulk relaxation time parameterization
    parameter_list.bulk_relaxation_type = (
        getParameter(input_file, "Bulk_relaxation_time_type", 0));
    parameter_list.bulk_relax_time_factor = (
        getParameter(input_file, "bulk_relax_time_factor", 1./14.55));

    parameter_list.T_dependent_bulk_to_s = (
        getParameter(input_file, "T_dependent_Bulk_to_S_ratio", 1));
    parameter_list.T_dependent_bulk_to_s = (
        getParameter(input_file, "T_dependent_zeta_over_s",
                     parameter_list.T_dependent_bulk_to_s));
    // "T_dependent_Bulk_to_S_ratio=2",
    // bulk viscosity is parametrized as with "A", "G" and "Tc" as
    // "A*(1/(1+((T-Tc)/G)^2)"
    parameter_list.bulk_2_normalisation = (
        getParameter(input_file, "bulk_viscosity_2_normalisation", 0.33));
    parameter_list.bulk_2_width_in_GeV = (
        getParameter(input_file, "bulk_viscosity_2_width_in_GeV", 0.08));
    parameter_list.bulk_2_peak_in_GeV = (
        getParameter(input_file, "bulk_viscosity_2_peak_in_GeV", 0.18));

    // "T_dependent_Bulk_to_S_ratio==3",
    parameter_list.bulk_3_max = (
        getParameter(input_file, "bulk_viscosity_3_zeta_over_s_max", 0.1));
    parameter_list.bulk_3_width_in_GeV = (
        getParameter(input_file,
                     "bulk_viscosity_3_zeta_over_s_width_in_GeV", 0.05));
    parameter_list.bulk_3_T_peak_in_GeV = (
        getParameter(input_file,
                     "bulk_viscosity_3_zeta_over_s_T_peak_in_GeV", 0.18));
    parameter_list.bulk_3_lambda_asymm = (
        getParameter(input_file,
                     "bulk_viscosity_3_zeta_over_s_lambda_asymm", 0.0));

    // "T_dependent_Bulk_to_S_ratio==10",
    parameter_list.bulk_10_max = (
        getParameter(input_file, "bulk_viscosity_10_max", 0.0));
    parameter_list.bulk_10_max_muB0p4 = (
        getParameter(input_file, "bulk_viscosity_10_max_muB0p4",
                     parameter_list.bulk_10_max));
    parameter_list.bulk_10_max_muB0p2 = (
        getParameter(input_file, "bulk_viscosity_10_max_muB0p2",
                     (parameter_list.bulk_10_max
                      + parameter_list.bulk_10_max_muB0p4)/2.));
    parameter_list.bulk_10_width_high = (
        getParameter(input_file, "bulk_viscosity_10_width_high", 0.1));  // GeV
    parameter_list.bulk_10_width_low = (
        getParameter(input_file, "bulk_viscosity_10_width_low", 0.015)); // GeV
    parameter_list.bulk_10_Tpeak = (
        getParameter(input_file, "bulk_viscosity_10_T_peak", 0.17));     // GeV
    parameter_list.bulk_10_Tpeak_muBcurv = (
        getParameter(input_file, "bulk_viscosity_10_T_peak_muBcurv", 0.0));

    // net baryon diffusion: kappa coefficient
    parameter_list.kappa_coefficient = (
        getParameter(input_file, "kappa_coefficient", 0.0));

    parameter_list.store_hydro_info_in_memory = (
        getParameter(input_file, "store_hydro_info_in_memory", 0));
    // output_evolution_data:
    // 2: output bulk information in binary format
    parameter_list.outputEvolutionData = (
        getParameter(input_file, "output_evolution_data", 0));

    // only works for output_evolution_data == 1
    parameter_list.outputBinaryEvolution = (
        getParameter(input_file, "outputBinaryEvolution", 0));
    parameter_list.output_hydro_params_header = (
        getParameter(input_file, "output_hydro_params_header", 0));

    parameter_list.output_movie_flag = (
            getParameter(input_file, "output_movie_flag", 0));
    parameter_list.output_outofequilibriumsize = (
            getParameter(input_file, "output_outofequilibriumsize", 0));
    parameter_list.output_vorticity = (
        getParameter(input_file, "output_vorticity", 0));
    parameter_list.output_hydro_debug_info = (
        getParameter(input_file, "output_hydro_debug_info", 0));

    // The evolution is outputted every
    // "output_evolution_every_N_timesteps" timesteps
    parameter_list.output_evolution_every_N_timesteps = (
        getParameter(input_file, "output_evolution_every_N_timesteps", 1));
    parameter_list.output_evolution_every_N_x = (
        getParameter(input_file, "output_evolution_every_N_x", 1));
    parameter_list.output_evolution_every_N_y = (
        parameter_list.output_evolution_every_N_x);
    parameter_list.output_evolution_every_N_eta = (
        getParameter(input_file, "output_evolution_every_N_eta", 1));
    parameter_list.output_evolution_T_cut = (
        getParameter(input_file, "output_evolution_T_cut", 0.105));
    parameter_list.output_evolution_e_cut = (
        getParameter(input_file, "output_evolution_e_cut", 0.15));  // GeV/fm^3

    //////////////////////////////////////////////////////////////////////////
    // Freeze-out surface parameters
    //////////////////////////////////////////////////////////////////////////
    // Do_FreezeOut_Yes_1_No_0
    // set to 0 to bypass freeze out surface finder
    parameter_list.doFreezeOut = (
        getParameter(input_file, "Do_FreezeOut_Yes_1_No_0", 1));
    parameter_list.doFreezeOut_lowtemp = (
        getParameter(input_file, "Do_FreezeOut_lowtemp", 1));
    // freeze_out_method:
    // 2: Schenke's more complex method
    // 4: Cornelius
    parameter_list.freezeOutMethod = (
            getParameter(input_file, "freeze_out_method", 4));

    // use_eps_for_freeze_out:
    // 0: freeze out at constant temperature T_freeze
    // 1: freeze out at constant energy density epsilon_freeze
    // if set in input input_file, overide above defaults
    parameter_list.useEpsFO = (
            getParameter(input_file, "use_eps_for_freeze_out", 1));
    parameter_list.freeze_eps_flag = (
            getParameter(input_file, "freeze_eps_flag", 0));

    // T_freeze: freeze out temperature
    // only used with use_eps_for_freeze_out = 0
    if (parameter_list.useEpsFO == 0) {
        parameter_list.TFO = getParameter(input_file, "T_freeze", 0.12);
        // only one freeze-out is allowed
        parameter_list.N_freeze_out = 1;
    } else {
        // epsilon_freeze: freeze-out energy density in GeV/fm^3
        // only used with use_eps_for_freeze_out = 1
        // allow multiple parameter names for backward compitibility
        parameter_list.epsilonFreeze = (
            getParameter(input_file, "eps_switch", 0.18));      // GeV/fm^3
        parameter_list.epsilonFreeze = (
            getParameter(input_file, "epsilon_freeze",
                         parameter_list.epsilonFreeze));        // GeV/fm^3
        parameter_list.epsilonFreeze = (
            getParameter(input_file, "eps_freeze_max",
                         parameter_list.epsilonFreeze));        // GeV/fm^3
        parameter_list.N_freeze_out = (
            getParameter(input_file, "N_freeze_out", 1));
    }
    if (parameter_list.N_freeze_out > 1) {
        parameter_list.eps_freeze_max = (
            getParameter(input_file, "eps_freeze_max", 0.18));
        parameter_list.eps_freeze_min = (
            getParameter(input_file, "eps_freeze_min", 0.18));
    } else {
        parameter_list.eps_freeze_min = parameter_list.epsilonFreeze;
        parameter_list.eps_freeze_max = parameter_list.epsilonFreeze;
    }

    //! Maximum starting time for freeze-out surface
    parameter_list.freezeOutTauStartMax = (
            getParameter(input_file, "freeze_out_tau_start_max", 2.));

    string temp_freeze_list_filename = "eps_freeze_list_s95p_v1.dat";
    tempinput = Util::StringFind4(input_file, "freeze_list_filename");
    if (tempinput != "empty")
        temp_freeze_list_filename.assign(tempinput);
    parameter_list.freeze_list_filename.assign(temp_freeze_list_filename);

    int temp_freeze_surface_binary = (
            getParameter(input_file, "freeze_surface_in_binary", 1));
    if (temp_freeze_surface_binary == 0) {
        parameter_list.freeze_surface_in_binary = false;
    } else {
        parameter_list.freeze_surface_in_binary = true;
    }

    // average_surface_over_this_many_time_steps:
    // Only save every N timesteps for finding freeze out surface
    parameter_list.facTau = (
            getParameter(input_file,
                         "average_surface_over_this_many_time_steps", 1));
    parameter_list.fac_x = getParameter(input_file, "freeze_Ncell_x_step", 1);
    parameter_list.fac_y = parameter_list.fac_x;
    parameter_list.fac_eta = (
            getParameter(input_file, "freeze_Ncell_eta_step", 1));


    //////////////////////////////////////////////////////////////////////////
    // Cooper-Frye (mode 3)
    //////////////////////////////////////////////////////////////////////////
    // particle_spectrum_to_compute:
    // 0: Do all up to number_of_particles_to_include
    // any natural number: Do the particle with this (internal) ID
    parameter_list.particleSpectrumNumber = (
            getParameter(input_file, "particle_spectrum_to_compute", 0));

    // number_of_particles_to_include:
    // This determines up to which particle in the list spectra
    // should be computed (mode=3) or resonances should be included (mode=4)
    // current maximum = 319
    parameter_list.NumberOfParticlesToInclude = (
        getParameter(input_file, "number_of_particles_to_include", 2));

    // for calculation of spectra:
    // max_pseudorapidity:
    // spectra calculated from zero to this pseudorapidity in +eta and -eta
    parameter_list.max_pseudorapidity = (
        getParameter(input_file, "max_pseudorapidity", 5.0));
    parameter_list.pseudo_steps = getParameter(input_file, "pseudo_steps", 51);

    // steps in azimuthal angle in calculation of spectra
    parameter_list.phi_steps = getParameter(input_file, "phi_steps", 48);

    // spectra calculated from min_pt to max_pt transverse momentum in GeV
    parameter_list.min_pt = getParameter(input_file, "min_pt", 0.0);
    parameter_list.max_pt = getParameter(input_file, "max_pt", 3.0);
    // steps in transverse momentum in calculation of spectra
    parameter_list.pt_steps = getParameter(input_file, "pt_steps", 30);

    // pseudofreeze
    // Calculate spectra at fixed,
    // equally-spaced grid in pseudorapidity, pt, and phi
    parameter_list.pseudofreeze = getParameter(input_file, "pseudofreeze", 1);

    // Include_deltaf:
    // Looks like 0 sets delta_f=0, 1 uses standard quadratic ansatz,
    // and 2 is supposed to use p^(2-alpha)
    parameter_list.include_deltaf = (
        getParameter(input_file, "Include_deltaf", 1));
    parameter_list.include_deltaf_bulk = (
        getParameter(input_file, "Include_deltaf_bulk", 0));
    parameter_list.include_deltaf_qmu = (
        getParameter(input_file, "Include_deltaf_qmu", 0));
    parameter_list.deltaf_14moments = (
        getParameter(input_file, "deltaf_14moments", 0));

    //////////////////////////////////////////////////////////////////////////
    // spectra vn analysis for mode 13 and 14
    //////////////////////////////////////////////////////////////////////////
    parameter_list.dNdy_y_min = getParameter(input_file, "dNdy_y_min", -0.5);
    parameter_list.dNdy_y_max = getParameter(input_file, "dNdy_y_max", 0.5);
    parameter_list.dNdy_eta_min = (
        getParameter(input_file, "dNdy_eta_min", -2.0));
    parameter_list.dNdy_eta_max = (
        getParameter(input_file, "dNdy_eta_max", 2.0));
    parameter_list.dNdy_nrap = getParameter(input_file, "dNdy_nrap", 30);

    parameter_list.dNdyptdpt_y_min = (
        getParameter(input_file, "dNdyptdpt_y_min", -0.5));
    parameter_list.dNdyptdpt_y_max = (
        getParameter(input_file, "dNdyptdpt_y_max", 0.5));
    parameter_list.dNdyptdpt_eta_min = (
        getParameter(input_file, "dNdyptdpt_eta_min", -0.5));
    parameter_list.dNdyptdpt_eta_max = (
        getParameter(input_file, "dNdyptdpt_eta_max", 0.5));

    music_message.info("Done read_in_parameters.");
    check_parameters(parameter_list, input_file);

    return parameter_list;
}


void set_parameter(InitData &parameter_list, std::string parameter_name,
                   double value) {
    if (parameter_name == "MUSIC_mode")
        parameter_list.mode = static_cast<int>(value);

    if (parameter_name == "Initial_time_tau_0")
        parameter_list.tau0 = value;

    if (parameter_name == "output_evolution_data")
        parameter_list.outputEvolutionData = static_cast<int>(value);

    if (parameter_name == "output_movie_flag")
        parameter_list.output_movie_flag = static_cast<int>(value);

    if (parameter_name == "store_hydro_info_in_memory")
        parameter_list.store_hydro_info_in_memory = static_cast<int>(value);

    if (parameter_name == "Viscosity_Flag_Yes_1_No_0")
        parameter_list.viscosity_flag = static_cast<int>(value);

    if (parameter_name == "Include_Shear_Visc_Yes_1_No_0")
        parameter_list.turn_on_shear = static_cast<int>(value);

    if (parameter_name == "Shear_to_S_ratio")
        parameter_list.shear_to_s = value;

    if (parameter_name == "T_freeze")
        parameter_list.TFO = value;

    if (parameter_name == "Include_Bulk_Visc_Yes_1_No_0")
        parameter_list.turn_on_bulk = static_cast<int>(value);

    if (parameter_name == "Include_second_order_terms")
        parameter_list.include_second_order_terms = static_cast<int>(value);

    if (parameter_name == "T_dependent_Shear_to_S_ratio")
        parameter_list.T_dependent_shear_to_s = static_cast<int>(value);

    if (parameter_name == "shear_viscosity_2_min")
        parameter_list.shear_2_min = value;
    if (parameter_name == "shear_viscosity_slope")
        parameter_list.shear_2_slope = value;
    if (parameter_name == "shear_viscosity_curv")
        parameter_list.shear_2_curv = value;

    if (parameter_name == "shear_viscosity_3_T_kink_in_GeV")
        parameter_list.shear_3_T_kink_in_GeV = value;
    if (parameter_name == "shear_viscosity_3_low_T_slope_in_GeV")
        parameter_list.shear_3_low_T_slope_in_GeV = value;
    if (parameter_name == "shear_viscosity_3_high_T_slope_in_GeV")
        parameter_list.shear_3_high_T_slope_in_GeV = value;
    if (parameter_name == "shear_viscosity_3_at_kink")
        parameter_list.shear_3_at_kink = value;

    if (parameter_name == "T_dependent_Bulk_to_S_ratio")
        parameter_list.T_dependent_bulk_to_s = static_cast<int>(value);

    if (parameter_name == "bulk_viscosity_2_normalisation")
        parameter_list.bulk_2_normalisation = value;
    if (parameter_name == "bulk_viscosity_2_peak_in_GeV")
        parameter_list.bulk_2_peak_in_GeV = value;
    if (parameter_name == "bulk_viscosity_2_width_in_GeV")
        parameter_list.bulk_2_width_in_GeV = value;

    if (parameter_name == "bulk_viscosity_3_max")
        parameter_list.bulk_3_max = value;
    if (parameter_name == "bulk_viscosity_3_width_in_GeV")
        parameter_list.bulk_3_width_in_GeV = value;
    if (parameter_name == "bulk_viscosity_3_T_peak_in_GeV")
        parameter_list.bulk_3_T_peak_in_GeV = value;
    if (parameter_name == "bulk_viscosity_3_lambda_asymm")
        parameter_list.bulk_3_lambda_asymm = value;
}


void check_parameters(InitData &parameter_list, std::string input_file) {
    music_message.info("Checking input parameter list ... ");

    if (parameter_list.Initial_profile < 0) {
        music_message << "Initial profile" << parameter_list.Initial_profile
                      << "not defined";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.initializeEntropy > 1
            || parameter_list.initializeEntropy < 0) {
        music_message.error("Must initialize with entropy or energy");
        exit(1);
    }

    if (parameter_list.useEpsFO == 0) {
        if (parameter_list.TFO < 0.0 || parameter_list.TFO > 0.2) {
            music_message << "T_freeze = " << parameter_list.TFO
                          << " is not physical!";
            music_message.flush("error");
            exit(1);
        }
    } else {
        if (parameter_list.epsilonFreeze <= 0) {
            music_message.error(
                    "Freeze out energy density must be greater than zero");
            exit(1);
        }
    }

    if (parameter_list.useEpsFO > 1 || parameter_list.useEpsFO < 0) {
        music_message << "did not set either freeze out energy density "
                      << "or temperature, or invalid option for "
                      << "use_eps_for_freeze_out:"
                      << parameter_list.useEpsFO;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.useEpsFO == 0 && parameter_list.turn_on_rhob == 1) {
        music_message << "freeze-out surface set by temperature is not "
                      << "support yet. reset use_eps_for_freeze_out to 1.";
        music_message.flush("warning");
        parameter_list.useEpsFO = 1;
    }

    if ((parameter_list.whichEOS > 19 && parameter_list.whichEOS != 91)
        || parameter_list.whichEOS < 0) {
        music_message << "EOS_to_use unspecified or invalid option: "
                      << parameter_list.whichEOS;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.whichEOS > 1 && parameter_list.whichEOS < 7
            && parameter_list.NumberOfParticlesToInclude > 320) {
        music_message << "Invalid option for number_of_particles_to_include:"
                      << parameter_list.NumberOfParticlesToInclude;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.freezeOutMethod != 4) {
        music_message << "Invalid option for freeze_out_method: "
                      << parameter_list.freezeOutMethod;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.facTau <= 0) {
        music_message << "average_surface_over_this_many_time_steps <= 0: "
                      << parameter_list.facTau;
        exit(1);
    }

    double freeze_dtau = parameter_list.facTau*parameter_list.delta_tau;
    if (freeze_dtau > 1.) {
        music_message << "freeze-out time setp is too large! "
                      << "freeze_dtau = " << freeze_dtau
                      << ", hydro_dtau = " << parameter_list.delta_tau
                      << ", average_surface_over_this_many_time_steps = "
                      << parameter_list.facTau;
        music_message.flush("warning");
    }

    if (parameter_list.fac_x <= 0) {
        music_message << "freeze out every x step <= 0: "
                      << parameter_list.fac_x;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.fac_eta <= 0) {
        music_message << "freeze out every eta step <= 0: "
                      << parameter_list.fac_eta;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.nx != parameter_list.ny) {
        music_message << "Grid size in x is not equal to grid size in y!";
        music_message.flush("warning");
    }

    if (parameter_list.neta < 2 && !parameter_list.boost_invariant) {
        music_message << "Grid size in eta = " << parameter_list.neta 
                      << "is too small for a (3+1)-d run! "
                      << "Please increase Grid_size_in_eta to be "
                      << "larger than 2 at least!";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.boost_invariant && parameter_list.neta > 1) {
        music_message << "Grid size in eta is set to "
                      << parameter_list.neta << " for a (2+1)-d simulation! "
                      << "This is redundant! Reset neta to 1! ";
        music_message.flush("warning");
        parameter_list.neta = 1;
    }

    if (parameter_list.boost_invariant) {
        music_message << "For a (2+1)-d simulation, "
                      << "reset deta = 0.1 and eta_size = 0.0";
        music_message.flush("info");
        parameter_list.delta_eta = 0.1;
        parameter_list.eta_size = 0.0;
    }

    if (parameter_list.delta_tau/parameter_list.delta_x
            > parameter_list.dtaudxRatio) {
        music_message << "Warning: Delta_Tau = " << parameter_list.delta_tau
                      << " maybe too large! ";
        music_message.flush("warning");

        bool reset_dtau_use_CFL_condition = true;
        int temp_CFL_condition = 1;
        string tempinput = Util::StringFind4(
                            input_file, "reset_dtau_use_CFL_condition");
        if (tempinput != "empty")
            istringstream(tempinput) >> temp_CFL_condition;
        if (temp_CFL_condition == 0)
            reset_dtau_use_CFL_condition = false;
        parameter_list.resetDtau = reset_dtau_use_CFL_condition;

        if (reset_dtau_use_CFL_condition) {
            music_message.info("reset dtau using CFL condition.");
            double dtau_CFL = std::min(
                std::min(parameter_list.delta_x*parameter_list.dtaudxRatio,
                         parameter_list.delta_y*parameter_list.dtaudxRatio),
                         parameter_list.tau0*parameter_list.delta_eta
                         *parameter_list.dtaudxRatio);
            parameter_list.delta_tau = dtau_CFL;
            parameter_list.nt = static_cast<int>(
                parameter_list.tau_size/(parameter_list.delta_tau) + 0.5);
            music_message << "read_in_parameters: Time step size = "
                          << parameter_list.delta_tau;
            music_message.flush("info");
            music_message << "read_in_parameters: "
                          << "Number of time steps required = "
                          << parameter_list.nt;
            music_message.flush("info");
        } else {
            exit(1);
        }
    }

    if (parameter_list.min_pt > parameter_list.max_pt) {
        music_message << "min_pt = " << parameter_list.min_pt << " > "
                      << "max_pt = " << parameter_list.max_pt;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.phi_steps < 20) {
        music_message << "phi_steps = " << parameter_list.phi_steps
                      << " is too small for computing vn, "
                      << "please increase to at least 40!";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.rk_order > 2 || parameter_list.rk_order < 0) {
        music_message << "Invalid option for Runge_Kutta_order: "
                      << parameter_list.rk_order;
        music_message.flush("error");
        exit(1);
    }
    if (parameter_list.rk_order != 2) {
        music_message << "Runge-Kutta order = " << parameter_list.rk_order;
        music_message.flush("info");
    }

    if (parameter_list.minmod_theta < 1.
            || parameter_list.minmod_theta > 2.) {
        music_message << "minmod = " << parameter_list.minmod_theta
                      << " is out of allowed range [1., 2]";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.turn_on_shear == 0 && parameter_list.shear_to_s > 0) {
        music_message << "non-zero eta/s = " << parameter_list.shear_to_s
                      << " is set with "
                      << "Include_Shear_Visc = "
                      << parameter_list.turn_on_shear
                      << ". Please check you want to run ideal hydro!";
        music_message.flush("warning");
    }

    if (parameter_list.turn_on_shear == 0
            && parameter_list.include_deltaf == 1) {
        music_message << "hydro with zero shear viscosity does not need "
                      << "shear delta f in Cooper-Frye! ";
        music_message << "input Include_deltaf = "
                      << parameter_list.include_deltaf
                      << ". Now rewrite it to 0!";
        music_message.flush("warning");
        parameter_list.include_deltaf = 0;
    }

    if (parameter_list.turn_on_bulk == 0
            && parameter_list.include_deltaf_bulk == 1) {
        music_message << "hydro with zero bulk viscosity does not need "
                      << "bulk delta f in Cooper-Frye!";
        music_message << "input Include_deltaf_bulk = "
                      << parameter_list.include_deltaf_bulk
                      << ". Now rewrite it to 0!";
        music_message.flush("warning");
        parameter_list.include_deltaf_bulk = 0;
    }

    if (parameter_list.output_evolution_every_N_timesteps <= 0) {
        music_message.error("output_evolution_every_N_timesteps < 0!");
        exit(1);
    }

    if (parameter_list.output_evolution_every_N_x <= 0) {
        music_message.error("output_evolution_every_N_x < 0!");
        exit(1);
    }

    if (parameter_list.output_evolution_every_N_y <= 0) {
        music_message.error("output_evolution_every_N_y < 0!");
        exit(1);
    }

    if (parameter_list.output_evolution_every_N_eta <= 0) {
        music_message.error("output_evolution_every_N_eta < 0!");
        exit(1);
    }

    if (parameter_list.dNdy_y_min > parameter_list.dNdy_y_max) {
        music_message << "dNdy_y_min = " << parameter_list.dNdy_y_min << " < " 
                      << "dNdy_y_max = " << parameter_list.dNdy_y_max << "!";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.dNdy_eta_min > parameter_list.dNdy_eta_max) {
        music_message << "dNdy_eta_min = " << parameter_list.dNdy_eta_min
                      << " < " 
                      << "dNdy_eta_max = " << parameter_list.dNdy_eta_max
                      << "!";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.dNdyptdpt_y_min > parameter_list.dNdyptdpt_y_max) {
        music_message << "dNdyptdpt_y_min = "
                      << parameter_list.dNdyptdpt_y_min << " < " 
                      << "dNdyptdpt_y_max = "
                      << parameter_list.dNdyptdpt_y_max << "!";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list.dNdyptdpt_eta_min
            > parameter_list.dNdyptdpt_eta_max) {
        music_message << "dNdyptdpt_eta_min = "
                      << parameter_list.dNdyptdpt_eta_min << " < " 
                      << "dNdyptdpt_eta_max = "
                      << parameter_list.dNdyptdpt_eta_max << "!";
        music_message.flush("error");
        exit(1);
    }

    music_message << "Finished checking input parameter list. "
                  << "Everything looks reasonable so far "
                  << emoji::success() << emoji::thumbup()
                  << emoji::beer() << emoji::beer()
                  << emoji::beerclinking() << emoji::beerclinking();
    music_message.flush("info");
}

}
