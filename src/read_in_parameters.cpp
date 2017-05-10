
#include <iostream>
#include <cstring>
#include "./read_in_parameters.h"

using namespace std;

ReadInParameters::ReadInParameters() {
    util = new Util();
    emoji_face = new emoji();
}

ReadInParameters::~ReadInParameters() {
    delete util;
    delete emoji_face;
}


void ReadInParameters::read_in_parameters(InitData *parameter_list,
                                          string input_file_in) {
    input_file = input_file_in;
    // this function reads in parameters
    int m, n;
    string tempinput;

    // echo_level controls the mount of
    // warning message output during the evolution
    double temp_echo_level = 9;
    tempinput = util->StringFind4(input_file, "echo_level");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_echo_level;
    parameter_list->echo_level = temp_echo_level;

    
    // Initial_profile: 
    int tempInitial_profile = 1;
    tempinput = util->StringFind4(input_file, "Initial_profile");
    if (tempinput != "empty") istringstream(tempinput) >> tempInitial_profile;
    parameter_list->Initial_profile = tempInitial_profile;

    // Initial_profile: 
    int temp_string_dump_mode = 1;
    tempinput = util->StringFind4(input_file, "string_dump_mode");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_string_dump_mode;
    parameter_list->string_dump_mode = temp_string_dump_mode;
    
    // hydro source
    double temp_string_quench_factor = 1.;
    tempinput = util->StringFind4(input_file, "string_quench_factor");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_string_quench_factor;
    parameter_list->string_quench_factor = temp_string_quench_factor;
    
    // hydro source
    double temp_parton_quench_factor = 1.;
    tempinput = util->StringFind4(input_file, "parton_quench_factor");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_parton_quench_factor;
    parameter_list->parton_quench_factor = temp_parton_quench_factor;

    // boost-invariant
    int temp_boost_invariant = 0;
    tempinput = util->StringFind4(input_file, "boost_invariant");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_boost_invariant;
    if (temp_boost_invariant == 0) {
        parameter_list->boost_invariant = false;
    } else {
        parameter_list->boost_invariant = true;
    }
    
    int temp_output_initial_profile = 0;
    tempinput = util->StringFind4(input_file,
                                  "output_initial_density_profiles");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_output_initial_profile;
    parameter_list->output_initial_density_profiles =
                                            temp_output_initial_profile;

    //Select the profile to use in eta for the energy/entropy initialisation
    //1 for Hirano's central plateau + Gaussian decay
    //2 for a Woods-Saxon proinput_file
    int tempinitial_eta_profile = 1;
    tempinput = util->StringFind4(input_file, "initial_eta_profile");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempinitial_eta_profile;
    parameter_list->initial_eta_profile = tempinitial_eta_profile;
    
    // eta envelope function parameter for rhob
    int temp_rhob_flag = 1;
    tempinput = util->StringFind4(input_file, "initial_eta_rhob_profile");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_rhob_flag;
    parameter_list->initial_eta_rhob_profile = temp_rhob_flag;
    
    int temp_check_eos = 0;
    tempinput = util->StringFind4(input_file, "check_eos");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_check_eos;
    parameter_list->check_eos = temp_check_eos;
    
    //initialize_with_entropy:
    //0: scale with energy density
    //1: scale with entropy density
    int tempinitializeEntropy = 0;
    tempinput = util->StringFind4(input_file, "initialize_with_entropy");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempinitializeEntropy;
    parameter_list->initializeEntropy = tempinitializeEntropy;
    
    // T_freeze: freeze out temperature
    // only used with use_eps_for_freeze_out = 0
    double tempTFO = 0.12;
    int tfoset = 0;
    tempinput = util->StringFind4(input_file, "T_freeze");
    if (tempinput != "empty") {
        istringstream(tempinput) >> tempTFO;
        // if only freeze out temperature is set, freeze out by temperature
        parameter_list->useEpsFO = 0;
        tfoset = 1;
    }
    parameter_list->TFO = tempTFO;
    
    // epsilon_freeze: freeze-out energy density in GeV/fm^3
    // only used with use_eps_for_freeze_out = 1
    double tempepsilonFreeze = 0.12;
    tempinput = util->StringFind4(input_file, "epsilon_freeze");
    if (tempinput != "empty") {
        istringstream(tempinput) >> tempepsilonFreeze;
        // if epsilon_freeze is set, freeze out by epsilon
        parameter_list->useEpsFO = 1;
    }
    parameter_list->epsilonFreeze = tempepsilonFreeze;
    
    int temp_N_freeze_out = 1;
    tempinput = util->StringFind4(input_file, "N_freeze_out");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_N_freeze_out;
    parameter_list->N_freeze_out = temp_N_freeze_out;

    //use_eps_for_freeze_out: 
    // 0: freeze out at constant temperature T_freeze
    // 1: freeze out at constant energy density epsilon_freeze
    // if set in input input_file, overide above defaults
    int tempuseEpsFO = parameter_list->useEpsFO;
    tempinput = util->StringFind4(input_file, "use_eps_for_freeze_out");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempuseEpsFO;
    parameter_list->useEpsFO = tempuseEpsFO;
    if (tfoset == 1 && parameter_list->useEpsFO == 1) {
        cerr << "T_freeze set but overridden"
             << " -- freezing out by energy density at " 
             << parameter_list->epsilonFreeze << " GeV/fm^3\n";
    }
    
    string temp_freeze_list_filename = "eps_freeze_list_s95p_v1.dat";
    tempinput = util->StringFind4(input_file, "freeze_list_filename");
    if (tempinput != "empty")
        temp_freeze_list_filename.assign(tempinput);
    parameter_list->freeze_list_filename.assign(temp_freeze_list_filename);
    
    double temp_eps_freeze_max = 0.18;
    tempinput = util->StringFind4(input_file, "eps_freeze_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_eps_freeze_max;
    parameter_list->eps_freeze_max = temp_eps_freeze_max;
    
    double temp_eps_freeze_min = 0.18;
    tempinput = util->StringFind4(input_file, "eps_freeze_min");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_eps_freeze_min;
    parameter_list->eps_freeze_min = temp_eps_freeze_min;
    
    int temp_freeze_eps_flag = 0;
    tempinput = util->StringFind4(input_file, "freeze_eps_flag");
    if (tempinput != "empty")
        istringstream (tempinput) >> temp_freeze_eps_flag;
    parameter_list->freeze_eps_flag = temp_freeze_eps_flag;

    //particle_spectrum_to_compute:
    // 0: Do all up to number_of_particles_to_include
    // any natural number: Do the particle with this (internal) ID
    int tempparticleSpectrumNumber = 0;
    tempinput = util->StringFind4(input_file, "particle_spectrum_to_compute");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempparticleSpectrumNumber;
    parameter_list->particleSpectrumNumber = tempparticleSpectrumNumber;
    
    // mode: 
    // 1: Does everything. Evolution. Computation of thermal spectra.
    //    Resonance decays. Observables.
    // 2: Evolution only.
    // 3: Compute all thermal spectra only.
    // 4: Resonance decays only.
    // 13: Compute observables from previously-computed thermal spectra
    // 14: Compute observables from post-decay spectra
    int tempmode = 1;
    tempinput = util->StringFind4(input_file, "mode");
    if (tempinput != "empty") {
        istringstream(tempinput) >> tempmode;
    } else {
        cerr << "Must specify mode. Exiting.\n";
        exit(1);
    }
    parameter_list->mode = tempmode;
    
    //EOS_to_use:
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
    // 10: finite muB EOS from A. Monnai (up to mu_B^4)
    // 11: finite muB EOS from Pasi
    // 12: finite muB EOS from A. Monnai (up to mu_B^6)
    int tempwhichEOS = 2;
    tempinput = util->StringFind4(input_file, "EOS_to_use");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempwhichEOS;
    parameter_list->whichEOS = tempwhichEOS;
    
    // number_of_particles_to_include:
    // This determines up to which particle in the list spectra
    // should be computed (mode=3) or resonances should be included (mode=4)
    // current maximum = 319
    int tempNumberOfParticlesToInclude = 2;
    tempinput = util->StringFind4(input_file,
                                  "number_of_particles_to_include");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempNumberOfParticlesToInclude;
    parameter_list->NumberOfParticlesToInclude = tempNumberOfParticlesToInclude;
    
    // freeze_out_method:
    // 2: Schenke's more complex method
    int tempfreezeOutMethod = 2;
    tempinput = util->StringFind4(input_file, "freeze_out_method");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempfreezeOutMethod;
    parameter_list->freezeOutMethod = tempfreezeOutMethod;

    // average_surface_over_this_many_time_steps:
    // Only save every N timesteps for finding freeze out surface
    int tempfacTau = 1;
    tempinput = util->StringFind4(input_file,
                                  "average_surface_over_this_many_time_steps");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempfacTau;
    parameter_list->facTau = tempfacTau;
    
    int tempfac_x = 1;
    tempinput = util->StringFind4(input_file, "freeze_Ncell_x_step");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempfac_x;
    parameter_list->fac_x = tempfac_x;
    parameter_list->fac_y = tempfac_x;
    
    int tempfac_eta = 1;
    tempinput = util->StringFind4(input_file, "freeze_Ncell_eta_step");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempfac_eta;
    parameter_list->fac_eta = tempfac_eta;
    
    // Grid_size_in_*
    // number of cells in x,y direction
    int tempnx = 10;
    tempinput = util->StringFind4(input_file, "Grid_size_in_x");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempnx;
    parameter_list->nx = tempnx;
    int tempny = 10;
    tempinput = util->StringFind4(input_file, "Grid_size_in_y");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempny;
    parameter_list->ny = tempny;
    
    // Grid_size_in_eta
    // number of cells in eta direction.
    // One cell is positioned at eta=0, 
    // half the cells are at negative eta,
    // the rest (one fewer) are at positive eta
    int tempneta = 1;
    tempinput = util->StringFind4(input_file, "Grid_size_in_eta");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempneta;
    parameter_list->neta = tempneta;
    
    // grid_size_in_fm:
    // total length of box in x,y direction in fm (minus delta_*)
    double tempx_size = 25.;
    tempinput = util->StringFind4(input_file, "X_grid_size_in_fm");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempx_size;
    parameter_list->x_size = tempx_size;
    double tempy_size = 25.;
    tempinput = util->StringFind4(input_file, "Y_grid_size_in_fm");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempy_size;
    parameter_list->y_size = tempy_size;
    
    
    // switch for baryon current propagation
    int tempturn_on_rhob = 0;
    tempinput = util->StringFind4(input_file, "Include_Rhob_Yes_1_No_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_rhob;
    parameter_list->turn_on_rhob = tempturn_on_rhob;
    if (parameter_list->turn_on_rhob == 1)
       parameter_list->alpha_max = 5;
    else
       parameter_list->alpha_max = 4;
    
    // Eta_grid_size:  total length of box in eta direction (minus delta_eta)
    // e.g., neta=8 and eta_size=8 has 8 cells that run from eta=-4 to eta=3
    double tempeta_size = 8.;
    tempinput = util->StringFind4(input_file, "Eta_grid_size");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempeta_size;
    parameter_list->eta_size = tempeta_size;
    
    // Total_evolution_time_tau
    // total evolution time in [fm]. in case of freeze_out_method = 2,3,
    // evolution will halt earlier if all cells are frozen out.
    double temptau_size = 50.;
    tempinput = util->StringFind4(input_file, "Total_evolution_time_tau");
    if (tempinput != "empty")
        istringstream(tempinput) >> temptau_size;
    parameter_list->tau_size = temptau_size;
    
    // Initial_time_tau_0:  in fm
    double temptau0 = 0.4;
    tempinput = util->StringFind4(input_file, "Initial_time_tau_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> temptau0;
    parameter_list->tau0 = temptau0;
    
    /* x-grid, for instance, runs from 0 to nx */
    parameter_list->delta_x =
            parameter_list->x_size/static_cast<double>(parameter_list->nx - 1); 
    parameter_list->delta_y =
            parameter_list->y_size/static_cast<double>(parameter_list->ny - 1); 
    parameter_list->delta_eta =
            parameter_list->eta_size/static_cast<double>(parameter_list->neta); 
    
    cerr << " DeltaX = " << parameter_list->delta_x << " fm" << endl;
    cerr << " DeltaY = " << parameter_list->delta_y << " fm" << endl;
    cerr << " DeltaETA = " << parameter_list->delta_eta << endl;
    
    // Delta_Tau: 
    // time step to use in [fm].
    double tempdelta_tau = 0.02;
    tempinput = util->StringFind4(input_file, "Delta_Tau");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempdelta_tau;
    parameter_list->delta_tau = tempdelta_tau;
    cerr << " DeltaTau = " << parameter_list->delta_tau << " fm" << endl;
    
    // output_evolution_data:  
    // 1: output bulk information at every grid point at every time step
    int tempoutputEvolutionData = 0;
    tempinput = util->StringFind4(input_file, "output_evolution_data");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutputEvolutionData;
    parameter_list->outputEvolutionData = tempoutputEvolutionData;
    
    int temp_output_movie_flag = 0;
    tempinput = util->StringFind4(input_file, "output_movie_flag");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_output_movie_flag;
    parameter_list->output_movie_flag = temp_output_movie_flag;
    
    parameter_list->nt = static_cast<int>(
            parameter_list->tau_size/(parameter_list->delta_tau) + 0.5);
    cout << "read_in_parameters: Time step size = "
         << parameter_list->delta_tau << endl;
    cout << "read_in_parameters: Number of time steps required = "
         << parameter_list->nt << endl;
    
    double temp_eta_0 = 3.0;
    tempinput = util->StringFind4(input_file, "eta_rhob_0");
    if (tempinput != "empty") istringstream (tempinput) >> temp_eta_0;
    parameter_list->eta_rhob_0 = temp_eta_0;
    double temp_eta_width = 1.0;
    tempinput = util->StringFind4(input_file, "eta_rhob_width");
    if (tempinput != "empty") istringstream (tempinput) >> temp_eta_width;
    parameter_list->eta_rhob_width = temp_eta_width;
    double temp_eta_plateau_height = 0.5;
    tempinput = util->StringFind4(input_file, "eta_rhob_plateau_height");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_eta_plateau_height;
    parameter_list->eta_rhob_plateau_height = temp_eta_plateau_height;
    double temp_eta_width_1 = 1.0;
    tempinput = util->StringFind4(input_file, "eta_rhob_width_1");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_eta_width_1;
    parameter_list->eta_rhob_width_1 = temp_eta_width_1;
    double temp_eta_width_2 = 1.0;
    tempinput = util->StringFind4(input_file, "eta_rhob_width_2");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_eta_width_2;
    parameter_list->eta_rhob_width_2 = temp_eta_width_2;
  
    // Eta_fall_off:
    // width of half-Gaussian on each side of a central pleateau in eta
    double tempeta_fall_off  = 0.4;
    tempinput = util->StringFind4(input_file, "Eta_fall_off");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempeta_fall_off ;
    parameter_list->eta_fall_off  = tempeta_fall_off;
    
    // Eta_plateau_size:
    // width of the flat region symmetrical around eta=0
    double tempeta_flat = 20.0;
    tempinput = util->StringFind4(input_file, "Eta_plateau_size");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempeta_flat;
    parameter_list->eta_flat = tempeta_flat;

    // s_factor:  for use with IP-Glasma initial conditions
    //   parameter_list->sFactor = util->DFind(input_file, "s_factor");
    double tempsFactor   = 20.;
    tempinput = util->StringFind4(input_file, "s_factor");
    if (tempinput != "empty") istringstream ( tempinput ) >> tempsFactor  ;
    parameter_list->sFactor   = tempsFactor;
    
    // for calculation of spectra:
    // max_pseudorapidity:
    // spectra calculated from zero to this pseudorapidity in +eta and -eta
    double tempmax_pseudorapidity = 5.0;
    tempinput = util->StringFind4(input_file, "max_pseudorapidity");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempmax_pseudorapidity;
    parameter_list->max_pseudorapidity = tempmax_pseudorapidity;
    
    // pseudo_steps:
    // steps in pseudorapidity in calculation of spectra
    int temppseudo_steps = 51;
    tempinput = util->StringFind4(input_file, "pseudo_steps");
    if (tempinput != "empty")
        istringstream(tempinput) >> temppseudo_steps;
    parameter_list->pseudo_steps = temppseudo_steps; 
    
    // phi_steps
    // steps in azimuthal angle in calculation of spectra
    int tempphi_steps = 48;
    tempinput = util->StringFind4(input_file, "phi_steps");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempphi_steps  ;
    parameter_list->phi_steps = tempphi_steps; 
    
    // min_pt: 
    // spectra calculated from this to max_pt transverse momentum in GeV
    double tempmin_pt   = 0.0;
    tempinput = util->StringFind4(input_file, "min_pt");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempmin_pt  ;
    parameter_list->min_pt = tempmin_pt;
      
    // max_pt: 
    // spectra calculated from min_pt to this transverse momentum in GeV
    double tempmax_pt   = 3.0;
    tempinput = util->StringFind4(input_file, "max_pt");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempmax_pt;
    parameter_list->max_pt = tempmax_pt;
    
    // pt_steps:
    // steps in transverse momentum in calculation of spectra
    int temppt_steps   = 60;
    tempinput = util->StringFind4(input_file, "pt_steps");
    if (tempinput != "empty")
        istringstream(tempinput) >> temppt_steps  ;
    parameter_list->pt_steps = temppt_steps;   
    
    // pseudofreeze
    // Calculate spectra at fixed,
    // equally-spaced grid in pseudorapidity, pt, and phi
    int temppseudofreeze = 1;
    tempinput = util->StringFind4(input_file, "pseudofreeze");
    if (tempinput != "empty")
        istringstream(tempinput) >> temppseudofreeze;
    parameter_list->pseudofreeze = temppseudofreeze;
    
    // Runge_Kutta_order:  must be 1 or 2
    int temprk_order = 1;
    tempinput = util->StringFind4(input_file, "Runge_Kutta_order");
    if (tempinput != "empty")
        istringstream(tempinput) >> temprk_order;
    parameter_list->rk_order = temprk_order;
    
    // reconstruction type
    int tempreconst_type = 0;
    tempinput = util->StringFind4(input_file, "reconst_type");
    if (tempinput != "empty")
        istringstream (tempinput) >> tempreconst_type;
    parameter_list->reconst_type = tempreconst_type;
    

    // Minmod_Theta: theta parameter in the min-mod like limiter
    double tempminmod_theta   = 1.8;
    tempinput = util->StringFind4(input_file, "Minmod_Theta");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempminmod_theta  ;
    parameter_list->minmod_theta = tempminmod_theta;
    
    // Viscosity_Flag_Yes_1_No_0:   set to 0 for ideal hydro
    int tempviscosity_flag = 1;
    tempinput = util->StringFind4(input_file, "Viscosity_Flag_Yes_1_No_0");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempviscosity_flag;
    parameter_list->viscosity_flag = tempviscosity_flag;
    
    // Include_Shear_Visc_Yes_1_No_0
    int tempturn_on_shear = 0;
    tempinput = util->StringFind4(input_file, "Include_Shear_Visc_Yes_1_No_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_shear;
    parameter_list->turn_on_shear = tempturn_on_shear;

    // T_dependent_Shear_to_S_ratio:
    // if 1, ignore constant eta/s
    // and use hard-coded T-dependent shear viscosity
    int tempT_dependent_shear_to_s = 0;
    tempinput = util->StringFind4(input_file, "T_dependent_Shear_to_S_ratio");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempT_dependent_shear_to_s;
    parameter_list->T_dependent_shear_to_s = tempT_dependent_shear_to_s;

    //Shear_to_S_ratio:  constant eta/s
    double tempshear_to_s = 0.08;
    tempinput = util->StringFind4(input_file, "Shear_to_S_ratio");
    if (tempinput != "empty") {
        istringstream(tempinput) >> tempshear_to_s;
    } else if (parameter_list->turn_on_shear == 1
                && parameter_list->T_dependent_shear_to_s == 0) {
        cerr << "please define Shear_to_S_ratio!" << endl;
        exit(1);
    }
    parameter_list->shear_to_s = tempshear_to_s;

    // Include_Bulk_Visc_Yes_1_No_0
    int tempturn_on_bulk = 0;
    tempinput = util->StringFind4(input_file, "Include_Bulk_Visc_Yes_1_No_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_bulk;
    parameter_list->turn_on_bulk = tempturn_on_bulk;
    
    // Include secord order terms
    int tempturn_on_second_order = 0;
    tempinput = util->StringFind4(input_file, "Include_second_order_terms");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_second_order;
    parameter_list->include_second_order_terms = tempturn_on_second_order;
    
    int tempturn_on_diff = 0;
    tempinput = util->StringFind4(input_file, "turn_on_baryon_diffusion");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_diff;
    parameter_list->turn_on_diff = tempturn_on_diff;
    
    // kappa coefficient
    double temp_kappa_coefficient = 0.0;
    tempinput = util->StringFind4(input_file, "kappa_coefficient");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_kappa_coefficient;
    parameter_list->kappa_coefficient = temp_kappa_coefficient;
    
    // Include_deltaf:  
    // Looks like 0 sets delta_f=0, 1 uses standard quadratic ansatz,
    // and 2 is supposed to use p^(2-alpha)
    int tempinclude_deltaf = 1;
    tempinput = util->StringFind4(input_file, "Include_deltaf");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempinclude_deltaf;
    parameter_list->include_deltaf = tempinclude_deltaf;
    
    int tempinclude_deltaf_bulk = 0;
    tempinput = util->StringFind4(input_file, "Include_deltaf_bulk");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempinclude_deltaf_bulk;
    parameter_list->include_deltaf_bulk = tempinclude_deltaf_bulk;
    
    int tempinclude_deltaf_qmu = 0;
    tempinput = util->StringFind4(input_file, "Include_deltaf_qmu");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempinclude_deltaf_qmu;
    parameter_list->include_deltaf_qmu = tempinclude_deltaf_qmu;
    
    int temp_deltaf_14moments = 0;
    tempinput = util->StringFind4(input_file, "deltaf_14moments");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_deltaf_14moments;
    parameter_list->deltaf_14moments = temp_deltaf_14moments;
    
    // Do_FreezeOut_Yes_1_No_0
    // set to 0 to bypass freeze out surface finder
    int tempdoFreezeOut = 1;
    tempinput = util->StringFind4(input_file, "Do_FreezeOut_Yes_1_No_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempdoFreezeOut;
    parameter_list->doFreezeOut = tempdoFreezeOut;
    
    int tempdoFreezeOut_lowtemp = 1;
    tempinput = util->StringFind4(input_file, "Do_FreezeOut_lowtemp");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempdoFreezeOut_lowtemp;
    parameter_list->doFreezeOut_lowtemp = tempdoFreezeOut_lowtemp;
    
    // Initial_Distribution_input_filename
    string tempinitName = "initial/initial_ed.dat";
    tempinput = util->StringFind4(input_file,
                                  "Initial_Distribution_input_filename");
    if (tempinput != "empty")
        tempinitName.assign(tempinput);
    parameter_list->initName.assign(tempinitName);
    
    // Initial_Distribution_Filename for rhob
    string tempinitName_rhob = "initial/initial_rhob.dat";
    tempinput = util->StringFind4(input_file,
                                  "Initial_Rhob_Distribution_Filename");
    if (tempinput != "empty")
        tempinitName_rhob.assign(tempinput);
    parameter_list->initName_rhob.assign(tempinitName_rhob);

    // Initial_Distribution_Filename for ux
    string tempinitName_ux = "initial/initial_ux.dat";
    tempinput = util->StringFind4(input_file,
                                  "Initial_ux_Distribution_Filename");
    if (tempinput != "empty")
        tempinitName_ux.assign(tempinput);
    parameter_list->initName_ux.assign(tempinitName_ux);
    // Initial_Distribution_Filename for uy
    string tempinitName_uy = "initial/initial_uy.dat";
    tempinput = util->StringFind4(input_file,
                                  "Initial_uy_Distribution_Filename");
    if (tempinput != "empty")
        tempinitName_uy.assign(tempinput);
    parameter_list->initName_uy.assign(tempinitName_uy);
    // Initial_Distribution_Filename for TA
    string tempinitName_TA = "initial/initial_TA.dat";
    tempinput = util->StringFind4(input_file,
                                  "Initial_TA_Distribution_Filename");
    if (tempinput != "empty")
        tempinitName_TA.assign(tempinput);
    parameter_list->initName_TA.assign(tempinitName_TA);
    // Initial_Distribution_Filename for TB
    string tempinitName_TB = "initial/initial_TB.dat";
    tempinput = util->StringFind4(input_file,
                                  "Initial_TB_Distribution_Filename");
    if (tempinput != "empty")
        tempinitName_TB.assign(tempinput);
    parameter_list->initName_TB.assign(tempinitName_TB);
    // Initial_Distribution_Filename for rhob TA
    string tempinitName_rhob_TA = "initial/initial_rhob_TA.dat";
    tempinput = util->StringFind4(input_file,
                                  "Initial_rhob_TA_Distribution_Filename");
    if (tempinput != "empty")
        tempinitName_rhob_TA.assign(tempinput);
    parameter_list->initName_rhob_TA.assign(tempinitName_rhob_TA);
    // Initial_Distribution_Filename for rhob TB
    string tempinitName_rhob_TB = "initial/initial_TB.dat";
    tempinput = util->StringFind4(input_file,
                                  "Initial_rhob_TB_Distribution_Filename");
    if (tempinput != "empty")
        tempinitName_rhob_TB.assign(tempinput);
    parameter_list->initName_rhob_TB.assign(tempinitName_rhob_TB);
    
    // Initial_Distribution_AMPT_filename for AMPT
    string tempinitName_AMPT = "initial/initial_AMPT.dat";
    tempinput = util->StringFind4(input_file,
                                  "Initial_Distribution_AMPT_filename");
    if (tempinput != "empty")
        tempinitName_AMPT.assign(tempinput);
    parameter_list->initName_AMPT.assign(tempinitName_AMPT);
    
    // compute beam rapidity according to the collision energy
    double temp_ecm = 200;
    tempinput = util->StringFind4(input_file, "ecm");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_ecm;
    parameter_list->ecm = temp_ecm;
    double y_beam = atanh(sqrt(1. - 1./pow(temp_ecm/2., 2.)));
    parameter_list->beam_rapidity = y_beam;
    
    // initialize the metric, mostly plus
    parameter_list->gmunu = util->mtx_malloc(4, 4);
    for (m=0; m<4; m++) {
        for (n=0; n<4; n++) {
            if (m == n)
                (parameter_list->gmunu)[m][n] = 1.0;
            else
                (parameter_list->gmunu)[m][n] = 0.0;
            if (m==0 && n==0)
                (parameter_list->gmunu)[m][n] *= -1.0;
        }
    }  /* m */
    
    int tempoutputBinaryEvolution = 0;
    tempinput = util->StringFind4(input_file, "outputBinaryEvolution");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutputBinaryEvolution;
    parameter_list->outputBinaryEvolution = tempoutputBinaryEvolution;

    //  Make MUSIC output additionnal hydro information
    //  0 for false (do not output), 1 for true
    bool tempoutput_hydro_debug_info = false;
    tempinput = util->StringFind4(input_file, "output_hydro_debug_info");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutput_hydro_debug_info;
    parameter_list->output_hydro_debug_info = tempoutput_hydro_debug_info;

    // The evolution is outputted every
    // "output_evolution_every_N_timesteps" timesteps
    int temp_evo_N_tau = 1;
    tempinput = util->StringFind4(input_file,
                                  "output_evolution_every_N_timesteps");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_evo_N_tau;
    parameter_list->output_evolution_every_N_timesteps = temp_evo_N_tau;
    
    int temp_evo_N_x = 1;
    tempinput = util->StringFind4(input_file, "output_evolution_every_N_x");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_evo_N_x;
    parameter_list->output_evolution_every_N_x = temp_evo_N_x;
    parameter_list->output_evolution_every_N_y = temp_evo_N_x;
    
    int temp_evo_N_eta = 1;
    tempinput = util->StringFind4(input_file, "output_evolution_every_N_eta");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_evo_N_eta;
    parameter_list->output_evolution_every_N_eta = temp_evo_N_eta;
    
    double temp_evo_T_cut = 0.105;  // GeV
    tempinput = util->StringFind4(input_file, "output_evolution_T_cut");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_evo_T_cut;
    parameter_list->output_evolution_T_cut = temp_evo_T_cut;
    
    // Make MUSIC output a C header input_file containing
    // informations about the hydro parameters used
    // 0 for false (do not output), 1 for true
    bool tempoutput_hydro_params_header = false;
    tempinput = util->StringFind4(input_file, "output_hydro_params_header");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutput_hydro_params_header;
    parameter_list->output_hydro_params_header = tempoutput_hydro_params_header;

    // initial parameters for mode 14
    double temp_dNdy_y_min = -0.5;
    tempinput = util->StringFind4(input_file, "dNdy_y_min");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdy_y_min;
    parameter_list->dNdy_y_min = temp_dNdy_y_min;
    
    double temp_dNdy_y_max = 0.5;
    tempinput = util->StringFind4(input_file, "dNdy_y_max");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdy_y_max;
    parameter_list->dNdy_y_max = temp_dNdy_y_max;
    
    double temp_dNdy_eta_min = -2.0;
    tempinput = util->StringFind4(input_file, "dNdy_eta_min");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdy_eta_min;
    parameter_list->dNdy_eta_min = temp_dNdy_eta_min;
    
    double temp_dNdy_eta_max = 2.0;
    tempinput = util->StringFind4(input_file, "dNdy_eta_max");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdy_eta_max;
    parameter_list->dNdy_eta_max = temp_dNdy_eta_max;
    
    int temp_dNdy_nrap = 30;
    tempinput = util->StringFind4(input_file, "dNdy_nrap");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdy_nrap;
    parameter_list->dNdy_nrap = temp_dNdy_nrap;
    
    double temp_dNdyptdpt_y_min = -0.5;
    tempinput = util->StringFind4(input_file, "dNdyptdpt_y_min");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdyptdpt_y_min;
    parameter_list->dNdyptdpt_y_min = temp_dNdyptdpt_y_min;
    
    double temp_dNdyptdpt_y_max = 0.5;
    tempinput = util->StringFind4(input_file, "dNdyptdpt_y_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdyptdpt_y_max;
    parameter_list->dNdyptdpt_y_max = temp_dNdyptdpt_y_max;
    
    double temp_dNdyptdpt_eta_min = -0.5;
    tempinput = util->StringFind4(input_file, "dNdyptdpt_eta_min");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdyptdpt_eta_min;
    parameter_list->dNdyptdpt_eta_min = temp_dNdyptdpt_eta_min;
    
    double temp_dNdyptdpt_eta_max = 0.5;
    tempinput = util->StringFind4(input_file, "dNdyptdpt_eta_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdyptdpt_eta_max;
    parameter_list->dNdyptdpt_eta_max = temp_dNdyptdpt_eta_max;

    cout << "Done read_in_parameters." << endl;
    check_parameters(parameter_list);
}

void ReadInParameters::check_parameters(InitData *parameter_list) {
    cout << "Checking input parameter list ... " << endl;

    if (parameter_list->Initial_profile < 0) {
        cerr << "Initial profile" << parameter_list->Initial_profile
             << "not defined\n";
        exit(1);
    }

    if (parameter_list->initial_eta_profile > 2
            || parameter_list->initial_eta_profile < 0) {
        cerr << "Initial eta profile" << parameter_list->initial_eta_profile
             << "not defined\n";
        exit(1);
    }

    if (parameter_list->initializeEntropy > 1
            || parameter_list->initializeEntropy < 0) {
        cerr << "Must initialize with entropy (initialize_with_entropy=1) "
              << "or energy (0)\n";
        exit(1);
    }

    if (parameter_list->TFO < 0.0 || parameter_list->TFO > 0.2) {
        cerr << "T_freeze = " << parameter_list->TFO << " is not physical!"
             << endl;
        exit(1);
    }
    
    if (parameter_list->epsilonFreeze <= 0) {
        cerr << "Freeze out energy density must be greater than zero\n";
        exit(1);
    }
    
    if (parameter_list->useEpsFO > 1 || parameter_list->useEpsFO < 0) {
        cerr << "Error: did not set either freeze out energy density "
             << "or temperature, or invalid option for use_eps_for_freeze_out:"
             << parameter_list->useEpsFO << endl;
        exit(1);
    }

    if (parameter_list->whichEOS > 12 || parameter_list->whichEOS < 0) {
        cerr << "EOS_to_use unspecified or invalid option: "
             << parameter_list->whichEOS << endl;
        exit(1);
    }

    if (parameter_list->NumberOfParticlesToInclude > 320) {
        cerr << "Invalid option for number_of_particles_to_include:"
             << parameter_list->NumberOfParticlesToInclude << endl;
        exit(1);
    }
    
    if (parameter_list->freezeOutMethod > 4) {
        cerr << "Invalid option for freeze_out_method: "
             << parameter_list->freezeOutMethod << endl;
        exit(1);
    }

    if (parameter_list->facTau <= 0) {
        cerr << "average_surface_over_this_many_time_steps <= 0: "
             << parameter_list->facTau << endl;
        exit(1);
    }

    double freeze_dtau = parameter_list->facTau*parameter_list->delta_tau;
    if (freeze_dtau > 1.) {
        cout << "[Warning]: freeze-out time setp is too large! "
             << "freeze_dtau = " << freeze_dtau
             << ", hydro_dtau = " << parameter_list->delta_tau
             << ", average_surface_over_this_many_time_steps = "
             << parameter_list->facTau << endl;
    }
    
    if (parameter_list->fac_x <= 0) {
        cerr << "freeze out every x step <= 0: "
             << parameter_list->fac_x << endl;
        exit(1);
    }
    
    if (parameter_list->fac_eta <= 0) {
        cerr << "freeze out every eta step <= 0: "
             << parameter_list->fac_eta << endl;
        exit(1);
    }

    if (parameter_list->nx != parameter_list->ny) {
        cout << "[Warning]: Grid size in x is not equal to grid size in y!"
             << endl;
    }

    if (parameter_list->neta < 32 && !parameter_list->boost_invariant) {
        cerr << "Grid size in eta = " << parameter_list->neta 
             << "is too small for a (3+1)-d run! "
             << "Please increase Grid_size_in_eta to be larger than 32 "
             << "at least!" << endl;
        exit(1);
    }

    if (parameter_list->boost_invariant && parameter_list->neta > 1) {
        cout << "[Warning]: Grid size in eta is set to "
             << parameter_list->neta << " for a (2+1)-d simulation! "
             << "This is redundant! Reset neta to 1!"
             << endl;
        parameter_list->neta = 1;
    }
    
    if (parameter_list->delta_tau > 0.1) {
        cout << "Warning: Delta_Tau = " << parameter_list->delta_tau
             << " maybe too large! "
             << "Please choose a dtau < 0.1 fm." << endl;

        bool reset_dtau_use_CFL_condition = true;
        int temp_CFL_condition = 1;
        string tempinput = util->StringFind4(input_file,
                                      "reset_dtau_use_CFL_condition");
        if (tempinput != "empty")
            istringstream(tempinput) >> temp_CFL_condition;
        if (temp_CFL_condition == 0)
            reset_dtau_use_CFL_condition = false;

        if (reset_dtau_use_CFL_condition) {
            cout << "reset dtau using CFL condition." << endl;
            double dtau_CFL = mini(
                    parameter_list->delta_x/10.0,
                    parameter_list->tau0*parameter_list->delta_eta/10.0);
            parameter_list->delta_tau = dtau_CFL;
            parameter_list->nt = static_cast<int>(
                parameter_list->tau_size/(parameter_list->delta_tau) + 0.5);
            cout << "read_in_parameters: Time step size = "
                 << parameter_list->delta_tau << endl;
            cout << "read_in_parameters: Number of time steps required = "
                 << parameter_list->nt << endl;
        } else {
            exit(1);
        }
    }

    if (parameter_list->min_pt > parameter_list->max_pt) {
        cerr << "min_pt = " << parameter_list->min_pt << " > "
             << "max_pt = " << parameter_list->max_pt << endl;
        exit(1);
    }

    if (parameter_list->phi_steps < 20) {
        cerr << "phi_steps = " << parameter_list->phi_steps
             << " is too small for computing vn, please increase to at least "
             << " 40!" << endl;
        exit(1);
    }
    
    if (parameter_list->rk_order > 2 || parameter_list->rk_order < 0) {
        cerr << "Invalid option for Runge_Kutta_order: "
             << parameter_list->rk_order << endl;
        exit(1);
    }
    if (parameter_list->rk_order != 2) {
        cout << "Runge-Kutta order = " << parameter_list->rk_order << endl;
    }
    
    if (parameter_list->reconst_type > 2
            && parameter_list->reconst_type != 99) {
        cerr << "unrecognized reconst_type: "
             << parameter_list->reconst_type << endl;
        exit(1);
    }

    if (parameter_list->minmod_theta < 1.
            || parameter_list->minmod_theta > 2.) {
        cerr << "minmod = " << parameter_list->minmod_theta << " is out of "
             << "allowed range [1., 2]" << endl;
        exit(1);
    }

    if (parameter_list->turn_on_shear == 0 && parameter_list->shear_to_s > 0) {
        cout << "[Warning]: non-zero eta/s = " << parameter_list->shear_to_s
             << " is set with "
             << "Include_Shear_Visc = " << parameter_list->turn_on_shear
             << ". Please check you want to run ideal hydro!" << endl;
    }

    if (parameter_list->turn_on_shear == 0
            && parameter_list->include_deltaf == 1) {
        cout << "[Warning]: hydro with zero shear viscosity does not need "
             << "shear delta f in Cooper-Fyre!" << endl;
        cout << "input Include_deltaf = " << parameter_list->include_deltaf
             << ". Now rewrite it to 0!" << endl;
        parameter_list->include_deltaf = 0;
    }

    if (parameter_list->turn_on_bulk == 0
            && parameter_list->include_deltaf_bulk == 1) {
        cout << "[Warning]: hydro with zero bulk viscosity does not need "
             << "bulk delta f in Cooper-Fyre!" << endl;
        cout << "input Include_deltaf_bulk = "
             << parameter_list->include_deltaf_bulk
             << ". Now rewrite it to 0!" << endl;
        parameter_list->include_deltaf_bulk = 0;
    }

    if (parameter_list->output_evolution_every_N_timesteps <= 0) {
        cerr << "output_evolution_every_N_timesteps < 0!" << endl;
        exit(1);
    }
    
    if (parameter_list->output_evolution_every_N_x <= 0) {
        cerr << "output_evolution_every_N_x < 0!" << endl;
        exit(1);
    }
    
    if (parameter_list->output_evolution_every_N_y <= 0) {
        cerr << "output_evolution_every_N_y < 0!" << endl;
        exit(1);
    }
    
    if (parameter_list->output_evolution_every_N_eta <= 0) {
        cerr << "output_evolution_every_N_eta < 0!" << endl;
        exit(1);
    }

    if (parameter_list->dNdy_y_min > parameter_list->dNdy_y_max) {
        cerr << "dNdy_y_min = " << parameter_list->dNdy_y_min << " < " 
             << "dNdy_y_max = " << parameter_list->dNdy_y_max << "!" << endl;
        exit(1);
    }
    
    if (parameter_list->dNdy_eta_min > parameter_list->dNdy_eta_max) {
        cerr << "dNdy_eta_min = " << parameter_list->dNdy_eta_min << " < " 
             << "dNdy_eta_max = " << parameter_list->dNdy_eta_max << "!"
             << endl;
        exit(1);
    }
    
    if (parameter_list->dNdyptdpt_y_min > parameter_list->dNdyptdpt_y_max) {
        cerr << "dNdyptdpt_y_min = " << parameter_list->dNdyptdpt_y_min
             << " < " 
             << "dNdyptdpt_y_max = " << parameter_list->dNdyptdpt_y_max
             << "!" << endl;
        exit(1);
    }
    
    if (parameter_list->dNdyptdpt_eta_min
            > parameter_list->dNdyptdpt_eta_max) {
        cerr << "dNdyptdpt_eta_min = " << parameter_list->dNdyptdpt_eta_min
             << " < " 
             << "dNdyptdpt_eta_max = " << parameter_list->dNdyptdpt_eta_max
             << "!" << endl;
        exit(1);
    }

    cout << "Finished checking input parameter list. "
         << "Everything looks reasonable so far " << emoji_face->success()
         << endl;
}
