// Original copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
// Massively cleaned up and improved by Chun Shen 2015-2016
#include <stdio.h>
#include <sys/stat.h>

#include "./util.h"
#include "./grid.h"
#include "./data.h"
#include "./init.h"
#include "./eos.h"
#include "./evolve.h"

using namespace std;

void ReadInData3(InitData *DATA, string file);

// main program
int main(int argc, char *argv[]) {
    string input_file;
    static InitData DATA;
    // you have the option to give a second command line option,
    // which is an integer to be added to the random seed from the current time
    // because on the cluster it will happen that two jobs start at exactly
    // the same time, this makes sure that they dont run with excatly the same
    // seed
    string sseed;
    int seed = 0;
    if (argc > 2) {
        sseed = argv[2];
        seed = atoi(sseed.c_str());
    }
    seed *= 10000;
    DATA.seed = seed;

    if (argc > 1)
        input_file = *(argv+1);
    else
        input_file = "";

    ReadInData3(&DATA, input_file);
    //cout << sizeof(Grid *) << endl;
    //exit(1);

    EOS *eos = new EOS(&DATA);

    Grid ***arena;
    if (DATA.mode == 1 || DATA.mode == 2) {
        // clean all the surface files
        system("rm surface.dat surface?.dat surface??.dat 2> /dev/null");
        Init *init = new Init(eos);
        init->InitArena(&DATA, &arena);
        
        Evolve *evolve = new Evolve(eos, &DATA);
        evolve->EvolveIt(&DATA, arena);
    }
}  /* main */

void ReadInData3(InitData *DATA, string file) {
    // this function reads in parameters for MUSIC
    // this is an improved version that supports comments
    // in the input parameter file
    int m, n;
    Util *util;
    util = new Util();
    string tempinput;

    // Impact_parameter: in fm, >=0
    double tempb = 0.0;
    tempinput = util->StringFind4(file, "Impact_parameter");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempb;
    DATA->b = tempb;
    if (DATA->b < 0)  {
        cerr << "Impact parameter must be greater than zero\n";
        exit(1);
    }

    // echo_level controls the mount of warning message output
    // during the evolution
    double temp_echo_level = 9;
    tempinput = util->StringFind4(file, "echo_level");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_echo_level;
    DATA->echo_level = temp_echo_level;

    // Target, Projectile:  any name as defined in known_nuclei.dat
    string tempTarget = "Pb";
    tempinput = util->StringFind4(file, "Target");
    if (tempinput != "empty")
        tempTarget.assign(tempinput);
    DATA->Target.assign(tempTarget);
    string tempProjectile = "Pb";
    tempinput = util->StringFind4(file, "Projectile");
    if(tempinput != "empty")
        tempProjectile.assign(tempinput);
    DATA->Projectile.assign(tempProjectile);
  
    // SigmaNN: nucleon-nucleon cross section in mb
    double tempSigmaNN = 60.;
    tempinput = util->StringFind4(file, "SigmaNN");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempSigmaNN;
    DATA->SigmaNN = tempSigmaNN;
    if (DATA->SigmaNN < 0) {
        cerr << "NN cross section must be greater than zero (in mb)\n";
        exit(1);
    }
  
    // Initial_profile: 
    //   0: Gubser flow test
    //   8: read in from file for 2d IP-Glasma initial conditions
    //  11: Read in initial profiles for energy density and rhob
    //      in the transverse plane from files
    int tempInitial_profile = 1;
    tempinput = util->StringFind4(file, "Initial_profile");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempInitial_profile;
    DATA->Initial_profile = tempInitial_profile;
    if (DATA->Initial_profile < 0) {
        cerr << "Initial profile " << DATA->Initial_profile
             << " not defined\n";
        exit(1);
    }

    // Select the profile to use in eta for the energy/entropy initialisation
    // 1 for Hirano's central plateau + Gaussian decay
    // 2 for a Woods-Saxon profile
    int tempinitial_eta_profile = 1;
    tempinput = util->StringFind4(file, "initial_eta_profile");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempinitial_eta_profile;
    DATA->initial_eta_profile = tempinitial_eta_profile;

    // initialize_with_entropy:
    // 0: scale with energy density
    // 1: scale with entropy density
    int tempinitializeEntropy = 0;
    tempinput = util->StringFind4(file, "initialize_with_entropy");
    if(tempinput != "empty")
        istringstream(tempinput) >> tempinitializeEntropy;
    DATA->initializeEntropy = tempinitializeEntropy;
    if (DATA->initializeEntropy > 1 || DATA->initializeEntropy < 0) {
        cerr << "Must initialize with entropy (initialize_with_entropy=1)"
             << "or energy (0)\n";
        exit(1);
    }

    // Must set either freeze out energy density or temperature,
    // otherwise generate error message below.
    DATA->useEpsFO = 2;
    // T_freeze:  freeze out temperature
    // only used with use_eps_for_freeze_out = 0
    double tempTFO = 0.12;
    tempinput = util->StringFind4(file, "T_freeze");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempTFO;
    DATA->TFO = tempTFO;
    // if only freeze out temperature is set, freeze out by temperature
    if (tempinput != "empty")
        DATA->useEpsFO = 0;
    int tfoset = 0;
    if (tempinput != "empty")
        tfoset = 1;
    // epsilon_freeze: freeze-out energy density in GeV/fm^3
    // only used with use_eps_for_freeze_out = 1
    double tempepsilonFreeze = 0.12;
    tempinput = util->StringFind4(file, "epsilon_freeze");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempepsilonFreeze;
    DATA->epsilonFreeze = tempepsilonFreeze;
    if (DATA->epsilonFreeze <= 0) {
        cerr << "Freeze out energy density must be greater than zero\n";
        exit(1);
    }
    if (tempinput != "empty")
        DATA->useEpsFO = 1;  // if epsilon_freeze is set, freeze out by epsilon

    // use_eps_for_freeze_out:
    // 0: freeze out at constant temperature T_freeze
    // 1: freeze out at constant energy density epsilon_freeze
    // if set in input file, overide above defaults
    int tempuseEpsFO = DATA->useEpsFO;
    tempinput = util->StringFind4(file, "use_eps_for_freeze_out");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempuseEpsFO;
    DATA->useEpsFO = tempuseEpsFO;
    if (DATA->useEpsFO > 1 || DATA->useEpsFO < 0) {
        cerr << "Error: did not set either freeze out energy density "
             << "or temperature, or invalid option for use_eps_for_freeze_out:"
             << DATA->useEpsFO << endl;
        exit(1);
    }
    if (tfoset == 1 && DATA->useEpsFO == 1)
        cerr << "T_freeze set but overridden -- "
             << "freezing out by energy density at "
             << DATA->epsilonFreeze << " GeV/fm^3\n";

    //particle_spectrum_to_compute:
    // 0: Do all up to number_of_particles_to_include
    // any natural number: Do the particle with this (internal) ID
    int tempparticleSpectrumNumber = 0;
    tempinput = util->StringFind4(file, "particle_spectrum_to_compute");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempparticleSpectrumNumber;
    DATA->particleSpectrumNumber = tempparticleSpectrumNumber;
  
    // mode: 
    // 1: Does everything. Evolution. Computation of thermal spectra.
    //    Resonance decays. Observables.
    //    Only compatible with freeze_out_method=3 and pseudofreeze=1
    // 2: Evolution only.
    // 3: Compute all thermal spectra only.
    // 4: Resonance decays only.
    // 5: Resonance decays for just one specific particle
    //    (only for testing - this will miss the complete chain of decays).
    // 6: only combine charged hadrons 
    //    - can be used for any postprocessing with the stored results
    // 13: Compute observables from previously-computed thermal spectra
    // 14: Compute observables from post-decay spectra
    int tempmode = 1;
    tempinput = util->StringFind4(file, "mode");
    if (tempinput != "empty") {
        istringstream ( tempinput ) >> tempmode;
    } else {
      cerr << "Must specify mode. Exiting.\n";
      exit(1);
    }
    DATA->mode = tempmode;
    
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
    // 10: lattice EOS at finite mu_B from A. Monnai
    int tempwhichEOS = 9;
    tempinput = util->StringFind4(file, "EOS_to_use");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempwhichEOS;
    DATA->whichEOS = tempwhichEOS;
    if(DATA->whichEOS>20 || DATA->whichEOS<0) 
    {
      cerr << "EOS_to_use unspecified or invalid option: "
           << DATA->whichEOS << endl;
      exit(1);
    }
  
    // number_of_particles_to_include:
    // This determines up to which particle in the list spectra 
    // should be computed (mode=3) or resonances should be included (mode=4)
    // current maximum = 319
    int tempNumberOfParticlesToInclude = 2;
    tempinput = util->StringFind4(file, "number_of_particles_to_include");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempNumberOfParticlesToInclude;
    DATA->NumberOfParticlesToInclude = tempNumberOfParticlesToInclude;
    if (DATA->NumberOfParticlesToInclude > 320) {
      cerr << "Invalid option for number_of_particles_to_include:"
           << DATA->NumberOfParticlesToInclude << endl;
      exit(1);
    }
    
    // freeze_out_method:
    // 1: Hirano's simplified method
    // 2: Schenke's more complex method
    // 3: Luzum's simple method
    // 4: Cornelius 
    int tempfreezeOutMethod = 4;
    tempinput = util->StringFind4(file, "freeze_out_method");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempfreezeOutMethod;
    DATA->freezeOutMethod = tempfreezeOutMethod;
    if (DATA->freezeOutMethod > 4) {
      cerr << "Invalid option for freeze_out_method:"
           << DATA->freezeOutMethod << endl;
      exit(1);
    }

    int temp_freeze_eps_flag = 0;
    tempinput = util->StringFind4(file, "freeze_eps_flag");
    if (tempinput != "empty") istringstream (tempinput) >> temp_freeze_eps_flag;
    DATA->freeze_eps_flag = temp_freeze_eps_flag;
    
    string temp_freeze_list_filename = "eps_freeze_list_s95p_v1.dat";
    tempinput = util->StringFind4(file, "freeze_list_filename");
    if (tempinput != "empty") temp_freeze_list_filename.assign(tempinput);
    DATA->freeze_list_filename.assign(temp_freeze_list_filename);

    int temp_N_freeze_out = 1;
    tempinput = util->StringFind4(file, "N_freeze_out");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_N_freeze_out;
    DATA->N_freeze_out = temp_N_freeze_out;
    
    double temp_eps_freeze_max = 0.18;
    tempinput = util->StringFind4(file, "eps_freeze_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_eps_freeze_max;
    DATA->eps_freeze_max = temp_eps_freeze_max;
    
    double temp_eps_freeze_min = 0.18;
    tempinput = util->StringFind4(file, "eps_freeze_min");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_eps_freeze_min;
    DATA->eps_freeze_min = temp_eps_freeze_min;
    
    // average_surface_over_this_many_time_steps:
    // Only save every N timesteps for finding freeze out surface
    int tempfacTau = 10;
    tempinput = util->StringFind4(file,
                                  "average_surface_over_this_many_time_steps");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempfacTau;
    DATA->facTau = tempfacTau;
    
    int tempfac_x = 1;
    tempinput = util->StringFind4(file, "Ncell_skip_x");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempfac_x;
    DATA->fac_x = tempfac_x;
    
    int tempfac_y = 1;
    tempinput = util->StringFind4(file, "Ncell_skip_y");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempfac_y;
    DATA->fac_y = tempfac_y;
    
    
    //  Grid_size_in_*
    // number of cells in x,y direction
    int tempnx = 10;
    tempinput = util->StringFind4(file, "Grid_size_in_x");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempnx;
    DATA->nx = tempnx;
    int tempny = 10;
    tempinput = util->StringFind4(file, "Grid_size_in_y");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempny;
    DATA->ny = tempny;
    
    // Grid_size_in_eta
    // number of cells in eta direction.
    // Must have at least 4 cells per processor.
    // Must be an even number.
    // One cell is positioned at eta=0, 
    // half the cells are at negative eta,
    // the rest (one fewer) are at positive eta
    int tempneta = 4;
    tempinput = util->StringFind4(file, "Grid_size_in_eta");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempneta;
    DATA->neta = tempneta;
  
    // *_grid_size_in_fm: 
    // total length of box in x,y direction in fm (minus delta_*)
    double tempx_size = 25.;
    tempinput = util->StringFind4(file, "X_grid_size_in_fm");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempx_size;
    DATA->x_size = tempx_size;
    double tempy_size = 25.;
    tempinput = util->StringFind4(file, "Y_grid_size_in_fm");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempy_size;
    DATA->y_size = tempy_size;
    
    // Eta_grid_size:  total length of box in eta direction (minus delta_eta)
    // e.g., neta=8 and eta_size=8 has 8 cells that run from eta=-4 to eta=3
    double tempeta_size = 8.;
    tempinput = util->StringFind4(file, "Eta_grid_size");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempeta_size;
    DATA->eta_size = tempeta_size;
    
    // Total_evolution_time_tau
    // total evolution time in [fm]. 
    // in case of freeze_out_method = 2,3, 
    // evolution will halt earlier if all cells are frozen out.
    double temptau_size = 50.;
    tempinput = util->StringFind4(file, "Total_evolution_time_tau");
    if (tempinput != "empty")
        istringstream(tempinput) >> temptau_size;
    DATA->tau_size = temptau_size;
    
    // Initial_time_tau_0:  in fm
    double temptau0 = 0.4;
    tempinput = util->StringFind4(file, "Initial_time_tau_0");
    if(tempinput != "empty")
        istringstream(tempinput) >> temptau0;
    DATA->tau0 = temptau0;
    
    /* x-grid, for instance, runs from 0 to nx */
    DATA->delta_x = DATA->x_size/(static_cast<double>(DATA->nx) - 1.);
    DATA->delta_y = DATA->y_size/(static_cast<double>(DATA->ny) - 1.);
    DATA->delta_eta = DATA->eta_size/(static_cast<double>(DATA->neta));
    
    cout << "DeltaX = " << DATA->delta_x << " fm." << endl;
    cout << "DeltaY = " << DATA->delta_y << " fm." << endl;
    cout << "DeltaETA = " << DATA->delta_eta << endl;
    
    // Delta_Tau: 
    // time step to use in [fm].
    // If a too large value is given, it will automatically be reduced to the 
    // maximal acceptable value according to the CFL condition.
    /* CFL condition : delta_tau < min(delta_x/10, tau0 delta_eta/10) */
    int tempUseCFL_condition = 0;
    tempinput = util->StringFind4(file, "UseCFL_condition");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempUseCFL_condition;

    double tempdelta_tau = 8.;
    tempinput = util->StringFind4(file, "Delta_Tau");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempdelta_tau;
    DATA->delta_tau = tempdelta_tau;
 
    if (tempUseCFL_condition == 1) {
        double tempf = mini(DATA->delta_x/10.0,
                            DATA->tau0*(DATA->delta_eta/10.0));
        if (tempf < DATA->delta_tau)
            DATA->delta_tau = tempf;
    }
    cout << " DeltaTau = " << DATA->delta_tau << " fm." << endl;
    
    DATA->nt = (int) (floor(DATA->tau_size/(DATA->delta_tau) + 0.5));
    cout << "ReadInData: Time step size = " << DATA->delta_tau << endl;
    cout << "ReadInData: Number of time steps required = " << DATA->nt << endl;
    
    // output_evolution_data:  
    // 1:  output bulk information at every grid point at every time step
    int tempoutputEvolutionData = 0;
    tempinput = util->StringFind4(file, "output_evolution_data");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutputEvolutionData;
    DATA->outputEvolutionData = tempoutputEvolutionData;
    
    // Eta_fall_off:
    // width of half-Gaussian on each side of a central pleateau in eta
    double tempeta_fall_off  = 0.4;
    tempinput = util->StringFind4(file, "Eta_fall_off");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempeta_fall_off ;
    DATA->eta_fall_off  = tempeta_fall_off;
  
    // Eta_plateau_size:
    // width of the flat region symmetrical around eta=0
    double tempeta_flat   = 5.9;
    tempinput = util->StringFind4(file, "Eta_plateau_size");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempeta_flat  ;
    DATA->eta_flat = tempeta_flat;

    // eta envelope function parameter for rhob
    int temp_rhob_flag = 1;
    tempinput = util->StringFind4(file, "initial_eta_rhob_profile");
    if (tempinput != "empty") istringstream(tempinput) >> temp_rhob_flag;
    DATA->initial_eta_rhob_profile = temp_rhob_flag;
    double temp_eta_0 = 3.0;
    tempinput = util->StringFind4(file, "eta_rhob_0");
    if (tempinput != "empty") istringstream (tempinput) >> temp_eta_0;
    DATA->eta_rhob_0 = temp_eta_0;
    double temp_eta_width = 1.0;
    tempinput = util->StringFind4(file, "eta_rhob_width");
    if (tempinput != "empty") istringstream (tempinput) >> temp_eta_width;
    DATA->eta_rhob_width = temp_eta_width;
    double temp_eta_plateau_height = 0.5;
    tempinput = util->StringFind4(file, "eta_rhob_plateau_height");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_eta_plateau_height;
    DATA->eta_rhob_plateau_height = temp_eta_plateau_height;
    double temp_eta_width_1 = 1.0;
    tempinput = util->StringFind4(file, "eta_rhob_width_1");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_eta_width_1;
    DATA->eta_rhob_width_1 = temp_eta_width_1;
    double temp_eta_width_2 = 1.0;
    tempinput = util->StringFind4(file, "eta_rhob_width_2");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_eta_width_2;
    DATA->eta_rhob_width_2 = temp_eta_width_2;

  // s_factor
  double tempsFactor = 20.;
    tempinput = util->StringFind4(file, "s_factor");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempsFactor  ;
    DATA->sFactor = tempsFactor;
    
    // for calculation of spectra: maximal_rapidity:
    // spectra calculated from zero to this rapidity in +y and -y
    double tempymax   = 4.8;
    tempinput = util->StringFind4(file, "maximal_rapidity");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempymax;
    DATA->ymax = tempymax;
    
    // delta_y:
    // step size in rapidity in calculation of spectra
    double tempdeltaY   = 0.1;
    tempinput = util->StringFind4(file, "delta_y");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempdeltaY;
    DATA->deltaY = tempdeltaY;
    
    // max_pseudorapidity:
    // spectra calculated from zero to this pseudorapidity in +eta and -eta
    double tempmax_pseudorapidity   = 2.5;
    tempinput = util->StringFind4(file, "max_pseudorapidity");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempmax_pseudorapidity  ;
    DATA->max_pseudorapidity = tempmax_pseudorapidity;
    
     // pseudo_steps:
    // steps in pseudorapidity in calculation of spectra
    int temppseudo_steps = 26;
    tempinput = util->StringFind4(file, "pseudo_steps");
    if (tempinput != "empty")
        istringstream(tempinput) >> temppseudo_steps;
    DATA->pseudo_steps = temppseudo_steps; 
    
    // phi_steps
    // steps in azimuthal angle in calculation of spectra
    int tempphi_steps = 48;
    tempinput = util->StringFind4(file, "phi_steps");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempphi_steps;
    DATA->phi_steps = tempphi_steps; 
  
    // min_pt:
    // spectra calculated from this to max_pt transverse momentum in GeV
    double tempmin_pt = 0.0;
    tempinput = util->StringFind4(file, "min_pt");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempmin_pt  ;
    DATA->min_pt = tempmin_pt;
      
    // max_pt:
    // spectra calculated from min_pt to this transverse momentum in GeV
    double tempmax_pt   = 3.0;
    tempinput = util->StringFind4(file, "max_pt");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempmax_pt  ;
    DATA->max_pt   = tempmax_pt;
    
     // pt_steps:
    // steps in transverse momentum in calculation of spectra
    int temppt_steps   = 60;
    tempinput = util->StringFind4(file, "pt_steps");
    if(tempinput != "empty") istringstream ( tempinput ) >> temppt_steps  ;
    DATA->pt_steps   = temppt_steps;   
    
    // pseudofreeze
    // Calculate spectra at fixed, 
    // equally-spaced grid in pseudorapidity, pt, and phi
    int temppseudofreeze = 1;
    tempinput = util->StringFind4(file, "pseudofreeze");
    if (tempinput != "empty")
        istringstream(tempinput) >> temppseudofreeze;
    DATA->pseudofreeze = temppseudofreeze;
    
    // switch for baryon current propagation
    int tempturn_on_rhob = 0;
    tempinput = util->StringFind4(file, "Include_Rhob_Yes_1_No_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_rhob;
    DATA->turn_on_rhob = tempturn_on_rhob;
    if (DATA->turn_on_rhob == 1)
       DATA->alpha_max = 5;
    else
       DATA->alpha_max = 4;
    int tempturn_on_diff = 0;
    tempinput = util->StringFind4(file, "turn_on_baryon_diffusion");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_diff;
    DATA->turn_on_diff = tempturn_on_diff;
  
    // Runge_Kutta_order:  must be 1 or 2
    int temprk_order = 1;
    tempinput = util->StringFind4(file, "Runge_Kutta_order");
    if (tempinput != "empty")
        istringstream(tempinput) >> temprk_order;
    DATA->rk_order = temprk_order;
    if (DATA->rk_order>2 || DATA->rk_order <0) {
        cerr << "Invalid option for Runge_Kutta_order: "
             << DATA->rk_order << endl;
        exit(1);
    }
    cout << "Runge-Kutta order = " << DATA->rk_order << endl;

    // reconstruction type
    int tempreconst_type = 0;
    tempinput = util->StringFind4(file, "reconst_type");
    if (tempinput != "empty")
        istringstream (tempinput) >> tempreconst_type;
    DATA->reconst_type = tempreconst_type;
    cout << "reconst type = " << DATA->reconst_type << endl;
  
  
    // Minmod_Theta: theta parameter in the min-mod like limiter
    double tempminmod_theta = 1.8;
    tempinput = util->StringFind4(file, "Minmod_Theta");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempminmod_theta  ;
    DATA->minmod_theta = tempminmod_theta;
    cout << "minmod theta = " << DATA->minmod_theta << endl;
  
  
    // Viscosity_Flag_Yes_1_No_0: set to 0 for ideal hydro
    int tempviscosity_flag = 1;
    tempinput = util->StringFind4(file, "Viscosity_Flag_Yes_1_No_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempviscosity_flag;
    DATA->viscosity_flag = tempviscosity_flag;
  
    // Bulk_to_S_ratio:  constant zeta/s
    double tempbulk_to_s = 0.0;
    tempinput = util->StringFind4(file, "Bulk_to_S_ratio");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempbulk_to_s;
    DATA->bulk_to_s = tempbulk_to_s;
  
    // Include_Bulk_Visc_Yes_1_No_0
    int tempturn_on_bulk = 0;
    tempinput = util->StringFind4(file, "Include_Bulk_Visc_Yes_1_No_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_bulk;
    DATA->turn_on_bulk = tempturn_on_bulk;
  
    // Shear_relaxation_time_tau_pi:  in fm?
    double temptau_pi = 0.01;
    tempinput = util->StringFind4(file, "Shear_relaxation_time_tau_pi");
    if (tempinput != "empty")
        istringstream(tempinput) >> temptau_pi  ;
    DATA->tau_pi = temptau_pi;
  
  
  // Bulk_relaxation_time_tau_b_pi:  in fm?
  double temptau_b_pi = 0.6;
  tempinput = util->StringFind4(file, "Bulk_relaxation_time_tau_b_pi");
  if(tempinput != "empty") istringstream ( tempinput ) >> temptau_b_pi  ;
  DATA->tau_b_pi = temptau_b_pi;
  
  //Shear_to_S_ratio:  constant eta/s
  double tempshear_to_s = 0.6;
  tempinput = util->StringFind4(file, "Shear_to_S_ratio");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempshear_to_s  ;
  DATA->shear_to_s   = tempshear_to_s;
    
  // Include_Shear_Visc_Yes_1_No_0
  int tempturn_on_shear = 0;
  tempinput = util->StringFind4(file, "Include_Shear_Visc_Yes_1_No_0");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempturn_on_shear;
  DATA->turn_on_shear = tempturn_on_shear;
  
    // T_dependent_Shear_to_S_ratio:
    // if 1, use hard-coded T-dependent shear viscosity
    int tempT_dependent_shear_to_s = 1;
    tempinput = util->StringFind4(file, "T_dependent_Shear_to_S_ratio");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempT_dependent_shear_to_s;
    DATA->T_dependent_shear_to_s = tempT_dependent_shear_to_s;
  
  // Include_deltaf:  
  // Looks like 0 sets delta_f=0, 1 uses standard quadratic ansatz,
  // and 2 is supposed to use p^(2-alpha)
  int tempinclude_deltaf = 1;
  tempinput = util->StringFind4(file, "Include_deltaf");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempinclude_deltaf;
  DATA->include_deltaf = tempinclude_deltaf;
  
  int tempinclude_deltaf_qmu = 0;
  tempinput = util->StringFind4(file, "Include_deltaf_qmu");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempinclude_deltaf_qmu;
  DATA->include_deltaf_qmu = tempinclude_deltaf_qmu;
  
  int temp_deltaf_14moments = 0;
  tempinput = util->StringFind4(file, "deltaf_14moments");
  if(tempinput != "empty") istringstream ( tempinput ) >> temp_deltaf_14moments;
  DATA->deltaf_14moments = temp_deltaf_14moments;
  
  int tempinclude_deltaf_bulk = 0;
  tempinput = util->StringFind4(file, "Include_deltaf_bulk");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempinclude_deltaf_bulk;
  DATA->include_deltaf_bulk = tempinclude_deltaf_bulk;
  
  // Do_FreezeOut_Yes_1_No_0
  // set to 0 to bypass freeze out surface finder
  int tempdoFreezeOut = 1;
  tempinput = util->StringFind4(file, "Do_FreezeOut_Yes_1_No_0");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempdoFreezeOut;
  DATA->doFreezeOut = tempdoFreezeOut;
  int tempdoFreezeOut_lowtemp = 1;
  tempinput = util->StringFind4(file, "Do_FreezeOut_lowtemp");
  if (tempinput != "empty")
      istringstream(tempinput) >> tempdoFreezeOut_lowtemp;
  DATA->doFreezeOut_lowtemp = tempdoFreezeOut_lowtemp;
  
  // Initial_Distribution_Filename
  string tempinitName = "initial/initial_ed.dat";
  tempinput = util->StringFind4(file, "Initial_Distribution_Filename");
  if(tempinput != "empty") tempinitName.assign(tempinput);
  DATA->initName.assign(tempinitName);
  // Initial_Distribution_Filename for rhob
  string tempinitName_rhob = "initial/initial_rhob.dat";
  tempinput = util->StringFind4(file, "Initial_Rhob_Distribution_Filename");
  if(tempinput != "empty") tempinitName_rhob.assign(tempinput);
  DATA->initName_rhob.assign(tempinitName_rhob);
  // Initial_Distribution_Filename for ux
  string tempinitName_ux = "initial/initial_ux.dat";
  tempinput = util->StringFind4(file, "Initial_ux_Distribution_Filename");
  if(tempinput != "empty") tempinitName_ux.assign(tempinput);
  DATA->initName_ux.assign(tempinitName_ux);
  // Initial_Distribution_Filename for uy
  string tempinitName_uy = "initial/initial_uy.dat";
  tempinput = util->StringFind4(file, "Initial_uy_Distribution_Filename");
  if(tempinput != "empty") tempinitName_uy.assign(tempinput);
  DATA->initName_uy.assign(tempinitName_uy);
  // Initial_Distribution_Filename for TA
  string tempinitName_TA = "initial/initial_TA.dat";
  tempinput = util->StringFind4(file, "Initial_TA_Distribution_Filename");
  if(tempinput != "empty") tempinitName_TA.assign(tempinput);
  DATA->initName_TA.assign(tempinitName_TA);
  // Initial_Distribution_Filename for TB
  string tempinitName_TB = "initial/initial_TB.dat";
  tempinput = util->StringFind4(file, "Initial_TB_Distribution_Filename");
  if(tempinput != "empty") tempinitName_TB.assign(tempinput);
  DATA->initName_TB.assign(tempinitName_TB);
  
    // compute beam rapidity according to the collision energy
    double temp_ecm = 200;
    tempinput = util->StringFind4(file, "ecm");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_ecm;
    DATA->ecm = temp_ecm;
    double y_beam = atanh(sqrt(1. - 1./pow(temp_ecm/2., 2.)));
    DATA->beam_rapidity = y_beam;

    /* initialize the metric, mostly plus */
    DATA->gmunu = util->mtx_malloc(4, 4);
    for (m = 0; m < 4; m++) {
        for (n = 0; n < 4; n++) {
            if (m == n)
                (DATA->gmunu)[m][n] = 1.0;
            else
                (DATA->gmunu)[m][n] = 0.0;
            if (m == 0 && n == 0)
                (DATA->gmunu)[m][n] *= -1.0;
        }
    }  /* m */

    int tempoutputBinaryEvolution = 0;
    tempinput = util->StringFind4(file, "outputBinaryEvolution");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutputBinaryEvolution;
    DATA->outputBinaryEvolution = tempoutputBinaryEvolution;
    // End MARTINI parameters

    // Set to 1 if initial condition is boost-invariant
    // for compatibility with code from MUSIC light fork
    int tempboost_invariant = 0;
    tempinput = util->StringFind4(file, "boost_invariant");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempboost_invariant;
    DATA->boost_invariant = tempboost_invariant;
    DATA->boostInvariant = DATA->boost_invariant;

    // Make MUSIC output additionnal hydro information
    // 0 for false (do not output), 1 for true
    bool tempoutput_hydro_debug_info = true;
    tempinput = util->StringFind4(file, "output_hydro_debug_info");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutput_hydro_debug_info;
    DATA->output_hydro_debug_info = tempoutput_hydro_debug_info;

    // The evolution is outputted every "output_evolution_every_N_timesteps"
    // timesteps Can't be modified from the input file for now, for safety.
    int temp_evo_N_tau = 1;
    tempinput = util->StringFind4(file, "output_evolution_every_N_timesteps");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_evo_N_tau;
    DATA->output_evolution_every_N_timesteps = temp_evo_N_tau;

    int temp_evo_N_x = 1;
    tempinput = util->StringFind4(file, "output_evolution_every_N_x");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_evo_N_x;
    DATA->output_evolution_every_N_x = temp_evo_N_x;

    int temp_evo_N_y = 1;
    tempinput = util->StringFind4(file, "output_evolution_every_N_y");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_evo_N_y;
    DATA->output_evolution_every_N_y = temp_evo_N_y;

    int temp_evo_N_eta = 1;
    tempinput = util->StringFind4(file, "output_evolution_every_N_eta");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_evo_N_eta;
    DATA->output_evolution_every_N_eta = temp_evo_N_eta;

    // Make MUSIC output a C header file containing informations about
    // the hydro parameters used
    // 0 for false (do not output), 1 for true
    bool tempoutput_hydro_params_header = false;
    tempinput = util->StringFind4(file, "output_hydro_params_header");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutput_hydro_params_header;
    DATA->output_hydro_params_header = tempoutput_hydro_params_header;

    // QuestRevert_rho_shear_max: QuestRevert has condition
    // rho_shear > rho_shear_max
    double tempQuestRevert_rho_shear_max = 0.1;
    tempinput = util->StringFind4(file, "QuestRevert_rho_shear_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempQuestRevert_rho_shear_max;
    DATA->QuestRevert_rho_shear_max = tempQuestRevert_rho_shear_max;

    // QuestRevert_rho_bulk_max: QuestRevert has condition
    // rho_bulk > rho_bulk_max
    double tempQuestRevert_rho_bulk_max = 0.1;
    tempinput = util->StringFind4(file, "QuestRevert_rho_bulk_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempQuestRevert_rho_bulk_max;
    DATA->QuestRevert_rho_bulk_max = tempQuestRevert_rho_bulk_max;

    // QuestRevert_rho_q_max: QuestRevert has condition rho_q > rho_q_max
    double tempQuestRevert_rho_q_max = 0.1;
    tempinput = util->StringFind4(file, "QuestRevert_rho_q_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempQuestRevert_rho_q_max;
    DATA->QuestRevert_rho_q_max = tempQuestRevert_rho_q_max;

    // QuestRevert_factor: defines aggressiveness of QuestRevert
    double tempQuestRevert_factor = 0.;
    tempinput = util->StringFind4(file, "QuestRevert_factor");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempQuestRevert_factor;
    DATA->QuestRevert_factor = tempQuestRevert_factor;

    // QuestRevert_prefactor: defines aggressiveness of QuestRevert
    double tempQuestRevert_prefactor = 300.;
    tempinput = util->StringFind4(file, "QuestRevert_prefactor");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempQuestRevert_prefactor;
    DATA->QuestRevert_prefactor = tempQuestRevert_prefactor;

    // QuestRevert_eps_factor:
    // factor \sim tanh ( eps / eps_F0 * QuestRevert_eps_factor)
    double tempQuestRevert_eps_factor = 1.22149;
    tempinput = util->StringFind4(file, "QuestRevert_eps_factor");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempQuestRevert_eps_factor;
    DATA->QuestRevert_eps_factor = tempQuestRevert_eps_factor;

    // QuestRevert_epsilon_min: QuestRevert is applied if epsilon > epsilon_min
    double tempQuestRevert_epsilon_min = 39061.;
    tempinput = util->StringFind4(file, "QuestRevert_epsilon_min");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempQuestRevert_epsilon_min;
    DATA->QuestRevert_epsilon_min = tempQuestRevert_epsilon_min;

    // initial parameters for mode 14
    double temp_dNdy_y_min = -0.5;
    tempinput = util->StringFind4(file, "dNdy_y_min");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdy_y_min;
    DATA->dNdy_y_min = temp_dNdy_y_min;

    double temp_dNdy_y_max = 0.5;
    tempinput = util->StringFind4(file, "dNdy_y_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdy_y_max;
    DATA->dNdy_y_max = temp_dNdy_y_max;

    double temp_dNdy_eta_min = -2.0;
    tempinput = util->StringFind4(file, "dNdy_eta_min");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdy_eta_min;
    DATA->dNdy_eta_min = temp_dNdy_eta_min;
    double temp_dNdy_eta_max = 2.0;
    tempinput = util->StringFind4(file, "dNdy_eta_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdy_eta_max;
    DATA->dNdy_eta_max = temp_dNdy_eta_max;

    int temp_dNdy_nrap = 30;
    tempinput = util->StringFind4(file, "dNdy_nrap");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdy_nrap;
    DATA->dNdy_nrap = temp_dNdy_nrap;

    double temp_dNdyptdpt_y_min = -0.5;
    tempinput = util->StringFind4(file, "dNdyptdpt_y_min");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdyptdpt_y_min;
    DATA->dNdyptdpt_y_min = temp_dNdyptdpt_y_min;

    double temp_dNdyptdpt_y_max = 0.5;
    tempinput = util->StringFind4(file, "dNdyptdpt_y_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdyptdpt_y_max;
    DATA->dNdyptdpt_y_max = temp_dNdyptdpt_y_max;

    double temp_dNdyptdpt_eta_min = -0.5;
    tempinput = util->StringFind4(file, "dNdyptdpt_eta_min");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdyptdpt_eta_min;
    DATA->dNdyptdpt_eta_min = temp_dNdyptdpt_eta_min;

    double temp_dNdyptdpt_eta_max = 0.5;
    tempinput = util->StringFind4(file, "dNdyptdpt_eta_max");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_dNdyptdpt_eta_max;
    DATA->dNdyptdpt_eta_max = temp_dNdyptdpt_eta_max;

    cout << "Done ReadInData3." << endl;
    delete util;
}   /* ReadInData3 */
