#include "mpi.h"
#include <stdio.h>
#include "util.h"
#include "grid.h"
#include "data.h"
#include "init.h"
#include "eos.h"
#include "freeze.h"
#include "evolve.h"
#include <sys/stat.h>// for mkdir
// #include "int.h"

using namespace std;
void ReadInData(InitData *DATA, string file);
void ReadInData2(InitData *DATA, string file);


// main program
int main(int argc, char *argv[])
{
  int rank;
  int size;
  Grid ***arena; 
  Grid ***Lneighbor; 
  Grid ***Rneighbor;
 
  string input_file;
  static InitData DATA;
//   int flag;
//   double  tau;


  // you have the option to give a second command line option, which is an integer to be added to the random seed from the current time
  // because on the cluster it will happen that two jobs start at exactly the same time, this makes sure that they dont run with excatly the same seed
  string sseed;
  int seed = 0;
  if (argc>2)
    {
      sseed = argv[2];
      seed = atoi(sseed.c_str()); 
    }
  seed *= 10000;

  DATA.seed = seed;

  if(argc >1)   input_file = *(argv+1);
  else input_file = "";

  
  ReadInData2(&DATA, input_file);

  // initialize MPI
  if(DATA.mode!=8)
    {
      MPI::Init(argc, argv);

      rank = MPI::COMM_WORLD.Get_rank(); //number of current processor
      size = MPI::COMM_WORLD.Get_size(); //total number of processors

      DATA.rank = rank;
      DATA.size = size;

      cout << "This is processor " << rank << "/" << size << ": READY." << endl;
    }
  else
    {  
      rank=0;
      size=1;
    }

  for(int i=1; i<50; i++)
    {
      cout << "#" << endl;
    }

  if (DATA.neta%size!=0 && DATA.mode < 3)
    {
      cerr << " Number of cells in eta direction " << DATA.neta << " is not a multiple of the number of processors " << size 
	   << ". Exiting." << endl;
      if(DATA.mode!=8)
	{
	  MPI::Finalize();
	}
      exit(1);
    }

  if (DATA.neta/size<4 && DATA.mode < 3)
    {
      cerr << " Number of cells in eta direction per processor " << DATA.neta/size << " is less than 4. Exiting." << endl;
      if(DATA.mode!=8)
	{
	  MPI::Finalize();
	}
      exit(1);
    }

  // reduce the lattice size on each processor (this is slicing it in 'size' pieces)
  if (size>0)
    DATA.neta=DATA.neta/size;
  
  EOS *eos;
  eos = new EOS;
  
  if (DATA.whichEOS==0)
    {
      if (rank == 0) 
        cout << "Using the ideal gas EOS" << endl;
        eos->init_eos0();
    }
  else if (DATA.whichEOS==1)
    {
      if (rank == 0)
	cout << "Using EOS-Q from AZHYDRO" << endl;
      eos->init_eos();
    }
  else if (DATA.whichEOS==2)
    {
      if (rank == 0)
	cout << "Using lattice EOS from Huovinen/Petreczky" << endl;
      eos->init_eos2();
    }
  else if (DATA.whichEOS==3)
    {
      if (rank == 0)
	cout << "Using lattice EOS from Huovinen/Petreczky with partial chemical equilibrium (PCE) chem. f.o. at 150 MeV" << endl;
      eos->init_eos3(1);
    }
  else if (DATA.whichEOS==4)
    {
      if (rank == 0)
	cout << "Using lattice EOS from Huovinen/Petreczky with partial chemical equilibrium (PCE) chem. f.o. at 155 MeV" << endl;
      eos->init_eos3(2);
    }
  else if (DATA.whichEOS==5)
    {
      if (rank == 0)
	cout << "Using lattice EOS from Huovinen/Petreczky with partial chemical equilibrium (PCE) chem. f.o. at 160 MeV" << endl;
      eos->init_eos3(3);
    }
  else if (DATA.whichEOS==6)
    {
      if (rank == 0)
	cout << "Using lattice EOS from Huovinen/Petreczky with partial chemical equilibrium (PCE) chem. f.o. at 165 MeV" << endl;
      eos->init_eos3(4);
    }
  else 
    {
      if (rank == 0)
	cout << "No EOS for whichEOS = " << DATA.whichEOS << ". Use EOS_to_use = 0 (ideal gas) 1 (AZHYDRO EOS-Q), 2 (s95p-v1), 3 (s95p-PCE150-v1), 4 (s95p-PCE155-v1), 5 (s95p-PCE160-v1), 6 (s95p-PCE165-v1)" << endl;
      exit(1);
    }
      
  if (DATA.mode == 1 || DATA.mode == 2)
    {
      system("rm surface.dat surface?.dat surface??.dat 2> /dev/null");
      Glauber *glauber;
      glauber = new Glauber;
      
      Init *init;
      init = new Init(eos,glauber);
      
      Evolve *evolve;
      evolve = new Evolve(eos);
      
      cout << "init Glauber" << endl;
      glauber->initGlauber(DATA.SigmaNN, DATA.Target, DATA.Projectile, DATA.b, DATA.LexusImax, size, rank);
      cout << "size=" << size << endl;
      cout << "rank=" << rank << endl;
      
      init->InitArena(&DATA, &arena, &Lneighbor, &Rneighbor, size, rank);
     
      size = DATA.size;
      rank = DATA.rank;
 
      cout << "size=" << size << endl;
      cout << "rank=" << rank << endl;
      evolve->EvolveIt(&DATA, arena, Lneighbor, Rneighbor, size, rank); 

      MPI_Barrier(MPI_COMM_WORLD);
      
//       tau = DATA.tau0 + DATA.tau_size; 
     
      if (DATA.output_hydro_debug_info) {
        FILE *t4_file;
        const char* t4_name = "avgT.dat";
        t4_file = fopen(t4_name, "a");
        fprintf(t4_file,"time average: %f GeV and %f GeV", DATA.avgT/DATA.nSteps*hbarc, DATA.avgT2/DATA.nSteps2*hbarc );
        fclose(t4_file);
      }
      
      cout << "average plasma T=" << DATA.avgT/DATA.nSteps*hbarc << " GeV." << endl; 
      cout << "average plasma+mixed T=" << DATA.avgT2/DATA.nSteps2*hbarc << " GeV." << endl; 
      //    delete evolve;
      //delete glauber;
      //delete init;
    }
  
  if (DATA.mode == 1 || DATA.mode == 3 || DATA.mode == 4 || DATA.mode >= 5)
    {
      mkdir("./outputs", 0755);
      //if (rank == 0) int ret = system("rm yptphiSpectra*");
      //  freeze-out
      Freeze *freeze;
      freeze = new Freeze();
      if(DATA.pseudofreeze) freeze->CooperFrye_pseudo(DATA.particleSpectrumNumber, DATA.mode, &DATA, eos, size, rank);
      else  freeze->CooperFrye(DATA.particleSpectrumNumber, DATA.mode, &DATA, eos, size, rank);
      delete freeze;
    }
  if(DATA.mode!=8)
    {
      MPI::Finalize();
    }
}/* main */
  

void ReadInData(InitData *DATA, string file)
{
  int m, n;
  double tempf;
  Util *util;
  util = new Util();
  
//   DATA->Projectile = util->char_malloc(20);
//   DATA->Target = util->char_malloc(20);
//   DATA->initName = util->char_malloc(80);
  
  DATA->LexusImax = util->IFind(file, "Lexus_Imax");
  DATA->b = util->DFind(file, "Impact_parameter");
  DATA->Target = util->StringFind(file, "Target");
  DATA->Projectile = util->StringFind(file, "Projectile");
  DATA->SigmaNN = util->DFind(file, "SigmaNN");
  DATA->Initial_profile = util->IFind(file, "Initial_profile");
  DATA->initializeEntropy = util->IFind(file, "initialize_with_entropy");
  DATA->useEpsFO = util->IFind(file, "use_eps_for_freeze_out");
  DATA->epsilonFreeze = util->DFind(file, "epsilon_freeze");
  DATA->TFO = util->DFind(file, "T_freeze");
  DATA->particleSpectrumNumber = util->IFind(file, "particle_spectrum_to_compute");
  DATA->mode = util->IFind(file, "mode"); // 1: do everything; 2: do hydro evolution only; 3: do calculation of thermal spectra only;
  // 4: do resonance decays only
  DATA->whichEOS = util->IFind(file, "EOS_to_use");
  DATA->NumberOfParticlesToInclude = util->IFind(file, "number_of_particles_to_include");
  DATA->freezeOutMethod = util->IFind(file, "freeze_out_method");
  DATA->rhoB0 = util->DFind(file, "rho_b_max");
  DATA->hard = util->DFind(file, "binary_collision_scaling_fraction");
  DATA->facTau = util->IFind(file, "average_surface_over_this_many_time_steps"); 
  
  DATA->nx = util->IFind(file, "Grid_size_in_x");
  DATA->ny = util->IFind(file, "Grid_size_in_y");
  DATA->neta = util->IFind(file, "Grid_size_in_eta");
  
  DATA->x_size = util->DFind(file, "X_grid_size_in_fm");
  DATA->y_size = util->DFind(file, "Y_grid_size_in_fm");
  DATA->eta_size = util->DFind(file, "Eta_grid_size");
  DATA->tau_size = util->DFind(file, "Total_evolution_time_tau");
  DATA->tau0 = util->DFind(file, "Initial_time_tau_0");
  
  /* x-grid, for instance, runs from 0 to nx */
  DATA->delta_x = (DATA->x_size)/( ((double) DATA->nx) ); 
  DATA->delta_y = (DATA->y_size)/( ((double) DATA->ny) ); 
  DATA->delta_eta = (DATA->eta_size)/( ((double) DATA->neta) ); 
  
  cerr << " DeltaX=" << DATA->delta_x << endl;
  cerr << " DeltaY=" << DATA->delta_y << endl;
  cerr << " DeltaETA=" << DATA->delta_eta << endl;
  
  /* CFL condition : delta_tau < min(delta_x/10, tau0 delta_eta/10) */
  
  DATA->delta_tau = util->DFind(file, "Delta_Tau");
  tempf = mini(DATA->delta_x/10.0, (DATA->tau0)*(DATA->delta_eta/10.0));
  if(tempf < DATA->delta_tau) DATA->delta_tau = tempf;
  
  DATA->rotateBy45degrees = util->IFind(file, "rotate_by_45_degrees");
  DATA->outputEvolutionData = util->IFind(file, "output_evolution_data");
  
  DATA->nt = (int) (floor(DATA->tau_size/(DATA->delta_tau) + 0.5));
  cout << "ReadInData: Time step size = " << DATA->delta_tau << endl;
  cout << "ReadInData: Number of time steps required = " << DATA->nt << endl;
  
  DATA->epsilon0 = util->DFind(file, "Maximum_energy_density");
  DATA->eta_fall_off = util->DFind(file, "Eta_fall_off");
  DATA->eta_flat = util->DFind(file, "Eta_plateau_size");
  DATA->R_A = util->DFind(file, "Initial_radius_size_in_fm");
  DATA->a_A = util->DFind(file, "Initial_radial_fall_off_in_fm");
  
  DATA->a_short = util->DFind(file, "Initial_asym_short_axis_in_fm");
  DATA->a_long  = util->DFind(file, "Initial_asym_long_axis_in_fm");
  
  DATA->sFactor = util->DFind(file, "s_factor");
  
  // for calculation of spectra:
  DATA->ymax = util->DFind(file, "maximal_rapidity");
  DATA->deltaY = util->DFind(file, "delta_y");
  
  if(DATA->turn_on_rhob == 1)
    {
      DATA->alpha_max = 5;
    }
  else
    {
      DATA->alpha_max = 4;
    }
  
  DATA->rk_order = util->IFind(file, "Runge_Kutta_order");
  
  if( (DATA->rk_order != 1) )
    {
      cout << "Runge-Kutta order = " << DATA->rk_order << endl;
    }
  
  // in case of using Initial_Distribution 3, these are the limits between which to sample the impact parameter
  DATA->bmin = util->DFind(file, "bmin");
  DATA->bmax = util->DFind(file, "bmax");
  
  DATA->includeJet = util->IFind(file,"Include_Jet");
  DATA->includeTrigger = util->IFind(file,"Include_Trigger");
  DATA->minmod_theta = util->DFind(file, "Minmod_Theta");
  
  DATA->local_y_max = util->DFind(file, "Maximum_Local_Rapidity");
  
  DATA->viscosity_flag = util->IFind(file, "Viscosity_Flag_Yes_1_No_0");
  
  DATA->verbose_flag = util->IFind(file, "Verbose_Flag_Yes_1_No_0");
  
  DATA->eps_limit  = util->DFind(file, "Minimum_Epsilon_Frac_For_Visc");
  
  DATA->turn_on_rhob = util->IFind(file, "Include_Rhob_Yes_1_No_0");
  
  DATA->turn_on_shear = util->IFind(file, "Include_Shear_Visc_Yes_1_No_0");
  
  DATA->turn_on_bulk = util->IFind(file, "Include_Bulk_Visc_Yes_1_No_0");
  
  DATA->zero_or_stop = util->IFind(file, "Vicous_Trouble_Zero=0_or_Stop=1:");
  
  DATA->tau_pi = util->DFind(file, "Shear_relaxation_time_tau_pi");
  DATA->tau_b_pi = util->DFind(file, "Bulk_relaxation_time_tau_b_pi");
  DATA->shear_to_s = util->DFind(file, "Shear_to_S_ratio");
  DATA->T_dependent_shear_to_s = util->IFind(file, "T_dependent_Shear_to_S_ratio");
  DATA->bulk_to_s = util->DFind(file, "Bulk_to_S_ratio");
  
  DATA->include_deltaf = util->IFind(file, "Include_deltaf");
  
  DATA->doFreezeOut = util->IFind(file, "Do_FreezeOut_Yes_1_No_0");
  DATA->sigma0 = util->DFind(file, "sigma_0");
  
  DATA->initName = util->StringFind(file,"Initial_Distribution_Filename");
  
  /* initialize the metric, mostly plus */
  
  DATA->avgT=0.;
  DATA->nSteps=0;
  DATA->avgT2=0.;
  DATA->nSteps2=0;
  DATA->plasmaEvolutionTime=0.;
  DATA->plasmaEvolutionTime2=0.;
  
  DATA->gmunu = util->mtx_malloc(4,4);
  for(m=0; m<4; m++)
    {
      for(n=0; n<4; n++)
	{
	  if(m == n) (DATA->gmunu)[m][n] = 1.0;
	  else (DATA->gmunu)[m][n] = 0.0;
	  
	  if(m==0 && n==0) (DATA->gmunu)[m][n] *= -1.0;
	}
    }/* m */
  
 //These were added for writing the coordinates of collisions to file for
 //MARTINI to sample. -CFY 11/2/2010
 DATA->Nbin_to_file = util->IFind(file, "Write_Nbin_to_File");

 //
 DATA->outputBinaryEvolution = util->IFind(file, "outputBinaryEvolution");

  cout << "Done ReadInData." << endl;
  delete util;
}/* ReadInData */




// Set default values for all runtime parameters, 
// then change them to match the values optionally placed in the input file.
// This way one can clean the clutterred input file of unused or unimportant
// parameters, and we hopefully no longer need to avoid adding parameters that can be changed
// at runtime.  Any old or new input file should still be functional, but 
// it should be checked that it runs as expected (i.e., that a new parameter
// hasn't changed an old default behavior, which I can't vouch for).  ML 04/2013
void ReadInData2(InitData *DATA, string file)
{
  int m, n;
  double tempf;
  Util *util;
  util = new Util();
  string tempinput;
  
//   DATA->Projectile = util->char_malloc(20);
//   DATA->Target = util->char_malloc(20);
//   DATA->initName = util->char_malloc(80);

  // Lexus_Imax: number of points to use to store the data of the thickness function
//   DATA->LexusImax = util->IFind(file, "Lexus_Imax");
  int tempLexusImax = 100;
  tempinput = util->StringFind3(file, "Lexus_Imax");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempLexusImax;
  DATA->LexusImax = tempLexusImax;
  

  // Impact_parameter: in fm, >=0
//   DATA->b = util->DFind(file, "Impact_parameter");
  double tempb = 0;
  tempinput = util->StringFind3(file, "Impact_parameter");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempb;
  DATA->b = tempb;
  if(DATA->b<0) 
  {
    cerr << "Impact parameter must be greater than zero\n";
    exit(1);
  }


  // Target, Projectile:  any name as defined in known_nuclei.dat
//   DATA->Target = util->StringFind(file, "Target");
  string tempTarget = "Pb";
  tempinput = util->StringFind3(file, "Target");
  if(tempinput != "empty") tempTarget.assign(tempinput);
  DATA->Target.assign(tempTarget);
//   DATA->Projectile = util->StringFind(file, "Projectile");
  string tempProjectile = "Pb";
  tempinput = util->StringFind3(file, "Projectile");
  if(tempinput != "empty") tempProjectile.assign(tempinput);
  DATA->Projectile.assign(tempProjectile);

  
  // SigmaNN: nucleon-nucleon cross section in mb
//   DATA->SigmaNN = util->DFind(file, "SigmaNN");
  double tempSigmaNN = 60.;
  tempinput = util->StringFind3(file, "SigmaNN");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempSigmaNN;
  DATA->SigmaNN = tempSigmaNN;
  if(DATA->SigmaNN<0) 
  {
    cerr << "NN cross section must be greater than zero (in mb)\n";
    exit(1);
  }
  
  
  // Initial_profile: 
//   0: Sangyong’s simple profile
//   1: Glauber model with variable binary collision scaling fraction and scaling of either entropy- or energy-density.
//   2: Test scenario for the freeze-out surface finder
//   3: Event-by-event Glauber MC
//   4: for testing the Glauber MC initial condition
//   5: Something like p+p
//   6: Read in initial profile from a file
//   DATA->Initial_profile = util->IFind(file, "Initial_profile");
  int tempInitial_profile = 1;
  tempinput = util->StringFind3(file, "Initial_profile");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempInitial_profile;
  DATA->Initial_profile = tempInitial_profile;
  if(DATA->Initial_profile>6 || DATA->Initial_profile<0) 
  {
    cerr << "Initial profile" << DATA->Initial_profile << "not defined\n";
    exit(1);
  }

  //Select the profile to use in eta for the energy/entropy initialisation
  //Warning: profile 2 is only available when smooth transverse initial conditions (DATA->Initial_profile == 2) are used
  //1 for Hirano's central plateau + Gaussian decay
  //2 for a Woods-Saxon profile
  int tempinitial_eta_profile = 1;
  tempinput = util->StringFind3(file, "initial_eta_profile");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempinitial_eta_profile;
  if (2 == DATA->Initial_profile) {
  	DATA->initial_eta_profile = tempinitial_eta_profile;  
  }
  else {
	if (tempinitial_eta_profile != 1) {
		cerr << "Initial eta profile " << tempinitial_eta_profile << " can only be used with transverse Initial_profile==2. Using the default eta profile.\n";
	}
	DATA->initial_eta_profile=1;
  }
  
  
  //initialize_with_entropy:
  //0: scale with energy density
  //1: scale with entropy density
//   DATA->initializeEntropy = util->IFind(file, "initialize_with_entropy");
  int tempinitializeEntropy = 0;
  tempinput = util->StringFind3(file, "initialize_with_entropy");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempinitializeEntropy;
  DATA->initializeEntropy = tempinitializeEntropy;
  if(DATA->initializeEntropy>1 || DATA->initializeEntropy<0) 
  {
    cerr << "Must initialize with entropy (initialize_with_entropy=1) or energy (0)\n";
    exit(1);
  }
  
  DATA->useEpsFO = 2; // Must set either freeze out energy density or temperature, otherwise generate error message below.
  
  // T_freeze:  freeze out temperature
  // only used with use_eps_for_freeze_out = 0
//   DATA->TFO = util->DFind(file, "T_freeze");
  double tempTFO = 0.12;
  tempinput = util->StringFind3(file, "T_freeze");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempTFO;
  DATA->TFO = tempTFO;
//   if(tempinput != "empty" && DATA->useEpsFO != 0) cerr << "T_freeze unused when use_eps_for_freeze_out set\n";
  if(tempinput != "empty") DATA->useEpsFO = 0;  // if only freeze out temperature is set, freeze out by temperature
  int tfoset = 0;
  if(tempinput != "empty") tfoset = 1;
  
  
  // epsilon_freeze: freeze-out energy density in GeV/fm^3
  // only used with use_eps_for_freeze_out = 1
//   DATA->epsilonFreeze = util->DFind(file, "epsilon_freeze");
  double tempepsilonFreeze = 0.12;
  tempinput = util->StringFind3(file, "epsilon_freeze");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempepsilonFreeze;
  DATA->epsilonFreeze = tempepsilonFreeze;
  if(DATA->epsilonFreeze<=0) 
  {
    cerr << "Freeze out energy density must be greater than zero\n";
    exit(1);
  }
  if(tempinput != "empty") DATA->useEpsFO = 1; // if epsilon_freeze is set, freeze out by epsilon
  

  //use_eps_for_freeze_out: 
  // 0: freeze out at constant temperature T_freeze
  // 1: freeze out at constant energy density epsilon_freeze
  // if set in input file, overide above defaults
//   DATA->useEpsFO = util->IFind(file, "use_eps_for_freeze_out");
  int tempuseEpsFO = DATA->useEpsFO;
  tempinput = util->StringFind3(file, "use_eps_for_freeze_out");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempuseEpsFO;
  DATA->useEpsFO = tempuseEpsFO;
  if(DATA->useEpsFO>1 || DATA->useEpsFO<0) 
  {
    cerr << "Error: did not set either freeze out energy density or temperature, or invalid option for use_eps_for_freeze_out:" << DATA->useEpsFO << endl;
    exit(1);
  }
  if(tfoset == 1 && DATA->useEpsFO == 1) 
    cerr << "T_freeze set but overridden -- freezing out by energy density at " 
    << DATA->epsilonFreeze << " GeV/fm^3\n";
  
  
  
  //particle_spectrum_to_compute:
  // 0: Do all up to number_of_particles_to_include
  // any natural number: Do the particle with this (internal) ID
//   DATA->particleSpectrumNumber = util->IFind(file, "particle_spectrum_to_compute");
  int tempparticleSpectrumNumber = 0;
  tempinput = util->StringFind3(file, "particle_spectrum_to_compute");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempparticleSpectrumNumber;
  DATA->particleSpectrumNumber = tempparticleSpectrumNumber;
  
  
  // mode: 
  // 1: Does everything. Evolution. Computation of thermal spectra. Resonance decays. Observables.  Only compatible with freeze_out_method=3 and pseudofreeze=1
  // 2: Evolution only.
  // 3: Compute all thermal spectra only.
  // 4: Resonance decays only.
  // 5: Resonance decays for just one specific particle (only for testing - this will miss the complete chain of decays).
  // 6: only combine charged hadrons - can be used for any postprocessing with the stored results */
//   DATA->mode = util->IFind(file, "mode"); // 1: do everything; 2: do hydro evolution only; 3: do calculation of thermal spectra only;
  // 13: Compute observables from previously-computed thermal spectra
  // 14: Compute observables from post-decay spectra
  int tempmode = 1;
  tempinput = util->StringFind3(file, "mode");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempmode;
  else
  {
    cerr << "Must specify mode. Exiting.\n";
    exit(1);
  }
//   if (tempmode == 1) 
//   {
//     cerr << "mode=1 not currently functional.  Run each step separately with mode>=2\n";
//     exit(1);
//   }
  DATA->mode = tempmode;
  
  
  
  //EOS_to_use:
  // 0: ideal gas
  // 1: EOS-Q from azhydro
  // 2: lattice EOS from Huovinen and Petreczky
  // 3: lattice EOS from Huovinen and Petreczky with partial chemical equilibrium (PCE) at 150 MeV (see https://wiki.bnl.gov/TECHQM/index.php/QCD_Equation_of_State)
  // 4: PCE EOS with chemical freeze out at 155 MeV
  // 5: PCE EOS at 160 MeV
  // 6: PCE EOS at 165 MeV
//   DATA->whichEOS = util->IFind(file, "EOS_to_use");
  int tempwhichEOS = 9;
  tempinput = util->StringFind3(file, "EOS_to_use");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempwhichEOS;
  DATA->whichEOS = tempwhichEOS;
  if(DATA->whichEOS>6 || DATA->whichEOS<0) 
  {
    cerr << "EOS_to_use unspecified or invalid option:" << DATA->whichEOS << endl;
    exit(1);
  }
  
  
  //number_of_particles_to_include:  This determines up to which particle in the list spectra should be computed (mode=3) or resonances should be included (mode=4)
  // current maximum = 319
//   DATA->NumberOfParticlesToInclude = util->IFind(file, "number_of_particles_to_include");
  int tempNumberOfParticlesToInclude = 2;
  tempinput = util->StringFind3(file, "number_of_particles_to_include");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempNumberOfParticlesToInclude;
  DATA->NumberOfParticlesToInclude = tempNumberOfParticlesToInclude;
  if(DATA->NumberOfParticlesToInclude>320) 
  {
    cerr << "Invalid option for number_of_particles_to_include:" << DATA->NumberOfParticlesToInclude << endl;
    exit(1);
  }
  
  
  
  
  // freeze_out_method:
  // 1: Hirano's simplified method
  // 2: Schenke's more complex method
  // 3: Luzum's simple method
//   DATA->freezeOutMethod = util->IFind(file, "freeze_out_method");
  int tempfreezeOutMethod = 3;//  This default allows for new users to run MUSIC without warnings.  Should set to 2 in input file for production use.
  tempinput = util->StringFind3(file, "freeze_out_method");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempfreezeOutMethod;
  DATA->freezeOutMethod = tempfreezeOutMethod;
  if(DATA->freezeOutMethod>3) 
  {
    cerr << "Invalid option for freeze_out_method:" << DATA->freezeOutMethod << endl;
    exit(1);
  }
  
  
  // max_delta_eta:  maximum size of freeze out surface segment in eta direction.
  // Even when hydro variables vary slowly in eta (e.g., in a boost-invariant solution),
  // the Cosh(y-eta) factor can vary within a step in eta of delta_eta is too large.
  // The freeze out surface will be subdivided into identical slices in eta with size
  // less than max_delta_eta.
  // Only used with freeze_out_method = 3.
  double tempmax_delta_eta = 5.;  // if this is larger than delta_eta, it does nothing
  tempinput = util->StringFind3(file, "max_delta_eta");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempmax_delta_eta;
  DATA->max_delta_eta = tempmax_delta_eta;
  if(tempinput != "empty" && DATA->freezeOutMethod != 3) cerr << "max_delta_eta unused when freeze_out_method != 3\n";
  
  
  // max_delta_eta2:  maximum size of freeze out surface segment in eta direction.
  // Alternate to max_delta_eta above 
  // (Though both can be used in combination.  Only has effect if max_delta_eta2 < max_delta_eta):
  // Integrate each segment over eta directly in Cooper-Frye calculation, rather than 
  // writing extra lines to surface.dat that need to be read in. Faster this way.
  double tempmax_delta_eta2 = 5.;
  tempinput = util->StringFind3(file, "max_delta_eta2");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempmax_delta_eta2;
  DATA->max_delta_eta2 = tempmax_delta_eta2;
  if(tempinput != "empty" && DATA->freezeOutMethod != 3) cerr << "max_delta_eta2 unused when freeze_out_method != 3\n";

  
  
  
  // rho_b_max: maximum baryon density for zero impact parameter. The shape of the ρB distribution is the same as that for the energy/entropy density distribution
  // any number > 0.  (don’t use 0, however if the EOS does not support finite baryon chemical potential, it does not matter what you put)
//   DATA->rhoB0 = util->DFind(file, "rho_b_max");
  double temprhoB0 = 0.000000001;
  tempinput = util->StringFind3(file, "rho_b_max");
  if(tempinput != "empty") istringstream ( tempinput ) >> temprhoB0;
  DATA->rhoB0 = temprhoB0;
  
  
  
  
  
  // binary_collision_scaling_fraction: (1-f)Nwn + f Nbc
//   DATA->hard = util->DFind(file, "binary_collision_scaling_fraction");
  double temphard = 0.05;
  tempinput = util->StringFind3(file, "binary_collision_scaling_fraction");
  if(tempinput != "empty") istringstream ( tempinput ) >> temphard;
  DATA->hard = temphard;
  
  
  
  
  // average_surface_over_this_many_time_steps:
  // Only save every N timesteps for finding freeze out surface
//   DATA->facTau = util->IFind(file, "average_surface_over_this_many_time_steps"); 
  int tempfacTau = 10;
  tempinput = util->StringFind3(file, "average_surface_over_this_many_time_steps");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempfacTau;
  DATA->facTau = tempfacTau;
  
  
  
  
  //  Grid_size_in_*
  // number of cells in x,y direction
//   DATA->nx = util->IFind(file, "Grid_size_in_x");
  int tempnx = 10;
  tempinput = util->StringFind3(file, "Grid_size_in_x");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempnx;
  DATA->nx = tempnx;
//   DATA->ny = util->IFind(file, "Grid_size_in_y");
  int tempny = 10;
  tempinput = util->StringFind3(file, "Grid_size_in_y");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempny;
  DATA->ny = tempny;
  
  
  
  // Grid_size_in_eta
  // number of cells in eta direction.
  // Must have at least 4 cells per processor.
  // Must be an even number.
  // One cell is positioned at eta=0, 
  // half the cells are at negative eta,
  // the rest (one fewer) are at positive eta
//   DATA->neta = util->IFind(file, "Grid_size_in_eta");
  int tempneta = 4;
  tempinput = util->StringFind3(file, "Grid_size_in_eta");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempneta;
  DATA->neta = tempneta;
  
  
  // *_grid_size_in_fm:  total length of box in x,y direction in fm (minus delta_*)
//   DATA->x_size = util->DFind(file, "X_grid_size_in_fm");
//   DATA->y_size = util->DFind(file, "Y_grid_size_in_fm");
  double tempx_size = 25.;
  tempinput = util->StringFind3(file, "X_grid_size_in_fm");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempx_size;
  DATA->x_size = tempx_size;
  double tempy_size = 25.;
  tempinput = util->StringFind3(file, "Y_grid_size_in_fm");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempy_size;
  DATA->y_size = tempy_size;
  
  
  
  // Eta_grid_size:  total length of box in eta direction (minus delta_eta)
  // e.g., neta=8 and eta_size=8 has 8 cells that run from eta=-4 to eta=3
//   DATA->eta_size = util->DFind(file, "Eta_grid_size");
  double tempeta_size = 8.;
  tempinput = util->StringFind3(file, "Eta_grid_size");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempeta_size;
  DATA->eta_size = tempeta_size;
  
  
  // Total_evolution_time_tau
  // otal evolution time in [fm]. in case of freeze_out_method = 2,3, evolution will halt erlier if all cells are frozen out.
//   DATA->tau_size = util->DFind(file, "Total_evolution_time_tau");
  double temptau_size = 30.;
  tempinput = util->StringFind3(file, "Total_evolution_time_tau");
  if(tempinput != "empty") istringstream ( tempinput ) >> temptau_size;
  DATA->tau_size = temptau_size;
  
  
  
  // Initial_time_tau_0:  in fm
//   DATA->tau0 = util->DFind(file, "Initial_time_tau_0");
  double temptau0 = 0.4;
  tempinput = util->StringFind3(file, "Initial_time_tau_0");
  if(tempinput != "empty") istringstream ( tempinput ) >> temptau0;
  DATA->tau0 = temptau0;
  
  
  
  /* x-grid, for instance, runs from 0 to nx */
  DATA->delta_x = (DATA->x_size)/( ((double) DATA->nx) ); 
  DATA->delta_y = (DATA->y_size)/( ((double) DATA->ny) ); 
  DATA->delta_eta = (DATA->eta_size)/( ((double) DATA->neta) ); 
  
  cerr << " DeltaX=" << DATA->delta_x << endl;
  cerr << " DeltaY=" << DATA->delta_y << endl;
  cerr << " DeltaETA=" << DATA->delta_eta << endl;
  
  
  // Delta_Tau: 
  // time step to use in [fm]. If a too large value is given, it will automatically be reduced to the maximal acceptable value according to the CFL condition.
  /* CFL condition : delta_tau < min(delta_x/10, tau0 delta_eta/10) */
//   DATA->delta_tau = util->DFind(file, "Delta_Tau");
  double tempdelta_tau = 8.;
  tempinput = util->StringFind3(file, "Delta_Tau");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempdelta_tau;
  DATA->delta_tau = tempdelta_tau;

  tempf = mini(DATA->delta_x/10.0, (DATA->tau0)*(DATA->delta_eta/10.0));
  if(tempf < DATA->delta_tau) DATA->delta_tau = tempf;
  
  cerr << " DeltaTau=" << DATA->delta_tau << endl;
  
  
  
  // rotate_by_45_degrees: rotates initial condition by Pi/4
  // only used if Initial_profile==1
//   DATA->rotateBy45degrees = util->IFind(file, "rotate_by_45_degrees");
  int temprotateBy45degrees = 0;
  tempinput = util->StringFind3(file, "rotate_by_45_degrees");
  if(tempinput != "empty") istringstream ( tempinput ) >> temprotateBy45degrees;
  DATA->rotateBy45degrees = temprotateBy45degrees;
  
  
  
  // output_evolution_data:  
  // 1:  output bulk information at every grid point at every time step
//   DATA->outputEvolutionData = util->IFind(file, "output_evolution_data");
  int tempoutputEvolutionData = 0;
  tempinput = util->StringFind3(file, "output_evolution_data");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempoutputEvolutionData;
  DATA->outputEvolutionData = tempoutputEvolutionData;
  
  DATA->nt = (int) (floor(DATA->tau_size/(DATA->delta_tau) + 0.5));
  cout << "ReadInData: Time step size = " << DATA->delta_tau << endl;
  cout << "ReadInData: Number of time steps required = " << DATA->nt << endl;
  
  
/*  // Maximum_energy_density:
  // determines the maximum energy density at zero impact parameter given in [GeV/fm3] (for initialize_with_entropy = 0)
or the maximum entropy density at zero impact parameter given in [1/fm3]
(for initialize_with_entropy = 1)
//   DATA->epsilon0 = util->DFind(file, "Maximum_energy_density");*/
  double tempepsilon0 = 21.67;
  tempinput = util->StringFind3(file, "Maximum_energy_density");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempepsilon0;
  DATA->epsilon0 = tempepsilon0;
  
  
  
  // Eta_fall_off:
  // width of half-Gaussian on each side of a central pleateau in eta
//   DATA->eta_fall_off = util->DFind(file, "Eta_fall_off");
  double tempeta_fall_off  = 0.4;
  tempinput = util->StringFind3(file, "Eta_fall_off");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempeta_fall_off ;
  DATA->eta_fall_off  = tempeta_fall_off;
  
  
  // Eta_plateau_size:
  // width of the flat region symmetrical around eta=0
//   DATA->eta_flat = util->DFind(file, "Eta_plateau_size");
  double tempeta_flat   = 5.9;
  tempinput = util->StringFind3(file, "Eta_plateau_size");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempeta_flat  ;
  DATA->eta_flat   = tempeta_flat;
  
  
  //Initial_radius_size_in_fm
  // initial radial size in [fm] for Initial_profile = 0
//   DATA->R_A = util->DFind(file, "Initial_radius_size_in_fm");
  double tempR_A   = 2.6;
  tempinput = util->StringFind3(file, "Initial_radius_size_in_fm");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempR_A  ;
  DATA->R_A   = tempR_A;
  if(tempinput != "empty" && DATA->Initial_profile != 0) cerr << "Initial_radius_size_in_fm unused when Initial_profile != 0\n";

  
  // Initial_radial_fall_off_in_fm:
  // radial fall-off in [fm] for Initial_profile = 0
//   DATA->a_A = util->DFind(file, "Initial_radial_fall_off_in_fm");
  double tempa_A   = 0.5;
  tempinput = util->StringFind3(file, "Initial_radial_fall_off_in_fm");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempa_A  ;
  DATA->a_A   = tempa_A  ;
  if(tempinput != "empty" && DATA->Initial_profile != 0) cerr << "Initial_radial_fall_off_in_fm unused when Initial_profile != 0\n";

  
  // Initial_asym_long_axis_in_fm:
  // length of the long axis for Initial_profile = 0 in [fm]
//   DATA->a_long = util->DFind(file, "Initial_asym_long_axis_in_fm");
  double tempa_long   = 2.;
  tempinput = util->StringFind3(file, "Initial_asym_long_axis_in_fm");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempa_long  ;
  DATA->a_long   = tempa_long;
  if(tempinput != "empty" && DATA->Initial_profile != 0) cerr << "Initial_asym_long_axis_in_fm unused when Initial_profile != 0\n";

  
  
  
  // Initial_asym_short_axis_in_fm:
  // length of the short axis for Initial_profile = 0 in [fm]
//   DATA->a_short = util->DFind(file, "Initial_asym_short_axis_in_fm");
  double tempa_short   = 2.;
  tempinput = util->StringFind3(file, "Initial_asym_short_axis_in_fm");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempa_short  ;
  DATA->a_short   = tempa_short;
  if(tempinput != "empty" && DATA->Initial_profile != 0) cerr << "Initial_asym_short_axis_in_fm unused when Initial_profile != 0\n";


  
  // s_factor:  for use with IP-Glasma initial conditions
//   DATA->sFactor = util->DFind(file, "s_factor");
  double tempsFactor   = 20.;
  tempinput = util->StringFind3(file, "s_factor");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempsFactor  ;
  DATA->sFactor   = tempsFactor;
  
  
  // for calculation of spectra:
  
  // maximal_rapidity:  spectra calculated from zero to this rapidity in +y and -y
//   DATA->ymax = util->DFind(file, "maximal_rapidity");
  double tempymax   = 2.;
  tempinput = util->StringFind3(file, "maximal_rapidity");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempymax  ;
  DATA->ymax   = tempymax;
  
  
  // delta_y:
  // step size in rapidity in calculation of spectra
//   DATA->deltaY = util->DFind(file, "delta_y");
  double tempdeltaY   = 0.2;
  tempinput = util->StringFind3(file, "delta_y");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempdeltaY  ;
  DATA->deltaY   = tempdeltaY;
  
  
  // max_pseudorapidity:  spectra calculated from zero to this pseudorapidity in +eta and -eta
  double tempmax_pseudorapidity   = 2.4;
  tempinput = util->StringFind3(file, "max_pseudorapidity");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempmax_pseudorapidity  ;
  DATA->max_pseudorapidity   = tempmax_pseudorapidity;
  
   // pseudo_steps:
  // steps in pseudorapidity in calculation of spectra
  int temppseudo_steps   = 26;
  tempinput = util->StringFind3(file, "pseudo_steps");
  if(tempinput != "empty") istringstream ( tempinput ) >> temppseudo_steps  ;
  DATA->pseudo_steps   = temppseudo_steps; 
  
  // phi_steps
  // steps in azimuthal angle in calculation of spectra
  int tempphi_steps   = 30;
  tempinput = util->StringFind3(file, "phi_steps");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempphi_steps  ;
  DATA->phi_steps   = tempphi_steps; 
  
  // min_pt:  spectra calculated from this to max_pt transverse momentum in GeV
  double tempmin_pt   = 0.4;
  tempinput = util->StringFind3(file, "min_pt");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempmin_pt  ;
  DATA->min_pt   = tempmin_pt;
    
  // max_pt:  spectra calculated from min_pt to this transverse momentum in GeV
  double tempmax_pt   = 6;
  tempinput = util->StringFind3(file, "max_pt");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempmax_pt  ;
  DATA->max_pt   = tempmax_pt;
  
   // pt_steps:
  // steps in transverse momentum in calculation of spectra
  int temppt_steps   = 60;
  tempinput = util->StringFind3(file, "pt_steps");
  if(tempinput != "empty") istringstream ( tempinput ) >> temppt_steps  ;
  DATA->pt_steps   = temppt_steps;   
  
  // pseudofreeze
  // Calculate spectra at fixed, equally-spaced grid in pseudorapidity, pt, and phi
  int temppseudofreeze = 0;
  tempinput = util->StringFind3(file, "pseudofreeze");
  if(tempinput != "empty") istringstream ( tempinput ) >> temppseudofreeze;
  DATA->pseudofreeze = temppseudofreeze;
  if (DATA->mode == 1 && (DATA->pseudofreeze!=1 || DATA->freezeOutMethod!=3)) 
  {
    cerr << "mode=1 only works with freeze_out_method=3 and pseudofreeze=1.  Run each step separately with mode>=2\n";
    exit(1);
  }
  
  // switch for baryon current propagation
  int tempturn_on_rhob = 0;
  tempinput = util->StringFind3(file, "Include_Rhob_Yes_1_No_0");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempturn_on_rhob;
  DATA->turn_on_rhob = tempturn_on_rhob;
  if(DATA->turn_on_rhob == 1)
     DATA->alpha_max = 5;
  else
     DATA->alpha_max = 4;
  int tempturn_on_diff = 0;
  tempinput = util->StringFind3(file, "turn_on_baryon_diffusion");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempturn_on_diff;
  DATA->turn_on_diff = tempturn_on_diff;
  
  // Runge_Kutta_order:  must be 1 or 2
//   DATA->rk_order = util->IFind(file, "Runge_Kutta_order");
  int temprk_order = 1;
  tempinput = util->StringFind3(file, "Runge_Kutta_order");
  if(tempinput != "empty") istringstream ( tempinput ) >> temprk_order;
  DATA->rk_order = temprk_order;
  if(DATA->rk_order>2 || DATA->rk_order <0) 
  {
    cerr << "Invalid option for Runge_Kutta_order:" << DATA->freezeOutMethod << endl;
    exit(1);
  }
  
  if( (DATA->rk_order != 1) )
    {
      cout << "Runge-Kutta order = " << DATA->rk_order << endl;
    }
  
  // in case of using Initial_Distribution 3, these are the limits between which to sample the impact parameter
//   DATA->bmin = util->DFind(file, "bmin");
//   DATA->bmax = util->DFind(file, "bmax");
  double tempbmin   = 0.;
  tempinput = util->StringFind3(file, "bmin");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempbmin  ;
  DATA->bmin   = tempbmin;
    if(tempinput != "empty" && DATA->Initial_profile != 3) cerr << "bmin unused when Initial_profile != 3\n";

  
  double tempbmax   = 6.73;
  tempinput = util->StringFind3(file, "bmax");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempbmax  ;
  DATA->bmax   = tempbmax;
    if(tempinput != "empty" && DATA->Initial_profile != 3) cerr << "bmax unused when Initial_profile != 3\n";

  
  
  // Include_Jet: unused
//   DATA->includeJet = util->IFind(file,"Include_Jet");
  int tempincludeJet = 0;
  tempinput = util->StringFind3(file, "Include_Jet");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempincludeJet;
  DATA->includeJet = tempincludeJet;
  
  
  // Include_Trigger: unused
//   DATA->includeTrigger = util->IFind(file,"Include_Trigger");
  int tempincludeTrigger = 0;
  tempinput = util->StringFind3(file, "Include_Trigger");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempincludeTrigger;
  DATA->includeTrigger = tempincludeTrigger;
  
  // Minmod_Theta: theta parameter in the min-mod like limiter
//   DATA->minmod_theta = util->DFind(file, "Minmod_Theta");
  double tempminmod_theta   = 1.8;
  tempinput = util->StringFind3(file, "Minmod_Theta");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempminmod_theta  ;
  DATA->minmod_theta   = tempminmod_theta;
  
  
  // Maximum_Local_Rapidity: unknown.  used in reconst.cpp
//   DATA->local_y_max = util->DFind(file, "Maximum_Local_Rapidity");
  double templocal_y_max   = 20.;
  tempinput = util->StringFind3(file, "Maximum_Local_Rapidity");
  if(tempinput != "empty") istringstream ( tempinput ) >> templocal_y_max  ;
  DATA->local_y_max   = templocal_y_max;
  
  
  // Viscosity_Flag_Yes_1_No_0:   set to 0 for ideal hydro
//   DATA->viscosity_flag = util->IFind(file, "Viscosity_Flag_Yes_1_No_0");
  int tempviscosity_flag = 1;
  tempinput = util->StringFind3(file, "Viscosity_Flag_Yes_1_No_0");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempviscosity_flag;
  DATA->viscosity_flag = tempviscosity_flag;
  if(DATA->viscosity_flag == 1 && DATA->freezeOutMethod==1) 
    cerr << "freeze_out_method 1 does not currently work properly with viscosity.  Shear viscous tensor will be set to zero on freeze out surface.\n";
  
  
  // Verbose_Flag_Yes_1_No_0:  unused
//   DATA->verbose_flag = util->IFind(file, "Verbose_Flag_Yes_1_No_0");
  int tempverbose_flag = 0;
  tempinput = util->StringFind3(file, "Verbose_Flag_Yes_1_No_0");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempverbose_flag;
  DATA->verbose_flag = tempverbose_flag;
  
  
  // Minimum_Epsilon_Frac_For_Visc:  unused
//   DATA->eps_limit  = util->DFind(file, "Minimum_Epsilon_Frac_For_Visc");
  double tempeps_limit   = 0.00001;
  tempinput = util->StringFind3(file, "Minimum_Epsilon_Frac_For_Visc");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempeps_limit  ;
  DATA->eps_limit   = tempeps_limit;
  
  
  
  
  // Bulk_to_S_ratio:  constant zeta/s
//   DATA->bulk_to_s = util->DFind(file, "Bulk_to_S_ratio");
  double tempbulk_to_s   = 0.0;
  tempinput = util->StringFind3(file, "Bulk_to_S_ratio");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempbulk_to_s  ;
  DATA->bulk_to_s   = tempbulk_to_s;
  int tempturn_on_bulk = 0;
  if(tempinput != "empty") tempturn_on_bulk = 1;
  
  
  // Include_Bulk_Visc_Yes_1_No_0
//   DATA->turn_on_bulk = util->IFind(file, "Include_Bulk_Visc_Yes_1_No_0");
//   int tempturn_on_bulk = 0;
  tempinput = util->StringFind3(file, "Include_Bulk_Visc_Yes_1_No_0");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempturn_on_bulk;
  DATA->turn_on_bulk = tempturn_on_bulk;
  
  //Vicous_Trouble_Zero=0_or_Stop=1: unused
//   DATA->zero_or_stop = util->IFind(file, "Vicous_Trouble_Zero=0_or_Stop=1:");
  int tempzero_or_stop = 2;
  tempinput = util->StringFind3(file, "Vicous_Trouble_Zero=0_or_Stop=1");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempzero_or_stop;
  tempinput = util->StringFind3(file, "Viscous_Trouble_Zero=0_or_Stop=1");//correctly spelled version should work as well
  if(tempinput != "empty") istringstream ( tempinput ) >> tempzero_or_stop;
  DATA->zero_or_stop = tempzero_or_stop;
  
  
  // Shear_relaxation_time_tau_pi:  in fm?
//   DATA->tau_pi = util->DFind(file, "Shear_relaxation_time_tau_pi");
  double temptau_pi   = 0.01;
  tempinput = util->StringFind3(file, "Shear_relaxation_time_tau_pi");
  if(tempinput != "empty") istringstream ( tempinput ) >> temptau_pi  ;
  DATA->tau_pi   = temptau_pi;
  
  
  // Bulk_relaxation_time_tau_b_pi:  in fm?
//   DATA->tau_b_pi = util->DFind(file, "Bulk_relaxation_time_tau_b_pi");
  double temptau_b_pi   = 0.6;
  tempinput = util->StringFind3(file, "Bulk_relaxation_time_tau_b_pi");
  if(tempinput != "empty") istringstream ( tempinput ) >> temptau_b_pi  ;
  DATA->tau_b_pi   = temptau_b_pi;
  
  //Shear_to_S_ratio:  constant eta/s
//   DATA->shear_to_s = util->DFind(file, "Shear_to_S_ratio");
  double tempshear_to_s   = 0.6;
  tempinput = util->StringFind3(file, "Shear_to_S_ratio");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempshear_to_s  ;
  DATA->shear_to_s   = tempshear_to_s;
  int tempturn_on_shear = 0;
  if(tempinput != "empty") tempturn_on_shear = 1;
    
  // Include_Shear_Visc_Yes_1_No_0
//   DATA->turn_on_shear = util->IFind(file, "Include_Shear_Visc_Yes_1_No_0");
//   int tempturn_on_shear = 0;
  tempinput = util->StringFind3(file, "Include_Shear_Visc_Yes_1_No_0");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempturn_on_shear;
  DATA->turn_on_shear = tempturn_on_shear;
  
  // T_dependent_Shear_to_S_ratio:  if 1, ignore constant eta/s and use hard-coded T-dependent shear viscosity
//   DATA->T_dependent_shear_to_s = util->IFind(file, "T_dependent_Shear_to_S_ratio");
  int tempT_dependent_shear_to_s = 1;
  tempinput = util->StringFind3(file, "T_dependent_Shear_to_S_ratio");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempT_dependent_shear_to_s;
  DATA->T_dependent_shear_to_s = tempT_dependent_shear_to_s;
  
  // Include_deltaf:  
  // Looks like 0 sets delta_f=0, 1 uses standard quadratic ansatz,
  // and 2 is supposed to use p^(2-alpha)
//   DATA->include_deltaf = util->IFind(file, "Include_deltaf");
  int tempinclude_deltaf = 1;
  tempinput = util->StringFind3(file, "Include_deltaf");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempinclude_deltaf;
  DATA->include_deltaf = tempinclude_deltaf;
  
  
  // Do_FreezeOut_Yes_1_No_0
  // set to 0 to bypass freeze out surface finder
//   DATA->doFreezeOut = util->IFind(file, "Do_FreezeOut_Yes_1_No_0");
  int tempdoFreezeOut = 1;
  tempinput = util->StringFind3(file, "Do_FreezeOut_Yes_1_No_0");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempdoFreezeOut;
  DATA->doFreezeOut = tempdoFreezeOut;
  
  
  // sigma_0:  unknown
//   DATA->sigma0 = util->DFind(file, "sigma_0");
  double tempsigma0   = 0.4;
  tempinput = util->StringFind3(file, "sigma_0");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempsigma0  ;
  DATA->sigma0   = tempsigma0;
  
  
  // Initial_Distribution_Filename:  unknown
//   DATA->initName = util->StringFind(file,"Initial_Distribution_Filename");
  string tempinitName = "filename";
  tempinput = util->StringFind3(file, "Initial_Distribution_Filename");
  if(tempinput != "empty") tempinitName.assign(tempinput);
  DATA->initName.assign(tempinitName);
  
  
  
  /* initialize the metric, mostly plus */
  
  DATA->avgT=0.;
  DATA->nSteps=0;
  DATA->avgT2=0.;
  DATA->nSteps2=0;
  DATA->plasmaEvolutionTime=0.;
  DATA->plasmaEvolutionTime2=0.;
  
  DATA->gmunu = util->mtx_malloc(4,4);
  for(m=0; m<4; m++)
    {
      for(n=0; n<4; n++)
	{
	  if(m == n) (DATA->gmunu)[m][n] = 1.0;
	  else (DATA->gmunu)[m][n] = 0.0;
	  
	  if(m==0 && n==0) (DATA->gmunu)[m][n] *= -1.0;
	}
    }/* m */
  
 //These were added for writing the coordinates of collisions to file for
 //MARTINI to sample. -CFY 11/2/2010
 
//  DATA->Nbin_to_file = util->IFind(file, "Write_Nbin_to_File");
  int tempNbin_to_file = 0;
  tempinput = util->StringFind3(file, "Write_Nbin_to_File");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempNbin_to_file;
  DATA->Nbin_to_file = tempNbin_to_file;
 

 //
//  DATA->outputBinaryEvolution = util->IFind(file, "outputBinaryEvolution");
  int tempoutputBinaryEvolution = 0;
  tempinput = util->StringFind3(file, "outputBinaryEvolution");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempoutputBinaryEvolution;
  DATA->outputBinaryEvolution = tempoutputBinaryEvolution;
 
  // End MARTINI parameters


  // Set to 1 if initial condition is boost-invariant 
  int tempboost_invariant = 0;
  tempinput = util->StringFind3(file, "boost_invariant");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempboost_invariant;
  DATA->boost_invariant = tempboost_invariant;
  DATA->boostInvariant = DATA->boost_invariant; // for compatibility with code from MUSIC light fork

  //Check if the freeze-out surface ever reaches the edge of the grid in the x or y direction
  //Warning: this check is only made when freeze-out method 3 is used! (not implemented for the other freeze-out methods)
  //set to 0 to disable
  //set to 1 to just output an error message when this happens
  //set to 2 to stop the evolution and exit when this happens
  int tempcheck_FO3_at_boundary_xy = 2;
  tempinput = util->StringFind3(file, "check_FO3_at_boundary_xy");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempcheck_FO3_at_boundary_xy;
  DATA->check_FO3_at_boundary_xy = tempcheck_FO3_at_boundary_xy;
  if(0 < DATA->check_FO3_at_boundary_xy && 3 != DATA->freezeOutMethod) cerr << "check_FO3_at_boundary_xy only used when freeze_out_method = 3.  Ignoring.\n";

  //Check if the freeze-out surface ever reaches the edge of the grid in the eta direction
  //Warning: this check is only made when freeze-out method 3 is used! (not implemented for the other freeze-out methods)
  //set to 0 to disable
  //set to 1 to just output an error message when this happens
  //set to 2 to stop the evolution and exit when this happens
  int tempcheck_FO3_at_boundary_eta = 0;
  tempinput = util->StringFind3(file, "check_FO3_at_boundary_eta");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempcheck_FO3_at_boundary_eta;
  if (DATA->boost_invariant > 0) tempcheck_FO3_at_boundary_eta = 0;
  DATA->check_FO3_at_boundary_eta = tempcheck_FO3_at_boundary_eta;
  if(0 <  DATA->check_FO3_at_boundary_eta && 3 != DATA->freezeOutMethod) cerr << "check_FO3_at_boundary_eta only used when freeze_out_method = 3.  Ignoring.\n";


  //Make MUSIC output additionnal hydro information
  //0 for false (do not output), 1 for true
  bool tempoutput_hydro_debug_info = true;
  tempinput = util->StringFind3(file, "output_hydro_debug_info");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempoutput_hydro_debug_info;
  DATA->output_hydro_debug_info = tempoutput_hydro_debug_info;

  //The evolution is outputted every "output_evolution_every_N_timesteps" timesteps
  //Can't be modified from the input file for now, for safety.
  DATA->output_evolution_every_N_timesteps=10;
  
  //Make MUSIC output a C header file containing informations about the hydro parameters used
  //0 for false (do not output), 1 for true
  bool tempoutput_hydro_params_header = false;
  tempinput = util->StringFind3(file, "output_hydro_params_header");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempoutput_hydro_params_header;
  DATA->output_hydro_params_header = tempoutput_hydro_params_header;

  // QuestRevert_rho_shear_max: QuestRevert has condition rho_shear > rho_shear_max
  double tempQuestRevert_rho_shear_max = 0.1;
  tempinput = util->StringFind3(file, "QuestRevert_rho_shear_max");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempQuestRevert_rho_shear_max;
  DATA->QuestRevert_rho_shear_max = tempQuestRevert_rho_shear_max;

  // QuestRevert_rho_bulk_max: QuestRevert has condition rho_bulk > rho_bulk_max
  double tempQuestRevert_rho_bulk_max = 0.1;
  tempinput = util->StringFind3(file, "QuestRevert_rho_bulk_max");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempQuestRevert_rho_bulk_max;
  DATA->QuestRevert_rho_bulk_max = tempQuestRevert_rho_bulk_max;

  // QuestRevert_factor: defines aggressiveness of QuestRevert
  double tempQuestRevert_factor = 0.;
  tempinput = util->StringFind3(file, "QuestRevert_factor");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempQuestRevert_factor;
  DATA->QuestRevert_factor = tempQuestRevert_factor;

  // QuestRevert_prefactor: defines aggressiveness of QuestRevert
  double tempQuestRevert_prefactor = 300.;
  tempinput = util->StringFind3(file, "QuestRevert_prefactor");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempQuestRevert_prefactor;
  DATA->QuestRevert_prefactor = tempQuestRevert_prefactor;

  // QuestRevert_eps_factor: factor \sim tanh ( eps / eps_F0 * QuestRevert_eps_factor)
  double tempQuestRevert_eps_factor = 1.22149;
  tempinput = util->StringFind3(file, "QuestRevert_eps_factor");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempQuestRevert_eps_factor;
  DATA->QuestRevert_eps_factor = tempQuestRevert_eps_factor;

  // QuestRevert_epsilon_min: QuestRevert is applied if epsilon > epsilon_min
  double tempQuestRevert_epsilon_min = 39061.;
  tempinput = util->StringFind3(file, "QuestRevert_epsilon_min");
  if(tempinput != "empty") istringstream ( tempinput ) >> tempQuestRevert_epsilon_min;
  DATA->QuestRevert_epsilon_min = tempQuestRevert_epsilon_min;

  cout << "Done ReadInData2." << endl;
  delete util;
}/* ReadInData2 */
