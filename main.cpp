#include "mpi.h"
#include <stdio.h>
#include "util.h"
#include "grid.h"
#include "data.h"
#include "init.h"
#include "eos.h"
#include "freeze.h"
#include "evolve.h"

using namespace std;
void ReadInData(InitData *DATA, string file);


// main program
main(int argc, char *argv[])
{
  int rank;
  int size;
  Grid ***arena; 
  Grid ***Lneighbor; 
  Grid ***Rneighbor;
 
  string input_file;
  static InitData DATA;
  int ieta, ix, iy, flag;
  double f, tau;


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

  input_file = *(argv+1);


  ReadInData(&DATA, input_file);

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
      flag =  evolve->EvolveIt(&DATA, arena, Lneighbor, Rneighbor, size, rank); 
      
      tau = DATA.tau0 + DATA.tau_size; 
      
      FILE *t4_file;
      char* t4_name = "avgT.dat";
      t4_file = fopen(t4_name, "a");
      fprintf(t4_file,"time average: %f GeV and %f GeV", DATA.avgT/DATA.nSteps*hbarc, DATA.avgT2/DATA.nSteps2*hbarc );
      fclose(t4_file);
      
      cout << "average plasma T=" << DATA.avgT/DATA.nSteps*hbarc << " GeV." << endl; 
      cout << "average plasma+mixed T=" << DATA.avgT2/DATA.nSteps2*hbarc << " GeV." << endl; 
      //    delete evolve;
      //delete glauber;
      //delete init;
    }
  
  if (DATA.mode == 1 || DATA.mode == 3 || DATA.mode == 4 || DATA.mode >= 5)
    {
      //if (rank == 0) int ret = system("rm yptphiSpectra*");
      //  freeze-out
      Freeze *freeze;
      freeze = new Freeze();
      freeze->CooperFrye(DATA.particleSpectrumNumber, DATA.mode, &DATA, eos, size, rank);
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

