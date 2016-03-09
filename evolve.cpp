#include "evolve.h"
#include "util.h"
#include "data.h"
#include "grid.h"
#include "eos.h"
#include "reconst.h"
#include "advance.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

Evolve::Evolve(EOS *eosIn, int gsl_seed)
{
  eos = new EOS;
  eos = eosIn;
  grid = new Grid;
  reconst = new Reconst(eosIn, grid);
  util = new Util;
  advance = new Advance(eosIn, grid);
  u_derivative = new U_derivative();

  gsl_ran = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (gsl_ran, gsl_seed);
}

// destructor
Evolve::~Evolve()
{
  delete eos;
  delete reconst;
  delete grid;
  delete util;
  delete advance;
  delete u_derivative;
}

int Evolve::EvolveIt(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank)
{
  //cerr << "OK starting evolve..." << endl;
/* implement Kurganov-Tadmor */
 int ix, iy, ieta, nx, ny, neta, it, itmax, rk_flag, flag, cent_eta;
 double dt, tau0, tau, x;

 stringstream numIn;
 numIn << DATA->seed;
 FILE *pc_file;
 //char* pc_name = "cellsWithProblems.dat";
 //pc_file = fopen(pc_name, "w");
 string pc_name = "cellsWithProblems_"+numIn.str()+".dat";
 pc_file = fopen(pc_name.c_str(), "w");
 fprintf(pc_file,"tau, x, y, eta \n");
 fclose(pc_file);
 FILE *v2_file;
 //char* v2_name = "aniso.dat";
 //v2_file = fopen(v2_name, "w");
 string v2_name = "aniso_"+numIn.str()+".dat";
 v2_file = fopen(v2_name.c_str(), "w");
 fprintf(v2_file,"tau [fm], momentum anisotropy\n");
 fclose(v2_file);
 FILE *ecc_file;
 //char* ecc_name = "eccentricity.dat";
 //ecc_file = fopen(ecc_name, "w");
 string ecc_name = "eccentricity_"+numIn.str()+".dat";
 ecc_file = fopen(ecc_name.c_str(), "w");
 fprintf(ecc_file,"tau [fm], spatial anisotropy\n");
 fclose(ecc_file);
 FILE *t_file;
 //char* t_name = "tauf.dat";
 //t_file = fopen(t_name, "w");
 string t_name = "tauf_"+numIn.str()+".dat";
 t_file = fopen(t_name.c_str(), "w");
 fprintf(t_file,"x_f [fm], tau_f [fm], Sx[fm^3], Stau [fm^3] \n");
 fclose(t_file);
 FILE *t2_file;
 //char* t2_name = "taufx.dat";
 //t2_file = fopen(t2_name, "w");
 string t2_name = "taufx_"+numIn.str()+".dat";
 t2_file = fopen(t2_name.c_str(), "w");
 fprintf(t2_file,"x_f [fm], tau_f [fm], Sx[fm^3], Stau [fm^3] \n");
 fclose(t2_file);
 FILE *t3_file;
 //char* t3_name = "taufy.dat";
 //t3_file = fopen(t3_name, "w");
 string t3_name = "taufy_"+numIn.str()+".dat";
 t3_file = fopen(t3_name.c_str(), "w");
 fprintf(t3_file,"y_f [fm], tau_f [fm], Sy[fm^3], Stau [fm^3] \n");
 fclose(t3_file);

 stringstream rankIn;
 rankIn << rank;

 //char *buf;
 //buf = util->char_malloc(40);
 FILE *s_file;
 //char* s_name;
 //s_name = util->char_malloc(100);
 
 //sprintf (buf, "%d", rank);
 //strcat(s_name, "surface");
 //strcat(s_name,buf);
 //strcat(s_name, ".dat");
 
 string s_name = "surface"+rankIn.str()+"_"+numIn.str()+".dat";
 //s_file = fopen(s_name, "w");
 s_file = fopen(s_name.c_str(), "w");
 fclose(s_file);

 FILE *out_file;
 //char* out_name = "evolution.dat";
 //out_file = fopen(out_name, "w");
 string out_name = "evolution_"+numIn.str()+".dat";
 out_file = fopen(out_name.c_str(), "w");
 fprintf(out_file,"");
 fclose(out_file);

 FILE *out2_file;
 //char* out2_name = "entropy.dat";
 //out2_file = fopen(out2_name, "w");
 string out2_name = "entropy_"+numIn.str()+".dat";
 out2_file = fopen(out2_name.c_str(), "w");
 fprintf(out2_file,"");
 fclose(out2_file);

 FILE *cout_file;
 //char* cout_name = "contourPlot.dat";
 //cout_file = fopen(cout_name, "w");
 string cout_name = "contourPlot_"+numIn.str()+".dat";
 cout_file = fopen(cout_name.c_str(), "w");
 fprintf(cout_file,"");
 fclose(cout_file);

 facTau = DATA->facTau;

 fprintf(stderr, "Starting Evolve on rank %d...\n", rank);
 
 itmax = DATA->nt;
 tau0 = DATA->tau0;
 dt = DATA->delta_tau;
 weirdCases=0;
 SUM = 0.;
 SUM2 = 0.;
 warnings = 0;
 cells = 0;
       
 
 for(it=0; it<=itmax; it++)
   {
   //cout << "OK at it = " << it << endl;
   tau = tau0 + dt*it;
   //fprintf(stderr, "Starting time step %d/%d on rank %d.\n", it, itmax, rank);
   if(it==0) 
     {
       // storePreviousT(tau, DATA, arena);
       storePreviousEpsilon2(tau, DATA, arena);
       storePreviousW(tau, DATA, arena);
     }

   //storePreviousEpsilon(tau, DATA, arena);
/*    //for testing */
/*    FindFreezeOutSurface(tau, DATA, arena); */
/*    sleep(1); */
/*    FindFreezeOutSurface2(tau, DATA, arena); */
/*    exit(1); */
   //for testing


   if(it%10==0 && it>=0) 
     {
       grid->PrintEtaEpsilon(arena, DATA, tau, size, rank);
       grid->PrintxEpsilon(arena, DATA, tau, size, rank);
       grid->ComputeEccentricity(DATA, arena, tau);
       grid->ComputeAnisotropy(DATA, arena, eos, tau);
       grid->OutputXY(arena, DATA, eos, tau, size, rank);
     }

   //cerr << "OK after PrintEta..." << endl; fflush(stdout);
   //cout << DATA->outputEvolutionData << endl; fflush(stdout);

   if(it%10==0 && DATA->outputEvolutionData) 
   //if(DATA->outputEvolutionData) 
     {
       //cout << "OK starting outputting..." << endl; fflush(stdout);
       //cout << DATA->outputEvolutionData << endl; fflush(stdout);

       grid->OutputEvolutionDataXYZ(arena, DATA, eos, tau, size, rank); 
       grid->OutputEvolutionDataXYEta(arena, DATA, eos, tau, size, rank); 
       //grid->OutputPlotDataXYZ(arena, DATA, eos, tau, size, rank); 
       if(DATA->fluctuatingHydroFlag==2 && it%10==0){
	 grid->OutputFluctuatingPlotDataXYEta(arena, DATA, eos, tau, size, rank);
       }
       // this produces potentially huge outputs so beware
       //cerr << "OK stopping outputting." << endl;
     }
 

   /* execute rk steps */
  
   //cerr << "Finished outputting successfully, starting AdvanceRK..." << endl;

   flag = AdvanceRK(tau, DATA, arena, Lneighbor, Rneighbor, size, rank);

   //cerr << "Advanced successfully, starting UpdateArena..." << endl;

   int itPassed = it;
   UpdateArena(tau, DATA, arena, itPassed);
   
   //cerr << "Arena updated..." << endl;

   //check energy conservation
   //grid->ComputeEnergyConservation(DATA, arena, tau);
   //storePreviousT(tau, DATA, arena);

   //determine freeze-out surface
  
  if(DATA->doFreezeOut == 1)
   {
    if (it%facTau==0 && it>0) 
      {
	if (DATA->freezeOutMethod == 1)
	  FindFreezeOutSurface(tau, DATA, arena, size, rank);
	else if ((DATA->freezeOutMethod == 2 || DATA->freezeOutMethod == 4) && DATA->Initial_profile != 7)
	  FindFreezeOutSurface2(tau, DATA, arena, size, rank);
	else if (DATA->freezeOutMethod == 5)
	  FindFreezeOutSurfaceOfHypercubes(tau, DATA, arena, size, rank);
	storePreviousEpsilon2(tau, DATA, arena);
	storePreviousW(tau, DATA, arena);
      } 
   }/* do freeze-out determination */
    
    if (rank == 0) fprintf(stderr, "Done time step %d/%d.\n", it, itmax);
    
   }/* it */ 

 //Isochronal freeze-out:
 if(DATA->doFreezeOut == 1 && DATA->freezeOutMethod == 3){
   FindIsochronalFreezeOutSurface(DATA->tau0 + DATA->tau_size, DATA, arena, size, rank);
 }

 if(rank == 0)
 {
  grid->PrintAxy2(DATA, arena, tau);
 }
 
 fprintf(stderr,"SUM=%f\n", SUM);
 return 1; /* successful */

}/* Evolve */




// faster one, needs to be called only every facTau time steps
void Evolve::storePreviousEpsilon2(double tau, InitData *DATA, Grid ***arena)
{
  int ix, iy, ieta, nx, ny, neta;
  double x, y, eta, tau0;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;
  tau0 = DATA->tau0;

  for(ix=0; ix<=nx; ix++)
    {
      for(iy=0; iy<=ny; iy++)
	{
	  for(ieta=0; ieta<neta; ieta++)
	    {
	      //cout << ix << " " << iy << " " << ieta << " " << " " << arena[ix][iy][ieta].u[0][1] << endl;
	      arena[ix][iy][ieta].epsilon_prev = arena[ix][iy][ieta].epsilon;
	      //arena[ix][iy][ieta].epsilon_prev[facTau-1]=arena[ix][iy][ieta].epsilon;
	      arena[ix][iy][ieta].u_prev[0] = arena[ix][iy][ieta].u[0][0];
	      arena[ix][iy][ieta].u_prev[1] = arena[ix][iy][ieta].u[0][1];
	      arena[ix][iy][ieta].u_prev[2] = arena[ix][iy][ieta].u[0][2];
	      arena[ix][iy][ieta].u_prev[3] = arena[ix][iy][ieta].u[0][3];
	      arena[ix][iy][ieta].deltaP_prev = arena[ix][iy][ieta].deltaP[0];
	      arena[ix][iy][ieta].deltaU_prev[0] = arena[ix][iy][ieta].deltaU[0][0];
	      arena[ix][iy][ieta].deltaU_prev[1] = arena[ix][iy][ieta].deltaU[0][1];
	      arena[ix][iy][ieta].deltaU_prev[2] = arena[ix][iy][ieta].deltaU[0][2];
	      arena[ix][iy][ieta].deltaU_prev[3] = arena[ix][iy][ieta].deltaU[0][3];
	      arena[ix][iy][ieta].deltaW_prev[0][0] = arena[ix][iy][ieta].deltaW[0][0][0];
	      arena[ix][iy][ieta].deltaW_prev[0][1] = arena[ix][iy][ieta].deltaW[0][0][1];
	      arena[ix][iy][ieta].deltaW_prev[0][2] = arena[ix][iy][ieta].deltaW[0][0][2];
	      arena[ix][iy][ieta].deltaW_prev[0][3] = arena[ix][iy][ieta].deltaW[0][0][3];
	      arena[ix][iy][ieta].deltaW_prev[1][1] = arena[ix][iy][ieta].deltaW[0][1][1];
	      arena[ix][iy][ieta].deltaW_prev[1][2] = arena[ix][iy][ieta].deltaW[0][1][2];
	      arena[ix][iy][ieta].deltaW_prev[1][3] = arena[ix][iy][ieta].deltaW[0][1][3];
	      arena[ix][iy][ieta].deltaW_prev[2][2] = arena[ix][iy][ieta].deltaW[0][2][2];
	      arena[ix][iy][ieta].deltaW_prev[2][3] = arena[ix][iy][ieta].deltaW[0][2][3];
	      arena[ix][iy][ieta].deltaW_prev[3][3] = arena[ix][iy][ieta].deltaW[0][3][3];
	      //arena[ix][iy][ieta].u_prev[facTau-1][0]=arena[ix][iy][ieta].u[0][0];
	      //arena[ix][iy][ieta].u_prev[facTau-1][1]=arena[ix][iy][ieta].u[0][1];
	      //arena[ix][iy][ieta].u_prev[facTau-1][2]=arena[ix][iy][ieta].u[0][2];
	      //arena[ix][iy][ieta].u_prev[facTau-1][3]=arena[ix][iy][ieta].u[0][3];
	    }
	}
    }
}

// faster one, needs to be called only every facTau time steps
void Evolve::storePreviousW(double tau, InitData *DATA, Grid ***arena)
{
  int ix, iy, ieta, nx, ny, neta;
  double x, y, eta, tau0;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;
  tau0 = DATA->tau0;

  for(ix=0; ix<=nx; ix++)
    {
      for(iy=0; iy<=ny; iy++)
	{
	  for(ieta=0; ieta<neta; ieta++)
	    {
	      arena[ix][iy][ieta].W_prev[0][0]=arena[ix][iy][ieta].Wmunu[0][0][0];
	      arena[ix][iy][ieta].W_prev[0][1]=arena[ix][iy][ieta].Wmunu[0][0][1];
	      arena[ix][iy][ieta].W_prev[0][2]=arena[ix][iy][ieta].Wmunu[0][0][2];
	      arena[ix][iy][ieta].W_prev[0][3]=arena[ix][iy][ieta].Wmunu[0][0][3];
	      arena[ix][iy][ieta].W_prev[1][1]=arena[ix][iy][ieta].Wmunu[0][1][1];
	      arena[ix][iy][ieta].W_prev[1][2]=arena[ix][iy][ieta].Wmunu[0][1][2];
	      arena[ix][iy][ieta].W_prev[1][3]=arena[ix][iy][ieta].Wmunu[0][1][3];
	      arena[ix][iy][ieta].W_prev[2][2]=arena[ix][iy][ieta].Wmunu[0][2][2];
	      arena[ix][iy][ieta].W_prev[2][3]=arena[ix][iy][ieta].Wmunu[0][2][3];
	      arena[ix][iy][ieta].W_prev[3][3]=arena[ix][iy][ieta].Wmunu[0][3][3];
	    }
	}
    }
}

void Evolve::storePreviousT(double tau, InitData *DATA, Grid ***arena)
{
  int ix, iy, ieta, nx, ny, neta;
  double x, y, eta, tau0;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;
  tau0 = DATA->tau0;

  for(ix=0; ix<=nx; ix++)
    {
      for(iy=0; iy<=ny; iy++)
	{
	  for(ieta=0; ieta<neta; ieta++)
	    {
	      arena[ix][iy][ieta].prev_T00=arena[ix][iy][ieta].TJb[0][0][0];
	      arena[ix][iy][ieta].prev_T33=arena[ix][iy][ieta].TJb[0][3][3];
	    }
	}
    }
}

int Evolve::UpdateArena(double tau, InitData *DATA, Grid ***arena, int it)
{
 int rk_flag, ix, iy, ieta, nx, ny, neta, flag, alpha, mu, rk_order, nu;
 double tempd;
 Grid *grid_pt;
 
 nx = DATA->nx;
 ny = DATA->ny;
 neta = DATA->neta-1;
 rk_order = DATA->rk_order;

 //The width needed for the sampling values for S. Because UpdateArena is called after AdvanceRK, using the same value for tau, 
 //dtau*(it+1)+tau_0 is the appropriate time to use:
 double sigmaV;
 if(DATA->Initial_profile != 7){
   sigmaV = sqrt(DATA->delta_tau * DATA->delta_x * DATA->delta_y * (tau + 1.5*DATA->delta_tau) * DATA->delta_eta * (double)(DATA->fluctuatingGridFactorTau) );
 }
 else{
   sigmaV = sqrt(DATA->delta_tau * DATA->delta_x * DATA->delta_x * DATA->delta_y );
 }

 //cerr << "sigmaV = " << sigmaV << endl;

 for(ix=0; ix<=nx; ix++)
  {
   for(iy=0; iy<=ny; iy++)
    {
     for(ieta=0; ieta<=neta; ieta++)
      {
       arena[ix][iy][ieta].prev_epsilon = arena[ix][iy][ieta].epsilon;
       arena[ix][iy][ieta].prev_p = arena[ix][iy][ieta].p;

       arena[ix][iy][ieta].p = arena[ix][iy][ieta].p_t;
       arena[ix][iy][ieta].epsilon = arena[ix][iy][ieta].epsilon_t;
       arena[ix][iy][ieta].rhob = arena[ix][iy][ieta].rhob_t;
       arena[ix][iy][ieta].T = arena[ix][iy][ieta].T_t;

       if(DATA->fluctuatingHydroFlag == 2){
	 arena[ix][iy][ieta].prev_deltaP = arena[ix][iy][ieta].deltaP[0];
	 arena[ix][iy][ieta].deltaP[0] = arena[ix][iy][ieta].deltaP[rk_order];
       }

       // /* this was the previous previous value */
       //arena[ix][iy][ieta].pprev_pi_b[0] = arena[ix][iy][ieta].prev_pi_b[0];
       ///* this was the previous value */
       //arena[ix][iy][ieta].prev_pi_b[0] = arena[ix][iy][ieta].pi_b[0];
       ///* this is the new value */
       //arena[ix][iy][ieta].pi_b[0] = arena[ix][iy][ieta].pi_b[rk_order];
      
       for(mu=0; mu<4; mu++)
        {
	 /* this was the previous previous value */
	 arena[ix][iy][ieta].pprev_u[0][mu] = 
	               arena[ix][iy][ieta].prev_u[0][mu]; 
	 
	 /* this was the previous value */
	 arena[ix][iy][ieta].prev_u[0][mu] = 
	               arena[ix][iy][ieta].u[0][mu]; 
	 
	 /* this is the new value */
	 arena[ix][iy][ieta].u[0][mu] = 
	               arena[ix][iy][ieta].u[rk_order][mu]; 
	 
	 if(DATA->fluctuatingHydroFlag == 1){
	   //This was the previous value:
	   arena[ix][iy][ieta].prev_Xi0[mu] = arena[ix][iy][ieta].Xi[0][0][mu];
	 }

	 if(DATA->fluctuatingHydroFlag == 2){
	   arena[ix][iy][ieta].prev_deltaU[mu] = arena[ix][iy][ieta].deltaU[0][mu];
	   arena[ix][iy][ieta].deltaU[0][mu] = arena[ix][iy][ieta].deltaU[rk_order][mu];
	 }

	 for(alpha=0; alpha<5; alpha++){
	   if(DATA->fluctuatingHydroFlag == 1 && alpha < 4){
	     //This is the new value:
	     arena[ix][iy][ieta].Xi[0][alpha][mu] = arena[ix][iy][ieta].Xi[rk_order][alpha][mu];
	   }

	   if(DATA->fluctuatingHydroFlag == 2 && alpha < 4){
	     arena[ix][iy][ieta].prev_deltaW[alpha][mu] = arena[ix][iy][ieta].deltaW[0][alpha][mu];
	     arena[ix][iy][ieta].deltaW[0][alpha][mu] = arena[ix][iy][ieta].deltaW[rk_order][alpha][mu];
	     arena[ix][iy][ieta].deltaT[0][alpha][mu] = arena[ix][iy][ieta].deltaT[rk_order][alpha][mu];
	     //arena[ix][iy][ieta].prev_Xi[alpha][mu] = arena[ix][iy][ieta].Xi[0][alpha][mu];
	     //arena[ix][iy][ieta].Xi[0][alpha][mu] = arena[ix][iy][ieta].Xi[rk_order][alpha][mu];
	   }

	   /* this is the new value */
	   arena[ix][iy][ieta].TJb[0][alpha][mu] = 
	     arena[ix][iy][ieta].TJb[rk_order][alpha][mu]; 
	   
	   /* this was the previous previous value */
	   arena[ix][iy][ieta].pprevWmunu[0][alpha][mu] = 
	     arena[ix][iy][ieta].prevWmunu[0][alpha][mu]; 
	   
	   /* this was the previous value */
	   arena[ix][iy][ieta].prevWmunu[0][alpha][mu] = 
	     arena[ix][iy][ieta].Wmunu[0][alpha][mu]; 
	   
	   /* this is the new value */
	   arena[ix][iy][ieta].Wmunu[0][alpha][mu] = 
	     arena[ix][iy][ieta].Wmunu[rk_order][alpha][mu]; 
	   
	   //if(isnan(arena[ix][iy][ieta].Wmunu[rk_order][alpha][mu]))
	   //    cout << "updateArena Wmunu[" << ix << "][" << iy << "][" << ieta << "],[" 
	   // << 0 << "][alpha=" << alpha << "][mu=" << "]=" << arena[ix][iy][ieta].Wmunu[0][alpha][mu] << endl;
	   
	   ///* this was the previous previous value */
	   //arena[ix][iy][ieta].pprevPimunu[0][alpha][mu] = 
	   //                arena[ix][iy][ieta].prevPimunu[0][alpha][mu]; 
	   
	   ///* this was the previous value */
	   //arena[ix][iy][ieta].prevPimunu[0][alpha][mu] = 
	   //                arena[ix][iy][ieta].Pimunu[0][alpha][mu]; 
	   
	   ///* this is the new value */
	   //arena[ix][iy][ieta].Pimunu[0][alpha][mu] = 
	   //          arena[ix][iy][ieta].Pimunu[rk_order][alpha][mu]; 
	 }
	}/* mu, alpha */
       
       //Update the sampling of the Wiener process. It should only be called every fluctuatingGridFactorTau'th timestep:
       //cout << "At (" << ix << ", " << iy << ", " << ieta << "):" << endl;
       if(it%(DATA->fluctuatingGridFactorTau) == 0 && DATA->fluctuatingHydroFlag != 0){
	 arena[ix][iy][ieta].updateS(gsl_ran, sigmaV, DATA);
	 //cout << "it = " << it << endl;
       }
       
       //Should not be commented out with new line for Xi above:
       ////Finally, determine Xi[0][mu][nu]:
       //arena[ix][iy][ieta].getXi(0, util, eos, DATA);
      }
    }
  }/* ix, iy, ieta */

 //cerr << "UpdateArena ended successfully..." << endl;

}/* UpdateArena */


int Evolve::AdvanceRK(double tau, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank)
{
 int rk_flag, ix, iy, ieta, nx, ny, neta, flag;

for(rk_flag = 0; rk_flag < DATA->rk_order; rk_flag++)
{
  // cout << "1 AdvanceRK Wmunu=" << (Lneighbor[1][1][0]).Wmunu[rk_flag][1][1] << endl;
  // advance->MPISendReceive(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);
  //cout << "2 AdvanceRK Wmunu=" << (Lneighbor[1][1][0]).Wmunu[rk_flag][1][1] << endl;
       
  advance->MPISendReceiveT(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag); 
  advance->MPISendReceiveW(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);
  
  //cerr << "OK up to here" << endl;

  if(DATA->fluctuatingHydroFlag == 1){
    advance->MPISendReceiveXi(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);
  }
  if(DATA->fluctuatingHydroFlag == 2){
    advance->MPISendReceiveXi2(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);
  }

  //cerr << "OK sending and receiving Xi" << endl;

  flag = u_derivative->MakedU(tau, DATA, arena, Lneighbor, Rneighbor, rk_flag, size, rank); 
  if(flag == 0) return 0;
  
  //cerr << "OK making dU" << endl;

  // advance->MPISendReceive(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);

  //cerr << "OK starting AdvanceIt..." << endl;
  flag = advance->AdvanceIt(tau, DATA, arena, Lneighbor, Rneighbor, rk_flag, size, rank);
  //cerr << "...OK finishing AdvanceIt." << endl;
  
  //AdvanceIt(tau, DATA, arena, rk_flag, size, rank); 
 }/* loop over rk_flag */
 
 if(flag == 0) return 0;
 
 return 1; /* successful */

}/* AdvanceRK1 *///      cout << "going through Advance rk_flag=" << rk_flag << " rank " << rank << endl;
      

void Evolve::FindFreezeOutSurface(double tau, InitData *DATA, Grid ***arena, int size, int rank)
{	
  double epsFO=DATA->epsilonFreeze/hbarc;
  double PFO=eos->p_func(epsFO, 0.);
  double sFO=eos->s_func(epsFO, PFO, 0., DATA);

  stringstream numIn;
  numIn << DATA->seed;
  stringstream rankIn;
  rankIn << rank;
  FILE *t_file;
  //char* t_name = "tauf.dat";
  //t_file = fopen(t_name, "a");
  string t_name = "tauf_"+numIn.str()+".dat";
  t_file = fopen(t_name.c_str(), "a");
  //char *buf;
  //buf = util->char_malloc(40);
  FILE *s_file;
  //char* s_name;
  string s_name = "surface"+rankIn.str()+"_"+numIn.str()+".dat";
  //s_name = util->char_malloc(100);
  //sprintf (buf, "%d", rank);
  //strcat(s_name, "surface");
  //strcat(s_name,buf);
  //strcat(s_name, ".dat");

  //s_file = fopen(s_name, "a");
  s_file = fopen(s_name.c_str(), "a");
  //fprintf(t_file,"x [fm], tau_f(x) [fm] \n");
  int ix, iy, ieta, nx, ny, neta;
  double x, y, eta;
  //double epsFO=DATA->epsilonFreeze/hbarc;
  double tauf, xf, yf, etaf;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;
  double FULLSU[4];
  int fac, intersect, intersectx, intersecty, intersecteta, intersecttau;
  double DX, DY, DETA, DTAU, SIG;
  double rhob, utau, ux, uy, ueta, TFO, muB;
  double utauX1, utauX2, utauX3, utauX4, utauY1, utauY2, utau1, utau2;
  double rhobX1, rhobX2, rhobX3, rhobX4, rhobY1, rhobY2, rhob1, rhob2;
  double uxX1, uxX2, uxX3, uxX4, uxY1, uxY2, ux1, ux2;
  double uyX1, uyX2, uyX3, uyX4, uyY1, uyY2, uy1, uy2;
  double uetaX1, uetaX2, uetaX3, uetaX4, uetaY1, uetaY2, ueta1, ueta2;
  double xfrac, yfrac, etafrac, taufrac;
  int shown;
  
  fac=1;
  facTau = DATA->facTau;
  DX=fac*DATA->delta_x;
  DY=fac*DATA->delta_y;
  DETA=fac*DATA->delta_eta;
  DTAU=facTau*DATA->delta_tau;
  
  fprintf(stderr,"DTAU=%f\n", DTAU);
  fprintf(stderr,"DX=%f\n", DX);
  fprintf(stderr,"DY=%f\n", DY);
  fprintf(stderr,"DETA=%f\n", DETA);
  shown = 0;

  int maxEta;
  int sizeOfData = (nx+1)*(ny+1);
  int position;
  double *package;
  double Rneighbor_eps[nx+1][ny+1];

  package = new double[sizeOfData];

  // receive from the right / send to the left
  int from = rank+1;
  int to = rank-1;

  // get cells from neighboring processors
  if ( rank != 0 )
    {
      //      cout << " sending to the left on rank " << rank << endl;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      //if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
	      //cout << "arena[ix][iy][0].TJb[i][alpha][0]=" << arena[ix][iy][0].TJb[i][alpha][0] << endl;
	      position = ix + (nx*iy);
	      package[position] = arena[ix][iy][0].epsilon;
	    }
	}
      MPI::COMM_WORLD.Send(package,sizeOfData,MPI::DOUBLE,to,1);
      //cout << " done sending to the left on rank " << rank << endl;
    }
  // receiving and unwrapping the package
  if ( rank != size-1 )
    {  
      MPI::COMM_WORLD.Recv(package,sizeOfData,MPI::DOUBLE,from,1);
      //cout << " receiving from the right on rank " << rank << endl;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      position = ix + (nx*iy);
	      //if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
	      //cout << "Rneighbor[ix][iy][0].TJb[i][alpha][0]=" << package[position] << endl;
	      Rneighbor_eps[ix][iy] = package[position];
	    }
	}
      //      cout << " done receiving from the right on rank " << rank << endl;
    }

  if (rank == size-1) maxEta = neta-fac;
  else maxEta = neta;

  for(ix=0; ix<nx; ix+=fac)
    {
      x = ix*(DATA->delta_x) - (DATA->x_size/2.0); 
      for(iy=0; iy<ny; iy+=fac)
	{
	  y = iy*(DATA->delta_y) - (DATA->y_size/2.0);
	  for(ieta=0; ieta<maxEta; ieta+=fac)
	    {
	      eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
	      intersect=0;
	      intersectx=0;
	      intersecty=0;
	      intersecteta=0;
	      intersecttau=0;
	      FULLSU[0]=0.;
	      FULLSU[1]=0.;
	      FULLSU[2]=0.;
	      FULLSU[3]=0.;

	      xf = x;
	      yf = y;
	      etaf = eta;
	      tauf = tau-DTAU;

	      if((arena[ix+fac][iy][ieta].epsilon-epsFO)*(arena[ix][iy][ieta].epsilon-epsFO)<0.)
		{
		  if (arena[ix+fac][iy][ieta].epsilon>epsFO)
		    SIG=-1.;
		  else SIG=1.;
		  intersect = 1;
		  intersectx = 1;
		  xf = x + DX * (arena[ix][iy][ieta].epsilon-epsFO)/(arena[ix][iy][ieta].epsilon-arena[ix+fac][iy][ieta].epsilon);
		  FULLSU[1]+=SIG*DTAU*DY*DETA;
		}
	      if((arena[ix][iy+fac][ieta].epsilon-epsFO)*(arena[ix][iy][ieta].epsilon-epsFO)<0.)
		{
		  if (arena[ix][iy+fac][ieta].epsilon>epsFO)
		    SIG=-1.;
		  else SIG=1.;
		  intersect = 1;
		  intersecty = 1;
		  yf = y + DY * (arena[ix][iy][ieta].epsilon-epsFO)/(arena[ix][iy][ieta].epsilon-arena[ix][iy+fac][ieta].epsilon);
		  FULLSU[2]+=SIG*DTAU*DX*DETA;
		}
	      //if((arena[ix][iy][ieta].epsilon-epsFO)*(arena[ix][iy][ieta].epsilon_prev[facTau-1]-epsFO)<0.)
	       if((arena[ix][iy][ieta].epsilon-epsFO)*(arena[ix][iy][ieta].epsilon_prev-epsFO)<0.)
		{
		  if (arena[ix][iy][ieta].epsilon>epsFO)
		    SIG=-1.;
		  else SIG=1.;
		  intersect = 1;
		  intersecttau = 1;
		  //tauf = tau - DTAU + DTAU * (arena[ix][iy][ieta].epsilon_prev[facTau-1]-epsFO)
		  //  /(arena[ix][iy][ieta].epsilon_prev[facTau-1]-arena[ix][iy][ieta].epsilon);
		  tauf = tau - DTAU + DTAU * (arena[ix][iy][ieta].epsilon_prev-epsFO)
		    /(arena[ix][iy][ieta].epsilon_prev-arena[ix][iy][ieta].epsilon);
		  FULLSU[0]+=SIG*DX*DY*DETA;
		}
	       if (ieta<neta-fac)
		 {
		   if((arena[ix][iy][ieta+fac].epsilon-epsFO)*(arena[ix][iy][ieta].epsilon-epsFO)<0.)
		     {
		       if (arena[ix][iy][ieta+fac].epsilon>epsFO)
			 SIG=-1.;
		       else SIG=1.;
		       intersect = 1;
		       intersecteta = 1;
		       etaf = eta + DETA * (arena[ix][iy][ieta].epsilon-epsFO)/(arena[ix][iy][ieta].epsilon-arena[ix][iy][ieta+fac].epsilon);
		       FULLSU[3]+=SIG*DTAU*DX*DY;
		     }
		 }
	       else
		 {
		   if((Rneighbor_eps[ix][iy]-epsFO)*(arena[ix][iy][ieta].epsilon-epsFO)<0.)
		     {
		       if (Rneighbor_eps[ix][iy]>epsFO)
			 SIG=-1.;
		       else SIG=1.;
		       intersect = 1;
		       intersecteta = 1;
		       etaf = eta + DETA * (arena[ix][iy][ieta].epsilon-epsFO)/(arena[ix][iy][ieta].epsilon-Rneighbor_eps[ix][iy]);
		       FULLSU[3]+=SIG*DTAU*DX*DY;
		     }
		 }
	       
	    
	       
	      if (intersect == 0 ) continue;

	      SUM += fabs(FULLSU[0])+fabs(FULLSU[1])+fabs(FULLSU[2])+fabs(FULLSU[3]);

	      xfrac = (xf-x)/(DX);
	      yfrac = (yf-y)/(DY);
	      etafrac = (etaf-eta)/(DETA);
	      taufrac = (tauf-(tau-DTAU))/(DTAU);

/* 	      fprintf(stderr,"before interpolation %d %d %d\n", ix, iy, ieta); */
	      //    fprintf(stderr,"%f %f %f %f\n", taufrac, xfrac, yfrac, etafrac);
	      //fprintf(stderr,"%f %f\n", tau, tauf);
	  

  
              // 4d interpolation:
	      //utau:
	      utauX1 = ((1.-xfrac)*arena[ix][iy][ieta].u[0][0] + xfrac*arena[ix+fac][iy][ieta].u[0][0]);
	      utauX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u[0][0] + xfrac*arena[ix+fac][iy+fac][ieta].u[0][0]);
	      utauX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u[0][0] + xfrac*arena[ix+fac][iy][ieta+fac].u[0][0]);
	      utauX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u[0][0] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u[0][0]);
	      
	      utauY1 = ((1.-yfrac)*utauX1+yfrac*utauX2);
	      utauY2 = ((1.-yfrac)*utauX3+yfrac*utauX4);
	      
	      utau1 = ((1.-etafrac)*utauY1+etafrac*utauY2);
	      
	      //utauX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][0]);
	      //utauX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][0]);
	      //utauX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][0]);
	      //utauX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][0]);
	      utauX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[0] + xfrac*arena[ix+fac][iy][ieta].u_prev[0]);
	      utauX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[0] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[0]);
	      utauX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[0] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[0]);
	      utauX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[0] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[0]);
	      	      
	      utauY1 = ((1.-yfrac)*utauX1+yfrac*utauX2);
	      utauY2 = ((1.-yfrac)*utauX3+yfrac*utauX4);
	      
	      utau2 = ((1.-etafrac)*utauY1+etafrac*utauY2);
	      
	      utau = (1.-taufrac)*utau2+taufrac*utau1;
	      
	      //ux:
	      uxX1 = ((1.-xfrac)*arena[ix][iy][ieta].u[0][1] + xfrac*arena[ix+fac][iy][ieta].u[0][1]);
	      uxX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u[0][1] + xfrac*arena[ix+fac][iy+fac][ieta].u[0][1]);
	      uxX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u[0][1] + xfrac*arena[ix+fac][iy][ieta+fac].u[0][1]);
	      uxX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u[0][1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u[0][1]);
	      
	      uxY1 = ((1.-yfrac)*uxX1+yfrac*uxX2);
	      uxY2 = ((1.-yfrac)*uxX3+yfrac*uxX4);
	      
	      ux1 = ((1.-etafrac)*uxY1+etafrac*uxY2);
	      
	      //uxX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][1]);
	      //uxX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][1]);
	      //uxX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][1]);
	      //uxX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][1]);
	      uxX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[1] + xfrac*arena[ix+fac][iy][ieta].u_prev[1]);
	      uxX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[1] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[1]);
	      uxX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[1] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[1]);
	      uxX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[1]);
	      
	      uxY1 = ((1.-yfrac)*uxX1+yfrac*uxX2);
	      uxY2 = ((1.-yfrac)*uxX3+yfrac*uxX4);
	      
	      ux2 = ((1.-etafrac)*uxY1+etafrac*uxY2);
	      
	      ux = (1.-taufrac)*ux2+taufrac*ux1;
	      
	      //uy:
	      uyX1 = ((1.-xfrac)*arena[ix][iy][ieta].u[0][2] + xfrac*arena[ix+fac][iy][ieta].u[0][2]);
	      uyX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u[0][2] + xfrac*arena[ix+fac][iy+fac][ieta].u[0][2]);
	      uyX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u[0][2] + xfrac*arena[ix+fac][iy][ieta+fac].u[0][2]);
	      uyX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u[0][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u[0][2]);
	      
	      uyY1 = ((1.-yfrac)*uyX1+yfrac*uyX2);
	      uyY2 = ((1.-yfrac)*uyX3+yfrac*uyX4);
	      
	      uy1 = ((1.-etafrac)*uyY1+etafrac*uyY2);
	      
	      //uyX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][2]);
	      //uyX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][2]);
	      //uyX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][2]);
	      //uyX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][2]);
	      uyX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[2] + xfrac*arena[ix+fac][iy][ieta].u_prev[2]);
	      uyX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[2] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[2]);
	      uyX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[2] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[2]);
	      uyX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[2]);
	      
	      uyY1 = ((1.-yfrac)*uyX1+yfrac*uyX2);
	      uyY2 = ((1.-yfrac)*uyX3+yfrac*uyX4);
	      
	      uy2 = ((1.-etafrac)*uyY1+etafrac*uyY2);
	      
	      uy = (1.-taufrac)*uy2+taufrac*uy1;
	      
	      
	      //ueta:
	      uetaX1 = ((1.-xfrac)*arena[ix][iy][ieta].u[0][3] + xfrac*arena[ix+fac][iy][ieta].u[0][3]);
	      uetaX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u[0][3] + xfrac*arena[ix+fac][iy+fac][ieta].u[0][3]);
	      uetaX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u[0][3] + xfrac*arena[ix+fac][iy][ieta+fac].u[0][3]);
	      uetaX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u[0][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u[0][3]);
	      
	      uetaY1 = ((1.-yfrac)*uetaX1+yfrac*uetaX2);
	      uetaY2 = ((1.-yfrac)*uetaX3+yfrac*uetaX4);
	      
	      ueta1 = ((1.-etafrac)*uetaY1+etafrac*uetaY2);
	      
	      //uetaX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][3]);
	      //uetaX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][3]);
	      //uetaX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][3]);
	      //uetaX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][3]);
	      uetaX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[3] + xfrac*arena[ix+fac][iy][ieta].u_prev[3]);
	      uetaX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[3] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[3]);
	      uetaX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[3] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[3]);
	      uetaX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[3]);
	      uetaY1 = ((1.-yfrac)*uetaX1+yfrac*uetaX2);
	      uetaY2 = ((1.-yfrac)*uetaX3+yfrac*uetaX4);
	      
	      ueta2 = ((1.-etafrac)*uetaY1+etafrac*uetaY2);
	      
	      ueta = (1.-taufrac)*ueta2+taufrac*ueta1;
	      
	      //CFY: For now, no interpolation of Wmunu:
	      double W00, W01, W02, W03, W11, W12, W13, W22, W23, W33;
	      if(DATA->fluctuatingHydroFlag == 2){
		W00 = arena[ix][iy][ieta].Wmunu[0][0][0] + arena[ix][iy][ieta].deltaT[0][0][0]
		  + arena[ix][iy][ieta].deltaW[0][0][0];
		W01 = arena[ix][iy][ieta].Wmunu[0][0][1] + arena[ix][iy][ieta].deltaT[0][0][1]
		  + arena[ix][iy][ieta].deltaW[0][0][1];
		W02 = arena[ix][iy][ieta].Wmunu[0][0][2] + arena[ix][iy][ieta].deltaT[0][0][2]
		  + arena[ix][iy][ieta].deltaW[0][0][2];
		W03 = arena[ix][iy][ieta].Wmunu[0][0][3] + arena[ix][iy][ieta].deltaT[0][0][3]
		  + arena[ix][iy][ieta].deltaW[0][0][3];
		W11 = arena[ix][iy][ieta].Wmunu[0][1][1] + arena[ix][iy][ieta].deltaT[0][1][1]
		  + arena[ix][iy][ieta].deltaW[0][1][1];
		W12 = arena[ix][iy][ieta].Wmunu[0][1][2] + arena[ix][iy][ieta].deltaT[0][1][2]
		  + arena[ix][iy][ieta].deltaW[0][1][2];
		W13 = arena[ix][iy][ieta].Wmunu[0][1][3] + arena[ix][iy][ieta].deltaT[0][1][3]
		  + arena[ix][iy][ieta].deltaW[0][1][3];
		W22 = arena[ix][iy][ieta].Wmunu[0][2][2] + arena[ix][iy][ieta].deltaT[0][2][2]
		  + arena[ix][iy][ieta].deltaW[0][2][2];
		W23 = arena[ix][iy][ieta].Wmunu[0][2][3] + arena[ix][iy][ieta].deltaT[0][2][3]
		  + arena[ix][iy][ieta].deltaW[0][2][3];
		W33 = arena[ix][iy][ieta].Wmunu[0][3][3] + arena[ix][iy][ieta].deltaT[0][3][3]
		  + arena[ix][iy][ieta].deltaW[0][3][3];
	      }
	      else{
		W00 = arena[ix][iy][ieta].Wmunu[0][0][0];
		W01 = arena[ix][iy][ieta].Wmunu[0][0][1];
		W02 = arena[ix][iy][ieta].Wmunu[0][0][2];
		W03 = arena[ix][iy][ieta].Wmunu[0][0][3];
		W11 = arena[ix][iy][ieta].Wmunu[0][1][1];
		W12 = arena[ix][iy][ieta].Wmunu[0][1][2];
		W13 = arena[ix][iy][ieta].Wmunu[0][1][3];
		W22 = arena[ix][iy][ieta].Wmunu[0][2][2];
		W23 = arena[ix][iy][ieta].Wmunu[0][2][3];
		W33 = arena[ix][iy][ieta].Wmunu[0][3][3];
	      }
	      
	      //rhob: 3D interpolation for now (no tau interpolation)
	      rhobX1 = ((1.-xfrac)*arena[ix][iy][ieta].rhob + xfrac*arena[ix+fac][iy][ieta].rhob);
	      rhobX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].rhob + xfrac*arena[ix+fac][iy+fac][ieta].rhob);
	      rhobX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].rhob + xfrac*arena[ix+fac][iy][ieta+fac].rhob);
	      rhobX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].rhob + xfrac*arena[ix+fac][iy+fac][ieta+fac].rhob);
	      
	      rhobY1 = ((1.-yfrac)*rhobX1+yfrac*rhobX2);
	      rhobY2 = ((1.-yfrac)*rhobX3+yfrac*rhobX4);
	      
	      rhob = ((1.-etafrac)*rhobY1+etafrac*rhobY2);
	      
	      if (DATA->whichEOS==1)
		{
		  TFO = eos->interpolate(epsFO, rhob, 0);
		  muB = eos->interpolate(epsFO, rhob, 1);
		}
	      else if (DATA->whichEOS==2)
		{
		  TFO = eos->interpolate2(epsFO, rhob, 1);
		  muB = 0.0;
		}

 
	   /*    if (xf==x) */
/* 		xf = x+DX/2.; */
/* 	      if (yf==y) */
/* 		yf = y+DY/2.; */
/* 	      if (etaf==eta) */
/* 		etaf = eta+DETA/2.; */
	      
	      if(shown==0)
		{
		  fprintf(stderr,"DTAU=%f\n", DTAU);
		  fprintf(stderr,"DX=%f\n", DX);
		  fprintf(stderr,"DY=%f\n", DY);
		  fprintf(stderr,"DETA=%f\n", DETA);
		  fprintf(stderr,"Volume=%f\n", SUM);
		  shown=1;
		}
		
	      if (intersecttau)
		{
		  //fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
		  //	  tauf, x, y, eta, FULLSU[0], 0., 0., 0.,
		  //	  utau, ux, uy, ueta, epsFO, TFO, muB);
		  fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
			  tauf, x, y, eta, FULLSU[0], 0., 0., 0.,
			  utau, ux, uy, ueta, epsFO, TFO, muB, sFO,
			  W00, W01, W02, W03,
			  W11, W12, W13,
			  W22, W23,
			  W33);
		  if(fabs(x)<0.05 && (fabs(y)<0.05))
		    {
		      fprintf(t_file,"%f %f %f %f %f %f %f %f\n",tauf, x, y, eta, FULLSU[0],0.,0.,0.);
		    }
		}
	      if (intersectx)
		{
		  fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
			  tau, xf, y, eta, 0., FULLSU[1], 0., 0.,
			  utau, ux, uy, ueta, epsFO, TFO, muB, sFO,
                          W00, W01, W02, W03,
                          W11, W12, W13,
                          W22, W23,
                          W33);

		  if(fabs(x)<0.05 && (fabs(y)<0.05))
		    {
		      fprintf(t_file,"%f %f %f %f %f %f %f %f\n",tau, xf, y, eta, 0., FULLSU[1],0.,0.);
		    }
		}
	      if (intersecty)
		{
		  fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
			  tau, x, yf, eta, 0., 0., FULLSU[2], 0.,
			  utau, ux, uy, ueta, epsFO, TFO, muB,  sFO,
                          W00, W01, W02, W03,
                          W11, W12, W13,
                          W22, W23,
                          W33);

		  if(fabs(x)<0.05 && (fabs(y)<0.05))
		    {
		      fprintf(t_file,"%f %f %f %f %f %f %f %f\n",tau, x, yf, eta, 0.,0.,FULLSU[2],0.);
		    }
		}
	      if (intersecteta)
		{
		  fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
			  tau, x, y, etaf, 0., 0., 0., FULLSU[3],
			  utau, ux, uy, ueta, epsFO, TFO, muB,  sFO,
                          W00, W01, W02, W03,
                          W11, W12, W13,
                          W22, W23,
                          W33);

		  if(fabs(x)<0.05 && (fabs(y)<0.05))
		    {
		      fprintf(t_file,"%f %f %f %f %f %f %f %f\n",tau, x, y, etaf, 0.,0.,0., FULLSU[3]);
		    }
		}
	      
	      //	      fprintf(stderr,"after interpolation\n");
	      //fprintf(stderr,"ix=%d, iy=%d, ieta=%d\n",ix,iy,ieta);
	      //fprintf(stderr,"xf=%f, yf=%f, etaf=%f\n",xf,yf,etaf);
	    }
	}
    }
  fprintf(stderr,"Volume=%f\n", SUM);
  delete []package;
  fclose(t_file);
  fclose(s_file);
}

void Evolve::FindIsochronalFreezeOutSurface(double tau, InitData *DATA, Grid ***arena, int size, int rank){
  stringstream numIn;
  numIn << DATA->seed;
  stringstream rankIn;
  rankIn << rank;

  FILE *t_file;
  string t_name = "tauf_"+numIn.str()+".dat";
  t_file = fopen(t_name.c_str(), "a");

  FILE *s_file;
  string s_name = "surface"+rankIn.str()+"_"+numIn.str()+".dat";
  s_file = fopen(s_name.c_str(), "a");

  FILE *ds_file;
  string ds_name = "dsurface"+rankIn.str()+"_"+numIn.str()+".dat";
  ds_file = fopen(ds_name.c_str(), "a");

  //Simply output all cell's surface elements:
  for(int ix=0; ix<DATA->nx; ix++){
    for(int iy=0; iy<DATA->ny; iy++){
      for(int ieta=0; ieta<DATA->neta; ieta++){
	//cout << "eta = " << DATA->delta_eta*((double)(ieta+rank*DATA->neta))-0.5*(double)DATA->eta_size << endl;
	if(arena[ix][iy][ieta].T > 0.01){
	  double entrDens = eos->s_func(arena[ix][iy][ieta].epsilon, arena[ix][iy][ieta].p, arena[ix][iy][ieta].rhob, DATA);
	  if(DATA->fluctuatingHydroFlag == 1){
	    fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
		    tau, DATA->delta_x*(double)ix-0.5*DATA->x_size, DATA->delta_y*(double)iy-0.5*DATA->y_size, 
		    DATA->delta_eta*((double)(ieta+rank*DATA->neta))-0.5*(double)DATA->eta_size, 
		    (DATA->delta_x)*(DATA->delta_y)*(DATA->delta_eta), 0., 0., 0.,
		    arena[ix][iy][ieta].u[0][0], arena[ix][iy][ieta].u[0][1], arena[ix][iy][ieta].u[0][2], arena[ix][iy][ieta].u[0][3], 
		    arena[ix][iy][ieta].epsilon, arena[ix][iy][ieta].T, arena[ix][iy][ieta].mu, entrDens,
		    //arena[ix][iy][ieta].Wmunu[0][0][0],
		    //arena[ix][iy][ieta].Wmunu[0][0][1],
		    //arena[ix][iy][ieta].Wmunu[0][0][2],
		    //arena[ix][iy][ieta].Wmunu[0][0][3],
		    //arena[ix][iy][ieta].Wmunu[0][1][1],
		    //arena[ix][iy][ieta].Wmunu[0][1][2],
		    //arena[ix][iy][ieta].Wmunu[0][1][3],
		    //arena[ix][iy][ieta].Wmunu[0][2][2],
		    //arena[ix][iy][ieta].Wmunu[0][2][3],
		    //arena[ix][iy][ieta].Wmunu[0][3][3]);
		    arena[ix][iy][ieta].Wmunu[0][0][0]+arena[ix][iy][ieta].Xi[0][0][0], 
		    arena[ix][iy][ieta].Wmunu[0][0][1]+arena[ix][iy][ieta].Xi[0][0][1], 
		    arena[ix][iy][ieta].Wmunu[0][0][2]+arena[ix][iy][ieta].Xi[0][0][2], 
		    arena[ix][iy][ieta].Wmunu[0][0][3]+arena[ix][iy][ieta].Xi[0][0][3], 
		    arena[ix][iy][ieta].Wmunu[0][1][1]+arena[ix][iy][ieta].Xi[0][1][1], 
		    arena[ix][iy][ieta].Wmunu[0][1][2]+arena[ix][iy][ieta].Xi[0][1][2], 
		    arena[ix][iy][ieta].Wmunu[0][1][3]+arena[ix][iy][ieta].Xi[0][1][3], 
		    arena[ix][iy][ieta].Wmunu[0][2][2]+arena[ix][iy][ieta].Xi[0][2][2], 
		    arena[ix][iy][ieta].Wmunu[0][2][3]+arena[ix][iy][ieta].Xi[0][2][3], 
		    arena[ix][iy][ieta].Wmunu[0][3][3]+arena[ix][iy][ieta].Xi[0][3][3]);
	  }
	  else if(DATA->fluctuatingHydroFlag == 2){
	    /*
	    fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
		    tau, DATA->delta_x*((double)ix - 0.5*(double)DATA->nx), DATA->delta_y*((double)iy - 0.5*(double)DATA->ny), 
		    DATA->delta_eta*((double)(ieta+rank*DATA->neta))-0.5*(double)DATA->eta_size, 
		    (DATA->delta_x)*(DATA->delta_y)*(DATA->delta_eta), 0., 0., 0.,
		    arena[ix][iy][ieta].u[0][0], arena[ix][iy][ieta].u[0][1], arena[ix][iy][ieta].u[0][2], arena[ix][iy][ieta].u[0][3], 
		    arena[ix][iy][ieta].epsilon, arena[ix][iy][ieta].T, arena[ix][iy][ieta].mu, entrDens,
		    arena[ix][iy][ieta].Wmunu[0][0][0] + arena[ix][iy][ieta].deltaT[0][0][0] + arena[ix][iy][ieta].deltaW[0][0][0], 
		    arena[ix][iy][ieta].Wmunu[0][0][1] + arena[ix][iy][ieta].deltaT[0][0][1] + arena[ix][iy][ieta].deltaW[0][0][1], 
		    arena[ix][iy][ieta].Wmunu[0][0][2] + arena[ix][iy][ieta].deltaT[0][0][2] + arena[ix][iy][ieta].deltaW[0][0][2], 
		    arena[ix][iy][ieta].Wmunu[0][0][3] + arena[ix][iy][ieta].deltaT[0][0][3] + arena[ix][iy][ieta].deltaW[0][0][3], 
		    arena[ix][iy][ieta].Wmunu[0][1][1] + arena[ix][iy][ieta].deltaT[0][1][1] + arena[ix][iy][ieta].deltaW[0][1][1], 
		    arena[ix][iy][ieta].Wmunu[0][1][2] + arena[ix][iy][ieta].deltaT[0][1][2] + arena[ix][iy][ieta].deltaW[0][1][2], 
		    arena[ix][iy][ieta].Wmunu[0][1][3] + arena[ix][iy][ieta].deltaT[0][1][3] + arena[ix][iy][ieta].deltaW[0][1][3], 
		    arena[ix][iy][ieta].Wmunu[0][2][2] + arena[ix][iy][ieta].deltaT[0][2][2] + arena[ix][iy][ieta].deltaW[0][2][2], 
		    arena[ix][iy][ieta].Wmunu[0][2][3] + arena[ix][iy][ieta].deltaT[0][2][3] + arena[ix][iy][ieta].deltaW[0][2][3], 
		    arena[ix][iy][ieta].Wmunu[0][3][3] + arena[ix][iy][ieta].deltaT[0][3][3] + arena[ix][iy][ieta].deltaW[0][3][3]);
	    */
	    fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
		    tau, DATA->delta_x*((double)ix - 0.5*(double)DATA->nx), DATA->delta_y*((double)iy - 0.5*(double)DATA->ny), 
		    DATA->delta_eta*((double)(ieta+rank*DATA->neta))-0.5*(double)DATA->eta_size, 
		    (DATA->delta_x)*(DATA->delta_y)*(DATA->delta_eta), 0., 0., 0.,
		    arena[ix][iy][ieta].u[0][0], arena[ix][iy][ieta].u[0][1], arena[ix][iy][ieta].u[0][2], arena[ix][iy][ieta].u[0][3], 
		    arena[ix][iy][ieta].epsilon, arena[ix][iy][ieta].T, arena[ix][iy][ieta].mu, entrDens,
		    arena[ix][iy][ieta].Wmunu[0][0][0], arena[ix][iy][ieta].Wmunu[0][0][1], arena[ix][iy][ieta].Wmunu[0][0][2], arena[ix][iy][ieta].Wmunu[0][0][3],
		    arena[ix][iy][ieta].Wmunu[0][1][1], arena[ix][iy][ieta].Wmunu[0][1][2], arena[ix][iy][ieta].Wmunu[0][1][3],
		    arena[ix][iy][ieta].Wmunu[0][2][2], arena[ix][iy][ieta].Wmunu[0][2][3],
		    arena[ix][iy][ieta].Wmunu[0][3][3]);
	    double dT = arena[ix][iy][ieta].deltaP[0]/eos->s_func(arena[ix][iy][ieta].epsilon, 
								  eos->p_func(arena[ix][iy][ieta].epsilon, arena[ix][iy][ieta].rhob), 
								  arena[ix][iy][ieta].rhob, DATA);
	    double de = arena[ix][iy][ieta].deltaP[0]/eos->get_dpOverde2(arena[ix][iy][ieta].epsilon, arena[ix][iy][ieta].rhob);
	    fprintf(ds_file, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e",  dT, 
		    arena[ix][iy][ieta].deltaU[0][0], arena[ix][iy][ieta].deltaU[0][1], arena[ix][iy][ieta].deltaU[0][2], arena[ix][iy][ieta].deltaU[0][3], de, 
		    arena[ix][iy][ieta].deltaP[0],
		    arena[ix][iy][ieta].deltaW[0][0][0], arena[ix][iy][ieta].deltaW[0][0][1], arena[ix][iy][ieta].deltaW[0][0][2], arena[ix][iy][ieta].deltaW[0][0][3],
		    arena[ix][iy][ieta].deltaW[0][1][1], arena[ix][iy][ieta].deltaW[0][1][2], arena[ix][iy][ieta].deltaW[0][1][3],
		    arena[ix][iy][ieta].deltaW[0][2][2], arena[ix][iy][ieta].deltaW[0][2][3],
		    arena[ix][iy][ieta].deltaW[0][3][3]);
	  }
	  else{
	    fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
		    tau, DATA->delta_x*((double)ix - 0.5*(double)DATA->nx), DATA->delta_y*((double)iy - 0.5*(double)DATA->ny), 
		    DATA->delta_eta*((double)(ieta+rank*DATA->neta))-0.5*(double)DATA->eta_size, 
		    (DATA->delta_x)*(DATA->delta_y)*(DATA->delta_eta), 0., 0., 0.,
		    arena[ix][iy][ieta].u[0][0], arena[ix][iy][ieta].u[0][1], arena[ix][iy][ieta].u[0][2], arena[ix][iy][ieta].u[0][3], 
		    arena[ix][iy][ieta].epsilon, arena[ix][iy][ieta].T, arena[ix][iy][ieta].mu, entrDens,
		    arena[ix][iy][ieta].Wmunu[0][0][0], 
		    arena[ix][iy][ieta].Wmunu[0][0][1], 
		    arena[ix][iy][ieta].Wmunu[0][0][2], 
		    arena[ix][iy][ieta].Wmunu[0][0][3], 
		    arena[ix][iy][ieta].Wmunu[0][1][1], 
		    arena[ix][iy][ieta].Wmunu[0][1][2], 
		    arena[ix][iy][ieta].Wmunu[0][1][3], 
		    arena[ix][iy][ieta].Wmunu[0][2][2], 
		    arena[ix][iy][ieta].Wmunu[0][2][3], 
		    arena[ix][iy][ieta].Wmunu[0][3][3]);

	  }
	}
      }
    }
  }

  return;

}

void Evolve::FindFreezeOutSurface2(double tau, InitData *DATA, Grid ***arena, int size, int rank)
{	

  stringstream numIn;
  numIn << DATA->seed;
  stringstream rankIn;
  rankIn << rank;
  FILE *t_file;
  //char* t_name = "tauf.dat";
  //t_file = fopen(t_name, "a");
  string t_name = "tauf_"+numIn.str()+".dat";
  t_file = fopen(t_name.c_str(), "a");
  FILE *t2_file;
  //char* t2_name = "taufx.dat";
  //t2_file = fopen(t2_name, "a");
  string t2_name = "taufx_"+numIn.str()+".dat";
  t2_file = fopen(t2_name.c_str(), "a");
  FILE *t3_file;
  //char* t3_name = "taufy.dat";
  //t3_file = fopen(t3_name, "a");
  string t3_name = "taufy_"+numIn.str()+".dat";
  t3_file = fopen(t3_name.c_str(), "a");
  //char *buf;
  //buf = util->char_malloc(40);
  FILE *s_file;
  //char* s_name;
  //s_name = util->char_malloc(100);
  string s_name = "surface"+rankIn.str()+"_"+numIn.str()+".dat";
  
  //sprintf (buf, "%d", rank);
  //strcat(s_name, "surface");
  //strcat(s_name,buf);
  //strcat(s_name, ".dat");

  //s_file = fopen(s_name, "a");
  s_file = fopen(s_name.c_str(), "a");

  //C-style sucks balls, switch to C++ style I/O:
  //FILE *ds_file;
  //string ds_name = "dsurface"+rankIn.str()+"_"+numIn.str()+".dat";
  //ds_file = fopen(ds_name.c_str(), "a");
  ofstream ds_file;
  string ds_name = "dsurface"+rankIn.str()+"_"+numIn.str()+".dat";
  ds_file.open(ds_name.c_str(), fstream::app);

  int ix, iy, ieta, nx, ny, neta;
  double x, y, eta;
  double epsFO=DATA->epsilonFreeze/hbarc;
  double tauf, xf, yf, etaf;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;
  double cube[16];
  double EK, EL, DEK, DEL, ADEK, ADEL, ELowerSum, EHigherSum;
  double cuts[32][4]; // a 4d hypercube has (2^(n-1)*n=32) edges.
  double VLower0, VHigher0, VLower1, VHigher1, VLower2, VHigher2, VLower3, VHigher3;
  double V0, V1, V2, V3, VD0, VD1, VD2, VD3;
  double VMID[4], AD[32][4], BD[4], CD[4], DD[32][4], SU[32][4], FULLSU[4];
  int NSurfaces;
  int iEdge[32];
  int prevEdge[32];
  int neighbors[32][6];
  int neighborsDone[32][6];
  int additionalNeighbors[32][6];
  int edge1[32];
  int edge2[32];
  int is, intersect, IE, JE, IS, KE, i, j, k, MINPTS, JMIN, IPTS, NSE, NSM, M, M2;
  double DX, DY, DETA, DTAU, APU, SIG;
  int countEdges, skip;
  int previousEdges[32], usedEdges1[32],usedEdges2[32],usedEdges3[32],usedEdges4[32], countAdditionalEdges;
  int fac,l, COUNTER, additionalEdges, temp, m, m2;
  int group[32][3];
  int tries, tries2, tries3;
  int shift;
  int addCOUNTER;
  int ISID[32];
  int countSingleEdges;
  int singleConnections[32][2];
  int singleConnectionsUsed[32];
  int numberConnectionIsUsed[32][32]; // should be 2 in the end
  double Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxeta, Wyy, Wyeta, Wetaeta;
  double Wtautau1, Wtaux1, Wtauy1, Wtaueta1, Wxx1, Wxy1, Wxeta1, Wyy1, Wyeta1, Wetaeta1;
  double Wtautau2, Wtaux2, Wtauy2, Wtaueta2, Wxx2, Wxy2, Wxeta2, Wyy2, Wyeta2, Wetaeta2;
  double WX1, WX2, WX3, WX4, WY1, WY2;
  double rhob, utau, ux, uy, ueta, TFO, muB;
  double utauX1, utauX2, utauX3, utauX4, utauY1, utauY2, utau1, utau2;
  double rhobX1, rhobX2, rhobX3, rhobX4, rhobY1, rhobY2, rhob1, rhob2;
  double uxX1, uxX2, uxX3, uxX4, uxY1, uxY2, ux1, ux2;
  double uyX1, uyX2, uyX3, uyX4, uyY1, uyY2, uy1, uy2;
  double uetaX1, uetaX2, uetaX3, uetaX4, uetaY1, uetaY2, ueta1, ueta2;
  double xfrac, yfrac, etafrac, taufrac;
  double dWtautau, dWtaux, dWtauy, dWtaueta, dWxx, dWxy, dWxeta, dWyy, dWyeta, dWetaeta;
  double dWtautau1, dWtaux1, dWtauy1, dWtaueta1, dWxx1, dWxy1, dWxeta1, dWyy1, dWyeta1, dWetaeta1;
  double dWtautau2, dWtaux2, dWtauy2, dWtaueta2, dWxx2, dWxy2, dWxeta2, dWyy2, dWyeta2, dWetaeta2;
  double dWX1, dWX2, dWX3, dWX4, dWY1, dWY2;
  double dutau, dux, duy, dueta, dp;
  double dutauX1, dutauX2, dutauX3, dutauX4, dutauY1, dutauY2, dutau1, dutau2;
  double duxX1, duxX2, duxX3, duxX4, duxY1, duxY2, dux1, dux2;
  double duyX1, duyX2, duyX3, duyX4, duyY1, duyY2, duy1, duy2;
  double duetaX1, duetaX2, duetaX3, duetaX4, duetaY1, duetaY2, dueta1, dueta2;
  // who is direct neighbor (1) or neighbor across the plane (2), more distant neighbor (3), (4)
  int intersections;
  int maxEta;
  double sFO, P;
  int const IBIT[32][32] = 
    { {0,1,2,1,1,1,3,3,2,3,4,3,2,3,4,3,3,3,0,0,4,0,0,0,1,1,3,3,3,3,0,0},
      {1,0,1,2,3,1,1,3,3,2,3,4,3,2,3,4,3,3,0,0,0,4,0,0,3,1,1,3,0,3,3,0},
      {2,1,0,1,3,3,1,1,4,3,2,3,4,3,2,3,0,0,3,3,0,0,4,0,3,3,1,1,0,0,3,3},
      {1,2,1,0,1,3,3,1,3,4,3,2,3,4,3,2,3,0,0,3,0,4,0,0,1,3,3,1,3,0,0,3},
      {1,3,3,1,0,2,4,2,1,3,3,1,3,0,0,3,2,4,0,4,3,0,0,3,1,3,0,3,1,3,0,3},
      {1,1,3,3,2,0,2,4,1,1,3,3,3,3,0,0,4,2,4,0,3,3,0,0,3,1,3,0,3,1,3,0},
      {3,1,1,3,4,2,0,2,3,1,1,3,0,3,3,0,0,4,2,4,0,3,3,0,0,3,1,3,0,3,1,3},
      {3,3,1,1,2,4,2,0,3,3,1,1,0,0,3,3,4,0,4,2,0,0,3,3,3,0,3,1,3,0,3,1},
      {2,3,4,3,1,1,3,3,0,1,2,1,4,0,0,0,3,3,0,0,2,3,4,3,3,3,0,0,1,1,3,3},
      {3,2,3,4,3,1,1,3,1,0,1,2,0,4,0,0,0,3,3,0,3,2,3,4,0,3,3,0,3,1,1,3},
      {4,3,2,3,3,3,1,1,2,1,0,1,0,0,4,0,0,0,3,3,4,3,2,3,0,0,3,3,3,3,1,1},
      {3,4,3,2,1,3,3,1,1,2,1,0,0,0,0,4,3,0,0,3,3,4,3,2,3,0,0,3,1,3,3,1},
      {2,3,4,3,3,3,0,0,4,0,0,0,0,1,2,1,1,1,3,3,2,3,4,3,1,1,3,3,3,3,0,0},
      {3,2,3,4,0,3,3,0,0,4,0,0,1,0,1,2,3,1,1,3,3,2,3,4,3,1,1,3,0,3,3,0},
      {4,3,2,3,0,0,3,3,0,0,4,0,2,1,0,1,3,3,1,1,4,3,2,3,3,3,1,1,0,0,3,3},
      {3,4,3,2,3,0,0,3,0,0,0,4,1,2,1,0,1,3,3,1,3,4,3,2,1,3,3,1,3,0,0,3},
      {3,0,0,3,2,4,0,4,3,0,0,3,1,3,3,1,0,2,4,2,1,3,3,1,1,3,0,3,1,3,0,3},
      {3,3,0,0,4,2,4,0,3,3,0,0,1,1,3,3,2,0,2,4,1,1,3,3,3,1,3,0,3,1,3,0},
      {0,3,3,0,0,4,2,4,0,3,3,0,3,1,1,3,4,2,0,2,3,1,1,3,0,3,1,3,0,3,1,3},
      {0,0,3,3,4,0,4,2,0,0,3,3,3,3,1,1,2,4,2,0,3,3,1,1,3,0,3,1,3,0,3,1},
      {4,0,0,0,3,3,0,0,2,3,4,3,2,3,4,3,1,1,3,3,0,1,2,1,3,3,0,0,1,1,3,3},
      {0,4,0,0,0,3,3,0,3,2,3,4,3,2,3,4,3,1,1,3,1,0,1,2,0,3,3,0,3,1,1,3},
      {0,0,4,0,0,0,3,3,4,3,2,3,4,3,2,3,3,3,1,1,2,1,0,1,0,0,3,3,3,3,1,1},
      {0,0,0,4,3,0,0,3,3,4,3,2,3,4,3,2,1,3,3,1,1,2,1,0,3,0,0,3,1,3,3,1},
      {1,3,3,1,1,3,0,3,3,0,0,3,1,3,3,1,1,3,0,3,3,0,0,3,0,2,4,2,2,4,0,4},
      {1,1,3,3,3,1,3,0,3,3,0,0,1,1,3,3,3,1,3,0,3,3,0,0,2,0,2,4,4,2,4,0},
      {3,1,1,3,0,3,1,3,0,3,3,0,3,1,1,3,0,3,1,3,0,3,3,0,4,2,0,2,0,4,2,4},
      {3,3,1,1,3,0,3,1,0,0,3,3,3,3,1,1,3,0,3,1,0,0,3,3,2,4,2,0,4,0,4,2},
      {3,0,0,3,1,3,0,3,1,3,3,1,3,0,0,3,1,3,0,3,1,3,3,1,2,4,0,4,0,2,4,2},
      {3,3,0,0,3,1,3,0,1,1,3,3,3,3,0,0,3,1,3,0,1,1,3,3,4,2,4,0,2,0,2,4},
      {0,3,3,0,0,3,1,3,3,1,1,3,0,3,3,0,0,3,1,3,3,1,1,3,0,4,2,4,4,2,0,2},
      {0,0,3,3,3,0,3,1,3,3,1,1,0,0,3,3,3,0,3,1,3,3,1,1,4,0,4,2,2,4,2,0},
    };
  
  fac=1; // keep fac=1 in this version for now (parallel needs it to be)
  facTau = DATA->facTau;
  DX=fac*DATA->delta_x;
  DY=fac*DATA->delta_y;
  DETA=fac*DATA->delta_eta;
  DTAU=facTau*DATA->delta_tau;
  intersections=0;
  
  //  int ix, iy, ieta, nx, ny, neta, i, alpha, iflag;
  int sizeOfData = (nx+1)*(ny+1);
  int position;
  double *package;
  double *package_prev;
  double *packageutau;
  double *packageux;
  double *packageuy;
  double *packageueta;
  double *packagerhob;

  double *package_dp;
  double *package_dp_prev;
  double *packagedutau;
  double *packagedux;
  double *packageduy;
  double *packagedueta;

  double *packageutau_prev;
  double *packageux_prev;
  double *packageuy_prev;
  double *packageueta_prev;
  double *packagerhob_prev;
  
  double *packagedutau_prev;
  double *packagedux_prev;
  double *packageduy_prev;
  double *packagedueta_prev;
  
  double **Rneighbor_eps;
  double **Rneighbor_eps_prev;
  double **Rneighbor_utau;
  double **Rneighbor_ux;
  double **Rneighbor_uy;
  double **Rneighbor_ueta;
  double **Rneighbor_rhob;
 
  double **Rneighbor_dp;
  double **Rneighbor_dp_prev;
  double **Rneighbor_dutau;
  double **Rneighbor_dux;
  double **Rneighbor_duy;
  double **Rneighbor_dueta;

  double **Rneighbor_utau_prev;
  double **Rneighbor_ux_prev;
  double **Rneighbor_uy_prev;
  double **Rneighbor_ueta_prev;
  double **Rneighbor_rhob_prev;

  double **Rneighbor_dutau_prev;
  double **Rneighbor_dux_prev;
  double **Rneighbor_duy_prev;
  double **Rneighbor_dueta_prev;

  double *packageWtautau;
  double *packageWtaux;
  double *packageWtauy;
  double *packageWtaueta;
  double *packageWxx;
  double *packageWxy;
  double *packageWxeta;
  double *packageWyy;
  double *packageWyeta;
  double *packageWetaeta;

  double *packagedWtautau;
  double *packagedWtaux;
  double *packagedWtauy;
  double *packagedWtaueta;
  double *packagedWxx;
  double *packagedWxy;
  double *packagedWxeta;
  double *packagedWyy;
  double *packagedWyeta;
  double *packagedWetaeta;

  double *packageWtautau_prev;
  double *packageWtaux_prev;
  double *packageWtauy_prev;
  double *packageWtaueta_prev;
  double *packageWxx_prev;
  double *packageWxy_prev;
  double *packageWxeta_prev;
  double *packageWyy_prev;
  double *packageWyeta_prev;
  double *packageWetaeta_prev;

  double *packagedWtautau_prev;
  double *packagedWtaux_prev;
  double *packagedWtauy_prev;
  double *packagedWtaueta_prev;
  double *packagedWxx_prev;
  double *packagedWxy_prev;
  double *packagedWxeta_prev;
  double *packagedWyy_prev;
  double *packagedWyeta_prev;
  double *packagedWetaeta_prev;

  double **Rneighbor_Wtautau_prev;
  double **Rneighbor_Wtaux_prev;
  double **Rneighbor_Wtauy_prev;
  double **Rneighbor_Wtaueta_prev;
  double **Rneighbor_Wxx_prev;
  double **Rneighbor_Wxy_prev;
  double **Rneighbor_Wxeta_prev;
  double **Rneighbor_Wyy_prev;
  double **Rneighbor_Wyeta_prev;
  double **Rneighbor_Wetaeta_prev;

  double **Rneighbor_dWtautau_prev;
  double **Rneighbor_dWtaux_prev;
  double **Rneighbor_dWtauy_prev;
  double **Rneighbor_dWtaueta_prev;
  double **Rneighbor_dWxx_prev;
  double **Rneighbor_dWxy_prev;
  double **Rneighbor_dWxeta_prev;
  double **Rneighbor_dWyy_prev;
  double **Rneighbor_dWyeta_prev;
  double **Rneighbor_dWetaeta_prev;

  double **Rneighbor_Wtautau;
  double **Rneighbor_Wtaux;
  double **Rneighbor_Wtauy;
  double **Rneighbor_Wtaueta;
  double **Rneighbor_Wxx;
  double **Rneighbor_Wxy;
  double **Rneighbor_Wxeta;
  double **Rneighbor_Wyy;
  double **Rneighbor_Wyeta;
  double **Rneighbor_Wetaeta;

  double **Rneighbor_dWtautau;
  double **Rneighbor_dWtaux;
  double **Rneighbor_dWtauy;
  double **Rneighbor_dWtaueta;
  double **Rneighbor_dWxx;
  double **Rneighbor_dWxy;
  double **Rneighbor_dWxeta;
  double **Rneighbor_dWyy;
  double **Rneighbor_dWyeta;
  double **Rneighbor_dWetaeta;

  package = new double[sizeOfData];
  package_prev = new double[sizeOfData];
  packageutau = new double[sizeOfData];
  packageux = new double[sizeOfData];
  packageuy = new double[sizeOfData];
  packageueta = new double[sizeOfData];
  packagerhob = new double[sizeOfData];

  packageWtautau = new double[sizeOfData];
  packageWtaux = new double[sizeOfData];
  packageWtauy = new double[sizeOfData];
  packageWtaueta = new double[sizeOfData];
  packageWxx = new double[sizeOfData];
  packageWxy = new double[sizeOfData];
  packageWxeta = new double[sizeOfData];
  packageWyy = new double[sizeOfData];
  packageWyeta = new double[sizeOfData];
  packageWetaeta = new double[sizeOfData];
  
  packageWtautau_prev = new double[sizeOfData];
  packageWtaux_prev = new double[sizeOfData];
  packageWtauy_prev = new double[sizeOfData];
  packageWtaueta_prev = new double[sizeOfData];
  packageWxx_prev = new double[sizeOfData];
  packageWxy_prev = new double[sizeOfData];
  packageWxeta_prev = new double[sizeOfData];
  packageWyy_prev = new double[sizeOfData];
  packageWyeta_prev = new double[sizeOfData];
  packageWetaeta_prev = new double[sizeOfData];
  
  packageutau_prev = new double[sizeOfData];
  packageux_prev = new double[sizeOfData];
  packageuy_prev = new double[sizeOfData];
  packageueta_prev = new double[sizeOfData];
  packagerhob_prev = new double[sizeOfData];
  
  Rneighbor_eps = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_eps_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_utau = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_ux = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_uy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_ueta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_rhob = util->mtx_malloc(nx+1,ny+1);
  
  Rneighbor_Wtautau = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wtaux = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wtauy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wtaueta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wxx = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wxy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wxeta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wyy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wyeta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wetaeta = util->mtx_malloc(nx+1,ny+1);
  
  Rneighbor_Wtautau_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wtaux_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wtauy_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wtaueta_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wxx_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wxy_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wxeta_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wyy_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wyeta_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wetaeta_prev = util->mtx_malloc(nx+1,ny+1);
  
  Rneighbor_utau_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_ux_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_uy_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_ueta_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_rhob_prev = util->mtx_malloc(nx+1,ny+1);
  
  if(DATA->fluctuatingHydroFlag == 2){

    package_dp = new double[sizeOfData];
    package_dp_prev = new double[sizeOfData];
    packagedutau = new double[sizeOfData];
    packagedux = new double[sizeOfData];
    packageduy = new double[sizeOfData];
    packagedueta = new double[sizeOfData];
    
    packagedWtautau = new double[sizeOfData];
    packagedWtaux = new double[sizeOfData];
    packagedWtauy = new double[sizeOfData];
    packagedWtaueta = new double[sizeOfData];
    packagedWxx = new double[sizeOfData];
    packagedWxy = new double[sizeOfData];
    packagedWxeta = new double[sizeOfData];
    packagedWyy = new double[sizeOfData];
    packagedWyeta = new double[sizeOfData];
    packagedWetaeta = new double[sizeOfData];
  
    packagedWtautau_prev = new double[sizeOfData];
    packagedWtaux_prev = new double[sizeOfData];
    packagedWtauy_prev = new double[sizeOfData];
    packagedWtaueta_prev = new double[sizeOfData];
    packagedWxx_prev = new double[sizeOfData];
    packagedWxy_prev = new double[sizeOfData];
    packagedWxeta_prev = new double[sizeOfData];
    packagedWyy_prev = new double[sizeOfData];
    packagedWyeta_prev = new double[sizeOfData];
    packagedWetaeta_prev = new double[sizeOfData];
  
    packagedutau_prev = new double[sizeOfData];
    packagedux_prev = new double[sizeOfData];
    packageduy_prev = new double[sizeOfData];
    packagedueta_prev = new double[sizeOfData];

    Rneighbor_dp = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dp_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dutau = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dux = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_duy = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dueta = util->mtx_malloc(nx+1,ny+1);
    
    Rneighbor_dWtautau = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWtaux = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWtauy = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWtaueta = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWxx = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWxy = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWxeta = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWyy = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWyeta = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWetaeta = util->mtx_malloc(nx+1,ny+1);
    
    Rneighbor_dWtautau_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWtaux_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWtauy_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWtaueta_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWxx_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWxy_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWxeta_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWyy_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWyeta_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dWetaeta_prev = util->mtx_malloc(nx+1,ny+1);
    
    Rneighbor_dutau_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dux_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_duy_prev = util->mtx_malloc(nx+1,ny+1);
    Rneighbor_dueta_prev = util->mtx_malloc(nx+1,ny+1);
  }
  // receive from the right / send to the left
  int from = rank+1;
  int to = rank-1;

  // get cells from neighboring processors
  if ( rank != 0 )
    {
      //      cout << " sending to the left on rank " << rank << endl;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      //if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
	      //cout << "arena[ix][iy][0].TJb[i][alpha][0]=" << arena[ix][iy][0].TJb[i][alpha][0] << endl;
	      position = ix + (nx*iy);
	      package[position] = arena[ix][iy][0].epsilon;
	      package_prev[position] = arena[ix][iy][0].epsilon_prev;
	      //package_prev[position] = arena[ix][iy][0].epsilon_prev[facTau-1];
	      packageutau[position] = arena[ix][iy][0].u[0][0];
	      packageux[position] = arena[ix][iy][0].u[0][1];
	      packageuy[position] = arena[ix][iy][0].u[0][2];
	      packageueta[position] = arena[ix][iy][0].u[0][3];
	      
	      packageWtautau[position] = arena[ix][iy][0].Wmunu[0][0][0];
	      packageWtaux[position] = arena[ix][iy][0].Wmunu[0][0][1];
	      packageWtauy[position] = arena[ix][iy][0].Wmunu[0][0][2];
	      packageWtaueta[position] = arena[ix][iy][0].Wmunu[0][0][3];
	      packageWxx[position] = arena[ix][iy][0].Wmunu[0][1][1];
	      packageWxy[position] = arena[ix][iy][0].Wmunu[0][1][2];
	      packageWxeta[position] = arena[ix][iy][0].Wmunu[0][1][3];
	      packageWyy[position] = arena[ix][iy][0].Wmunu[0][2][2];
	      packageWyeta[position] = arena[ix][iy][0].Wmunu[0][2][3];
	      packageWetaeta[position] = arena[ix][iy][0].Wmunu[0][3][3];
      
	      packageWtautau_prev[position] = arena[ix][iy][0].W_prev[0][0];
	      packageWtaux_prev[position] = arena[ix][iy][0].W_prev[0][1];
	      packageWtauy_prev[position] = arena[ix][iy][0].W_prev[0][2];
	      packageWtaueta_prev[position] = arena[ix][iy][0].W_prev[0][3];
	      packageWxx_prev[position] = arena[ix][iy][0].W_prev[1][1];
	      packageWxy_prev[position] = arena[ix][iy][0].W_prev[1][2];
	      packageWxeta_prev[position] = arena[ix][iy][0].W_prev[1][3];
	      packageWyy_prev[position] = arena[ix][iy][0].W_prev[2][2];
	      packageWyeta_prev[position] = arena[ix][iy][0].W_prev[2][3];
	      packageWetaeta_prev[position] = arena[ix][iy][0].W_prev[3][3];
      
	      if(DATA->fluctuatingHydroFlag == 2){
		package_dp[position] = arena[ix][iy][0].deltaP[0];
		package_dp_prev[position] = arena[ix][iy][0].deltaP_prev;
		//package_prev[position] = arena[ix][iy][0].epsilon_prev[facTau-1];
		packagedutau[position] = arena[ix][iy][0].deltaU[0][0];
		packagedux[position] = arena[ix][iy][0].deltaU[0][1];
		packageduy[position] = arena[ix][iy][0].deltaU[0][2];
		packagedueta[position] = arena[ix][iy][0].deltaU[0][3];
	      
		packagedWtautau[position] = arena[ix][iy][0].deltaW[0][0][0];
		packagedWtaux[position] = arena[ix][iy][0].deltaW[0][0][1];
		packagedWtauy[position] = arena[ix][iy][0].deltaW[0][0][2];
		packagedWtaueta[position] = arena[ix][iy][0].deltaW[0][0][3];
		packagedWxx[position] = arena[ix][iy][0].deltaW[0][1][1];
		packagedWxy[position] = arena[ix][iy][0].deltaW[0][1][2];
		packagedWxeta[position] = arena[ix][iy][0].deltaW[0][1][3];
		packagedWyy[position] = arena[ix][iy][0].deltaW[0][2][2];
		packagedWyeta[position] = arena[ix][iy][0].deltaW[0][2][3];
		packagedWetaeta[position] = arena[ix][iy][0].deltaW[0][3][3];
      
		packagedWtautau_prev[position] = arena[ix][iy][0].deltaW_prev[0][0];
		packagedWtaux_prev[position] = arena[ix][iy][0].deltaW_prev[0][1];
		packagedWtauy_prev[position] = arena[ix][iy][0].deltaW_prev[0][2];
		packagedWtaueta_prev[position] = arena[ix][iy][0].deltaW_prev[0][3];
		packagedWxx_prev[position] = arena[ix][iy][0].deltaW_prev[1][1];
		packagedWxy_prev[position] = arena[ix][iy][0].deltaW_prev[1][2];
		packagedWxeta_prev[position] = arena[ix][iy][0].deltaW_prev[1][3];
		packagedWyy_prev[position] = arena[ix][iy][0].deltaW_prev[2][2];
		packagedWyeta_prev[position] = arena[ix][iy][0].deltaW_prev[2][3];
		packagedWetaeta_prev[position] = arena[ix][iy][0].deltaW_prev[3][3];
	      }      
	      if(DATA->freezeOutMethod == 4){
		packageWtautau[position] += arena[ix][iy][0].deltaT[0][0][0] + arena[ix][iy][0].deltaW[0][0][0];
		packageWtaux[position] +=   arena[ix][iy][0].deltaT[0][0][1] + arena[ix][iy][0].deltaW[0][0][1];
		packageWtauy[position] +=   arena[ix][iy][0].deltaT[0][0][2] + arena[ix][iy][0].deltaW[0][0][2];
		packageWtaueta[position] += arena[ix][iy][0].deltaT[0][0][3] + arena[ix][iy][0].deltaW[0][0][3];
		packageWxx[position] +=     arena[ix][iy][0].deltaT[0][1][1] + arena[ix][iy][0].deltaW[0][1][1];
		packageWxy[position] +=     arena[ix][iy][0].deltaT[0][1][2] + arena[ix][iy][0].deltaW[0][1][2];
		packageWxeta[position] +=   arena[ix][iy][0].deltaT[0][1][3] + arena[ix][iy][0].deltaW[0][1][3];
		packageWyy[position] +=     arena[ix][iy][0].deltaT[0][2][2] + arena[ix][iy][0].deltaW[0][2][2];
		packageWyeta[position] +=   arena[ix][iy][0].deltaT[0][2][3] + arena[ix][iy][0].deltaW[0][2][3];
		packageWetaeta[position] += arena[ix][iy][0].deltaT[0][3][3] + arena[ix][iy][0].deltaW[0][3][3];

		packageWtautau_prev[position] += arena[ix][iy][0].getPrevDeltaT(0, 0, eos) + arena[ix][iy][0].prev_deltaW[0][0];
		packageWtaux_prev[position] +=   arena[ix][iy][0].getPrevDeltaT(0, 1, eos) + arena[ix][iy][0].prev_deltaW[0][1];
		packageWtauy_prev[position] +=   arena[ix][iy][0].getPrevDeltaT(0, 2, eos) + arena[ix][iy][0].prev_deltaW[0][2];
		packageWtaueta_prev[position] += arena[ix][iy][0].getPrevDeltaT(0, 3, eos) + arena[ix][iy][0].prev_deltaW[0][3];
		packageWxx_prev[position] +=     arena[ix][iy][0].getPrevDeltaT(1, 1, eos) + arena[ix][iy][0].prev_deltaW[1][1];
		packageWxy_prev[position] +=     arena[ix][iy][0].getPrevDeltaT(1, 2, eos) + arena[ix][iy][0].prev_deltaW[1][2];
		packageWxeta_prev[position] +=   arena[ix][iy][0].getPrevDeltaT(1, 3, eos) + arena[ix][iy][0].prev_deltaW[1][3];
		packageWyy_prev[position] +=     arena[ix][iy][0].getPrevDeltaT(2, 2, eos) + arena[ix][iy][0].prev_deltaW[2][2];
		packageWyeta_prev[position] +=   arena[ix][iy][0].getPrevDeltaT(2, 3, eos) + arena[ix][iy][0].prev_deltaW[2][3];
		packageWetaeta_prev[position] += arena[ix][iy][0].getPrevDeltaT(3, 3, eos) + arena[ix][iy][0].prev_deltaW[3][3];
	      }

	      packagerhob[position] = arena[ix][iy][0].rhob;
	      //packageutau_prev[position] = arena[ix][iy][0].u_prev[facTau-1][0];
	      //packageux_prev[position] = arena[ix][iy][0].u_prev[facTau-1][1];
	      //packageuy_prev[position] = arena[ix][iy][0].u_prev[facTau-1][2];
	      //packageueta_prev[position] = arena[ix][iy][0].u_prev[facTau-1][3];
	      packageutau_prev[position] = arena[ix][iy][0].u_prev[0];
	      packageux_prev[position] = arena[ix][iy][0].u_prev[1];
	      packageuy_prev[position] = arena[ix][iy][0].u_prev[2];
	      packageueta_prev[position] = arena[ix][iy][0].u_prev[3];

	      if(DATA->fluctuatingHydroFlag == 2){
		packagedutau_prev[position] = arena[ix][iy][0].deltaU_prev[0];
		packagedux_prev[position] = arena[ix][iy][0].deltaU_prev[1];
		packageduy_prev[position] = arena[ix][iy][0].deltaU_prev[2];
		packagedueta_prev[position] = arena[ix][iy][0].deltaU_prev[3];
	      }
	    }
	}
      MPI::COMM_WORLD.Send(package,sizeOfData,MPI::DOUBLE,to,1);
      MPI::COMM_WORLD.Send(package_prev,sizeOfData,MPI::DOUBLE,to,2);
      MPI::COMM_WORLD.Send(packageutau,sizeOfData,MPI::DOUBLE,to,3);
      MPI::COMM_WORLD.Send(packageux,sizeOfData,MPI::DOUBLE,to,4);
      MPI::COMM_WORLD.Send(packageuy,sizeOfData,MPI::DOUBLE,to,5);
      MPI::COMM_WORLD.Send(packageueta,sizeOfData,MPI::DOUBLE,to,6);
      MPI::COMM_WORLD.Send(packagerhob,sizeOfData,MPI::DOUBLE,to,7);
      MPI::COMM_WORLD.Send(packageutau_prev,sizeOfData,MPI::DOUBLE,to,8);
      MPI::COMM_WORLD.Send(packageux_prev,sizeOfData,MPI::DOUBLE,to,9);
      MPI::COMM_WORLD.Send(packageuy_prev,sizeOfData,MPI::DOUBLE,to,10);
      MPI::COMM_WORLD.Send(packageueta_prev,sizeOfData,MPI::DOUBLE,to,11);
      MPI::COMM_WORLD.Send(packageWtautau,sizeOfData,MPI::DOUBLE,to,12);
      MPI::COMM_WORLD.Send(packageWtaux,sizeOfData,MPI::DOUBLE,to,13);
      MPI::COMM_WORLD.Send(packageWtauy,sizeOfData,MPI::DOUBLE,to,14);
      MPI::COMM_WORLD.Send(packageWtaueta,sizeOfData,MPI::DOUBLE,to,15);
      MPI::COMM_WORLD.Send(packageWxx,sizeOfData,MPI::DOUBLE,to,16);
      MPI::COMM_WORLD.Send(packageWxy,sizeOfData,MPI::DOUBLE,to,17);
      MPI::COMM_WORLD.Send(packageWxeta,sizeOfData,MPI::DOUBLE,to,18);
      MPI::COMM_WORLD.Send(packageWyy,sizeOfData,MPI::DOUBLE,to,19);
      MPI::COMM_WORLD.Send(packageWyeta,sizeOfData,MPI::DOUBLE,to,20);
      MPI::COMM_WORLD.Send(packageWetaeta,sizeOfData,MPI::DOUBLE,to,21);
      MPI::COMM_WORLD.Send(packageWtautau_prev,sizeOfData,MPI::DOUBLE,to,22);
      MPI::COMM_WORLD.Send(packageWtaux_prev,sizeOfData,MPI::DOUBLE,to,23);
      MPI::COMM_WORLD.Send(packageWtauy_prev,sizeOfData,MPI::DOUBLE,to,24);
      MPI::COMM_WORLD.Send(packageWtaueta_prev,sizeOfData,MPI::DOUBLE,to,25);
      MPI::COMM_WORLD.Send(packageWxx_prev,sizeOfData,MPI::DOUBLE,to,26);
      MPI::COMM_WORLD.Send(packageWxy_prev,sizeOfData,MPI::DOUBLE,to,27);
      MPI::COMM_WORLD.Send(packageWxeta_prev,sizeOfData,MPI::DOUBLE,to,28);
      MPI::COMM_WORLD.Send(packageWyy_prev,sizeOfData,MPI::DOUBLE,to,29);
      MPI::COMM_WORLD.Send(packageWyeta_prev,sizeOfData,MPI::DOUBLE,to,30);
      MPI::COMM_WORLD.Send(packageWetaeta_prev,sizeOfData,MPI::DOUBLE,to,31);
      if(DATA->fluctuatingHydroFlag ==2){
	MPI::COMM_WORLD.Send(package_dp,sizeOfData,MPI::DOUBLE,to,32);
	MPI::COMM_WORLD.Send(package_dp_prev,sizeOfData,MPI::DOUBLE,to,33);
	MPI::COMM_WORLD.Send(packagedutau,sizeOfData,MPI::DOUBLE,to,34);
	MPI::COMM_WORLD.Send(packagedux,sizeOfData,MPI::DOUBLE,to,35);
	MPI::COMM_WORLD.Send(packageduy,sizeOfData,MPI::DOUBLE,to,36);
	MPI::COMM_WORLD.Send(packagedueta,sizeOfData,MPI::DOUBLE,to,37);
	MPI::COMM_WORLD.Send(packagedutau_prev,sizeOfData,MPI::DOUBLE,to,38);
	MPI::COMM_WORLD.Send(packagedux_prev,sizeOfData,MPI::DOUBLE,to,39);
	MPI::COMM_WORLD.Send(packageduy_prev,sizeOfData,MPI::DOUBLE,to,40);
	MPI::COMM_WORLD.Send(packagedueta_prev,sizeOfData,MPI::DOUBLE,to,41);
	MPI::COMM_WORLD.Send(packagedWtautau,sizeOfData,MPI::DOUBLE,to,42);
	MPI::COMM_WORLD.Send(packagedWtaux,sizeOfData,MPI::DOUBLE,to,43);
	MPI::COMM_WORLD.Send(packagedWtauy,sizeOfData,MPI::DOUBLE,to,44);
	MPI::COMM_WORLD.Send(packagedWtaueta,sizeOfData,MPI::DOUBLE,to,45);
	MPI::COMM_WORLD.Send(packagedWxx,sizeOfData,MPI::DOUBLE,to,46);
	MPI::COMM_WORLD.Send(packagedWxy,sizeOfData,MPI::DOUBLE,to,47);
	MPI::COMM_WORLD.Send(packagedWxeta,sizeOfData,MPI::DOUBLE,to,48);
	MPI::COMM_WORLD.Send(packagedWyy,sizeOfData,MPI::DOUBLE,to,49);
	MPI::COMM_WORLD.Send(packagedWyeta,sizeOfData,MPI::DOUBLE,to,50);
	MPI::COMM_WORLD.Send(packagedWetaeta,sizeOfData,MPI::DOUBLE,to,51);
	MPI::COMM_WORLD.Send(packagedWtautau_prev,sizeOfData,MPI::DOUBLE,to,52);
	MPI::COMM_WORLD.Send(packagedWtaux_prev,sizeOfData,MPI::DOUBLE,to,53);
	MPI::COMM_WORLD.Send(packagedWtauy_prev,sizeOfData,MPI::DOUBLE,to,54);
	MPI::COMM_WORLD.Send(packagedWtaueta_prev,sizeOfData,MPI::DOUBLE,to,55);
	MPI::COMM_WORLD.Send(packagedWxx_prev,sizeOfData,MPI::DOUBLE,to,56);
	MPI::COMM_WORLD.Send(packagedWxy_prev,sizeOfData,MPI::DOUBLE,to,57);
	MPI::COMM_WORLD.Send(packagedWxeta_prev,sizeOfData,MPI::DOUBLE,to,58);
	MPI::COMM_WORLD.Send(packagedWyy_prev,sizeOfData,MPI::DOUBLE,to,59);
	MPI::COMM_WORLD.Send(packagedWyeta_prev,sizeOfData,MPI::DOUBLE,to,60);
	MPI::COMM_WORLD.Send(packagedWetaeta_prev,sizeOfData,MPI::DOUBLE,to,61);
        //cout << " done sending to the left on rank " << rank << endl;
      }
    }
  // receiving and unwrapping the package
  if ( rank != size-1 )
    {  
      MPI::COMM_WORLD.Recv(package,sizeOfData,MPI::DOUBLE,from,1);
      MPI::COMM_WORLD.Recv(package_prev,sizeOfData,MPI::DOUBLE,from,2);
      MPI::COMM_WORLD.Recv(packageutau,sizeOfData,MPI::DOUBLE,from,3);
      MPI::COMM_WORLD.Recv(packageux,sizeOfData,MPI::DOUBLE,from,4);
      MPI::COMM_WORLD.Recv(packageuy,sizeOfData,MPI::DOUBLE,from,5);
      MPI::COMM_WORLD.Recv(packageueta,sizeOfData,MPI::DOUBLE,from,6);
      MPI::COMM_WORLD.Recv(packagerhob,sizeOfData,MPI::DOUBLE,from,7);
      MPI::COMM_WORLD.Recv(packageutau_prev,sizeOfData,MPI::DOUBLE,from,8);
      MPI::COMM_WORLD.Recv(packageux_prev,sizeOfData,MPI::DOUBLE,from,9);
      MPI::COMM_WORLD.Recv(packageuy_prev,sizeOfData,MPI::DOUBLE,from,10);
      MPI::COMM_WORLD.Recv(packageueta_prev,sizeOfData,MPI::DOUBLE,from,11);
      MPI::COMM_WORLD.Recv(packageWtautau,sizeOfData,MPI::DOUBLE,from,12);
      MPI::COMM_WORLD.Recv(packageWtaux,sizeOfData,MPI::DOUBLE,from,13);
      MPI::COMM_WORLD.Recv(packageWtauy,sizeOfData,MPI::DOUBLE,from,14);
      MPI::COMM_WORLD.Recv(packageWtaueta,sizeOfData,MPI::DOUBLE,from,15);
      MPI::COMM_WORLD.Recv(packageWxx,sizeOfData,MPI::DOUBLE,from,16);
      MPI::COMM_WORLD.Recv(packageWxy,sizeOfData,MPI::DOUBLE,from,17);
      MPI::COMM_WORLD.Recv(packageWxeta,sizeOfData,MPI::DOUBLE,from,18);
      MPI::COMM_WORLD.Recv(packageWyy,sizeOfData,MPI::DOUBLE,from,19);
      MPI::COMM_WORLD.Recv(packageWyeta,sizeOfData,MPI::DOUBLE,from,20);
      MPI::COMM_WORLD.Recv(packageWetaeta,sizeOfData,MPI::DOUBLE,from,21);
      MPI::COMM_WORLD.Recv(packageWtautau_prev,sizeOfData,MPI::DOUBLE,from,22);
      MPI::COMM_WORLD.Recv(packageWtaux_prev,sizeOfData,MPI::DOUBLE,from,23);
      MPI::COMM_WORLD.Recv(packageWtauy_prev,sizeOfData,MPI::DOUBLE,from,24);
      MPI::COMM_WORLD.Recv(packageWtaueta_prev,sizeOfData,MPI::DOUBLE,from,25);
      MPI::COMM_WORLD.Recv(packageWxx_prev,sizeOfData,MPI::DOUBLE,from,26);
      MPI::COMM_WORLD.Recv(packageWxy_prev,sizeOfData,MPI::DOUBLE,from,27);
      MPI::COMM_WORLD.Recv(packageWxeta_prev,sizeOfData,MPI::DOUBLE,from,28);
      MPI::COMM_WORLD.Recv(packageWyy_prev,sizeOfData,MPI::DOUBLE,from,29);
      MPI::COMM_WORLD.Recv(packageWyeta_prev,sizeOfData,MPI::DOUBLE,from,30);
      MPI::COMM_WORLD.Recv(packageWetaeta_prev,sizeOfData,MPI::DOUBLE,from,31);
      if(DATA->fluctuatingHydroFlag == 2){
	MPI::COMM_WORLD.Recv(package_dp,sizeOfData,MPI::DOUBLE,from,32);
	MPI::COMM_WORLD.Recv(package_dp_prev,sizeOfData,MPI::DOUBLE,from,33);
	MPI::COMM_WORLD.Recv(packagedutau,sizeOfData,MPI::DOUBLE,from,34);
	MPI::COMM_WORLD.Recv(packagedux,sizeOfData,MPI::DOUBLE,from,35);
	MPI::COMM_WORLD.Recv(packageduy,sizeOfData,MPI::DOUBLE,from,36);
	MPI::COMM_WORLD.Recv(packagedueta,sizeOfData,MPI::DOUBLE,from,37);
	MPI::COMM_WORLD.Recv(packagedutau_prev,sizeOfData,MPI::DOUBLE,from,38);
	MPI::COMM_WORLD.Recv(packagedux_prev,sizeOfData,MPI::DOUBLE,from,39);
	MPI::COMM_WORLD.Recv(packageduy_prev,sizeOfData,MPI::DOUBLE,from,40);
	MPI::COMM_WORLD.Recv(packagedueta_prev,sizeOfData,MPI::DOUBLE,from,41);
	MPI::COMM_WORLD.Recv(packagedWtautau,sizeOfData,MPI::DOUBLE,from,42);
	MPI::COMM_WORLD.Recv(packagedWtaux,sizeOfData,MPI::DOUBLE,from,43);
	MPI::COMM_WORLD.Recv(packagedWtauy,sizeOfData,MPI::DOUBLE,from,44);
	MPI::COMM_WORLD.Recv(packagedWtaueta,sizeOfData,MPI::DOUBLE,from,45);
	MPI::COMM_WORLD.Recv(packagedWxx,sizeOfData,MPI::DOUBLE,from,46);
	MPI::COMM_WORLD.Recv(packagedWxy,sizeOfData,MPI::DOUBLE,from,47);
	MPI::COMM_WORLD.Recv(packagedWxeta,sizeOfData,MPI::DOUBLE,from,48);
	MPI::COMM_WORLD.Recv(packagedWyy,sizeOfData,MPI::DOUBLE,from,49);
	MPI::COMM_WORLD.Recv(packagedWyeta,sizeOfData,MPI::DOUBLE,from,50);
	MPI::COMM_WORLD.Recv(packagedWetaeta,sizeOfData,MPI::DOUBLE,from,51);
	MPI::COMM_WORLD.Recv(packagedWtautau_prev,sizeOfData,MPI::DOUBLE,from,52);
	MPI::COMM_WORLD.Recv(packagedWtaux_prev,sizeOfData,MPI::DOUBLE,from,53);
	MPI::COMM_WORLD.Recv(packagedWtauy_prev,sizeOfData,MPI::DOUBLE,from,54);
	MPI::COMM_WORLD.Recv(packagedWtaueta_prev,sizeOfData,MPI::DOUBLE,from,55);
	MPI::COMM_WORLD.Recv(packagedWxx_prev,sizeOfData,MPI::DOUBLE,from,56);
	MPI::COMM_WORLD.Recv(packagedWxy_prev,sizeOfData,MPI::DOUBLE,from,57);
	MPI::COMM_WORLD.Recv(packagedWxeta_prev,sizeOfData,MPI::DOUBLE,from,58);
	MPI::COMM_WORLD.Recv(packagedWyy_prev,sizeOfData,MPI::DOUBLE,from,59);
	MPI::COMM_WORLD.Recv(packagedWyeta_prev,sizeOfData,MPI::DOUBLE,from,60);
	MPI::COMM_WORLD.Recv(packagedWetaeta_prev,sizeOfData,MPI::DOUBLE,from,61);
	//cout << " receiving from the right on rank " << rank << endl;
      }
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      position = ix + (nx*iy);
	      //if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
	      //cout << "Rneighbor[ix][iy][0].TJb[i][alpha][0]=" << package[position] << endl;
	      Rneighbor_eps[ix][iy] = package[position];
	      Rneighbor_eps_prev[ix][iy] = package_prev[position];
	      Rneighbor_utau[ix][iy] = packageutau[position];
	      Rneighbor_ux[ix][iy] = packageux[position];
	      Rneighbor_uy[ix][iy] = packageuy[position];
	      Rneighbor_ueta[ix][iy] = packageueta[position];
	      Rneighbor_rhob[ix][iy] = packagerhob[position];
	      Rneighbor_utau_prev[ix][iy] = packageutau_prev[position];
	      Rneighbor_ux_prev[ix][iy] = packageux_prev[position];
	      Rneighbor_uy_prev[ix][iy] = packageuy_prev[position];
	      Rneighbor_ueta_prev[ix][iy] = packageueta_prev[position];
	      Rneighbor_Wtautau[ix][iy] = packageWtautau[position];
	      Rneighbor_Wtaux[ix][iy] = packageWtaux[position];
	      Rneighbor_Wtauy[ix][iy] = packageWtauy[position];
	      Rneighbor_Wtaueta[ix][iy] = packageWtaueta[position];
	      Rneighbor_Wxx[ix][iy] = packageWxx[position];
	      Rneighbor_Wxy[ix][iy] = packageWxy[position];
	      Rneighbor_Wxeta[ix][iy] = packageWxeta[position];
	      Rneighbor_Wyy[ix][iy] = packageWyy[position];
	      Rneighbor_Wyeta[ix][iy] = packageWyeta[position];
	      Rneighbor_Wetaeta[ix][iy] = packageWetaeta[position];
	      Rneighbor_Wtautau_prev[ix][iy] = packageWtautau_prev[position];
	      Rneighbor_Wtaux_prev[ix][iy] = packageWtaux_prev[position];
	      Rneighbor_Wtauy_prev[ix][iy] = packageWtauy_prev[position];
	      Rneighbor_Wtaueta_prev[ix][iy] = packageWtaueta_prev[position];
	      Rneighbor_Wxx_prev[ix][iy] = packageWxx_prev[position];
	      Rneighbor_Wxy_prev[ix][iy] = packageWxy_prev[position];
	      Rneighbor_Wxeta_prev[ix][iy] = packageWxeta_prev[position];
	      Rneighbor_Wyy_prev[ix][iy] = packageWyy_prev[position];
	      Rneighbor_Wyeta_prev[ix][iy] = packageWyeta_prev[position];
	      Rneighbor_Wetaeta_prev[ix][iy] = packageWetaeta_prev[position];
	      if(DATA->fluctuatingHydroFlag == 2){
		Rneighbor_dp[ix][iy] = package_dp[position];
		Rneighbor_dp_prev[ix][iy] = package_dp_prev[position];
		Rneighbor_dutau[ix][iy] = packagedutau[position];
		Rneighbor_dux[ix][iy] = packagedux[position];
		Rneighbor_duy[ix][iy] = packageduy[position];
		Rneighbor_dueta[ix][iy] = packagedueta[position];
		Rneighbor_dutau_prev[ix][iy] = packagedutau_prev[position];
		Rneighbor_dux_prev[ix][iy] = packagedux_prev[position];
		Rneighbor_duy_prev[ix][iy] = packageduy_prev[position];
		Rneighbor_dueta_prev[ix][iy] = packagedueta_prev[position];
		Rneighbor_dWtautau[ix][iy] = packagedWtautau[position];
		Rneighbor_dWtaux[ix][iy] = packagedWtaux[position];
		Rneighbor_dWtauy[ix][iy] = packagedWtauy[position];
		Rneighbor_dWtaueta[ix][iy] = packagedWtaueta[position];
		Rneighbor_dWxx[ix][iy] = packagedWxx[position];
		Rneighbor_dWxy[ix][iy] = packagedWxy[position];
		Rneighbor_dWxeta[ix][iy] = packagedWxeta[position];
		Rneighbor_dWyy[ix][iy] = packagedWyy[position];
		Rneighbor_dWyeta[ix][iy] = packagedWyeta[position];
		Rneighbor_dWetaeta[ix][iy] = packagedWetaeta[position];
		Rneighbor_dWtautau_prev[ix][iy] = packagedWtautau_prev[position];
		Rneighbor_dWtaux_prev[ix][iy] = packagedWtaux_prev[position];
		Rneighbor_dWtauy_prev[ix][iy] = packagedWtauy_prev[position];
		Rneighbor_dWtaueta_prev[ix][iy] = packagedWtaueta_prev[position];
		Rneighbor_dWxx_prev[ix][iy] = packagedWxx_prev[position];
		Rneighbor_dWxy_prev[ix][iy] = packagedWxy_prev[position];
		Rneighbor_dWxeta_prev[ix][iy] = packagedWxeta_prev[position];
		Rneighbor_dWyy_prev[ix][iy] = packagedWyy_prev[position];
		Rneighbor_dWyeta_prev[ix][iy] = packagedWyeta_prev[position];
		Rneighbor_dWetaeta_prev[ix][iy] = packagedWetaeta_prev[position];
	      }
	    }
	}
      //cout << " done receiving from the right on rank " << rank << endl;
    }

  if (rank == size-1) maxEta = neta-fac;
  else maxEta = neta;
  
  for(ix=0; ix<=nx-fac; ix+=fac)
    {
      // fprintf(stderr,"IBIT[%d][%d]=%d\n",0,0,IBIT[0][0]);
      // fprintf(stderr,"IBIT[%d][%d]=%d\n",20,6,IBIT[20-1][6-1]);
      x = ix*(DATA->delta_x) - (DATA->x_size/2.0); 
      for(iy=0; iy<=ny-fac; iy+=fac)
	{
	  y = iy*(DATA->delta_y) - (DATA->y_size/2.0);
	  for(ieta=0; ieta<maxEta; ieta+=fac)
	    {
	      eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
	      //fprintf(stderr, "%d, %d, %d\n",ix,iy,ieta);
	      
	      if (ieta<neta-fac)
		{
		  intersect=1;
		  //if((arena[ix+fac][iy+fac][ieta+fac].epsilon-epsFO)*(arena[ix][iy][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
		  //  if((arena[ix+fac][iy][ieta].epsilon-epsFO)*(arena[ix][iy+fac][ieta+fac].epsilon_prev[facTau-1]-epsFO)>0.)
		  //    if((arena[ix][iy+fac][ieta].epsilon-epsFO)*(arena[ix+fac][iy][ieta+fac].epsilon_prev[facTau-1]-epsFO)>0.)
		  //	if((arena[ix][iy][ieta+fac].epsilon-epsFO)*(arena[ix+fac][iy+fac][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
		  //	  if((arena[ix+fac][iy+fac][ieta].epsilon-epsFO)*(arena[ix][iy][ieta+fac].epsilon_prev[facTau-1]-epsFO)>0.)
		  //	    if((arena[ix+fac][iy][ieta+fac].epsilon-epsFO)*(arena[ix][iy+fac][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
		  //	      if((arena[ix][iy+fac][ieta+fac].epsilon-epsFO)*(arena[ix+fac][iy][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
		  //		if((arena[ix][iy][ieta].epsilon-epsFO)*(arena[ix+fac][iy+fac][ieta+fac].epsilon_prev[facTau-1]-epsFO)>0.)
		  if((arena[ix+fac][iy+fac][ieta+fac].epsilon-epsFO)*(arena[ix][iy][ieta].epsilon_prev-epsFO)>0.)
		    if((arena[ix+fac][iy][ieta].epsilon-epsFO)*(arena[ix][iy+fac][ieta+fac].epsilon_prev-epsFO)>0.)
		      if((arena[ix][iy+fac][ieta].epsilon-epsFO)*(arena[ix+fac][iy][ieta+fac].epsilon_prev-epsFO)>0.)
			if((arena[ix][iy][ieta+fac].epsilon-epsFO)*(arena[ix+fac][iy+fac][ieta].epsilon_prev-epsFO)>0.)
			  if((arena[ix+fac][iy+fac][ieta].epsilon-epsFO)*(arena[ix][iy][ieta+fac].epsilon_prev-epsFO)>0.)
			    if((arena[ix+fac][iy][ieta+fac].epsilon-epsFO)*(arena[ix][iy+fac][ieta].epsilon_prev-epsFO)>0.)
			      if((arena[ix][iy+fac][ieta+fac].epsilon-epsFO)*(arena[ix+fac][iy][ieta].epsilon_prev-epsFO)>0.)
				if((arena[ix][iy][ieta].epsilon-epsFO)*(arena[ix+fac][iy+fac][ieta+fac].epsilon_prev-epsFO)>0.)
				  intersect=0;
		}
	      else // if this is the right most edge
		{
		  intersect=1;
		  if((Rneighbor_eps[ix+fac][iy+fac] - epsFO)*(arena[ix][iy][ieta].epsilon_prev-epsFO)>0.)
		    //if((Rneighbor_eps[ix+fac][iy+fac] - epsFO)*(arena[ix][iy][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
		    if((arena[ix+fac][iy][ieta].epsilon-epsFO)*(Rneighbor_eps_prev[ix][iy+fac]-epsFO)>0.)
		      if((arena[ix][iy+fac][ieta].epsilon-epsFO)*(Rneighbor_eps_prev[ix+fac][iy]-epsFO)>0.)
			//if((Rneighbor_eps[ix][iy]-epsFO)*(arena[ix+fac][iy+fac][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
			if((Rneighbor_eps[ix][iy]-epsFO)*(arena[ix+fac][iy+fac][ieta].epsilon_prev-epsFO)>0.)
			  if((arena[ix+fac][iy+fac][ieta].epsilon-epsFO)*(Rneighbor_eps_prev[ix][iy]-epsFO)>0.)
			    //if((Rneighbor_eps[ix+fac][iy]-epsFO)*(arena[ix][iy+fac][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
			    //if((Rneighbor_eps[ix][iy+fac]-epsFO)*(arena[ix+fac][iy][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
			    //	if((arena[ix][iy][ieta].epsilon-epsFO)*(Rneighbor_eps_prev[ix+fac][iy+fac]-epsFO)>0.)
			    if((Rneighbor_eps[ix+fac][iy]-epsFO)*(arena[ix][iy+fac][ieta].epsilon_prev-epsFO)>0.)
			      if((Rneighbor_eps[ix][iy+fac]-epsFO)*(arena[ix+fac][iy][ieta].epsilon_prev-epsFO)>0.)
				if((arena[ix][iy][ieta].epsilon-epsFO)*(Rneighbor_eps_prev[ix+fac][iy+fac]-epsFO)>0.)
				  intersect=0;
		}

	      //fprintf(stderr, "done checks\n");
	      if (intersect==0)
		{
		  //fprintf(stderr, "continue...\n");
		  continue;
		}
	      else		 
		{
		  //fprintf(stderr, "starting freeze out\n");
		  intersections++;
		  if (ieta<neta-fac)
		    {
		      //cube[0]=arena[ix][iy][ieta].epsilon_prev[facTau-1];
		      //cube[1]=arena[ix+fac][iy][ieta].epsilon_prev[facTau-1];
		      //cube[2]=arena[ix+fac][iy+fac][ieta].epsilon_prev[facTau-1];
		      //cube[3]=arena[ix][iy+fac][ieta].epsilon_prev[facTau-1];
		      cube[0]=arena[ix][iy][ieta].epsilon_prev;
		      cube[1]=arena[ix+fac][iy][ieta].epsilon_prev;
		      cube[2]=arena[ix+fac][iy+fac][ieta].epsilon_prev;
		      cube[3]=arena[ix][iy+fac][ieta].epsilon_prev;
		      cube[4]=arena[ix][iy][ieta].epsilon;
		      cube[5]=arena[ix+fac][iy][ieta].epsilon;
		      cube[6]=arena[ix+fac][iy+fac][ieta].epsilon;
		      cube[7]=arena[ix][iy+fac][ieta].epsilon;
		      //cube[8]=arena[ix][iy][ieta+fac].epsilon_prev[facTau-1];
		      //cube[9]=arena[ix+fac][iy][ieta+fac].epsilon_prev[facTau-1];
		      //cube[10]=arena[ix+fac][iy+fac][ieta+fac].epsilon_prev[facTau-1];
		      //cube[11]=arena[ix][iy+fac][ieta+fac].epsilon_prev[facTau-1];
		      cube[8]=arena[ix][iy][ieta+fac].epsilon_prev;
		      cube[9]=arena[ix+fac][iy][ieta+fac].epsilon_prev;
		      cube[10]=arena[ix+fac][iy+fac][ieta+fac].epsilon_prev;
		      cube[11]=arena[ix][iy+fac][ieta+fac].epsilon_prev;
		      cube[12]=arena[ix][iy][ieta+fac].epsilon;
		      cube[13]=arena[ix+fac][iy][ieta+fac].epsilon;
		      cube[14]=arena[ix+fac][iy+fac][ieta+fac].epsilon;
		      cube[15]=arena[ix][iy+fac][ieta+fac].epsilon;
		    }
		  else
		    {
		      //cube[0]=arena[ix][iy][ieta].epsilon_prev[facTau-1];
		      //cube[1]=arena[ix+fac][iy][ieta].epsilon_prev[facTau-1];
		      //cube[2]=arena[ix+fac][iy+fac][ieta].epsilon_prev[facTau-1];
		      //cube[3]=arena[ix][iy+fac][ieta].epsilon_prev[facTau-1];
		      cube[0]=arena[ix][iy][ieta].epsilon_prev;
		      cube[1]=arena[ix+fac][iy][ieta].epsilon_prev;
		      cube[2]=arena[ix+fac][iy+fac][ieta].epsilon_prev;
		      cube[3]=arena[ix][iy+fac][ieta].epsilon_prev;
		      cube[4]=arena[ix][iy][ieta].epsilon;
		      cube[5]=arena[ix+fac][iy][ieta].epsilon;
		      cube[6]=arena[ix+fac][iy+fac][ieta].epsilon;
		      cube[7]=arena[ix][iy+fac][ieta].epsilon;
		      cube[8]=Rneighbor_eps_prev[ix][iy];
		      cube[9]=Rneighbor_eps_prev[ix+fac][iy];
		      cube[10]=Rneighbor_eps_prev[ix+fac][iy+fac];
		      cube[11]=Rneighbor_eps_prev[ix][iy+fac];
		      cube[12]=Rneighbor_eps[ix][iy];
		      cube[13]=Rneighbor_eps[ix+fac][iy];
		      cube[14]=Rneighbor_eps[ix+fac][iy+fac];
		      cube[15]=Rneighbor_eps[ix][iy+fac];
		    }

		  /* 	  // set values for testing: */
/* 		  cube[0]=epsFO+0.01; // 0    0   0   0 */
/* 		  cube[1]=epsFO-0.01; // 0    dx  0   0 */
/* 		  cube[2]=epsFO-0.01; // 0    dx  dy  0 */
/* 		  cube[3]=epsFO+0.01; // 0    0   dy  0 */
/* 		  cube[4]=epsFO+0.01; // dtau 0   0   0 */
/* 		  cube[5]=epsFO-0.01; // dtau dx  0   0 */
/* 		  cube[6]=epsFO-0.01; // dtau dx  dy  0 */
/* 		  cube[7]=epsFO+0.01; // dtau 0   dy  0 */
/* 		  cube[8]=epsFO+0.01; // 0    0   0   deta */
/* 		  cube[9]=epsFO-0.01; // 0    dx  0   deta */
/* 		  cube[10]=epsFO-0.01;// 0    dx  dy  deta */
/* 		  cube[11]=epsFO+0.01;// 0    0   dy  deta */
/* 		  cube[12]=epsFO+0.01;// dtau 0   0   deta */
/* 		  cube[13]=epsFO-0.01;// dtau dx  0   deta */
/* 		  cube[14]=epsFO-0.01;// dtau dx  dy  deta */
/* 		  cube[15]=epsFO+0.01;// dtau 0   dy  deta */
		  
                  // set values for testing:
/* 		  cube[0]=epsFO+0.01; // 0    0   0   0 */
/* 		  cube[1]=epsFO-0.01; // 0    dx  0   0 */
/* 		  cube[2]=epsFO-0.01; // 0    dx  dy  0 */
/* 		  cube[3]=epsFO-0.01; // 0    0   dy  0 */
/* 		  cube[4]=epsFO+0.01; // dtau 0   0   0 */
/* 		  cube[5]=epsFO+0.01; // dtau dx  0   0 */
/* 		  cube[6]=epsFO+0.01; // dtau dx  dy  0 */
/* 		  cube[7]=epsFO+0.01; // dtau 0   dy  0 */
/* 		  cube[8]=epsFO+0.01; // 0    0   0   deta */
/* 		  cube[9]=epsFO-0.01; // 0    dx  0   deta */
/* 		  cube[10]=epsFO-0.01;// 0    dx  dy  deta */
/* 		  cube[11]=epsFO-0.01;// 0    0   dy  deta */
/* 		  cube[12]=epsFO+0.01;// dtau 0   0   deta */
/* 		  cube[13]=epsFO+0.01;// dtau dx  0   deta */
/* 		  cube[14]=epsFO+0.01;// dtau dx  dy  deta */
/* 		  cube[15]=epsFO+0.01;// dtau 0   dy  deta */
		  
		  NSurfaces = 0;
		  ELowerSum = 0.;
		  EHigherSum = 0.;
		  VLower0=0.;
		  VHigher0=0.;
		  VLower1=0.;
		  VHigher1=0.;
		  VLower2=0.;
		  VHigher2=0.;
		  VLower3=0.;
		  VHigher3=0.;


		  // edge 1:
		  EK = cube[0];
		  EL = cube[1];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEK = fabs(DEK);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 1;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = DEK/(EL-EK)*DX;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = 0.;
		    }
		  if ( DEK > 0 )
		    ELowerSum+=ADEK;
		  else
		    EHigherSum+=ADEK;

		  // edge 2: v for point 1
		  EK = cube[1];
		  EL = cube[2];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEK = fabs(DEK);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 2;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = DEK/(EL-EK)*DY;
		      cuts[NSurfaces-1][3] = 0.;
		    }
		  if ( DEK > 0 )
		    {
		      ELowerSum+=ADEK;
		      VLower1+=ADEK;
		    }
		  else
		    {
		      EHigherSum+=ADEK;
		      VHigher1+=ADEK;
		    }
		  
		  // edge 3: v for point 2
		  EK = cube[3];
		  EL = cube[2];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 3;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = DEK/(EL-EK)*DX;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = 0.;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower1+=ADEL;
		      VLower2+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher1+=ADEL;
		      VHigher2+=ADEL;
		    }

		  // edge 4: v for point 3
		  EK = cube[0];
		  EL = cube[3];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 4;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = DEK/(EL-EK)*DY;
		      cuts[NSurfaces-1][3] = 0.;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower2+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher2+=ADEL;
		    }
		
		  // edge 5: v for point 4
		  EK = cube[0];
		  EL = cube[4];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 5;
		      cuts[NSurfaces-1][0] = DEK/(EL-EK)*DTAU;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = 0.;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower0+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher0+=ADEL;
		    }


		  // edge 6: v for point 5
		  EK = cube[1];
		  EL = cube[5];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 6;
		      cuts[NSurfaces-1][0] = DEK/(EL-EK)*DTAU;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = 0.;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower0+=ADEL;
		      VLower1+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher0+=ADEL;
		      VHigher1+=ADEL;
		    }

		  // edge 7: v for point 6
		  EK = cube[2];
		  EL = cube[6];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 7;
		      cuts[NSurfaces-1][0] = DEK/(EL-EK)*DTAU;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = 0.;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower0+=ADEL;
		      VLower1+=ADEL;
		      VLower2+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher0+=ADEL;
		      VHigher1+=ADEL;
		      VHigher2+=ADEL;
		    }

		  // edge 8: v for point 7
		  EK = cube[3];
		  EL = cube[7];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 8;
		      cuts[NSurfaces-1][0] = DEK/(EL-EK)*DTAU;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = 0.;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower0+=ADEL;
		      VLower2+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher0+=ADEL;
		      VHigher2+=ADEL;
		    }

		  // edge 9:
		  EK = cube[4];
		  EL = cube[5];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 9;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = DEK/(EL-EK)*DX;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = 0.;
		    }

		  // edge 10:
		  EK = cube[5];
		  EL = cube[6];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 10;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = DEK/(EL-EK)*DY;
		      cuts[NSurfaces-1][3] = 0.;
		    }

		  // edge 11:
		  EK = cube[7];
		  EL = cube[6];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 11;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = DEK/(EL-EK)*DX;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = 0.;
		    }

		  // edge 12:
		  EK = cube[4];
		  EL = cube[7];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 12;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = DEK/(EL-EK)*DY;
		      cuts[NSurfaces-1][3] = 0.;
		    }

		  // edge 13: v for point 8
		  EK = cube[8];
		  EL = cube[9];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEK = fabs(DEK);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 13;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = DEK/(EL-EK)*DX;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = DETA;
		    }
		  if ( DEK > 0 )
		    {
		      VLower3+=ADEK;
		      ELowerSum+=ADEK;
		    }
		  else
		    {
		      VHigher3+=ADEK;
		      EHigherSum+=ADEK;
		    }

		  // edge 14: v for point 9
		  EK = cube[9];
		  EL = cube[10];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEK = fabs(DEK);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 14;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = DEK/(EL-EK)*DY;
		      cuts[NSurfaces-1][3] = DETA;
		    }
		  if ( DEK > 0 )
		    {
		      ELowerSum+=ADEK;
		      VLower1+=ADEK;
		      VLower3+=ADEK;
		    }
		  else
		    {
		      EHigherSum+=ADEK;
		      VHigher1+=ADEK;
		      VHigher3+=ADEK;
		    }
		  
		  // edge 15: v for point 10
		  EK = cube[11];
		  EL = cube[10];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 15;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = DEK/(EL-EK)*DX;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = DETA;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower1+=ADEL;
		      VLower2+=ADEL;
		      VLower3+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher1+=ADEL;
		      VHigher2+=ADEL;
		      VHigher3+=ADEL;
		    }

		  // edge 16: v for point 11
		  EK = cube[8];
		  EL = cube[11];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 16;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = DEK/(EL-EK)*DY;
		      cuts[NSurfaces-1][3] = DETA;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower2+=ADEL;
		      VLower3+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher2+=ADEL;
		      VHigher3+=ADEL;
		    }
		
		  // edge 17: v for point 12
		  EK = cube[8];
		  EL = cube[12];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 17;
		      cuts[NSurfaces-1][0] = DEK/(EL-EK)*DTAU;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = DETA;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower0+=ADEL;
		      VLower3+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher0+=ADEL;
		      VHigher3+=ADEL;
		    }


		  // edge 18: v for point 13
		  EK = cube[9];
		  EL = cube[13];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 18;
		      cuts[NSurfaces-1][0] = DEK/(EL-EK)*DTAU;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = DETA;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower0+=ADEL;
		      VLower1+=ADEL;
		      VLower3+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher0+=ADEL;
		      VHigher1+=ADEL;
		      VHigher3+=ADEL;
		    }

		  // edge 19: v for point 14
		  EK = cube[10];
		  EL = cube[14];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 19;
		      cuts[NSurfaces-1][0] = DEK/(EL-EK)*DTAU;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = DETA;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower0+=ADEL;
		      VLower1+=ADEL;
		      VLower2+=ADEL;
		      VLower3+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher0+=ADEL;
		      VHigher1+=ADEL;
		      VHigher2+=ADEL;
		      VHigher3+=ADEL;
		    }

		  // edge 20: v for point 15
		  EK = cube[11];
		  EL = cube[15];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 20;
		      cuts[NSurfaces-1][0] = DEK/(EL-EK)*DTAU;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = DETA;
		    }
		  if ( DEL > 0 )
		    {
		      ELowerSum+=ADEL;
		      VLower0+=ADEL;
		      VLower2+=ADEL;
		      VLower3+=ADEL;
		    }
		  else
		    {
		      EHigherSum+=ADEL;
		      VHigher0+=ADEL;
		      VHigher2+=ADEL;
		      VHigher3+=ADEL;
		    }

		  // edge 21:
		  EK = cube[12];
		  EL = cube[13];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 21;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = DEK/(EL-EK)*DX;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = DETA;
		    }

		  // edge 22:
		  EK = cube[13];
		  EL = cube[14];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 22;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = DEK/(EL-EK)*DY;
		      cuts[NSurfaces-1][3] = DETA;
		    }

		  // edge 23:
		  EK = cube[15];
		  EL = cube[14];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 23;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = DEK/(EL-EK)*DX;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = DETA;
		    }

		  // edge 24:
		  EK = cube[12];
		  EL = cube[15];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 24;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = DEK/(EL-EK)*DY;
		      cuts[NSurfaces-1][3] = DETA;
		    }

		  // edge 25:
		  EK = cube[0];
		  EL = cube[8];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 25;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = DEK/(EL-EK)*DETA;
		    }

		  // edge 26:
		  EK = cube[1];
		  EL = cube[9];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 26;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = DEK/(EL-EK)*DETA;
		    }

		  // edge 27:
		  EK = cube[2];
		  EL = cube[10];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 27;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = DEK/(EL-EK)*DETA;
		    }

		  // edge 28:
		  EK = cube[3];
		  EL = cube[11];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 28;
		      cuts[NSurfaces-1][0] = 0.;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = DEK/(EL-EK)*DETA;
		    }

		  // edge 29:
		  EK = cube[4];
		  EL = cube[12];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 29;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = DEK/(EL-EK)*DETA;
		    }

		  // edge 30:
		  EK = cube[5];
		  EL = cube[13];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 30;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = 0.;
		      cuts[NSurfaces-1][3] = DEK/(EL-EK)*DETA;
		    }

		  // edge 31:
		  EK = cube[6];
		  EL = cube[14];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 31;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = DX;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = DEK/(EL-EK)*DETA;
		    }

		  // edge 32:
		  EK = cube[7];
		  EL = cube[15];
		  DEK = epsFO-EK;
		  DEL = epsFO-EL;
		  ADEL = fabs(DEL);
		  if( DEK*DEL<0 )
		    {
		      NSurfaces += 1;
		      iEdge[NSurfaces-1] = 32;
		      cuts[NSurfaces-1][0] = DTAU;
		      cuts[NSurfaces-1][1] = 0.;
		      cuts[NSurfaces-1][2] = DY;
		      cuts[NSurfaces-1][3] = DEK/(EL-EK)*DETA;
		    }


//______________________got all intersection points_________________________

		  if (ELowerSum!=0)
		    { 
		      VLower0=VLower0/ELowerSum*DTAU;
		      VLower1=VLower1/ELowerSum*DX;
		      VLower2=VLower2/ELowerSum*DY;
		      VLower3=VLower3/ELowerSum*DETA;
		    }
		  if (EHigherSum!=0)
		    { 
		      VHigher0=VHigher0/EHigherSum*DTAU;
		      VHigher1=VHigher1/EHigherSum*DX;
		      VHigher2=VHigher2/EHigherSum*DY;
		      VHigher3=VHigher3/EHigherSum*DETA;
		    }

		  VD0 = VLower0-VHigher0;
		  VD1 = VLower1-VHigher1;
		  VD2 = VLower2-VHigher2;
		  VD3 = VLower3-VHigher3;
		
		  // compute the mean vector of intersection points:
		  
		  V0=0.0;
		  V1=0.0;
		  V2=0.0;
		  V3=0.0;
 
		  for (is=0; is<NSurfaces; is++)
		    {
		      V0+=cuts[is][0];
		      V1+=cuts[is][1];
		      V2+=cuts[is][2];
		      V3+=cuts[is][3];
		    }

		  VMID[0]=V0/NSurfaces;
		  VMID[1]=V1/NSurfaces;
		  VMID[2]=V2/NSurfaces;
		  VMID[3]=V3/NSurfaces;
 
		  xf = x + VMID[1];
		  yf = y + VMID[2];
		  etaf = eta + VMID[3];
		  tauf = tau - DTAU + VMID[0];
		  
		  
		  //fprintf(stderr, "NUMBER OF HYPERSURFACES = %d\n",NSurfaces);
		  if (NSurfaces>12 || NSurfaces%2>0) fprintf(stderr, "********************NUMBER OF HYPERSURFACES = %d\n",NSurfaces);
		  
		  // find the neighbors
		  tries3=0;
		  COUNTER=0;
		  while(((NSurfaces==10 && (COUNTER!=16 && COUNTER!=15)) || (NSurfaces==12 && (COUNTER!=19 && COUNTER!=20 && COUNTER!=21))
			||(NSurfaces==4 && (COUNTER!=4)) || (NSurfaces==6 && (COUNTER!=8))
			 ||(NSurfaces==8 && (COUNTER!=12))) && tries3<2*NSurfaces)
		    {
		 
		      for (i=0; i<NSurfaces; i++)
			{
			  for (j=0; j<6; j++)
			    {
			      additionalNeighbors[i][j]=0;
			      neighborsDone[i][j]=0;
			    }
			}
		     
		      for(is=0; is<32; is++)
			for(i=0; i<32; i++)
			  numberConnectionIsUsed[is][i] = -1;
	
		      //if (tries3>0) sleep(1);
		      tries3+=1;
		      //fprintf(stderr,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++tries3=%d\n",tries3);
		      for (m=0; m<NSurfaces; m++)
			{
			  i=abs((m+tries3)%NSurfaces);
			  IE = iEdge[i];
			  IS = 0;
			  for (m2=0; m2<NSurfaces; m2++)
			    {
			      if (tries3<10) j=m2;//abs((m2+tries)%NSurfaces);
			      else 
				j=abs(m2+2*tries3)%NSurfaces;
			      JE = iEdge[j];
			      if ( IBIT[IE-1][JE-1]==1 ) 
				{
				  IS = IS+1;
				  ISID[IS-1] = j;
				  neighbors[i][IS-1]=JE; // determines the neighbors of ith edge (which is iEdge[i])
				}
			    }
			  for (m2=0; m2<NSurfaces; m2++)
			    {
			      j=abs((m2+3*tries3)%NSurfaces);
			      JE = iEdge[j];
			      if ( IS<6 && IBIT[IE-1][JE-1]==2 ) 
				{
				  IS = IS+1;
				  ISID[IS-1] = j;
				  neighbors[i][IS-1]=JE; // determines the neighbors of ith edge (which is iEdge[i])
				}
			    }
			  if (IS!=3) 
			    {
			      fprintf(stderr, "WARNING: total number of neighbors=%d\n",IS);
			      //sleep(1);
			      if (tries3>10) break;
			    }
			}
		      
		      
		      FULLSU[0]=0.;
		      FULLSU[1]=0.;
		      FULLSU[2]=0.;
		      FULLSU[3]=0.;
		      
	      
		      // do the first group of tetraedra that do not share any edges:
		      COUNTER=0;
		      countEdges=0;
		      for (is=0; is<NSurfaces; is++)
			{
			  skip=0;
			  if (is>0)
			    for (j=0; j<countEdges; j++)
			      {
				if (iEdge[is]==neighbors[previousEdges[j]][0] || iEdge[is]==neighbors[previousEdges[j]][1] 
				    || iEdge[is]==neighbors[previousEdges[j]][2])
				  {
				    //fprintf(stderr,"************skipping edge %d\n",iEdge[is]);
				    skip=1;
				    break;
				  }
			      }
			  if (skip==1) continue;
			  // fprintf(stderr,"is=%d, doing edge %d\n",is, iEdge[is]);
			  
			  previousEdges[countEdges]=is;
			  countEdges+=1;
			  IS=0;
			  for(k=0;k<4;k++)
			    {
			      BD[k] = VMID[k]-cuts[is][k];
			    }
			  for (j=0; j<3; j++)
			    {
			      for (i=0; i<NSurfaces; i++)
				{
				  if (neighbors[is][j]==iEdge[i])
				    {
				      neighborsDone[is][j]=1;
				      for (l=0; l<3; l++)
					{
					  if (neighbors[i][l] == iEdge[is]) 
					    {
					      neighborsDone[i][l]=1;
					    }
					}
				      IS = IS+1;
				      prevEdge[IS-1]=iEdge[i];
				      if (iEdge[is]>iEdge[i])
					numberConnectionIsUsed[iEdge[is]-1][iEdge[i]-1]+=3;
				      else
					numberConnectionIsUsed[iEdge[i]-1][iEdge[is]-1]+=3;
				      
				      
				      for(k=0;k<4;k++)
					{
					  AD[IS-1][k] = cuts[is][k]-cuts[i][k];
					}
				    }
				}
			    }
			  // update the number of connections for every used edge:
			  for (m=0; m<IS; m++)
			    {
			      for (m2=0; m2<m; m2++)
				{
				  if (prevEdge[m]>prevEdge[m2])
				    if (numberConnectionIsUsed[prevEdge[m]-1][prevEdge[m2]-1]==-1)
				      numberConnectionIsUsed[prevEdge[m]-1][prevEdge[m2]-1]+=2;
				    else
				      numberConnectionIsUsed[prevEdge[m]-1][prevEdge[m2]-1]+=1;
				  else
				    if (numberConnectionIsUsed[prevEdge[m2]-1][prevEdge[m]-1]==-1)
				      numberConnectionIsUsed[prevEdge[m2]-1][prevEdge[m]-1]+=2;
				    else
				      numberConnectionIsUsed[prevEdge[m2]-1][prevEdge[m]-1]+=1;
				}
			    }
			  
			  if (IS==3)
			    {
			      COUNTER+=3;
			      
			      SU[COUNTER-3][0] = AD[0][1]*AD[1][2]*BD[3] - AD[0][1]*AD[1][3]*BD[2] + AD[0][3]*AD[1][1]*BD[2]
				-  AD[0][3]*AD[1][2]*BD[1] + AD[0][2]*AD[1][3]*BD[1] - AD[0][2]*AD[1][1]*BD[3];
			      SU[COUNTER-3][1] = AD[0][2]*AD[1][0]*BD[3] - AD[0][0]*AD[1][2]*BD[3] + AD[0][0]*AD[1][3]*BD[2]
				-  AD[0][3]*AD[1][0]*BD[2] + AD[0][3]*AD[1][2]*BD[0] - AD[0][2]*AD[1][3]*BD[0];
			      SU[COUNTER-3][2] = AD[0][0]*AD[1][1]*BD[3] - AD[0][0]*AD[1][3]*BD[1] + AD[0][3]*AD[1][0]*BD[1]
				-  AD[0][3]*AD[1][1]*BD[0] + AD[0][1]*AD[1][3]*BD[0] - AD[0][1]*AD[1][0]*BD[3];
			      SU[COUNTER-3][3] = AD[0][2]*AD[1][1]*BD[0] - AD[0][1]*AD[1][2]*BD[0] + AD[0][0]*AD[1][2]*BD[1]
				-  AD[0][2]*AD[1][0]*BD[1] + AD[0][1]*AD[1][0]*BD[2] - AD[0][0]*AD[1][1]*BD[2];
			      
			      SU[COUNTER-2][0] = AD[0][1]*AD[2][2]*BD[3] - AD[0][1]*AD[2][3]*BD[2] + AD[0][3]*AD[2][1]*BD[2]
				-  AD[0][3]*AD[2][2]*BD[1] + AD[0][2]*AD[2][3]*BD[1] - AD[0][2]*AD[2][1]*BD[3];
			      SU[COUNTER-2][1] = AD[0][2]*AD[2][0]*BD[3] - AD[0][0]*AD[2][2]*BD[3] + AD[0][0]*AD[2][3]*BD[2]
				-  AD[0][3]*AD[2][0]*BD[2] + AD[0][3]*AD[2][2]*BD[0] - AD[0][2]*AD[2][3]*BD[0];
			      SU[COUNTER-2][2] = AD[0][0]*AD[2][1]*BD[3] - AD[0][0]*AD[2][3]*BD[1] + AD[0][3]*AD[2][0]*BD[1]
				-  AD[0][3]*AD[2][1]*BD[0] + AD[0][1]*AD[2][3]*BD[0] - AD[0][1]*AD[2][0]*BD[3];
			      SU[COUNTER-2][3] = AD[0][2]*AD[2][1]*BD[0] - AD[0][1]*AD[2][2]*BD[0] + AD[0][0]*AD[2][2]*BD[1]
				-  AD[0][2]*AD[2][0]*BD[1] + AD[0][1]*AD[2][0]*BD[2] - AD[0][0]*AD[2][1]*BD[2];
			      
			      SU[COUNTER-1][0] = AD[1][1]*AD[2][2]*BD[3] - AD[1][1]*AD[2][3]*BD[2] + AD[1][3]*AD[2][1]*BD[2]
				-  AD[1][3]*AD[2][2]*BD[1] + AD[1][2]*AD[2][3]*BD[1] - AD[1][2]*AD[2][1]*BD[3];
			      SU[COUNTER-1][1] = AD[1][2]*AD[2][0]*BD[3] - AD[1][0]*AD[2][2]*BD[3] + AD[1][0]*AD[2][3]*BD[2]
				-  AD[1][3]*AD[2][0]*BD[2] + AD[1][3]*AD[2][2]*BD[0] - AD[1][2]*AD[2][3]*BD[0];
			      SU[COUNTER-1][2] = AD[1][0]*AD[2][1]*BD[3] - AD[1][0]*AD[2][3]*BD[1] + AD[1][3]*AD[2][0]*BD[1]
				-  AD[1][3]*AD[2][1]*BD[0] + AD[1][1]*AD[2][3]*BD[0] - AD[1][1]*AD[2][0]*BD[3];
			      SU[COUNTER-1][3] = AD[1][2]*AD[2][1]*BD[0] - AD[1][1]*AD[2][2]*BD[0] + AD[1][0]*AD[2][2]*BD[1]
				-  AD[1][2]*AD[2][0]*BD[1] + AD[1][1]*AD[2][0]*BD[2] - AD[1][0]*AD[2][1]*BD[2];
			      
			      for (k=0; k<4;k++)
				{
				  SU[COUNTER-3][k]/=6.; // to get tetraedron's volume
				  SU[COUNTER-2][k]/=6.; // to get tetraedron's volume
				  SU[COUNTER-1][k]/=6.; // to get tetraedron's volume
				}
			      
			    }
			}
		      
		      for(j=0; j<32; j++)
			{
			  singleConnections[j][0]=0;
			  singleConnections[j][1]=0;
			  singleConnectionsUsed[j]=0;
			}

		      countSingleEdges=0;
		      for(j=1; j<=32; j++)
			for (k=1; k<j; k++)
			  {
			    if ( numberConnectionIsUsed[j-1][k-1]>-1 )
			      {
				// fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]);
				if (numberConnectionIsUsed[j-1][k-1]<2)
				  {
				    if (numberConnectionIsUsed[j-1][k-1]==0)
				      countSingleEdges+=2;
				    else
				      countSingleEdges+=1;
				    singleConnections[countSingleEdges-1][0]=j;
				    singleConnections[countSingleEdges-1][1]=k;
				    //fprintf(stderr,"-------------------------------------------\n");
				  }
				if (numberConnectionIsUsed[j-1][k-1]>2)
				  {
				    fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]);
				    fprintf(stderr,"Three connections is messed up. That can't happen.\n");
				    break;
				  }
			      }
			  }

		         
		      // find additional connections between edges of the previously defined tetraedra
		      additionalEdges=1; // just to get started in the while loop at all
		      // do the missed edges:
		      countSingleEdges=1; // just to get started in the while loop at all
		      shift=0;
		      //addCOUNTER=0;
		      tries2=0;
		   /*    for(j=1; j<=32; j++) */
/* 			for (k=1; k<j; k++) */
/* 			  { */
/* 			    if ( numberConnectionIsUsed[j-1][k-1]>-1 ) */
/* 			      { */
/* 				//fprintf(stderr,"BEFORE - connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]); */
/* 				if (numberConnectionIsUsed[j-1][k-1]<2) */
/* 				  { */
/* 				    fprintf(stderr,"-------------------------------------------\n"); */
/* 				  } */
/* 			      } */
/* 			  } */
		      
		      for (i=0; i<NSurfaces; i++)
			{
			  for (j=0; j<3; j++)
			    {
			      if (additionalNeighbors[i][j]==1) neighborsDone[i][j]=0;
			      additionalNeighbors[i][j]=0;
			    }
			}
		      //fprintf(stderr,"NUMBER OF SURFACES=%d\n",NSurfaces);
		      //fprintf(stderr,"COUNTER=%d\n",COUNTER);
		      //fprintf(stderr,"additionalEdges=%d\n",additionalEdges);
		      //fprintf(stderr,"countSingleEdges=%d\n",countSingleEdges);
		      //fprintf(stderr,"shift=%d\n",shift);
		      additionalEdges=0;
		      IS=0;
		      for (m=0; m<NSurfaces; m++)
			{
			  is = abs((m-shift)%NSurfaces);
			  for (j=0; j<3; j++)
			    {
			      //fprintf(stderr,"%d's neighbor %d done = %d\n",iEdge[is],neighbors[is][j],neighborsDone[is][j]);
			      if (neighborsDone[is][j] == 0)
				{
				  for (m2=0; m2<NSurfaces; m2++)
				    {
				      i = (tries3+shift+m2)%NSurfaces;
				      for (l=0; l<3; l++)
					{
					  if (neighborsDone[i][l] == 0 && neighbors[is][j]==iEdge[i] && neighbors[i][l]==iEdge[is])
					    {
					      additionalEdges+=1;
					      additionalNeighbors[is][j]=1;
					      additionalNeighbors[i][j]=1;
					      neighborsDone[is][j]=1;
					      neighborsDone[i][l]=1;
					      IS=IS+1;
					      edge1[IS-1]=iEdge[is];
					      edge2[IS-1]=iEdge[i];
					      //edge1Done[IS-1]=0;
					      //edge2Done[IS-1]=0;
					      //fprintf(stderr,"doing for edge %d and %d\n",iEdge[is],iEdge[i]);
					      //fprintf(stderr,"is=%d, i=%d, IS=%d\n",is,i,IS);
					      //fprintf(stderr,"edge1=%d, edge2=%d \n",edge1[IS-1],edge2[IS-1]);
					      
					      for(k=0;k<4;k++)
						{
						  DD[IS-1][k] = VMID[k]-cuts[is][k];
						  AD[IS-1][k] = cuts[is][k]-cuts[i][k];
						}
					    }
					}
				    }
				}
			    }
			}
		      
		      countAdditionalEdges=0;
		      //fprintf(stderr,"IS=%d\n",IS);
		      
		      
		      for (m=0; m<IS; m++)
			{
			  srand(time(NULL));
			  if (tries3<5) i=m;
			  else
			    i=abs(m-tries3*2)%IS;
			  for (m2=0; m2<IS; m2++)
			    {
			      if (tries3==0) j=m2;
			      else
				j=abs(m2+tries3)%IS;
			      //fprintf(stderr,"j=%d\n",j);
			      if ( additionalEdges>0 &&
				   (
				    ((edge1[i] == edge2[j]) || (edge2[i] == edge1[j]) || (edge1[i] == edge1[j])|| (edge2[i] == edge2[j]))))
				{
				  if (edge1[i] == edge1[j] && edge2[i] == edge2[j]) continue;
				  if (edge1[i] == edge2[j] && edge2[i] == edge1[j]) continue;
				  
				  usedEdges1[countAdditionalEdges]=edge1[i];
				  usedEdges2[countAdditionalEdges]=edge1[j];
				  usedEdges3[countAdditionalEdges]=edge2[i];
				  usedEdges4[countAdditionalEdges]=edge2[j];
				  
				  // sort them
				  if (usedEdges3[countAdditionalEdges]>usedEdges4[countAdditionalEdges])
				    {
				      temp = usedEdges3[countAdditionalEdges];
				      usedEdges3[countAdditionalEdges] = usedEdges4[countAdditionalEdges];
				      usedEdges4[countAdditionalEdges] = temp;
				    }
				  if (usedEdges2[countAdditionalEdges]>usedEdges3[countAdditionalEdges])
				    {
				      temp = usedEdges2[countAdditionalEdges];
				      usedEdges2[countAdditionalEdges] = usedEdges3[countAdditionalEdges];
				      usedEdges3[countAdditionalEdges] = temp;
				    }
				  if (usedEdges1[countAdditionalEdges]>usedEdges2[countAdditionalEdges])
				    {
				      temp = usedEdges1[countAdditionalEdges];
				      usedEdges1[countAdditionalEdges] = usedEdges2[countAdditionalEdges];
				      usedEdges2[countAdditionalEdges] = temp;
				    }
				  
				  // remove double one
				  if (usedEdges2[countAdditionalEdges] == usedEdges3[countAdditionalEdges])
				    usedEdges3[countAdditionalEdges] = usedEdges4[countAdditionalEdges];
				  if (usedEdges1[countAdditionalEdges] == usedEdges2[countAdditionalEdges])
				    {
				      usedEdges2[countAdditionalEdges] = usedEdges3[countAdditionalEdges];
				      usedEdges3[countAdditionalEdges] = usedEdges4[countAdditionalEdges];
				    }
				  
				 /*  fprintf(stderr,"usedEdges1[%d]=%d\n",countAdditionalEdges, usedEdges1[countAdditionalEdges]); */
/* 				  fprintf(stderr,"usedEdges2[%d]=%d\n",countAdditionalEdges, usedEdges2[countAdditionalEdges]); */
/* 				  fprintf(stderr,"usedEdges3[%d]=%d\n",countAdditionalEdges, usedEdges3[countAdditionalEdges]); */
/* 				  fprintf(stderr,"usedEdges4[%d]=%d\n",countAdditionalEdges, usedEdges4[countAdditionalEdges]); */
/* 				  fprintf(stderr,"edge1[%d]=%d\n",i, edge1[i]); */
/* 				  fprintf(stderr,"edge2[%d]=%d\n",i, edge2[i]); */
/* 				  fprintf(stderr,"edge1[%d]=%d\n",j, edge1[j]); */
/* 				  fprintf(stderr,"edge2[%d]=%d\n",j, edge2[j]); */
				  
				  if (countAdditionalEdges>0)
				    {
				      for (l=0; l<countAdditionalEdges; l++)
					{
					  skip = 0;
					  if ( usedEdges1[countAdditionalEdges] == usedEdges1[l]
					       && usedEdges2[countAdditionalEdges] == usedEdges2[l]
					       && usedEdges3[countAdditionalEdges] == usedEdges3[l]
					       )
					    {
					      //fprintf(stderr,"usedEdges1[%d]=%d\n",countAdditionalEdges, usedEdges1[countAdditionalEdges]);
					      //fprintf(stderr,"usedEdges2[%d]=%d\n",countAdditionalEdges, usedEdges2[countAdditionalEdges]);
					      //fprintf(stderr,"usedEdges3[%d]=%d\n",countAdditionalEdges, usedEdges3[countAdditionalEdges]);
					      // fprintf(stderr,"usedEdges1[%d]=%d\n",l, usedEdges1[l]);
					      //fprintf(stderr,"usedEdges2[%d]=%d\n",l, usedEdges2[l]);
					      //fprintf(stderr,"usedEdges3[%d]=%d\n",l, usedEdges3[l]);
					      skip=1;
					      //      fprintf(stderr,"skipping surface\n");
					      break;
					    }
					}
				      if (skip==1) continue;
				    }
				  
				  skip=0;
				  countAdditionalEdges+=1;
				  //fprintf(stderr,"BEFORE adding\n");
				      if (edge1[i]!=edge2[i]) 
					{
					  //fprintf(stderr,"1 adding 1 for %d and %d\n",edge1[i], edge2[i]);
					  if (edge1[i]>edge2[i])
					    {
					      if (numberConnectionIsUsed[edge1[i]-1][edge2[i]-1]==2) continue;
					      if (numberConnectionIsUsed[edge1[i]-1][edge2[i]-1]==-1)
						numberConnectionIsUsed[edge1[i]-1][edge2[i]-1]+=2;
					      else
						numberConnectionIsUsed[edge1[i]-1][edge2[i]-1]+=1;
					    }
					  else
					    {
					      if (numberConnectionIsUsed[edge2[i]-1][edge1[i]-1]==2) continue;
					      if (numberConnectionIsUsed[edge2[i]-1][edge1[i]-1]==-1)
						numberConnectionIsUsed[edge2[i]-1][edge1[i]-1]+=2;
					      else
						numberConnectionIsUsed[edge2[i]-1][edge1[i]-1]+=1;
					    }
					}
				      if (edge1[j]!=edge2[i])
					if (edge1[i]!=edge1[j]) 
					  {
					    //fprintf(stderr,"2 adding 1 for %d and %d\n",edge1[i], edge1[j]);
					    if(edge1[i]>edge1[j])
					      {
						if (numberConnectionIsUsed[edge1[i]-1][edge1[j]-1]==2) continue;
						if (numberConnectionIsUsed[edge1[i]-1][edge1[j]-1]==-1)
						  numberConnectionIsUsed[edge1[i]-1][edge1[j]-1]+=2;
						else
						  numberConnectionIsUsed[edge1[i]-1][edge1[j]-1]+=1;
					      }
					    else
					      {
						if (numberConnectionIsUsed[edge1[j]-1][edge1[i]-1]==2) continue;
						if (numberConnectionIsUsed[edge1[j]-1][edge1[i]-1]==-1)
						  numberConnectionIsUsed[edge1[j]-1][edge1[i]-1]+=2;
						else
						  numberConnectionIsUsed[edge1[j]-1][edge1[i]-1]+=1;
					      }  
					  }
				      
				      if (edge1[j]!=edge2[j])
					{			
					  //fprintf(stderr,"3 adding 1 for %d and %d\n",edge1[j], edge2[j]);
					  if(edge1[j]>edge2[j])
					    {
					      if (numberConnectionIsUsed[edge1[j]-1][edge2[j]-1]==2) continue;
					      if (numberConnectionIsUsed[edge1[j]-1][edge2[j]-1]==-1)
						numberConnectionIsUsed[edge1[j]-1][edge2[j]-1]+=2;
					      else
						numberConnectionIsUsed[edge1[j]-1][edge2[j]-1]+=1;
					    }
					  else
					    {
					      if (numberConnectionIsUsed[edge2[j]-1][edge1[j]-1]==2) continue;
					      if (numberConnectionIsUsed[edge2[j]-1][edge1[j]-1]==-1)
						numberConnectionIsUsed[edge2[j]-1][edge1[j]-1]+=2;
					      else
						numberConnectionIsUsed[edge2[j]-1][edge1[j]-1]+=1;
					    }
					}
				      
				      if (edge1[j]!=edge2[i])
					if (edge2[i]!=edge2[j]) 
					  {
					    //fprintf(stderr,"4 adding 1 for %d and %d\n",edge2[i], edge2[j]);
					    if(edge2[i]>edge2[j])
					      {
						if (numberConnectionIsUsed[edge2[i]-1][edge2[j]-1]==2) continue;
						if (numberConnectionIsUsed[edge2[i]-1][edge2[j]-1]==-1)
						  numberConnectionIsUsed[edge2[i]-1][edge2[j]-1]+=2;
						else
						  numberConnectionIsUsed[edge2[i]-1][edge2[j]-1]+=1;
					      }
					    else
					      {
						if (numberConnectionIsUsed[edge2[j]-1][edge2[i]-1]==2) continue;
						if (numberConnectionIsUsed[edge2[j]-1][edge2[i]-1]==-1)
						  numberConnectionIsUsed[edge2[j]-1][edge2[i]-1]+=2;
						else
						  numberConnectionIsUsed[edge2[j]-1][edge2[i]-1]+=1;
					      }
					  }
				      
				      if (edge2[j]!=edge1[j] && edge2[j]!=edge2[i] && edge1[i]!=edge1[j])
					if (edge1[i]!=edge2[j]) 
					  {
					    //fprintf(stderr,"5 adding 1 for %d and %d\n",edge1[i], edge2[j]);
					    if(edge1[i]>edge1[j])
					      {
						if (numberConnectionIsUsed[edge1[i]-1][edge2[j]-1]==2) continue;
						if (numberConnectionIsUsed[edge1[i]-1][edge2[j]-1]==-1)
						  numberConnectionIsUsed[edge1[i]-1][edge2[j]-1]+=2;
						else
						  numberConnectionIsUsed[edge1[i]-1][edge2[j]-1]+=1;
					      }
					    else
					      {
						if (numberConnectionIsUsed[edge2[j]-1][edge1[i]-1]==2) continue;
						if (numberConnectionIsUsed[edge2[j]-1][edge1[i]-1]==-1)
						  numberConnectionIsUsed[edge2[j]-1][edge1[i]-1]+=2;
						else
						  numberConnectionIsUsed[edge2[j]-1][edge1[i]-1]+=1;
					      }
					  }
				      
				      if (edge1[j]!=edge2[j] && edge2[i]!=edge1[i] && edge2[i]!=edge2[j] && edge1[j]!=edge1[i])
					if (edge1[j]!=edge2[i])
					  {
					    //fprintf(stderr,"6 adding 1 for %d and %d\n",edge1[j], edge2[i]);
					    if(edge1[j]>edge2[i])
					      {
						if (numberConnectionIsUsed[edge1[j]-1][edge2[i]-1]==2) continue;
						if (numberConnectionIsUsed[edge1[j]-1][edge2[i]-1]==-1)
						  numberConnectionIsUsed[edge1[j]-1][edge2[i]-1]+=2;
						else
						  numberConnectionIsUsed[edge1[j]-1][edge2[i]-1]+=1;
					      }
					    else
					      {				      
						if (numberConnectionIsUsed[edge2[i]-1][edge1[j]-1]==2) continue;
						if (numberConnectionIsUsed[edge2[i]-1][edge1[j]-1]==-1)
						  numberConnectionIsUsed[edge2[i]-1][edge1[j]-1]+=2;
						else
						  numberConnectionIsUsed[edge2[i]-1][edge1[j]-1]+=1;
					      }
					  }
				      
				      COUNTER+=1;
				      //addCOUNTER+=1;
				      //fprintf(stderr,"counter=%d\n",COUNTER);
				      //fprintf(stderr,"edge1[%d]=%d\n",i,edge1[i]);
				      //fprintf(stderr,"edge2[%d]=%d\n",i,edge2[i]);
				      //fprintf(stderr,"edge1[%d]=%d\n",j,edge1[j]);
				      //fprintf(stderr,"edge2[%d]=%d\n",j,edge2[j]);
				      SU[COUNTER-1][0] = AD[i][1]*AD[j][2]*DD[i][3] - AD[i][1]*AD[j][3]*DD[i][2] + AD[i][3]*AD[j][1]*DD[i][2]
					-  AD[i][3]*AD[j][2]*DD[i][1] + AD[i][2]*AD[j][3]*DD[i][1] - AD[i][2]*AD[j][1]*DD[i][3];
				      SU[COUNTER-1][1] = AD[i][2]*AD[j][0]*DD[i][3] - AD[i][0]*AD[j][2]*DD[i][3] + AD[i][0]*AD[j][3]*DD[i][2]
					-  AD[i][3]*AD[j][0]*DD[i][2] + AD[i][3]*AD[j][2]*DD[i][0] - AD[i][2]*AD[j][3]*DD[i][0];
				      SU[COUNTER-1][2] = AD[i][0]*AD[j][1]*DD[i][3] - AD[i][0]*AD[j][3]*DD[i][1] + AD[i][3]*AD[j][0]*DD[i][1]
					-  AD[i][3]*AD[j][1]*DD[i][0] + AD[i][1]*AD[j][3]*DD[i][0] - AD[i][1]*AD[j][0]*DD[i][3];
				      SU[COUNTER-1][3] = AD[i][2]*AD[j][1]*DD[i][0] - AD[i][1]*AD[j][2]*DD[i][0] + AD[i][0]*AD[j][2]*DD[i][1]
					-  AD[i][2]*AD[j][0]*DD[i][1] + AD[i][1]*AD[j][0]*DD[i][2] - AD[i][0]*AD[j][1]*DD[i][2];
				      
				      for (k=0; k<4; k++)
					{
					  SU[COUNTER-1][k]/=6.; // to get tetraedron's volume
					}
				}
			    }
			}
		      
/* 		      countSingleEdges=0; */
/* 		      for(j=1; j<=32; j++) */
/* 			for (k=1; k<j; k++) */
/* 			  { */
/* 			    if ( numberConnectionIsUsed[j-1][k-1]>-1 ) */
/* 			      { */
/* 				//fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]); */
/* 				if (numberConnectionIsUsed[j-1][k-1]<2) */
/* 				  { */
/* 				    if (numberConnectionIsUsed[j-1][k-1]==0) */
/* 				      countSingleEdges+=2; */
/* 				    else */
/* 				      countSingleEdges+=1; */
/* 				    singleConnections[countSingleEdges-1][0]=j; */
/* 				    singleConnections[countSingleEdges-1][1]=k; */
/* 				    //fprintf(stderr,"-------------------------------------------\n"); */
/* 				  } */
/* 				if (numberConnectionIsUsed[j-1][k-1]>2) */
/* 				  { */
/* 				    fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]); */
/* 				    fprintf(stderr,"Three connections is messed up. That can't happen.\n"); */
/* 				    exit(1); */
/* 				  } */
/* 			      } */
/* 			  } */


		      for(j=0; j<32; j++)
			{
			  singleConnections[j][0]=0;
			  singleConnections[j][1]=0;
			  singleConnectionsUsed[j]=0;
			}
			  
			  
			  countSingleEdges=0;
			  for(j=1; j<=32; j++)
			    for (k=1; k<j; k++)
			      {
				if ( numberConnectionIsUsed[j-1][k-1]>-1 )
				  {
				    //fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]);
				    if (numberConnectionIsUsed[j-1][k-1]<2)
				      {
					if (numberConnectionIsUsed[j-1][k-1]==0)
					  countSingleEdges+=2;
					else
					  countSingleEdges+=1;
					singleConnections[countSingleEdges-1][0]=j;
					singleConnections[countSingleEdges-1][1]=k;
					//fprintf(stderr,"-------------------------------------------\n");
				      }
				    if (numberConnectionIsUsed[j-1][k-1]>2)
				      {
					fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]);
					fprintf(stderr,"Three connections is messed up. That can't happen.\n");
					break;
				      }
				  }
			      }
			  
			  
			  /* 		      countSingleEdges=0; */
			  /* 		      for(j=1; j<=32; j++) */
			  /* 			for (k=1; k<j; k++) */
			  /* 			  { */
			  /* 			    if ( numberConnectionIsUsed[j-1][k-1]>0 ) */
			  /* 			      { */
			  /* 				if (numberConnectionIsUsed[j-1][k-1]<2) */
			  /* 				  { */
			  /* 				    countSingleEdges+=1; */
			  /* 				    singleConnections[countSingleEdges-1][0]=j; */
			  /* 				    singleConnections[countSingleEdges-1][1]=k; */
			  /* 				  } */
			  /* 			      } */
			  /* 			  } */
			  
			  shift+=1;
			  
			  
			  if (countSingleEdges>0)
			    {
			      if (1==1)//countSingleEdges%3==0) //if we have a number of single edges that is a multiple of 3
				{
				  IS=0;
				  for (j=0; j<countSingleEdges; j++)
				    {
				      if(NSurfaces==12)
					{
					  // fprintf(stderr,"singleConnections[%d][0]=%d\n",j,singleConnections[j][0]);
					  //fprintf(stderr,"singleConnections[%d][1]=%d\n",j,singleConnections[j][1]);
					}
				      for (i=0; i<countSingleEdges; i++)
					{
					  if (i!=j && singleConnectionsUsed[j]==0 && singleConnectionsUsed[i]==0 )
					    {
					      for ( k=0; k<countSingleEdges; k++ )
						{
						  if(k!=j && k!=i)
						    {
						      if ( ((singleConnections[j][0] == singleConnections[i][0]
							     && singleConnections[j][1] == singleConnections[k][0]
							     && singleConnections[k][1] == singleConnections[i][1])
							    ||
							    (singleConnections[j][0] == singleConnections[i][1]
							     && singleConnections[j][1] == singleConnections[k][0]
							     && singleConnections[k][1] == singleConnections[i][0])
							    ||
							    (singleConnections[j][0] == singleConnections[i][0]
							     && singleConnections[j][1] == singleConnections[k][1]
							     && singleConnections[k][0] == singleConnections[i][1])
							    ||
							    (singleConnections[j][0] == singleConnections[i][1]
							     && singleConnections[j][1] == singleConnections[k][1]
							     && singleConnections[k][0] == singleConnections[i][0])
							    ||
							    (singleConnections[j][1] == singleConnections[i][0]
							     && singleConnections[j][0] == singleConnections[k][0]
							     && singleConnections[k][1] == singleConnections[i][1])
							    ||
							    (singleConnections[j][1] == singleConnections[i][1]
							     && singleConnections[j][0] == singleConnections[k][0]
							     && singleConnections[k][1] == singleConnections[i][0])
							    ||
							    (singleConnections[j][1] == singleConnections[i][0]
							     && singleConnections[j][0] == singleConnections[k][1]
							     && singleConnections[k][0] == singleConnections[i][1])
							    ||
							    (singleConnections[j][1] == singleConnections[i][1]
							     && singleConnections[j][0] == singleConnections[k][1]
							     && singleConnections[k][0] == singleConnections[i][0])
							    )&& singleConnectionsUsed[k]==0 )
							{
						      singleConnectionsUsed[k]+=1;
						      singleConnectionsUsed[i]+=1;
						      singleConnectionsUsed[j]+=1;
						      group[IS][0]=singleConnections[j][0];
						      group[IS][1]=singleConnections[j][1];
						      if (singleConnections[i][0]!=singleConnections[j][0] 
							  && singleConnections[i][0]!=singleConnections[j][1] )
							group[IS][2]=singleConnections[i][0];
						      else
							group[IS][2]=singleConnections[i][1];
						      IS+=1;
							}
						    }
						}
					    }
					}
				    }
				  
				  if (IS!=countSingleEdges/3) 
				    {
				      srand(time(NULL));
				      m=0;
				      tries=0;
				      while(IS!=countSingleEdges/3 && tries<5) 
					{
					  tries+=1;
					  for(j=0; j<32; j++)
					    {
					      singleConnectionsUsed[j]=0;
					    }
					  //fprintf(stderr,"MISSED SOME TETRAEDRA. got %d out of %d\n",IS, countSingleEdges/3);
					  IS=0;
					  for (j=0; j<countSingleEdges; j++)
					    {
					      m=abs(rand()%countSingleEdges);
					      //fprintf(stderr,"singleConnections[%d][0]=%d\n",j,singleConnections[j][0]);
					      //fprintf(stderr,"singleConnections[%d][1]=%d\n",j,singleConnections[j][1]);
					      for (l=0; l<countSingleEdges; l++)
						{
						  i=(l+m)%countSingleEdges;
						  //RANDOM ORDERING		
						  if (i!=j && singleConnectionsUsed[j]==0 && singleConnectionsUsed[i]==0 )
						    {
						      for ( k=0; k<countSingleEdges; k++ )
							{
							  if(k!=j && k!=i)
							    {
							      if ( ((singleConnections[j][0] == singleConnections[i][0]
								     && singleConnections[j][1] == singleConnections[k][0]
								     && singleConnections[k][1] == singleConnections[i][1])
								    ||
								    (singleConnections[j][0] == singleConnections[i][1]
								     && singleConnections[j][1] == singleConnections[k][0]
								     && singleConnections[k][1] == singleConnections[i][0])
								    ||
								    (singleConnections[j][0] == singleConnections[i][0]
								     && singleConnections[j][1] == singleConnections[k][1]
								     && singleConnections[k][0] == singleConnections[i][1])
								    ||
								    (singleConnections[j][0] == singleConnections[i][1]
								     && singleConnections[j][1] == singleConnections[k][1]
								     && singleConnections[k][0] == singleConnections[i][0])
								    ||
								    (singleConnections[j][1] == singleConnections[i][0]
								     && singleConnections[j][0] == singleConnections[k][0]
								     && singleConnections[k][1] == singleConnections[i][1])
								    ||
								    (singleConnections[j][1] == singleConnections[i][1]
								     && singleConnections[j][0] == singleConnections[k][0]
								     && singleConnections[k][1] == singleConnections[i][0])
								    ||
								    (singleConnections[j][1] == singleConnections[i][0]
								     && singleConnections[j][0] == singleConnections[k][1]
								     && singleConnections[k][0] == singleConnections[i][1])
								    ||
								    (singleConnections[j][1] == singleConnections[i][1]
								     && singleConnections[j][0] == singleConnections[k][1]
								     && singleConnections[k][0] == singleConnections[i][0])
								    )&& singleConnectionsUsed[k]==0 )
								{
								  singleConnectionsUsed[k]+=1;
								  singleConnectionsUsed[i]+=1;
								  singleConnectionsUsed[j]+=1;
								  group[IS][0]=singleConnections[j][0];
								  group[IS][1]=singleConnections[j][1];
								  if (singleConnections[i][0]!=singleConnections[j][0] 
								      && singleConnections[i][0]!=singleConnections[j][1] )
								    group[IS][2]=singleConnections[i][0];
								  else
								    group[IS][2]=singleConnections[i][1];
								  IS+=1;
								  //fprintf(stderr,"IS=%d\n",IS);
								}
							    }
							}
						    }
						}
					    }
					/*   for (j=0; j<IS; j++) */
/* 					    { */
/* 					      fprintf(stderr,"group[%d]=(%d,%d,%d)\n",j,group[j][0],group[j][1],group[j][2]); */
/* 					    } */
					}
				      //sleep(50);
				    }
				}
			      
			      if (countSingleEdges>0)
				{
				  if (1==1)//countSingleEdges%3==0) //if we have a number of single edges that is a multiple of 3
				    {
				      for (j=0; j<IS; j++)
					{
					  skip=0;
					  for(is=0; is<NSurfaces; is++)
					    {
					      if (iEdge[is]==group[j][0])
						{
						  for(i=0; i<NSurfaces; i++)
						    {
						      if (iEdge[i]==group[j][1])
							{
							  if(iEdge[is]>iEdge[i])
							    numberConnectionIsUsed[iEdge[is]-1][iEdge[i]-1]+=1;
							  else
							    numberConnectionIsUsed[iEdge[i]-1][iEdge[is]-1]+=1;
							  if(iEdge[is]>group[j][2])
							    numberConnectionIsUsed[iEdge[is]-1][group[j][2]-1]+=1;
							  else
							    numberConnectionIsUsed[group[j][2]-1][iEdge[is]-1]+=1;
							  if(iEdge[i]>group[j][2])
							    numberConnectionIsUsed[iEdge[i]-1][group[j][2]-1]+=1;
							  else
							    numberConnectionIsUsed[group[j][2]-1][iEdge[i]-1]+=1;
							  
							  //fprintf(stderr,"IN THE LOOP 1 and j=%d\n",j);
							  //fprintf(stderr,"edges= %d, %d, %d\n",iEdge[is], iEdge[i], group[j][2]);
							  for(k=0;k<4;k++)
							    {
							      AD[0][k] = cuts[is][k]-cuts[i][k];
							      BD[k] = VMID[k]-cuts[is][k];
							    }
							  skip=1;
							  continue;
							}
						    }
						  for(i=0; i<NSurfaces; i++)
						    {
						      if (iEdge[i]==group[j][2])
							{
							  for(k=0;k<4;k++)
							    {
							      AD[1][k] = cuts[is][k]-cuts[i][k];
							    }
							  skip=1;
							  continue;
							}
						    }
						}
					      if ( skip==1 ) continue;
					    }
					  COUNTER+=1;
					  SU[COUNTER-1][0] = AD[0][1]*AD[1][2]*BD[3] - AD[0][1]*AD[1][3]*BD[2] + AD[0][3]*AD[1][1]*BD[2]
					    -  AD[0][3]*AD[1][2]*BD[1] + AD[0][2]*AD[1][3]*BD[1] - AD[0][2]*AD[1][1]*BD[3];
					  SU[COUNTER-1][1] = AD[0][2]*AD[1][0]*BD[3] - AD[0][0]*AD[1][2]*BD[3] + AD[0][0]*AD[1][3]*BD[2]
					    -  AD[0][3]*AD[1][0]*BD[2] + AD[0][3]*AD[1][2]*BD[0] - AD[0][2]*AD[1][3]*BD[0];
					  SU[COUNTER-1][2] = AD[0][0]*AD[1][1]*BD[3] - AD[0][0]*AD[1][3]*BD[1] + AD[0][3]*AD[1][0]*BD[1]
					    -  AD[0][3]*AD[1][1]*BD[0] + AD[0][1]*AD[1][3]*BD[0] - AD[0][1]*AD[1][0]*BD[3];
					  SU[COUNTER-1][3] = AD[0][2]*AD[1][1]*BD[0] - AD[0][1]*AD[1][2]*BD[0] + AD[0][0]*AD[1][2]*BD[1]
					    -  AD[0][2]*AD[1][0]*BD[1] + AD[0][1]*AD[1][0]*BD[2] - AD[0][0]*AD[1][1]*BD[2];
					  
					  for (k=0; k<4;k++)
					    {
					      SU[COUNTER-1][k]/=6.; // to get tetraedron's volume
					    }
					  
					}
				      
/* 				      for (j=0; j<IS; j++) */
/* 					{ */
/* 					  fprintf(stderr,"group[%d]=(%d,%d,%d)\n",j,group[j][0],group[j][1],group[j][2]); */
/* 					} */
				      // if it is not a multiple of 3 try something else...
				    }
				}
			    }

			  //check if there are still issues
		  
			  countSingleEdges=0;
			  for(j=1; j<=32; j++)
			    for (k=1; k<j; k++)
			      {
				if ( numberConnectionIsUsed[j-1][k-1]>-1 )
				  {
				    //fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]);
				    if (numberConnectionIsUsed[j-1][k-1]<2)
				      {
					if (numberConnectionIsUsed[j-1][k-1]==0)
					  countSingleEdges+=2;
					else
					  countSingleEdges+=1;
					singleConnections[countSingleEdges-1][0]=j;
					singleConnections[countSingleEdges-1][1]=k;
					//fprintf(stderr,"-------------------------------------------\n");
				      }
				    if (numberConnectionIsUsed[j-1][k-1]>2)
				      {
					fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]);
					fprintf(stderr,"Three connections is messed up. That can't happen.\n");
					break;
				      }
				  }
			      }
		




			  if (countSingleEdges>0)
			    {
			      IS=0;
			      for (i=0; i<NSurfaces; i++)
				for (j=0; j<NSurfaces; j++)
				  {
				    if (numberConnectionIsUsed[iEdge[i]-1][iEdge[j]-1]==1)
				      {
					IS+=1;
					//fprintf(stderr,"Missing ones: %d-%d \n",iEdge[i],iEdge[j]);
					edge1[IS-1]=iEdge[i];
					edge2[IS-1]=iEdge[j];
				      }
				  }
			      for (is=0; is<IS; is++)
				{
				  //				  fprintf(stderr,"edge1[%d]=%d edge2[%d]=%d \n",is,edge1[is],is,edge2[is]);
				  for (m=0; m<IS; m++)
				    {
				      if (edge1[is]==edge1[m])
					{
					  if(IBIT[edge2[is]-1][edge2[m]-1]==1 || IBIT[edge2[is]-1][edge2[m]-1]==2 
					     || IBIT[edge2[is]-1][edge2[m]-1]==3 || IBIT[edge2[is]-1][edge2[m]-1]==4)
					    {
					      //      fprintf(stderr,"add %d - %d \n",edge2[is],edge2[m]);
					      if (numberConnectionIsUsed[edge2[is]-1][edge2[m]-1]==-1)
						numberConnectionIsUsed[edge2[is]-1][edge2[m]-1]=0;
					    }
					}
				      if (edge1[is]==edge2[m])
					{
					  if(IBIT[edge2[is]-1][edge1[m]-1]==1 || IBIT[edge2[is]-1][edge1[m]-1]==2
					     || IBIT[edge2[is]-1][edge1[m]-1]==3 || IBIT[edge2[is]-1][edge1[m]-1]==4)
					    {
					      //fprintf(stderr,"add %d - %d \n",edge2[is],edge1[m]);
					      if (numberConnectionIsUsed[edge2[is]-1][edge1[m]-1]==-1)
						numberConnectionIsUsed[edge2[is]-1][edge1[m]-1]=0;
					    }
					}
				      if (edge2[is]==edge2[m])
					{
					  if(IBIT[edge1[is]-1][edge1[m]-1]==1 || IBIT[edge1[is]-1][edge1[m]-1]==2
					     || IBIT[edge1[is]-1][edge1[m]-1]==3 || IBIT[edge1[is]-1][edge1[m]-1]==4)
					    {
					      //fprintf(stderr,"add %d - %d \n",edge1[is],edge1[m]);
					      if (numberConnectionIsUsed[edge1[is]-1][edge1[m]-1]==-1)
						numberConnectionIsUsed[edge1[is]-1][edge1[m]-1]=0;
					    }
					}
				      
				    }
				}
			      
			      
			      
			      countSingleEdges=0;
			      for(j=1; j<=32; j++)
				for (k=1; k<j; k++)
				  {
				    if ( numberConnectionIsUsed[j-1][k-1]>-1 )
				      {
					if (numberConnectionIsUsed[j-1][k-1]<2)
					  {
					    if (numberConnectionIsUsed[j-1][k-1]==0)
					      {
						countSingleEdges+=2;
						singleConnections[countSingleEdges-1][0]=j;
						singleConnections[countSingleEdges-1][1]=k;
						singleConnectionsUsed[countSingleEdges-1]=-1;
					      }
					    else
					      {
						countSingleEdges+=1;
						singleConnections[countSingleEdges-1][0]=j;
						singleConnections[countSingleEdges-1][1]=k;
						singleConnectionsUsed[countSingleEdges-1]=0;
					      }
					    //  fprintf(stderr,"-------------------------------------------\n");
					  }
					if (numberConnectionIsUsed[j-1][k-1]>2)
					  {
					    fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]);
					    fprintf(stderr,"Three connections is messed up. That can't happen.\n");
					    break;
					  }
				      }
				  }
			      
			      //fprintf(stderr,"DO THE VERY LAST CHECK, single edges=%d, tries2=%d\n",countSingleEdges,tries2);
			      
			      if (countSingleEdges>0)
				{
				  if (countSingleEdges%3==0 || NSurfaces>8) //if we have a number of single edges that is a multiple of 3
				    {
				      IS=0;
				      for (j=0; j<countSingleEdges; j++)
					{
					  if(NSurfaces==12)
					    {
					      //	      fprintf(stderr,"singleConnections[%d][0]=%d\n",j,singleConnections[j][0]);
					      //fprintf(stderr,"singleConnections[%d][1]=%d\n",j,singleConnections[j][1]);
					    }
					  for (i=0; i<countSingleEdges; i++)
					    {
					      if (i!=j && singleConnectionsUsed[j]<1 && singleConnectionsUsed[i]<1 )
						{
						  for ( k=0; k<countSingleEdges; k++ )
						    {
						      if(k!=j && k!=i)
							{
							  if ( ((singleConnections[j][0] == singleConnections[i][0]
								 && singleConnections[j][1] == singleConnections[k][0]
								 && singleConnections[k][1] == singleConnections[i][1])
								||
								(singleConnections[j][0] == singleConnections[i][1]
								 && singleConnections[j][1] == singleConnections[k][0]
								 && singleConnections[k][1] == singleConnections[i][0])
								||
								(singleConnections[j][0] == singleConnections[i][0]
								 && singleConnections[j][1] == singleConnections[k][1]
								 && singleConnections[k][0] == singleConnections[i][1])
								||
								(singleConnections[j][0] == singleConnections[i][1]
								 && singleConnections[j][1] == singleConnections[k][1]
								 && singleConnections[k][0] == singleConnections[i][0])
								||
								(singleConnections[j][1] == singleConnections[i][0]
								 && singleConnections[j][0] == singleConnections[k][0]
								 && singleConnections[k][1] == singleConnections[i][1])
								||
								(singleConnections[j][1] == singleConnections[i][1]
								 && singleConnections[j][0] == singleConnections[k][0]
								 && singleConnections[k][1] == singleConnections[i][0])
								||
								(singleConnections[j][1] == singleConnections[i][0]
								 && singleConnections[j][0] == singleConnections[k][1]
								 && singleConnections[k][0] == singleConnections[i][1])
								||
								(singleConnections[j][1] == singleConnections[i][1]
								 && singleConnections[j][0] == singleConnections[k][1]
								 && singleConnections[k][0] == singleConnections[i][0])
								)&& singleConnectionsUsed[k]==0 )
							    {
							      singleConnectionsUsed[k]+=1;
							      singleConnectionsUsed[i]+=1;
							      singleConnectionsUsed[j]+=1;
							      group[IS][0]=singleConnections[j][0];
							      group[IS][1]=singleConnections[j][1];
							      if (singleConnections[i][0]!=singleConnections[j][0] 
								  && singleConnections[i][0]!=singleConnections[j][1] )
								group[IS][2]=singleConnections[i][0];
							      else
								group[IS][2]=singleConnections[i][1];
							      IS+=1;
							    }
							}
						    }
						}
					    }
					}
				      
				      if (IS!=countSingleEdges/3) 
					{
					  srand(time(NULL));
					  m=0;
					  tries=0;
					  while(IS!=countSingleEdges/3 && tries<5) 
					    {
					      tries+=1;
					      
					      
					      
					      for(j=1; j<=32; j++)
						for (k=1; k<j; k++)
						  {
						    if (numberConnectionIsUsed[j-1][k-1]==0)
						      {
							singleConnectionsUsed[countSingleEdges-1]=-1;
						      }
						    else
						      {
							singleConnectionsUsed[countSingleEdges-1]=0;
						      }
						  }
					      
					      
					      //fprintf(stderr,"MISSED SOME TETRAEDRA. got %d out of %d\n",IS, countSingleEdges/3);
					      IS=0;
					      for (j=0; j<countSingleEdges; j++)
						{
						  m=abs(rand()%countSingleEdges);
						  //  fprintf(stderr,"singleConnections[%d][0]=%d\n",j,singleConnections[j][0]);
						  //fprintf(stderr,"singleConnections[%d][1]=%d\n",j,singleConnections[j][1]);
						  for (l=0; l<countSingleEdges; l++)
						    {
						      i=(l+m)%countSingleEdges;
						      //RANDOM ORDERING		
						      if (i!=j && singleConnectionsUsed[j]<1 && singleConnectionsUsed[i]<1 )
							{
							  for ( k=0; k<countSingleEdges; k++ )
							    {
							      if(k!=j && k!=i)
								{
								  if ( ((singleConnections[j][0] == singleConnections[i][0]
									 && singleConnections[j][1] == singleConnections[k][0]
									 && singleConnections[k][1] == singleConnections[i][1])
									||
									(singleConnections[j][0] == singleConnections[i][1]
									 && singleConnections[j][1] == singleConnections[k][0]
									 && singleConnections[k][1] == singleConnections[i][0])
									||
									(singleConnections[j][0] == singleConnections[i][0]
									 && singleConnections[j][1] == singleConnections[k][1]
									 && singleConnections[k][0] == singleConnections[i][1])
									||
									(singleConnections[j][0] == singleConnections[i][1]
									 && singleConnections[j][1] == singleConnections[k][1]
									 && singleConnections[k][0] == singleConnections[i][0])
									||
									(singleConnections[j][1] == singleConnections[i][0]
									 && singleConnections[j][0] == singleConnections[k][0]
									 && singleConnections[k][1] == singleConnections[i][1])
									||
									(singleConnections[j][1] == singleConnections[i][1]
									 && singleConnections[j][0] == singleConnections[k][0]
									 && singleConnections[k][1] == singleConnections[i][0])
									||
									(singleConnections[j][1] == singleConnections[i][0]
									 && singleConnections[j][0] == singleConnections[k][1]
									 && singleConnections[k][0] == singleConnections[i][1])
									||
									(singleConnections[j][1] == singleConnections[i][1]
									 && singleConnections[j][0] == singleConnections[k][1]
									 && singleConnections[k][0] == singleConnections[i][0])
									)&& singleConnectionsUsed[k]==0 )
								    {
								      if (singleConnectionsUsed[k]<1)
									singleConnectionsUsed[k]+=1;
								      if (singleConnectionsUsed[i]<1)
									singleConnectionsUsed[i]+=1;
								      if (singleConnectionsUsed[j]<1)
									singleConnectionsUsed[j]+=1;
								      
								      group[IS][0]=singleConnections[j][0];
								      group[IS][1]=singleConnections[j][1];
								      if (singleConnections[i][0]!=singleConnections[j][0] 
									  && singleConnections[i][0]!=singleConnections[j][1] )
									group[IS][2]=singleConnections[i][0];
								      else
									group[IS][2]=singleConnections[i][1];
								      IS+=1;
								      //      fprintf(stderr,"IS=%d\n",IS);
								    }
								}
							    }
							}
						    }
						}
/* 					      for (j=0; j<IS; j++) */
/* 						{ */
/* 						  fprintf(stderr,"group[%d]=(%d,%d,%d)\n",j,group[j][0],group[j][1],group[j][2]); */
/* 						} */
					    }
					  //sleep(50);
					}
				    }
				}
			      
			      if (countSingleEdges>0)
				{
				  if (countSingleEdges%3==0 || NSurfaces>8) //if we have a number of single edges that is a multiple of 3
				    {
				      for (j=0; j<IS; j++)
					{
					  skip=0;
					  for(is=0; is<NSurfaces; is++)
					    {
					      if (iEdge[is]==group[j][0])
						{
						  for(i=0; i<NSurfaces; i++)
						    {
						      if (iEdge[i]==group[j][1])
							{
							  if(iEdge[is]>iEdge[i])
							    {
							      if (numberConnectionIsUsed[iEdge[is]-1][iEdge[i]-1]>1)
								{skip=1; continue;}
							      numberConnectionIsUsed[iEdge[is]-1][iEdge[i]-1]+=1;
							    }
							  else
							    {
							      if (numberConnectionIsUsed[iEdge[i]-1][iEdge[is]-1]>1)
								{skip=1; continue;}
							      numberConnectionIsUsed[iEdge[i]-1][iEdge[is]-1]+=1;
							    }
							  if(iEdge[is]>group[j][2])
							    {
							      if (numberConnectionIsUsed[iEdge[is]-1][group[j][2]-1]>1)
								{skip=1; continue;}
							      numberConnectionIsUsed[iEdge[is]-1][group[j][2]-1]+=1;
							    }						      
							  else
							    {
							      if (numberConnectionIsUsed[group[j][2]-1][iEdge[is]-1]>1)
								{skip=1; continue;}
							      numberConnectionIsUsed[group[j][2]-1][iEdge[is]-1]+=1;
							    }
							  if(iEdge[i]>group[j][2])
							    {
							      if (numberConnectionIsUsed[iEdge[i]-1][group[j][2]-1]>1)
								{skip=1; continue;}
							      numberConnectionIsUsed[iEdge[i]-1][group[j][2]-1]+=1;
							    }
							  else
							    {
							      if (numberConnectionIsUsed[group[j][2]-1][iEdge[i]-1]>1)
								{skip=1; continue;}
							      numberConnectionIsUsed[group[j][2]-1][iEdge[i]-1]+=1;
							    }
							  //fprintf(stderr,"IN THE LOOP 1 and j=%d\n",j);
							  //fprintf(stderr,"edges= %d, %d, %d\n",iEdge[is], iEdge[i], group[j][2]);
							  for(k=0;k<4;k++)
							    {
							      AD[0][k] = cuts[is][k]-cuts[i][k];
							      BD[k] = VMID[k]-cuts[is][k];
							    }
							  skip=1;
							  continue;
							}
						    }
						  for(i=0; i<NSurfaces; i++)
						    {
						      if (iEdge[i]==group[j][2])
							{
							  for(k=0;k<4;k++)
							    {
							      AD[1][k] = cuts[is][k]-cuts[i][k];
							    }
							  skip=1;
							  continue;
							}
						    }
						}
					      if ( skip==1 ) continue;
					    }
					  COUNTER+=1;
					  SU[COUNTER-1][0] = AD[0][1]*AD[1][2]*BD[3] - AD[0][1]*AD[1][3]*BD[2] + AD[0][3]*AD[1][1]*BD[2]
					    -  AD[0][3]*AD[1][2]*BD[1] + AD[0][2]*AD[1][3]*BD[1] - AD[0][2]*AD[1][1]*BD[3];
					  SU[COUNTER-1][1] = AD[0][2]*AD[1][0]*BD[3] - AD[0][0]*AD[1][2]*BD[3] + AD[0][0]*AD[1][3]*BD[2]
					    -  AD[0][3]*AD[1][0]*BD[2] + AD[0][3]*AD[1][2]*BD[0] - AD[0][2]*AD[1][3]*BD[0];
					  SU[COUNTER-1][2] = AD[0][0]*AD[1][1]*BD[3] - AD[0][0]*AD[1][3]*BD[1] + AD[0][3]*AD[1][0]*BD[1]
					    -  AD[0][3]*AD[1][1]*BD[0] + AD[0][1]*AD[1][3]*BD[0] - AD[0][1]*AD[1][0]*BD[3];
					  SU[COUNTER-1][3] = AD[0][2]*AD[1][1]*BD[0] - AD[0][1]*AD[1][2]*BD[0] + AD[0][0]*AD[1][2]*BD[1]
					    -  AD[0][2]*AD[1][0]*BD[1] + AD[0][1]*AD[1][0]*BD[2] - AD[0][0]*AD[1][1]*BD[2];
					  
					  for (k=0; k<4;k++)
					    {
					      SU[COUNTER-1][k]/=6.; // to get tetraedron's volume
					    }
					  
					}
				      
				/*       for (j=0; j<IS; j++) */
/* 					{ */
/* 					  fprintf(stderr,"group[%d]=(%d,%d,%d)\n",j,group[j][0],group[j][1],group[j][2]); */
/* 					} */
				      // if it is not a multiple of 3 try something else...
				    }
				}
			    }
		    }		  
		  
		  //--------------------------------------------------------------------------------------------------------------
		  // at this point all tetraedra should have been found and added
		  //--------------------------------------------------------------------------------------------------------------

		  
		  // make sure every vector points towards lower energy densities.
		  for (j=0; j<COUNTER; j++)
		    {
		      if ( SU[j][0]*VD0 < 0)
			SU[j][0]*=-1;
		      if ( SU[j][1]*VD1 < 0)
			SU[j][1]*=-1;
		      if ( SU[j][2]*VD2 < 0)
			SU[j][2]*=-1;
		      if ( SU[j][3]*VD3 < 0)
			SU[j][3]*=-1;
		      for (k=0; k<4; k++)
			{
			  FULLSU[k]+=SU[j][k];
			}
		    }
		  
		  //fprintf(stderr,"number of tetraedra = %d\n",COUNTER);
		  
		  countSingleEdges=0;
		  for(j=1; j<=32; j++)
		    for (k=1; k<j; k++)
		      {
			if ( numberConnectionIsUsed[j-1][k-1]>-1 )
			  {
			    //fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]);
			    if (numberConnectionIsUsed[j-1][k-1]<2)
			      {
				countSingleEdges+=1;
				//fprintf(stderr,"-------------------------------------------\n");
			      }
			    if (numberConnectionIsUsed[j-1][k-1]>2)
			      {
				fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]);
				fprintf(stderr,"Three connections is messed up. That can't happen.\n");
				break;
			      }
			  }
		      }
//		  fprintf(stderr,"number of single Edges = %d\n",countSingleEdges);
		  if (countSingleEdges%3!=0)
		    {
		      fprintf(stderr,"NUMBER OF SINGLE EDGES IS NOT A MULTIPLE OF 3, number=%d\n",countSingleEdges);
		      //fprintf(stderr,"NSurfaces=%d\n",NSurfaces);
		      //fprintf(stderr,"Tetraedra=%d\n",COUNTER);
		      //sleep(1);
		      weirdCases+=1;
		     /*  countSingleEdges=0; */
/* 		      for(j=1; j<=32; j++) */
/* 			for (k=1; k<j; k++) */
/* 			  { */
/* 			    if ( numberConnectionIsUsed[j-1][k-1]>0 ) */
/* 			      { */
/* 				fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]); */
/* 				if (numberConnectionIsUsed[j-1][k-1]<2) */
/* 				  { */
/* 				    countSingleEdges+=1; */
/* 				    fprintf(stderr,"-------------------------------------------\n"); */
/* 				  } */
/* 				if (numberConnectionIsUsed[j-1][k-1]>2) */
/* 				  { */
/* 				    fprintf(stderr,"connection between %d and %d is used %d times.\n", j, k, numberConnectionIsUsed[j-1][k-1]); */
/* 				    fprintf(stderr,"Three connections is messed up. That can't happen.\n"); */
/* 				    exit(1); */
/* 				  } */
/* 			      } */
/* 			  } */
		      
		      //if (NSurfaces<12) exit(1);
		      fprintf(stderr,"percent error=%f\n", 100*(warnings*1.0)/(cells*1.0));
		    }
		  

		  cells+=1;
		  SUM += sqrt(FULLSU[0]*FULLSU[0]+FULLSU[1]*FULLSU[1]+FULLSU[2]*FULLSU[2]+FULLSU[3]*FULLSU[3]);
		  

	 	  if ((NSurfaces==6 && COUNTER!=8) || (NSurfaces==8 && COUNTER!=12) 
		      || (NSurfaces==10 && (COUNTER!=16 && COUNTER!=15)) || (NSurfaces==12 && (COUNTER!=19 && COUNTER!=20 && COUNTER!=21)))
		    {
		      fprintf(stderr,"*************************************************** NSurfaces=%d but number of tetraedra=%d\n",NSurfaces,COUNTER);
		      fprintf(stderr,"tauf=%f, xf=%f, yf=%f, etaf=%f\n", tauf,xf,yf,etaf);
		      //if (abs(yf)<=0.1 && abs(etaf)<=0.1) 
		      //		      sleep(1);
		      // fprintf(t_file,"%f %f %f %f\n",xf,tauf,FULLSU[1],FULLSU[0]);
		      //if(NSurfaces<10) exit(1);
		      warnings+=1;
		      // sleep(1);
		      if (countSingleEdges%3!=0) 
			{
			  fprintf(stderr,"single edges=%d\n",countSingleEdges);
			  //sleep(1);
			  //break;
			}
		      fprintf(stderr,"percent error=%f\n", 100*(warnings*1.0)/(cells*1.0));
		      //continue; // don't add the flawed cells at all
		    }
	  
		  //fprintf(stderr,"Volume=%f\n", SUM);
		  //fprintf(stderr,"Cells done=%d\n", cells);
		  //fprintf(stderr,"Warnings=%d\n", warnings);
		  
		  xfrac = (xf-x)/(DX);
		  yfrac = (yf-y)/(DY);
		  etafrac = (etaf-eta)/(DETA);
		  taufrac = (tauf-(tau-DTAU))/(DTAU);

		  // 4d interpolation:
		  //utau:
		  utauX1 = ((1.-xfrac)*arena[ix][iy][ieta].u[0][0] + xfrac*arena[ix+fac][iy][ieta].u[0][0]);
		  utauX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u[0][0] + xfrac*arena[ix+fac][iy+fac][ieta].u[0][0]);
		  
		  if (ieta<neta-fac)
		    {
		      utauX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u[0][0] + xfrac*arena[ix+fac][iy][ieta+fac].u[0][0]);
		      utauX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u[0][0] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u[0][0]);
		    }
		  else
		    {
		      utauX3 = ((1.-xfrac)*Rneighbor_utau[ix][iy] + xfrac*Rneighbor_utau[ix+fac][iy]);
		      utauX4 = ((1.-xfrac)*Rneighbor_utau[ix][iy+fac] + xfrac*Rneighbor_utau[ix+fac][iy+fac]);
		      //cout << "utau=" << utauX3 << " " << utauX4 << endl;
		    }

		  utauY1 = ((1.-yfrac)*utauX1+yfrac*utauX2);
		  utauY2 = ((1.-yfrac)*utauX3+yfrac*utauX4);

		  utau1 = ((1.-etafrac)*utauY1+etafrac*utauY2);
		
		  //utauX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][0]);
		  //utauX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][0]);
		  utauX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[0] + xfrac*arena[ix+fac][iy][ieta].u_prev[0]);
		  utauX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[0] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[0]);
		  
		  if (ieta<neta-fac)
		    {
		      //utauX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][0] 
		      //	+ xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][0]);
		      //utauX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][0] 
		      //	+ xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][0]);
		      utauX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[0] 
				+ xfrac*arena[ix+fac][iy][ieta+fac].u_prev[0]);
		      utauX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[0] 
				+ xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[0]);
		    }
		  else
		    {
		      utauX3 = ((1.-xfrac)*Rneighbor_utau_prev[ix][iy] 
				+ xfrac*Rneighbor_utau_prev[ix+fac][iy]);
		      utauX4 = ((1.-xfrac)*Rneighbor_utau_prev[ix][iy+fac]
				+ xfrac*Rneighbor_utau_prev[ix+fac][iy+fac]);
		      //cout << "utau_prev=" << utauX3 << " " << utauX4 << endl;
		    }
		  
		  utauY1 = ((1.-yfrac)*utauX1+yfrac*utauX2);
		  utauY2 = ((1.-yfrac)*utauX3+yfrac*utauX4);

		  utau2 = ((1.-etafrac)*utauY1+etafrac*utauY2);
		  
		  utau = (1.-taufrac)*utau2+taufrac*utau1;

		  //ux:
		  uxX1 = ((1.-xfrac)*arena[ix][iy][ieta].u[0][1] + xfrac*arena[ix+fac][iy][ieta].u[0][1]);
		  uxX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u[0][1] + xfrac*arena[ix+fac][iy+fac][ieta].u[0][1]);
		
		  if (ieta<neta-fac)
		    {
		      uxX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u[0][1] + xfrac*arena[ix+fac][iy][ieta+fac].u[0][1]);
		      uxX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u[0][1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u[0][1]);
		    }
		  else
		    {
		      uxX3 = ((1.-xfrac)*Rneighbor_ux[ix][iy] + xfrac*Rneighbor_ux[ix+fac][iy]);
		      uxX4 = ((1.-xfrac)*Rneighbor_ux[ix][iy+fac] + xfrac*Rneighbor_ux[ix+fac][iy+fac]);
		      //cout << "ux=" << uxX3 << " " << uxX4 << endl;
		    }
		  
		  uxY1 = ((1.-yfrac)*uxX1+yfrac*uxX2);
		  uxY2 = ((1.-yfrac)*uxX3+yfrac*uxX4);

		  ux1 = ((1.-etafrac)*uxY1+etafrac*uxY2);
			  
		  //uxX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][1]);
		  //uxX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][1]);
		  uxX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[1] + xfrac*arena[ix+fac][iy][ieta].u_prev[1]);
		  uxX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[1] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[1]);
		  

		  if (ieta<neta-fac)
		    {
		      //uxX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][1]);
		      //uxX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][1]);
		      uxX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[1] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[1]);
		      uxX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[1]);
		    }
		  else
		    {
		      uxX3 = ((1.-xfrac)*Rneighbor_ux_prev[ix][iy] 
			      + xfrac*Rneighbor_ux_prev[ix+fac][iy]);
		      uxX4 = ((1.-xfrac)*Rneighbor_ux_prev[ix][iy+fac] 
			      + xfrac*Rneighbor_ux_prev[ix+fac][iy+fac]);
		      //cout << "ux_prev=" << uxX3 << " " << uxX4 << endl;
		    }

		  uxY1 = ((1.-yfrac)*uxX1+yfrac*uxX2);
		  uxY2 = ((1.-yfrac)*uxX3+yfrac*uxX4);

		  ux2 = ((1.-etafrac)*uxY1+etafrac*uxY2);

		  ux = (1.-taufrac)*ux2+taufrac*ux1;
		
		  //uy:
		  uyX1 = ((1.-xfrac)*arena[ix][iy][ieta].u[0][2] + xfrac*arena[ix+fac][iy][ieta].u[0][2]);
		  uyX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u[0][2] + xfrac*arena[ix+fac][iy+fac][ieta].u[0][2]);

		  if (ieta<neta-fac)
		    {
		      uyX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u[0][2] + xfrac*arena[ix+fac][iy][ieta+fac].u[0][2]);
		      uyX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u[0][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u[0][2]);
		    }
		  else
		    {
		      uyX3 = ((1.-xfrac)*Rneighbor_uy[ix][iy] + xfrac*Rneighbor_uy[ix+fac][iy]);
		      uyX4 = ((1.-xfrac)*Rneighbor_uy[ix][iy+fac] + xfrac*Rneighbor_uy[ix+fac][iy+fac]);
		      //cout << "uy=" << uyX3 << " " << uyX4 << endl;
		    }

		  uyY1 = ((1.-yfrac)*uyX1+yfrac*uyX2);
		  uyY2 = ((1.-yfrac)*uyX3+yfrac*uyX4);

		  uy1 = ((1.-etafrac)*uyY1+etafrac*uyY2);
			  
		  //uyX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][2]);
		  //uyX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][2]);
		  uyX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[2] + xfrac*arena[ix+fac][iy][ieta].u_prev[2]);
		  uyX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[2] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[2]);

		  if (ieta<neta-fac)
		    {
		      //uyX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][2]);
		      //uyX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][2]);
		      uyX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[2] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[2]);
		      uyX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[2]);
		    }
		  else
		    {
		      uyX3 = ((1.-xfrac)*Rneighbor_uy_prev[ix][iy] + xfrac*Rneighbor_uy_prev[ix+fac][iy]);
		      uyX4 = ((1.-xfrac)*Rneighbor_uy_prev[ix][iy+fac] + xfrac*Rneighbor_uy_prev[ix+fac][iy+fac]);
		      //cout << "uy_prev=" << uyX3 << " " << uyX4 << endl;
		    }

		  uyY1 = ((1.-yfrac)*uyX1+yfrac*uyX2);
		  uyY2 = ((1.-yfrac)*uyX3+yfrac*uyX4);

		  uy2 = ((1.-etafrac)*uyY1+etafrac*uyY2);

		  uy = (1.-taufrac)*uy2+taufrac*uy1;
		  

		  //ueta:
		  uetaX1 = ((1.-xfrac)*arena[ix][iy][ieta].u[0][3] + xfrac*arena[ix+fac][iy][ieta].u[0][3]);
		  uetaX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u[0][3] + xfrac*arena[ix+fac][iy+fac][ieta].u[0][3]);

		  if (ieta<neta-fac)
		    {
		      uetaX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u[0][3] + xfrac*arena[ix+fac][iy][ieta+fac].u[0][3]);
		      uetaX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u[0][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u[0][3]);
		    }
		  else
		    {
		      uetaX3 = ((1.-xfrac)*Rneighbor_ueta[ix][iy] + xfrac*Rneighbor_ueta[ix+fac][iy]);
		      uetaX4 = ((1.-xfrac)*Rneighbor_ueta[ix][iy+fac] + xfrac*Rneighbor_ueta[ix+fac][iy+fac]);
		      //cout << "ueta=" << uetaX3 << " " << uetaX4 << endl;
		    }
		  
		  uetaY1 = ((1.-yfrac)*uetaX1+yfrac*uetaX2);
		  uetaY2 = ((1.-yfrac)*uetaX3+yfrac*uetaX4);

		  ueta1 = ((1.-etafrac)*uetaY1+etafrac*uetaY2);
			  
		  //uetaX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][3]);
		  //uetaX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][3]);
		  uetaX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[3] + xfrac*arena[ix+fac][iy][ieta].u_prev[3]);
		  uetaX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[3] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[3]);

		  if (ieta<neta-fac)
		    {
		      //uetaX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][3]);
		      //uetaX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][3]);
		      uetaX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[3] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[3]);
		      uetaX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[3]);
		    }
		  else
		    {
		      uetaX3 = ((1.-xfrac)*Rneighbor_ueta_prev[ix][iy] 
				+ xfrac*Rneighbor_ueta_prev[ix+fac][iy]);
		      uetaX4 = ((1.-xfrac)*Rneighbor_ueta_prev[ix][iy+fac] 
				+ xfrac*Rneighbor_ueta_prev[ix+fac][iy+fac]);
		      //cout << "ueta_prev=" << uetaX3 << " " << uetaX4 << endl;
		    }
		  
		  uetaY1 = ((1.-yfrac)*uetaX1+yfrac*uetaX2);
		  uetaY2 = ((1.-yfrac)*uetaX3+yfrac*uetaX4);

		  ueta2 = ((1.-etafrac)*uetaY1+etafrac*uetaY2);

		  ueta = (1.-taufrac)*ueta2+taufrac*ueta1;
		  
		  
		  //rhob: 3D interpolation for now (no tau interpolation)
		  rhobX1 = ((1.-xfrac)*arena[ix][iy][ieta].rhob + xfrac*arena[ix+fac][iy][ieta].rhob);
		  rhobX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].rhob + xfrac*arena[ix+fac][iy+fac][ieta].rhob);

		  if (ieta<neta-fac)
		    {
		      rhobX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].rhob + xfrac*arena[ix+fac][iy][ieta+fac].rhob);
		      rhobX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].rhob + xfrac*arena[ix+fac][iy+fac][ieta+fac].rhob);
		    }
		  else
		    {
		      rhobX3 = ((1.-xfrac)*Rneighbor_rhob[ix][iy] + xfrac*Rneighbor_rhob[ix+fac][iy]);
		      rhobX4 = ((1.-xfrac)*Rneighbor_rhob[ix][iy+fac] + xfrac*Rneighbor_rhob[ix+fac][iy+fac]);
		    }
		  
		  rhobY1 = ((1.-yfrac)*rhobX1+yfrac*rhobX2);
		  rhobY2 = ((1.-yfrac)*rhobX3+yfrac*rhobX4);

		  rhob = ((1.-etafrac)*rhobY1+etafrac*rhobY2);
		
		  //Wtautau:
		  WX1 =     ((1.-xfrac)*arena[ix][iy][ieta].Wmunu[0][0][0] + xfrac*arena[ix+fac][iy][ieta].Wmunu[0][0][0]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].Wmunu[0][0][0] + xfrac*arena[ix+fac][iy+fac][ieta].Wmunu[0][0][0]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].deltaT[0][0][0] + arena[ix][iy][ieta].deltaW[0][0][0])
		             + xfrac*(arena[ix+fac][iy][ieta].deltaT[0][0][0] + arena[ix+fac][iy][ieta].deltaW[0][0][0]);
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].deltaT[0][0][0] + arena[ix][iy+fac][ieta].deltaW[0][0][0])
		              + xfrac*(arena[ix+fac][iy+fac][ieta].deltaT[0][0][0] + arena[ix+fac][iy+fac][ieta].deltaW[0][0][0]);
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].Wmunu[0][0][0] + xfrac*arena[ix+fac][iy][ieta+fac].Wmunu[0][0][0]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].Wmunu[0][0][0] + xfrac*arena[ix+fac][iy+fac][ieta+fac].Wmunu[0][0][0]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].deltaT[0][0][0] + arena[ix][iy][ieta+fac].deltaW[0][0][0])
			          + xfrac*(arena[ix+fac][iy][ieta+fac].deltaT[0][0][0] + arena[ix+fac][iy][ieta+fac].deltaW[0][0][0]);
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].deltaT[0][0][0] + arena[ix][iy+fac][ieta+fac].deltaW[0][0][0])
			          + xfrac*(arena[ix+fac][iy+fac][ieta+fac].deltaT[0][0][0] + arena[ix+fac][iy+fac][ieta+fac].deltaW[0][0][0]);
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wtautau[ix][iy] + xfrac*Rneighbor_Wtautau[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wtautau[ix][iy+fac] + xfrac*Rneighbor_Wtautau[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wtautau1 = ((1.-etafrac)*WY1+etafrac*WY2);
		  
		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].W_prev[0][0] + xfrac*arena[ix+fac][iy][ieta].W_prev[0][0]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].W_prev[0][0] + xfrac*arena[ix+fac][iy+fac][ieta].W_prev[0][0]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].prev_deltaW[0][0] + arena[ix][iy][ieta].getPrevDeltaT(0, 0, eos))
		      + xfrac*(arena[ix+fac][iy][ieta].prev_deltaW[0][0] + arena[ix+fac][iy][ieta].getPrevDeltaT(0, 0, eos));
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].prev_deltaW[0][0] + arena[ix][iy+fac][ieta].getPrevDeltaT(0, 0, eos))
		      + xfrac*(arena[ix+fac][iy+fac][ieta].prev_deltaW[0][0] + arena[ix+fac][iy+fac][ieta].getPrevDeltaT(0, 0, eos));
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].W_prev[0][0] + xfrac*arena[ix+fac][iy][ieta+fac].W_prev[0][0]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].W_prev[0][0] + xfrac*arena[ix+fac][iy+fac][ieta+fac].W_prev[0][0]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].prev_deltaW[0][0] + arena[ix][iy][ieta+fac].getPrevDeltaT(0, 0, eos))
			  + xfrac*(arena[ix+fac][iy][ieta+fac].prev_deltaW[0][0] + arena[ix+fac][iy][ieta+fac].getPrevDeltaT(0, 0, eos));
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].prev_deltaW[0][0] + arena[ix][iy+fac][ieta+fac].getPrevDeltaT(0, 0, eos))
			  + xfrac*(arena[ix+fac][iy+fac][ieta+fac].prev_deltaW[0][0] + arena[ix+fac][iy+fac][ieta+fac].getPrevDeltaT(0, 0, eos));
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wtautau_prev[ix][iy] + xfrac*Rneighbor_Wtautau_prev[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wtautau_prev[ix][iy+fac] + xfrac*Rneighbor_Wtautau_prev[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wtautau2 = ((1.-etafrac)*WY1+etafrac*WY2);

		  Wtautau = (1.-taufrac)*Wtautau2+taufrac*Wtautau1;

		  //Wtaux:

		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].Wmunu[0][0][1] + xfrac*arena[ix+fac][iy][ieta].Wmunu[0][0][1]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].Wmunu[0][0][1] + xfrac*arena[ix+fac][iy+fac][ieta].Wmunu[0][0][1]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].deltaT[0][0][1] + arena[ix][iy][ieta].deltaW[0][0][1])
		             + xfrac*(arena[ix+fac][iy][ieta].deltaT[0][0][1] + arena[ix+fac][iy][ieta].deltaW[0][0][1]);
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].deltaT[0][0][1] + arena[ix][iy+fac][ieta].deltaW[0][0][1])
		              + xfrac*(arena[ix+fac][iy+fac][ieta].deltaT[0][0][1] + arena[ix+fac][iy+fac][ieta].deltaW[0][0][1]);
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].Wmunu[0][0][1] + xfrac*arena[ix+fac][iy][ieta+fac].Wmunu[0][0][1]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].Wmunu[0][0][1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].Wmunu[0][0][1]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].deltaT[0][0][1] + arena[ix][iy][ieta+fac].deltaW[0][0][1])
			          + xfrac*(arena[ix+fac][iy][ieta+fac].deltaT[0][0][1] + arena[ix+fac][iy][ieta+fac].deltaW[0][0][1]);
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].deltaT[0][0][1] + arena[ix][iy+fac][ieta+fac].deltaW[0][0][1])
			          + xfrac*(arena[ix+fac][iy+fac][ieta+fac].deltaT[0][0][1] + arena[ix+fac][iy+fac][ieta+fac].deltaW[0][0][1]);
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wtaux[ix][iy] + xfrac*Rneighbor_Wtaux[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wtaux[ix][iy+fac] + xfrac*Rneighbor_Wtaux[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wtaux1 = ((1.-etafrac)*WY1+etafrac*WY2);
		
		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].W_prev[0][1] + xfrac*arena[ix+fac][iy][ieta].W_prev[0][1]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].W_prev[0][1] + xfrac*arena[ix+fac][iy+fac][ieta].W_prev[0][1]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].prev_deltaW[0][1] + arena[ix][iy][ieta].getPrevDeltaT(0, 1, eos))
		      + xfrac*(arena[ix+fac][iy][ieta].prev_deltaW[0][1] + arena[ix+fac][iy][ieta].getPrevDeltaT(0, 1, eos));
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].prev_deltaW[0][1] + arena[ix][iy+fac][ieta].getPrevDeltaT(0, 1, eos))
		      + xfrac*(arena[ix+fac][iy+fac][ieta].prev_deltaW[0][1] + arena[ix+fac][iy+fac][ieta].getPrevDeltaT(0, 1, eos));
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].W_prev[0][1] + xfrac*arena[ix+fac][iy][ieta+fac].W_prev[0][1]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].W_prev[0][1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].W_prev[0][1]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].prev_deltaW[0][1] + arena[ix][iy][ieta+fac].getPrevDeltaT(0, 1, eos))
			  + xfrac*(arena[ix+fac][iy][ieta+fac].prev_deltaW[0][1] + arena[ix+fac][iy][ieta+fac].getPrevDeltaT(0, 1, eos));
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].prev_deltaW[0][1] + arena[ix][iy+fac][ieta+fac].getPrevDeltaT(0, 1, eos))
			  + xfrac*(arena[ix+fac][iy+fac][ieta+fac].prev_deltaW[0][1] + arena[ix+fac][iy+fac][ieta+fac].getPrevDeltaT(0, 1, eos));
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wtaux_prev[ix][iy] + xfrac*Rneighbor_Wtaux_prev[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wtaux_prev[ix][iy+fac] + xfrac*Rneighbor_Wtaux_prev[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wtaux2 = ((1.-etafrac)*WY1+etafrac*WY2);
		
		  Wtaux = (1.-taufrac)*Wtaux2+taufrac*Wtaux1;

		  //Wtauy:
 
		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].Wmunu[0][0][2] + xfrac*arena[ix+fac][iy][ieta].Wmunu[0][0][2]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].Wmunu[0][0][2] + xfrac*arena[ix+fac][iy+fac][ieta].Wmunu[0][0][2]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].deltaT[0][0][2] + arena[ix][iy][ieta].deltaW[0][0][2])
		             + xfrac*(arena[ix+fac][iy][ieta].deltaT[0][0][2] + arena[ix+fac][iy][ieta].deltaW[0][0][2]);
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].deltaT[0][0][2] + arena[ix][iy+fac][ieta].deltaW[0][0][2])
		              + xfrac*(arena[ix+fac][iy+fac][ieta].deltaT[0][0][2] + arena[ix+fac][iy+fac][ieta].deltaW[0][0][2]);
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].Wmunu[0][0][2] + xfrac*arena[ix+fac][iy][ieta+fac].Wmunu[0][0][2]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].Wmunu[0][0][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].Wmunu[0][0][2]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].deltaT[0][0][2] + arena[ix][iy][ieta+fac].deltaW[0][0][2])
			          + xfrac*(arena[ix+fac][iy][ieta+fac].deltaT[0][0][2] + arena[ix+fac][iy][ieta+fac].deltaW[0][0][2]);
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].deltaT[0][0][2] + arena[ix][iy+fac][ieta+fac].deltaW[0][0][2])
			          + xfrac*(arena[ix+fac][iy+fac][ieta+fac].deltaT[0][0][2] + arena[ix+fac][iy+fac][ieta+fac].deltaW[0][0][2]);
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wtauy[ix][iy] + xfrac*Rneighbor_Wtauy[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wtauy[ix][iy+fac] + xfrac*Rneighbor_Wtauy[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wtauy1 = ((1.-etafrac)*WY1+etafrac*WY2);
		
		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].W_prev[0][2] + xfrac*arena[ix+fac][iy][ieta].W_prev[0][2]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].W_prev[0][2] + xfrac*arena[ix+fac][iy+fac][ieta].W_prev[0][2]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].prev_deltaW[0][2] + arena[ix][iy][ieta].getPrevDeltaT(0, 2, eos))
		      + xfrac*(arena[ix+fac][iy][ieta].prev_deltaW[0][2] + arena[ix+fac][iy][ieta].getPrevDeltaT(0, 2, eos));
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].prev_deltaW[0][2] + arena[ix][iy+fac][ieta].getPrevDeltaT(0, 2, eos))
		      + xfrac*(arena[ix+fac][iy+fac][ieta].prev_deltaW[0][2] + arena[ix+fac][iy+fac][ieta].getPrevDeltaT(0, 2, eos));
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].W_prev[0][2] + xfrac*arena[ix+fac][iy][ieta+fac].W_prev[0][2]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].W_prev[0][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].W_prev[0][2]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].prev_deltaW[0][2] + arena[ix][iy][ieta+fac].getPrevDeltaT(0, 2, eos))
			  + xfrac*(arena[ix+fac][iy][ieta+fac].prev_deltaW[0][2] + arena[ix+fac][iy][ieta+fac].getPrevDeltaT(0, 2, eos));
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].prev_deltaW[0][2] + arena[ix][iy+fac][ieta+fac].getPrevDeltaT(0, 2, eos))
			  + xfrac*(arena[ix+fac][iy+fac][ieta+fac].prev_deltaW[0][2] + arena[ix+fac][iy+fac][ieta+fac].getPrevDeltaT(0, 2, eos));
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wtauy_prev[ix][iy] + xfrac*Rneighbor_Wtauy_prev[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wtauy_prev[ix][iy+fac] + xfrac*Rneighbor_Wtauy_prev[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wtauy2 = ((1.-etafrac)*WY1+etafrac*WY2);

		  Wtauy = (1.-taufrac)*Wtauy2+taufrac*Wtauy1;
		
		  //Wtaueta

		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].Wmunu[0][0][3] + xfrac*arena[ix+fac][iy][ieta].Wmunu[0][0][3]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].Wmunu[0][0][3] + xfrac*arena[ix+fac][iy+fac][ieta].Wmunu[0][0][3]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].deltaT[0][0][3] + arena[ix][iy][ieta].deltaW[0][0][3])
		             + xfrac*(arena[ix+fac][iy][ieta].deltaT[0][0][3] + arena[ix+fac][iy][ieta].deltaW[0][0][3]);
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].deltaT[0][0][3] + arena[ix][iy+fac][ieta].deltaW[0][0][3])
		              + xfrac*(arena[ix+fac][iy+fac][ieta].deltaT[0][0][3] + arena[ix+fac][iy+fac][ieta].deltaW[0][0][3]);
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].Wmunu[0][0][3] + xfrac*arena[ix+fac][iy][ieta+fac].Wmunu[0][0][3]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].Wmunu[0][0][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].Wmunu[0][0][3]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].deltaT[0][0][3] + arena[ix][iy][ieta+fac].deltaW[0][0][3])
			          + xfrac*(arena[ix+fac][iy][ieta+fac].deltaT[0][0][3] + arena[ix+fac][iy][ieta+fac].deltaW[0][0][3]);
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].deltaT[0][0][3] + arena[ix][iy+fac][ieta+fac].deltaW[0][0][3])
			          + xfrac*(arena[ix+fac][iy+fac][ieta+fac].deltaT[0][0][3] + arena[ix+fac][iy+fac][ieta+fac].deltaW[0][0][3]);
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wtaueta[ix][iy] + xfrac*Rneighbor_Wtaueta[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wtaueta[ix][iy+fac] + xfrac*Rneighbor_Wtaueta[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wtaueta1 = ((1.-etafrac)*WY1+etafrac*WY2);

		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].W_prev[0][3] + xfrac*arena[ix+fac][iy][ieta].W_prev[0][3]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].W_prev[0][3] + xfrac*arena[ix+fac][iy+fac][ieta].W_prev[0][3]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].prev_deltaW[0][3] + arena[ix][iy][ieta].getPrevDeltaT(0, 3, eos))
		      + xfrac*(arena[ix+fac][iy][ieta].prev_deltaW[0][3] + arena[ix+fac][iy][ieta].getPrevDeltaT(0, 3, eos));
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].prev_deltaW[0][3] + arena[ix][iy+fac][ieta].getPrevDeltaT(0, 3, eos))
		      + xfrac*(arena[ix+fac][iy+fac][ieta].prev_deltaW[0][3] + arena[ix+fac][iy+fac][ieta].getPrevDeltaT(0, 3, eos));
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].W_prev[0][3] + xfrac*arena[ix+fac][iy][ieta+fac].W_prev[0][3]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].W_prev[0][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].W_prev[0][3]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].prev_deltaW[0][3] + arena[ix][iy][ieta+fac].getPrevDeltaT(0, 3, eos))
			  + xfrac*(arena[ix+fac][iy][ieta+fac].prev_deltaW[0][3] + arena[ix+fac][iy][ieta+fac].getPrevDeltaT(0, 3, eos));
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].prev_deltaW[0][3] + arena[ix][iy+fac][ieta+fac].getPrevDeltaT(0, 3, eos))
			  + xfrac*(arena[ix+fac][iy+fac][ieta+fac].prev_deltaW[0][3] + arena[ix+fac][iy+fac][ieta+fac].getPrevDeltaT(0, 3, eos));
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wtaueta_prev[ix][iy] + xfrac*Rneighbor_Wtaueta_prev[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wtaueta_prev[ix][iy+fac] + xfrac*Rneighbor_Wtaueta_prev[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wtaueta2 = ((1.-etafrac)*WY1+etafrac*WY2);

		  Wtaueta = (1.-taufrac)*Wtaueta2+taufrac*Wtaueta1;

		  //Wxx:

		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].Wmunu[0][1][1] + xfrac*arena[ix+fac][iy][ieta].Wmunu[0][1][1]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].Wmunu[0][1][1] + xfrac*arena[ix+fac][iy+fac][ieta].Wmunu[0][1][1]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].deltaT[0][1][1] + arena[ix][iy][ieta].deltaW[0][1][1])
		             + xfrac*(arena[ix+fac][iy][ieta].deltaT[0][1][1] + arena[ix+fac][iy][ieta].deltaW[0][1][1]);
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].deltaT[0][1][1] + arena[ix][iy+fac][ieta].deltaW[0][1][1])
		              + xfrac*(arena[ix+fac][iy+fac][ieta].deltaT[0][1][1] + arena[ix+fac][iy+fac][ieta].deltaW[0][1][1]);
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].Wmunu[0][1][1] + xfrac*arena[ix+fac][iy][ieta+fac].Wmunu[0][1][1]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].Wmunu[0][1][1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].Wmunu[0][1][1]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].deltaT[0][1][1] + arena[ix][iy][ieta+fac].deltaW[0][1][1])
			          + xfrac*(arena[ix+fac][iy][ieta+fac].deltaT[0][1][1] + arena[ix+fac][iy][ieta+fac].deltaW[0][1][1]);
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].deltaT[0][1][1] + arena[ix][iy+fac][ieta+fac].deltaW[0][1][1])
			          + xfrac*(arena[ix+fac][iy+fac][ieta+fac].deltaT[0][1][1] + arena[ix+fac][iy+fac][ieta+fac].deltaW[0][1][1]);
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wxx[ix][iy] + xfrac*Rneighbor_Wxx[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wxx[ix][iy+fac] + xfrac*Rneighbor_Wxx[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wxx1 = ((1.-etafrac)*WY1+etafrac*WY2);
		
		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].W_prev[1][1] + xfrac*arena[ix+fac][iy][ieta].W_prev[1][1]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].W_prev[1][1] + xfrac*arena[ix+fac][iy+fac][ieta].W_prev[1][1]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].prev_deltaW[1][1] + arena[ix][iy][ieta].getPrevDeltaT(1, 1, eos))
		      + xfrac*(arena[ix+fac][iy][ieta].prev_deltaW[1][1] + arena[ix+fac][iy][ieta].getPrevDeltaT(1, 1, eos));
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].prev_deltaW[1][1] + arena[ix][iy+fac][ieta].getPrevDeltaT(1, 1, eos))
		      + xfrac*(arena[ix+fac][iy+fac][ieta].prev_deltaW[1][1] + arena[ix+fac][iy+fac][ieta].getPrevDeltaT(1, 1, eos));
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].W_prev[1][1] + xfrac*arena[ix+fac][iy][ieta+fac].W_prev[1][1]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].W_prev[1][1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].W_prev[1][1]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].prev_deltaW[1][1] + arena[ix][iy][ieta+fac].getPrevDeltaT(1, 1, eos))
			  + xfrac*(arena[ix+fac][iy][ieta+fac].prev_deltaW[1][1] + arena[ix+fac][iy][ieta+fac].getPrevDeltaT(1, 1, eos));
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].prev_deltaW[1][1] + arena[ix][iy+fac][ieta+fac].getPrevDeltaT(1, 1, eos))
			  + xfrac*(arena[ix+fac][iy+fac][ieta+fac].prev_deltaW[1][1] + arena[ix+fac][iy+fac][ieta+fac].getPrevDeltaT(1, 1, eos));
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wxx_prev[ix][iy] + xfrac*Rneighbor_Wxx_prev[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wxx_prev[ix][iy+fac] + xfrac*Rneighbor_Wxx_prev[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wxx2 = ((1.-etafrac)*WY1+etafrac*WY2);
		  
		  Wxx = (1.-taufrac)*Wxx2+taufrac*Wxx1;
		  
		  //Wxy: 

		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].Wmunu[0][1][2] + xfrac*arena[ix+fac][iy][ieta].Wmunu[0][1][2]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].Wmunu[0][1][2] + xfrac*arena[ix+fac][iy+fac][ieta].Wmunu[0][1][2]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].deltaT[0][1][2] + arena[ix][iy][ieta].deltaW[0][1][2])
		             + xfrac*(arena[ix+fac][iy][ieta].deltaT[0][1][2] + arena[ix+fac][iy][ieta].deltaW[0][1][2]);
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].deltaT[0][1][2] + arena[ix][iy+fac][ieta].deltaW[0][1][2])
		              + xfrac*(arena[ix+fac][iy+fac][ieta].deltaT[0][1][2] + arena[ix+fac][iy+fac][ieta].deltaW[0][1][2]);
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].Wmunu[0][1][2] + xfrac*arena[ix+fac][iy][ieta+fac].Wmunu[0][1][2]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].Wmunu[0][1][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].Wmunu[0][1][2]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].deltaT[0][1][2] + arena[ix][iy][ieta+fac].deltaW[0][1][2])
			          + xfrac*(arena[ix+fac][iy][ieta+fac].deltaT[0][1][2] + arena[ix+fac][iy][ieta+fac].deltaW[0][1][2]);
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].deltaT[0][1][2] + arena[ix][iy+fac][ieta+fac].deltaW[0][1][2])
			          + xfrac*(arena[ix+fac][iy+fac][ieta+fac].deltaT[0][1][2] + arena[ix+fac][iy+fac][ieta+fac].deltaW[0][1][2]);
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wxy[ix][iy] + xfrac*Rneighbor_Wxy[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wxy[ix][iy+fac] + xfrac*Rneighbor_Wxy[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wxy1 = ((1.-etafrac)*WY1+etafrac*WY2);
		  
		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].W_prev[1][2] + xfrac*arena[ix+fac][iy][ieta].W_prev[1][2]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].W_prev[1][2] + xfrac*arena[ix+fac][iy+fac][ieta].W_prev[1][2]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].prev_deltaW[1][2] + arena[ix][iy][ieta].getPrevDeltaT(1, 2, eos))
		      + xfrac*(arena[ix+fac][iy][ieta].prev_deltaW[1][2] + arena[ix+fac][iy][ieta].getPrevDeltaT(1, 2, eos));
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].prev_deltaW[1][2] + arena[ix][iy+fac][ieta].getPrevDeltaT(1, 2, eos))
		      + xfrac*(arena[ix+fac][iy+fac][ieta].prev_deltaW[1][2] + arena[ix+fac][iy+fac][ieta].getPrevDeltaT(1, 2, eos));
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].W_prev[1][2] + xfrac*arena[ix+fac][iy][ieta+fac].W_prev[1][2]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].W_prev[1][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].W_prev[1][2]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].prev_deltaW[1][2] + arena[ix][iy][ieta+fac].getPrevDeltaT(1, 2, eos))
			  + xfrac*(arena[ix+fac][iy][ieta+fac].prev_deltaW[1][2] + arena[ix+fac][iy][ieta+fac].getPrevDeltaT(1, 2, eos));
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].prev_deltaW[1][2] + arena[ix][iy+fac][ieta+fac].getPrevDeltaT(1, 2, eos))
			  + xfrac*(arena[ix+fac][iy+fac][ieta+fac].prev_deltaW[1][2] + arena[ix+fac][iy+fac][ieta+fac].getPrevDeltaT(1, 2, eos));
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wxy_prev[ix][iy] + xfrac*Rneighbor_Wxy_prev[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wxy_prev[ix][iy+fac] + xfrac*Rneighbor_Wxy_prev[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wxy2 = ((1.-etafrac)*WY1+etafrac*WY2);
		  
		  Wxy = (1.-taufrac)*Wxy2+taufrac*Wxy1;
		  
		  //Wxeta:

		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].Wmunu[0][1][3] + xfrac*arena[ix+fac][iy][ieta].Wmunu[0][1][3]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].Wmunu[0][1][3] + xfrac*arena[ix+fac][iy+fac][ieta].Wmunu[0][1][3]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].deltaT[0][1][3] + arena[ix][iy][ieta].deltaW[0][1][3])
		             + xfrac*(arena[ix+fac][iy][ieta].deltaT[0][1][3] + arena[ix+fac][iy][ieta].deltaW[0][1][3]);
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].deltaT[0][1][3] + arena[ix][iy+fac][ieta].deltaW[0][1][3])
		              + xfrac*(arena[ix+fac][iy+fac][ieta].deltaT[0][1][3] + arena[ix+fac][iy+fac][ieta].deltaW[0][1][3]);
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].Wmunu[0][1][3] + xfrac*arena[ix+fac][iy][ieta+fac].Wmunu[0][1][3]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].Wmunu[0][1][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].Wmunu[0][1][3]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].deltaT[0][1][3] + arena[ix][iy][ieta+fac].deltaW[0][1][3])
			          + xfrac*(arena[ix+fac][iy][ieta+fac].deltaT[0][1][3] + arena[ix+fac][iy][ieta+fac].deltaW[0][1][3]);
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].deltaT[0][1][3] + arena[ix][iy+fac][ieta+fac].deltaW[0][1][3])
			          + xfrac*(arena[ix+fac][iy+fac][ieta+fac].deltaT[0][1][3] + arena[ix+fac][iy+fac][ieta+fac].deltaW[0][1][3]);
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wxeta[ix][iy] + xfrac*Rneighbor_Wxeta[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wxeta[ix][iy+fac] + xfrac*Rneighbor_Wxeta[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wxeta1 = ((1.-etafrac)*WY1+etafrac*WY2);
		  
		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].W_prev[1][3] + xfrac*arena[ix+fac][iy][ieta].W_prev[1][3]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].W_prev[1][3] + xfrac*arena[ix+fac][iy+fac][ieta].W_prev[1][3]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].prev_deltaW[1][3] + arena[ix][iy][ieta].getPrevDeltaT(1, 3, eos))
		      + xfrac*(arena[ix+fac][iy][ieta].prev_deltaW[1][3] + arena[ix+fac][iy][ieta].getPrevDeltaT(1, 3, eos));
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].prev_deltaW[1][3] + arena[ix][iy+fac][ieta].getPrevDeltaT(1, 3, eos))
		      + xfrac*(arena[ix+fac][iy+fac][ieta].prev_deltaW[1][3] + arena[ix+fac][iy+fac][ieta].getPrevDeltaT(1, 3, eos));
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].W_prev[1][3] + xfrac*arena[ix+fac][iy][ieta+fac].W_prev[1][3]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].W_prev[1][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].W_prev[1][3]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].prev_deltaW[1][3] + arena[ix][iy][ieta+fac].getPrevDeltaT(1, 3, eos))
			  + xfrac*(arena[ix+fac][iy][ieta+fac].prev_deltaW[1][3] + arena[ix+fac][iy][ieta+fac].getPrevDeltaT(1, 3, eos));
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].prev_deltaW[1][3] + arena[ix][iy+fac][ieta+fac].getPrevDeltaT(1, 3, eos))
			  + xfrac*(arena[ix+fac][iy+fac][ieta+fac].prev_deltaW[1][3] + arena[ix+fac][iy+fac][ieta+fac].getPrevDeltaT(1, 3, eos));
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wxeta_prev[ix][iy] + xfrac*Rneighbor_Wxeta_prev[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wxeta_prev[ix][iy+fac] + xfrac*Rneighbor_Wxeta_prev[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wxeta2 = ((1.-etafrac)*WY1+etafrac*WY2);

		  Wxeta = (1.-taufrac)*Wxeta2+taufrac*Wxeta1;
		
		  //Wyy:
		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].Wmunu[0][2][2] + xfrac*arena[ix+fac][iy][ieta].Wmunu[0][2][2]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].Wmunu[0][2][2] + xfrac*arena[ix+fac][iy+fac][ieta].Wmunu[0][2][2]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].deltaT[0][2][2] + arena[ix][iy][ieta].deltaW[0][2][2])
		             + xfrac*(arena[ix+fac][iy][ieta].deltaT[0][2][2] + arena[ix+fac][iy][ieta].deltaW[0][2][2]);
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].deltaT[0][2][2] + arena[ix][iy+fac][ieta].deltaW[0][2][2])
		              + xfrac*(arena[ix+fac][iy+fac][ieta].deltaT[0][2][2] + arena[ix+fac][iy+fac][ieta].deltaW[0][2][2]);
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].Wmunu[0][2][2] + xfrac*arena[ix+fac][iy][ieta+fac].Wmunu[0][2][2]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].Wmunu[0][2][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].Wmunu[0][2][2]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].deltaT[0][2][2] + arena[ix][iy][ieta+fac].deltaW[0][2][2])
			          + xfrac*(arena[ix+fac][iy][ieta+fac].deltaT[0][2][2] + arena[ix+fac][iy][ieta+fac].deltaW[0][2][2]);
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].deltaT[0][2][2] + arena[ix][iy+fac][ieta+fac].deltaW[0][2][2])
			          + xfrac*(arena[ix+fac][iy+fac][ieta+fac].deltaT[0][2][2] + arena[ix+fac][iy+fac][ieta+fac].deltaW[0][2][2]);
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wyy[ix][iy] + xfrac*Rneighbor_Wyy[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wyy[ix][iy+fac] + xfrac*Rneighbor_Wyy[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wyy1 = ((1.-etafrac)*WY1+etafrac*WY2);
		
		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].W_prev[2][2] + xfrac*arena[ix+fac][iy][ieta].W_prev[2][2]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].W_prev[2][2] + xfrac*arena[ix+fac][iy+fac][ieta].W_prev[2][2]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].prev_deltaW[2][2] + arena[ix][iy][ieta].getPrevDeltaT(2, 2, eos))
		      + xfrac*(arena[ix+fac][iy][ieta].prev_deltaW[2][2] + arena[ix+fac][iy][ieta].getPrevDeltaT(2, 2, eos));
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].prev_deltaW[2][2] + arena[ix][iy+fac][ieta].getPrevDeltaT(2, 2, eos))
		      + xfrac*(arena[ix+fac][iy+fac][ieta].prev_deltaW[2][2] + arena[ix+fac][iy+fac][ieta].getPrevDeltaT(2, 2, eos));
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].W_prev[2][2] + xfrac*arena[ix+fac][iy][ieta+fac].W_prev[2][2]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].W_prev[2][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].W_prev[2][2]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].prev_deltaW[2][2] + arena[ix][iy][ieta+fac].getPrevDeltaT(2, 2, eos))
			  + xfrac*(arena[ix+fac][iy][ieta+fac].prev_deltaW[2][2] + arena[ix+fac][iy][ieta+fac].getPrevDeltaT(2, 2, eos));
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].prev_deltaW[2][2] + arena[ix][iy+fac][ieta+fac].getPrevDeltaT(2, 2, eos))
			  + xfrac*(arena[ix+fac][iy+fac][ieta+fac].prev_deltaW[2][2] + arena[ix+fac][iy+fac][ieta+fac].getPrevDeltaT(2, 2, eos));
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wyy_prev[ix][iy] + xfrac*Rneighbor_Wyy_prev[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wyy_prev[ix][iy+fac] + xfrac*Rneighbor_Wyy_prev[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wyy2 = ((1.-etafrac)*WY1+etafrac*WY2);
		  
		  Wyy = (1.-taufrac)*Wyy2+taufrac*Wyy1;

		  //Wyeta:

		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].Wmunu[0][2][3] + xfrac*arena[ix+fac][iy][ieta].Wmunu[0][2][3]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].Wmunu[0][2][3] + xfrac*arena[ix+fac][iy+fac][ieta].Wmunu[0][2][3]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].deltaT[0][2][3] + arena[ix][iy][ieta].deltaW[0][2][3])
		             + xfrac*(arena[ix+fac][iy][ieta].deltaT[0][2][3] + arena[ix+fac][iy][ieta].deltaW[0][2][3]);
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].deltaT[0][2][3] + arena[ix][iy+fac][ieta].deltaW[0][2][3])
		              + xfrac*(arena[ix+fac][iy+fac][ieta].deltaT[0][2][3] + arena[ix+fac][iy+fac][ieta].deltaW[0][2][3]);
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].Wmunu[0][2][3] + xfrac*arena[ix+fac][iy][ieta+fac].Wmunu[0][2][3]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].Wmunu[0][2][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].Wmunu[0][2][3]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].deltaT[0][2][3] + arena[ix][iy][ieta+fac].deltaW[0][2][3])
			          + xfrac*(arena[ix+fac][iy][ieta+fac].deltaT[0][2][3] + arena[ix+fac][iy][ieta+fac].deltaW[0][2][3]);
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].deltaT[0][2][3] + arena[ix][iy+fac][ieta+fac].deltaW[0][2][3])
			          + xfrac*(arena[ix+fac][iy+fac][ieta+fac].deltaT[0][2][3] + arena[ix+fac][iy+fac][ieta+fac].deltaW[0][2][3]);
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wyeta[ix][iy] + xfrac*Rneighbor_Wyeta[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wyeta[ix][iy+fac] + xfrac*Rneighbor_Wyeta[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wyeta1 = ((1.-etafrac)*WY1+etafrac*WY2);

		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].W_prev[2][3] + xfrac*arena[ix+fac][iy][ieta].W_prev[2][3]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].W_prev[2][3] + xfrac*arena[ix+fac][iy+fac][ieta].W_prev[2][3]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].prev_deltaW[2][3] + arena[ix][iy][ieta].getPrevDeltaT(2, 3, eos))
		      + xfrac*(arena[ix+fac][iy][ieta].prev_deltaW[2][3] + arena[ix+fac][iy][ieta].getPrevDeltaT(2, 3, eos));
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].prev_deltaW[2][3] + arena[ix][iy+fac][ieta].getPrevDeltaT(2, 3, eos))
		      + xfrac*(arena[ix+fac][iy+fac][ieta].prev_deltaW[2][3] + arena[ix+fac][iy+fac][ieta].getPrevDeltaT(2, 3, eos));
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].W_prev[2][3] + xfrac*arena[ix+fac][iy][ieta+fac].W_prev[2][3]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].W_prev[2][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].W_prev[2][3]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].prev_deltaW[2][3] + arena[ix][iy][ieta+fac].getPrevDeltaT(2, 3, eos))
			  + xfrac*(arena[ix+fac][iy][ieta+fac].prev_deltaW[2][3] + arena[ix+fac][iy][ieta+fac].getPrevDeltaT(2, 3, eos));
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].prev_deltaW[2][3] + arena[ix][iy+fac][ieta+fac].getPrevDeltaT(2, 3, eos))
			  + xfrac*(arena[ix+fac][iy+fac][ieta+fac].prev_deltaW[2][3] + arena[ix+fac][iy+fac][ieta+fac].getPrevDeltaT(2, 3, eos));
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wyeta_prev[ix][iy] + xfrac*Rneighbor_Wyeta_prev[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wyeta_prev[ix][iy+fac] + xfrac*Rneighbor_Wyeta_prev[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wyeta2 = ((1.-etafrac)*WY1+etafrac*WY2);

		  Wyeta = (1.-taufrac)*Wyeta2+taufrac*Wyeta1;

		  //Wetaeta:

		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].Wmunu[0][3][3] + xfrac*arena[ix+fac][iy][ieta].Wmunu[0][3][3]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].Wmunu[0][3][3] + xfrac*arena[ix+fac][iy+fac][ieta].Wmunu[0][3][3]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].deltaT[0][3][3] + arena[ix][iy][ieta].deltaW[0][3][3])
		             + xfrac*(arena[ix+fac][iy][ieta].deltaT[0][3][3] + arena[ix+fac][iy][ieta].deltaW[0][3][3]);
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].deltaT[0][3][3] + arena[ix][iy+fac][ieta].deltaW[0][3][3])
		              + xfrac*(arena[ix+fac][iy+fac][ieta].deltaT[0][3][3] + arena[ix+fac][iy+fac][ieta].deltaW[0][3][3]);
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].Wmunu[0][3][3] + xfrac*arena[ix+fac][iy][ieta+fac].Wmunu[0][3][3]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].Wmunu[0][3][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].Wmunu[0][3][3]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].deltaT[0][3][3] + arena[ix][iy][ieta+fac].deltaW[0][3][3])
			          + xfrac*(arena[ix+fac][iy][ieta+fac].deltaT[0][3][3] + arena[ix+fac][iy][ieta+fac].deltaW[0][3][3]);
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].deltaT[0][3][3] + arena[ix][iy+fac][ieta+fac].deltaW[0][3][3])
			          + xfrac*(arena[ix+fac][iy+fac][ieta+fac].deltaT[0][3][3] + arena[ix+fac][iy+fac][ieta+fac].deltaW[0][3][3]);
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wetaeta[ix][iy] + xfrac*Rneighbor_Wetaeta[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wetaeta[ix][iy+fac] + xfrac*Rneighbor_Wetaeta[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wetaeta1 = ((1.-etafrac)*WY1+etafrac*WY2);
		  
		  WX1 = ((1.-xfrac)*arena[ix][iy][ieta].W_prev[3][3] + xfrac*arena[ix+fac][iy][ieta].W_prev[3][3]);
		  WX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].W_prev[3][3] + xfrac*arena[ix+fac][iy+fac][ieta].W_prev[3][3]);

		  if(DATA->freezeOutMethod==4){
		    WX1 += (1.-xfrac)*(arena[ix][iy][ieta].prev_deltaW[3][3] + arena[ix][iy][ieta].getPrevDeltaT(3, 3, eos))
		      + xfrac*(arena[ix+fac][iy][ieta].prev_deltaW[3][3] + arena[ix+fac][iy][ieta].getPrevDeltaT(3, 3, eos));
		    WX2 += (1.-xfrac)*(arena[ix][iy+fac][ieta].prev_deltaW[3][3] + arena[ix][iy+fac][ieta].getPrevDeltaT(3, 3, eos))
		      + xfrac*(arena[ix+fac][iy+fac][ieta].prev_deltaW[3][3] + arena[ix+fac][iy+fac][ieta].getPrevDeltaT(3, 3, eos));
		  }

		  if (ieta<neta-fac)
		    {
		      WX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].W_prev[3][3] + xfrac*arena[ix+fac][iy][ieta+fac].W_prev[3][3]);
		      WX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].W_prev[3][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].W_prev[3][3]);

		      if(DATA->freezeOutMethod==4){
			WX3 += (1.-xfrac)*(arena[ix][iy][ieta+fac].prev_deltaW[3][3] + arena[ix][iy][ieta+fac].getPrevDeltaT(3, 3, eos))
			  + xfrac*(arena[ix+fac][iy][ieta+fac].prev_deltaW[3][3] + arena[ix+fac][iy][ieta+fac].getPrevDeltaT(3, 3, eos));
			WX4 += (1.-xfrac)*(arena[ix][iy+fac][ieta+fac].prev_deltaW[3][3] + arena[ix][iy+fac][ieta+fac].getPrevDeltaT(3, 3, eos))
			  + xfrac*(arena[ix+fac][iy+fac][ieta+fac].prev_deltaW[3][3] + arena[ix+fac][iy+fac][ieta+fac].getPrevDeltaT(3, 3, eos));
		      }
		    }
		  else
		    {
		      WX3 = ((1.-xfrac)*Rneighbor_Wetaeta_prev[ix][iy] + xfrac*Rneighbor_Wetaeta_prev[ix+fac][iy]);
		      WX4 = ((1.-xfrac)*Rneighbor_Wetaeta_prev[ix][iy+fac] + xfrac*Rneighbor_Wetaeta_prev[ix+fac][iy+fac]);
		    }
		  
		  WY1 = ((1.-yfrac)*WX1+yfrac*WX2);
		  WY2 = ((1.-yfrac)*WX3+yfrac*WX4);

		  Wetaeta2 = ((1.-etafrac)*WY1+etafrac*WY2);
		
		  Wetaeta = (1.-taufrac)*Wetaeta2+taufrac*Wetaeta1;
		  
		  if (DATA->whichEOS==1)
		    {
		      TFO = eos->interpolate(epsFO, rhob, 0);
		      muB = eos->interpolate(epsFO, rhob, 1);
		    }
		  else if (DATA->whichEOS==2)
		    {
		      TFO = eos->interpolate2(epsFO, 0., 1);
		      muB = 0.0;
		    }
		  
		  double dT, de;
		  if(DATA->fluctuatingHydroFlag == 2){
		    dp = 0.;
		    dux = 0.;
		    duy = 0.;
		    dueta = 0.;
		    dWtautau = 0.;
		    dWtaux = 0.;
		    dWtauy = 0.;
		    dWtaueta = 0.;
		    dWxx = 0.;
		    dWxy = 0.;
		    dWxeta = 0.;
		    dWyy = 0.;
		    dWyeta = 0.;
		    dWetaeta = 0.;
		    
		    double xfracs[2]; xfracs[0] = 1. - xfrac; xfracs[1] = xfrac;
		    double yfracs[2]; yfracs[0] = 1. - yfrac; yfracs[1] = yfrac;
		    double etafracs[2]; etafracs[0] = 1. - etafrac; etafracs[1] = etafrac;
		    for(int jx = 0; jx<2; jx++){
		      for(int jy = 0; jy<2; jy++){
			if(ieta < neta-fac){
			  for(int jeta = 0; jeta<2; jeta++){
			    dp += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaP[0]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaP_prev;
			    dux += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaU[0][1]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaU_prev[1];
			    duy += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaU[0][2]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaU_prev[2];
			    dueta += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaU[0][3]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaU_prev[3];
			    dWtautau += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW[0][0][0]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW_prev[0][0];
			    dWtaux += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW[0][0][1]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW_prev[0][1];
			    dWtauy += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW[0][0][2]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW_prev[0][2];
			    dWtaueta += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW[0][0][3]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW_prev[0][3];
			    dWxx += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW[0][1][1]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW_prev[1][1];
			    dWxy += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW[0][1][2]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW_prev[1][2];
			    dWxeta += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW[0][1][3]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW_prev[1][3];
			    dWyy += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW[0][1][2]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW_prev[1][2];
			    dWyeta += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW[0][2][3]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW_prev[2][3];
			    dWetaeta += taufrac*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW[0][3][3]
			      + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[jeta]*arena[ix+jx][iy+jy][ieta+jeta].deltaW_prev[3][3];
			  }
			}
			else{
			  dp += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaP[0]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dp[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaP_prev
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dp_prev[ix+jx][iy+jy];
			  dux += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaU[0][1]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dux[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaU_prev[1]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dux_prev[ix+jx][iy+jy];
			  duy += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaU[0][2]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_duy[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaU_prev[2]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_duy_prev[ix+jx][iy+jy];
			  dueta += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaU[0][3]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dueta[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaU_prev[3]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dueta_prev[ix+jx][iy+jy];
			  dWtautau += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW[0][0][0]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWtautau[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW_prev[0][0]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWtautau_prev[ix+jx][iy+jy];
			  dWtaux += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW[0][0][1]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWtaux[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW_prev[0][1]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWtaux_prev[ix+jx][iy+jy];
			  dWtauy += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW[0][0][2]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWtauy[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW_prev[0][2]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWtauy_prev[ix+jx][iy+jy];
			  dWtaueta += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW[0][0][3]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWtaueta[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW_prev[0][3]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWtaueta_prev[ix+jx][iy+jy];
			  dWxx += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW[0][1][1]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWxx[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW_prev[1][1]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWxx_prev[ix+jx][iy+jy];
			  dWxy += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW[0][1][2]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWxy[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW_prev[1][2]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWxy_prev[ix+jx][iy+jy];
			  dWxeta += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW[0][1][3]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWxeta[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW_prev[1][3]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWxeta_prev[ix+jx][iy+jy];
			  dWyy += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW[0][2][2]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWyy[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW_prev[2][2]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWyy_prev[ix+jx][iy+jy];
			  dWyeta += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW[0][2][3]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWyeta[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW_prev[2][3]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWyeta_prev[ix+jx][iy+jy];
			  dWetaeta += taufrac*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW[0][3][3]
			    + taufrac*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWetaeta[ix+jx][iy+jy]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[0]*arena[ix+jx][iy+jy][ieta].deltaW_prev[3][3]
			    + (1. - taufrac)*xfracs[jx]*yfracs[jy]*etafracs[1]*Rneighbor_dWetaeta_prev[ix+jx][iy+jy];
			}
		      }
		    }
		  
		    dT = dp/eos->s_func(epsFO, eos->p_func(epsFO, rhob), rhob, DATA);
		    de = dp/eos->get_dpOverde2(epsFO, rhob);
		    dutau = (ux*dux + uy*duy + ueta*dueta)/utau;

		    dWtaux = (ux*dWxx + uy*dWxy + ueta*dWxeta - dutau*Wtaux + dux*Wxx + duy*Wxy + dueta*Wxeta)/utau;
		    dWtauy = (ux*dWxy + uy*dWyy + ueta*dWyeta - dutau*Wtauy + dux*Wxy + duy*Wyy + dueta*Wyeta)/utau;
		    dWtaueta = (ux*dWxeta + uy*dWyeta + ueta*dWetaeta - dutau*Wtaueta + dux*Wxeta + duy*Wyeta + dueta*Wetaeta)/utau;
		    dWtautau = (ux*dWtaux + uy*dWtauy + ueta*dWtaueta - dutau*Wtautau + dux*Wtaux + duy*Wtauy + dueta*Wtaueta)/utau;
		  }


		  P=eos->p_func(epsFO, rhob);
		  sFO=eos->s_func(epsFO, P, rhob, DATA);

		  if (fabs(FULLSU[0])>DX*DY*DETA+0.01)
		    {   
		      fprintf(stderr,"problem: volume in tau direction %f > DX*DY*DETA = %f\n",fabs(FULLSU[0]),DX*DY*DETA);
		      //FULLSU[0] = DX*DY*DETA*(FULLSU[0])/fabs(FULLSU[0]);
		    }
		  if (fabs(FULLSU[1])>DTAU*DY*DETA+0.01)
		    {
		      fprintf(stderr,"problem: volume in x direction %f > DTAU*DY*DETA = %f\n",fabs(FULLSU[1]),DTAU*DY*DETA);
		      //FULLSU[1] = DTAU*DY*DETA*(FULLSU[1])/fabs(FULLSU[1]);
		    }
		  if (fabs(FULLSU[2])>DX*DTAU*DETA+0.01)
		    {
		      fprintf(stderr,"problem: volume in y direction %f > DX*DTAU*DETA = %f\n",fabs(FULLSU[2]),DX*DTAU*DETA);
		      //FULLSU[2] = DX*DTAU*DETA*(FULLSU[2])/fabs(FULLSU[2]);
		    }
		  if (fabs(FULLSU[3])>DX*DY*DTAU+0.01)
		    {
		      fprintf(stderr,"problem: volume in eta direction %f > DX*DY*DTAU = %f\n",fabs(FULLSU[3]),DX*DY*DTAU);
		      //FULLSU[3] = DX*DY*DTAU*(FULLSU[3])/fabs(FULLSU[3]);
		    }
		  fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
			  tauf, xf, yf, etaf, FULLSU[0],FULLSU[1],FULLSU[2],FULLSU[3],
			  utau, ux, uy, ueta, epsFO, TFO, muB, sFO, Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxeta, Wyy, Wyeta, Wetaeta);
		  if(DATA->fluctuatingHydroFlag == 2){
		    ds_file << dT << " " << dutau << " " << dux << " " << duy << " " << dueta
			    << " " << de << " " << dp
			    << " " << dWtautau << " " << dWtaux << " " << dWtauy << " " << dWtaueta
			    << " " << dWxx << " " << dWxy << " " << dWxeta
			    << " " << dWyy << " " << dWyeta
			    << " " << dWetaeta << endl;
		    //
		    //fprintf(ds_file, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
		    //	    dT, dutau, dux, duy, dueta, de, dp,
		    //	    dWtautau, dWtaux, dWtauy, dWtaueta,
		    //	    dWxx, dWxy, dWxeta,
		    //	    dWyy, dWyeta,
		    //	    dWetaeta);
		  }
		  
		  if(fabs(x)<0.05 && (fabs(y)<0.05))
		    //if(fabs(x)<0.05 && (fabs(y)<0.05))
		    {
		      fprintf(stderr,"tau_F(%f,%f,%f)=%f\n",xf,yf,etaf,tauf);
		      fprintf(t_file,"%f %f %f %f %f %f %f %f\n",tauf, xf, yf, etaf, FULLSU[0],FULLSU[1],FULLSU[2],FULLSU[3]);
		      //fprintf(stderr,"ix=%d, iy=%d, ieta=%d\n",ix,iy,ieta);
		    }
		  if(fabs(eta)<0.05 && (fabs(y)<0.05))
		    {
		      fprintf(t2_file,"%f %f %f %f %f %f %f %f\n",tauf, xf, yf, etaf, FULLSU[0],FULLSU[1],FULLSU[2],FULLSU[3]);
		    }
		  if(fabs(eta)<0.05 && (fabs(x)<0.1))
		    {
		      fprintf(t3_file,"%f %f %f %f %f %f %f %f\n",tauf, xf, yf, etaf, FULLSU[0],FULLSU[1],FULLSU[2],FULLSU[3]);
		    }
		}
	    }
	}
    }
  
  int intersectionsArray[1];
  int allIntersectionsArray[1];
  intersectionsArray[0] = intersections;
  if (rank!=0) //send to rank 0
    {
      MPI::COMM_WORLD.Send(intersectionsArray,1,MPI::INT,0,1);
    }
  if (rank==0) // receive from all ranks >0
    {
      //fprintf(stderr,"percent error=%f\n", 100*(warnings*1.0)/(cells*1.0));
		  
      int allIntersections = intersections;
      for (from = 1; from<size; from ++)
	{
	  MPI::COMM_WORLD.Recv(intersectionsArray,1,MPI::INT,from,1);
	  allIntersections+=intersectionsArray[0];
	}
      for (from = 1; from<size; from ++)
	{
	  allIntersectionsArray[0]=allIntersections;
	  MPI::COMM_WORLD.Send(allIntersectionsArray,1,MPI::INT,from,2);
	}
      if (allIntersections==0 && DATA->Initial_profile != 6)
	{
	  cout << "All cells frozen out. Exiting." << endl;
	  MPI::Finalize();
 	  exit(1);
	}
    }
  if (rank!=0)
    {
      MPI::COMM_WORLD.Recv(allIntersectionsArray,1,MPI::INT,0,2);
      if (allIntersectionsArray[0]==0 && DATA->Initial_profile != 6)
	{
	  cout << "All cells frozen out. Exiting." << endl;
	  MPI::Finalize();
 	  exit(1);
	}
    }

  delete[] packageWtautau;
  delete[] packageWtaux;
  delete[] packageWtauy;
  delete[] packageWtaueta;
  delete[] packageWxx;
  delete[] packageWxy;
  delete[] packageWxeta;
  delete[] packageWyy;
  delete[] packageWyeta; 
  delete[] packageWetaeta;

  delete[] packageWtautau_prev;
  delete[] packageWtaux_prev;
  delete[] packageWtauy_prev;
  delete[] packageWtaueta_prev;
  delete[] packageWxx_prev;
  delete[] packageWxy_prev;
  delete[] packageWxeta_prev;
  delete[] packageWyy_prev;
  delete[] packageWyeta_prev; 
  delete[] packageWetaeta_prev;

  delete[] package;
  delete[] package_prev;
  delete[] packageutau;
  delete[] packageux;
  delete[] packageuy;
  delete[] packageueta;
  delete[] packagerhob;

  delete[] packageutau_prev;
  delete[] packageux_prev;
  delete[] packageuy_prev;
  delete[] packageueta_prev;
  delete[] packagerhob_prev;

  util->mtx_free(Rneighbor_eps,nx+1,ny+1);
  util->mtx_free(Rneighbor_eps_prev,nx+1,ny+1);
  util->mtx_free(Rneighbor_utau,nx+1,ny+1);
  util->mtx_free(Rneighbor_ux,nx+1,ny+1);
  util->mtx_free(Rneighbor_uy,nx+1,ny+1);
  util->mtx_free(Rneighbor_ueta,nx+1,ny+1);
  util->mtx_free(Rneighbor_rhob,nx+1,ny+1);
 
  util->mtx_free(Rneighbor_utau_prev,nx+1,ny+1);
  util->mtx_free(Rneighbor_ux_prev,nx+1,ny+1);
  util->mtx_free(Rneighbor_uy_prev,nx+1,ny+1);
  util->mtx_free(Rneighbor_ueta_prev,nx+1,ny+1);
  util->mtx_free(Rneighbor_rhob_prev,nx+1,ny+1);
 
 
  fclose(t2_file);
  fclose(t_file);
  fclose(s_file);
}

// Modified version of simplified freeze out  (M. Luzum, 04/2013)
// Every element of the freezeout surface is a rectangular cuboid
// located mid-way between cells/grid points.
int Evolve::FindFreezeOutSurfaceOfHypercubes(double tau, InitData *DATA, Grid ***arena, int size, int rank)
{

  cerr << "Starting FindFreezeOutSurfaceOfHypercubes..." << endl;

  stringstream numIn;
  numIn << DATA->seed;
  
  stringstream strs_name;
  strs_name << "surface" << rank << "_" << numIn.str() <<  ".dat";
  string s_name = strs_name.str();
  
  ofstream s_file;
  s_file.open(s_name.c_str() , ios::out | ios::app );
  
  stringstream dstrs_name;
  dstrs_name << "dsurface" << rank << "_" << numIn.str() << ".dat";
  string ds_name = dstrs_name.str();
  
  ofstream ds_file;
  ds_file.open(ds_name.c_str() , ios::out | ios::app );
  
  stringstream fr_name;
  fr_name << "freezeout" << rank << "_" << numIn.str() << ".dat";
  string freeze_name = fr_name.str();
  
  ofstream freeze_file;
  freeze_file.open(freeze_name.c_str() , ios::out | ios::app );
  //   }
  
  int frozen;
  frozen = 1;
  int allfrozen;
  allfrozen = 1;
  
  int ix, iy, ieta, nx, ny, neta;
  double x, y, eta;
  //   double epsFO=DATA->epsilonFreeze/hbarc;
  double tauf, xf, yf, etaf;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;
  double FULLSU[4];
  int fac;
  double DX, DY, DETA, DTAU, SIG;
  double iepsFO;
  double home, neighbor, FO;
  int nix, niy, nieta;
  double Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxeta, Wyy, Wyeta, Wetaeta;
  double rhob, utau, ux, uy, ueta, TFO, muB, eps_plus_p_over_T_FO;
  //Added 7/15/2014 by CFY, for the thermally fluctuating case:
  double dT, dutau, dux, duy, dueta, de, dp;
  double dWtautau, dWtaux, dWtauy, dWtaueta, dWxx, dWxy, dWxeta, dWyy, dWyeta, dWetaeta;
  double pi_b; // bulk
  //   int shown;
  
  fac=1; // Non-unity value does not currently work
  facTau = DATA->facTau;
  DX=fac*DATA->delta_x;
  DY=fac*DATA->delta_y;
  DETA=fac*DATA->delta_eta;
  DTAU=facTau*DATA->delta_tau;
  
  double maxDETA = DATA->delta_eta;// maximum size of cuboid in eta direction.  If DETA > maxDETA, split the surface into several identical sections spread across eta. -CFY For now, max_delta_eta == delta_eta.
  int subsections = floor(DETA/maxDETA) + 1;// subdivide the blocks into this many subdivisions in eta
 
  //  cout << "subsections = " << subsections << endl;
 
  //   fprintf(stderr,"DTAU=%f\n", DTAU);
  //   fprintf(stderr,"DX=%f\n", DX);
  //   fprintf(stderr,"DY=%f\n", DY);
  //   fprintf(stderr,"DETA=%f\n", DETA);
  //   shown = 0;
  
  int maxEta;
  
  //  MPI code copied from FindFreezeOutSurface2
  //*******************************************
  //  int ix, iy, ieta, nx, ny, neta, i, alpha, iflag;
  int sizeOfData = (nx+1)*(ny+1);
  int position;
  double *package;
  //  double *packageutau;
  double *packageux;
  double *packageuy;
  double *packageueta;
  double *packagerhob;
  double *packagePi_b;

  double **Rneighbor_eps;
  //  double **Rneighbor_utau;
  double **Rneighbor_ux;
  double **Rneighbor_uy;
  double **Rneighbor_ueta;
  double **Rneighbor_rhob;
  double **Rneighbor_Pi_b;
 
  //  double *packageWtautau;
  //  double *packageWtaux;
  //  double *packageWtauy;
  //  double *packageWtaueta;
  double *packageWxx;
  double *packageWxy;
  double *packageWxeta;
  double *packageWyy;
  double *packageWyeta;
  //  double *packageWetaeta;

  //  double **Rneighbor_Wtautau;
  //  double **Rneighbor_Wtaux;
  //  double **Rneighbor_Wtauy;
  //  double **Rneighbor_Wtaueta;
  double **Rneighbor_Wxx;
  double **Rneighbor_Wxy;
  double **Rneighbor_Wxeta;
  double **Rneighbor_Wyy;
  double **Rneighbor_Wyeta;
  //  double **Rneighbor_Wetaeta;
  //dT, dutau, dux, duy, dueta, de, dp:
  double *packagedT;
  double *packagedutau;
  double *packagedux;
  double *packageduy;
  double *packagedueta;
  double *packagede;
  double *packagedp;

  double **Rneighbor_dT;
  double **Rneighbor_dutau;
  double **Rneighbor_dux;
  double **Rneighbor_duy;
  double **Rneighbor_dueta;
  double **Rneighbor_de;
  double **Rneighbor_dp;

  //dWmunu:
  double *packagedWtautau;
  double *packagedWtaux;
  double *packagedWtauy;
  double *packagedWtaueta;
  double *packagedWxx;
  double *packagedWxy;
  double *packagedWxeta;
  double *packagedWyy;
  double *packagedWyeta;
  double *packagedWetaeta;

  double **Rneighbor_dWtautau;
  double **Rneighbor_dWtaux;
  double **Rneighbor_dWtauy;
  double **Rneighbor_dWtaueta;
  double **Rneighbor_dWxx;
  double **Rneighbor_dWxy;
  double **Rneighbor_dWxeta;
  double **Rneighbor_dWyy;
  double **Rneighbor_dWyeta;
  double **Rneighbor_dWetaeta;

  //cerr << "rank = " << rank << endl;

  package = new double[sizeOfData];
  //  packageutau = new double[sizeOfData];
  packageux = new double[sizeOfData];
  packageuy = new double[sizeOfData];
  packageueta = new double[sizeOfData];
  packagerhob = new double[sizeOfData];
  packagePi_b = new double[sizeOfData];

  //  packageWtautau = new double[sizeOfData];
  //  packageWtaux = new double[sizeOfData];
  //  packageWtauy = new double[sizeOfData];
  //  packageWtaueta = new double[sizeOfData];
  packageWxx = new double[sizeOfData];
  packageWxy = new double[sizeOfData];
  packageWxeta = new double[sizeOfData];
  packageWyy = new double[sizeOfData];
  packageWyeta = new double[sizeOfData];
  //  packageWetaeta = new double[sizeOfData];
  packagedT = new double[sizeOfData];
  packagedutau = new double[sizeOfData];
  packagedux = new double[sizeOfData];
  packageduy = new double[sizeOfData];
  packagedueta = new double[sizeOfData];
  packagede = new double[sizeOfData];
  packagedp = new double[sizeOfData];
  packagedWtautau = new double[sizeOfData];
  packagedWtaux = new double[sizeOfData];
  packagedWtauy = new double[sizeOfData];
  packagedWtaueta = new double[sizeOfData];
  packagedWxx = new double[sizeOfData];
  packagedWxy = new double[sizeOfData];
  packagedWxeta = new double[sizeOfData];
  packagedWyy = new double[sizeOfData];
  packagedWyeta = new double[sizeOfData];
  packagedWetaeta = new double[sizeOfData];
  
  Rneighbor_eps = util->mtx_malloc(nx+1,ny+1);
  //  Rneighbor_utau = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_ux = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_uy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_ueta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_rhob = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Pi_b = util->mtx_malloc(nx+1,ny+1);
 
  //  Rneighbor_Wtautau = util->mtx_malloc(nx+1,ny+1);
  //  Rneighbor_Wtaux = util->mtx_malloc(nx+1,ny+1);
  //  Rneighbor_Wtauy = util->mtx_malloc(nx+1,ny+1);
  //  Rneighbor_Wtaueta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wxx = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wxy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wxeta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wyy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_Wyeta = util->mtx_malloc(nx+1,ny+1);
  //  Rneighbor_Wetaeta = util->mtx_malloc(nx+1,ny+1);

  Rneighbor_dT = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dutau = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dux = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_duy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dueta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_de = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dp = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dWtautau = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dWtaux = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dWtauy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dWtaueta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dWxx = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dWxy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dWxeta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dWyy = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dWyeta = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_dWetaeta = util->mtx_malloc(nx+1,ny+1);
  
  // receive from the right / send to the left
  int from = rank+1;
  int to = rank-1;
    
  if(size > 1){
    
    cerr << "from = " << from << ", to = " << to << endl;
    
    // get cells from neighboring processors
    if ( rank != 0 )
      {
	//      cout << " sending to the left on rank " << rank << endl;
	for(ix=0; ix<=nx; ix++)
	  {
	    for(iy=0; iy<=ny; iy++)
	      {
		//if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
		//cout << "arena[ix][iy][0].TJb[i][alpha][0]=" << arena[ix][iy][0].TJb[i][alpha][0] << endl;
		position = ix + ((nx+1)*iy);
		package[position] = arena[ix][iy][0].epsilon;
		//      packageutau[position] = arena[ix][iy][0].u[0][0];
		packageux[position] = arena[ix][iy][0].u[0][1];
		packageuy[position] = arena[ix][iy][0].u[0][2];
		packageueta[position] = arena[ix][iy][0].u[0][3];
		if(DATA->turn_on_rhob) packagerhob[position] = arena[ix][iy][0].rhob;
		//if(DATA->turn_on_bulk) packagePi_b[position] = arena[ix][iy][0].pi_b[0];
		if(DATA->turn_on_shear) 
		  {
		    //      packageWtautau[position] = arena[ix][iy][0].Wmunu[0][0][0];
		    //      packageWtaux[position] = arena[ix][iy][0].Wmunu[0][0][1];
		    //      packageWtauy[position] = arena[ix][iy][0].Wmunu[0][0][2];
		    //      packageWtaueta[position] = arena[ix][iy][0].Wmunu[0][0][3];
		    packageWxx[position] = arena[ix][iy][0].Wmunu[0][1][1];
		    packageWxy[position] = arena[ix][iy][0].Wmunu[0][1][2];
		    packageWxeta[position] = arena[ix][iy][0].Wmunu[0][1][3];
		    packageWyy[position] = arena[ix][iy][0].Wmunu[0][2][2];
		    packageWyeta[position] = arena[ix][iy][0].Wmunu[0][2][3];
		    //      packageWetaeta[position] = arena[ix][iy][0].Wmunu[0][3][3];
		  }
		if(DATA->fluctuatingHydroFlag == 2){
		  if(arena[ix][iy][0].deltaP[0] == 0.){packagedT[position] = 0.;}
		  else{packagedT[position] = 
		      arena[ix][iy][0].deltaP[0]
		      /eos->s_func(arena[ix][iy][0].epsilon, 
				   eos->p_func(arena[ix][iy][0].epsilon, 
					       arena[ix][iy][0].rhob), 
				   arena[ix][iy][0].rhob, DATA);}
		  //cerr << "packagedT[" << position << "] = " << packagedT[position] << ", deltaP = " << arena[ix][iy][0].deltaP[0] << endl;
		  packagedutau[position] = arena[ix][iy][0].deltaU[0][0];
		  packagedux[position] = arena[ix][iy][0].deltaU[0][1];
		  packageduy[position] = arena[ix][iy][0].deltaU[0][2];
		  packagedueta[position] = arena[ix][iy][0].deltaU[0][3];
		  packagede[position] = arena[ix][iy][0].deltaP[0]/eos->get_dpOverde2(fmax(0.03, arena[ix][iy][0].epsilon), 
										      arena[ix][iy][0].rhob);
		  packagedp[position] = arena[ix][iy][0].deltaP[0];
		  packagedWtautau[position] = arena[ix][iy][0].deltaW[0][0][0];
		  packagedWtaux[position] = arena[ix][iy][0].deltaW[0][0][1];
		  packagedWtauy[position] = arena[ix][iy][0].deltaW[0][0][2];
		  packagedWtaueta[position] = arena[ix][iy][0].deltaW[0][0][3];
		  packagedWxx[position] = arena[ix][iy][0].deltaW[0][1][1];
		  packagedWxy[position] = arena[ix][iy][0].deltaW[0][1][2];
		  packagedWxeta[position] = arena[ix][iy][0].deltaW[0][1][3];
		  packagedWyy[position] = arena[ix][iy][0].deltaW[0][2][2];
		  packagedWyeta[position] = arena[ix][iy][0].deltaW[0][2][3];
		  packagedWetaeta[position] = arena[ix][iy][0].deltaW[0][3][3];
		}
		
		
	      }
	  }
	
	//cerr << "In (rank != 0), rank = " << rank << ", all OK before the first MPI call..." << endl;
	//cerr << "In (rank != 0), packagedp = " << packagedp[129*64] << endl;
	//cerr << "In (rank != 0), packagedp = " << packagedp[129*64] << endl;
	MPI::COMM_WORLD.Send(package,sizeOfData,MPI::DOUBLE,to,1);
	//      MPI::COMM_WORLD.Send(packageutau,sizeOfData,MPI::DOUBLE,to,3);
	MPI::COMM_WORLD.Send(packageux,sizeOfData,MPI::DOUBLE,to,2);
	MPI::COMM_WORLD.Send(packageuy,sizeOfData,MPI::DOUBLE,to,3);
	//cerr << "In (rank != 0), rank = " << rank << ", all OK after the third MPI call..." << endl;
	MPI::COMM_WORLD.Send(packageueta,sizeOfData,MPI::DOUBLE,to,4);
	if(DATA->turn_on_rhob) MPI::COMM_WORLD.Send(packagerhob,sizeOfData,MPI::DOUBLE,to,5);
	if(DATA->turn_on_bulk) MPI::COMM_WORLD.Send(packagePi_b,sizeOfData,MPI::DOUBLE,to,6);
	
	//cerr << "On rank " << rank << ", sizeof(packagedT) = " << sizeof(packagedT[0]) << endl;
	//cerr << "On rank " << rank << ", sizeof(packagedWxx) = " << sizeof(packageWxx[0]) << endl;
	
	if(DATA->turn_on_shear) 
	  {
	    //      MPI::COMM_WORLD.Send(packageWtautau,sizeOfData,MPI::DOUBLE,to,12);
	    //      MPI::COMM_WORLD.Send(packageWtaux,sizeOfData,MPI::DOUBLE,to,13);
	    //      MPI::COMM_WORLD.Send(packageWtauy,sizeOfData,MPI::DOUBLE,to,14);
	    //      MPI::COMM_WORLD.Send(packageWtaueta,sizeOfData,MPI::DOUBLE,to,15);
	    MPI::COMM_WORLD.Send(packageWxx,sizeOfData,MPI::DOUBLE,to,16);
	    MPI::COMM_WORLD.Send(packageWxy,sizeOfData,MPI::DOUBLE,to,17);
	    MPI::COMM_WORLD.Send(packageWxeta,sizeOfData,MPI::DOUBLE,to,18);
	    MPI::COMM_WORLD.Send(packageWyy,sizeOfData,MPI::DOUBLE,to,19);
	    MPI::COMM_WORLD.Send(packageWyeta,sizeOfData,MPI::DOUBLE,to,20);
	    //cerr << "Finished sending on rank " << rank << endl;
	    //     MPI::COMM_WORLD.Send(packageWetaeta,sizeOfData,MPI::DOUBLE,to,21);
	    if(DATA->fluctuatingHydroFlag == 2){
	      //cerr << "Send loop entered on rank " << rank << endl;
	      MPI::COMM_WORLD.Send(packagedT,sizeOfData,MPI::DOUBLE,to,21);
	      MPI::COMM_WORLD.Send(packagedutau,sizeOfData,MPI::DOUBLE,to,22);
	      MPI::COMM_WORLD.Send(packagedux, sizeOfData, MPI::DOUBLE, to, 23);
	      MPI::COMM_WORLD.Send(packageduy, sizeOfData, MPI::DOUBLE, to, 24);
	      MPI::COMM_WORLD.Send(packagedueta, sizeOfData, MPI::DOUBLE, to, 25);
	      MPI::COMM_WORLD.Send(packagede, sizeOfData, MPI::DOUBLE, to, 26);
	      MPI::COMM_WORLD.Send(packagedp, sizeOfData, MPI::DOUBLE, to, 27);
	      MPI::COMM_WORLD.Send(packagedWtautau, sizeOfData, MPI::DOUBLE, to, 28);
	      MPI::COMM_WORLD.Send(packagedWtaux, sizeOfData, MPI::DOUBLE, to, 29);
	      MPI::COMM_WORLD.Send(packagedWtauy, sizeOfData, MPI::DOUBLE, to, 30);
	      MPI::COMM_WORLD.Send(packagedWtaueta, sizeOfData, MPI::DOUBLE, to, 31);
	      MPI::COMM_WORLD.Send(packagedWxx, sizeOfData, MPI::DOUBLE, to, 32);
	      MPI::COMM_WORLD.Send(packagedWxy, sizeOfData, MPI::DOUBLE, to, 33);
	      MPI::COMM_WORLD.Send(packagedWxeta, sizeOfData, MPI::DOUBLE, to, 34);
	      MPI::COMM_WORLD.Send(packagedWyy, sizeOfData, MPI::DOUBLE, to, 35);
	      MPI::COMM_WORLD.Send(packagedWyeta, sizeOfData, MPI::DOUBLE, to, 36);
	      MPI::COMM_WORLD.Send(packagedWetaeta, sizeOfData, MPI::DOUBLE, to, 37);
	      //cerr << "Send loop finished!" << endl;
	    }
	  }
	//cerr << " done sending to the left on rank " << rank << endl;
      }
    //MPI::COMM_WORLD.Barrier();
    // receiving and unwrapping the package
    if ( rank != size-1 )
      {  
	MPI::COMM_WORLD.Recv(package,sizeOfData,MPI::DOUBLE,from,1);
	//      MPI::COMM_WORLD.Recv(packageutau,sizeOfData,MPI::DOUBLE,from,3);
	MPI::COMM_WORLD.Recv(packageux,sizeOfData,MPI::DOUBLE,from,2);
	MPI::COMM_WORLD.Recv(packageuy,sizeOfData,MPI::DOUBLE,from,3);
	MPI::COMM_WORLD.Recv(packageueta,sizeOfData,MPI::DOUBLE,from,4);
	if(DATA->turn_on_rhob) MPI::COMM_WORLD.Recv(packagerhob,sizeOfData,MPI::DOUBLE,from,5);
	if(DATA->turn_on_bulk) MPI::COMM_WORLD.Recv(packagePi_b,sizeOfData,MPI::DOUBLE,from,6);
	if(DATA->turn_on_shear) 
	  {
	    //      MPI::COMM_WORLD.Recv(packageWtautau,sizeOfData,MPI::DOUBLE,from,12);
	    //      MPI::COMM_WORLD.Recv(packageWtaux,sizeOfData,MPI::DOUBLE,from,13);
	    //      MPI::COMM_WORLD.Recv(packageWtauy,sizeOfData,MPI::DOUBLE,from,14);
	    //      MPI::COMM_WORLD.Recv(packageWtaueta,sizeOfData,MPI::DOUBLE,from,15);
	    MPI::COMM_WORLD.Recv(packageWxx,sizeOfData,MPI::DOUBLE,from,16);
	    MPI::COMM_WORLD.Recv(packageWxy,sizeOfData,MPI::DOUBLE,from,17);
	    MPI::COMM_WORLD.Recv(packageWxeta,sizeOfData,MPI::DOUBLE,from,18);
	    MPI::COMM_WORLD.Recv(packageWyy,sizeOfData,MPI::DOUBLE,from,19);
	    MPI::COMM_WORLD.Recv(packageWyeta,sizeOfData,MPI::DOUBLE,from,20);
	    //MPI::COMM_WORLD.Recv(packagedutau, sizeOfData, MPI::DOUBLE,from,22);
	    //      MPI::COMM_WORLD.Recv(packageWetaeta,sizeOfData,MPI::DOUBLE,from,21);
	    if(DATA->fluctuatingHydroFlag == 2){
	      //cerr << "Recv loop entered on rank " << rank << endl;
	      MPI::COMM_WORLD.Recv(packagedT, sizeOfData, MPI::DOUBLE, from, 21);
	      MPI::COMM_WORLD.Recv(packagedutau, sizeOfData, MPI::DOUBLE, from, 22);
	      MPI::COMM_WORLD.Recv(packagedux, sizeOfData, MPI::DOUBLE, from, 23);
	      MPI::COMM_WORLD.Recv(packageduy, sizeOfData, MPI::DOUBLE, from, 24);
	      MPI::COMM_WORLD.Recv(packagedueta, sizeOfData, MPI::DOUBLE, from, 25);
	      MPI::COMM_WORLD.Recv(packagede, sizeOfData, MPI::DOUBLE, from, 26);
	      MPI::COMM_WORLD.Recv(packagedp, sizeOfData, MPI::DOUBLE, from, 27);
	      MPI::COMM_WORLD.Recv(packagedWtautau, sizeOfData, MPI::DOUBLE, from, 28);
	      MPI::COMM_WORLD.Recv(packagedWtaux, sizeOfData, MPI::DOUBLE, from, 29);
	      MPI::COMM_WORLD.Recv(packagedWtauy, sizeOfData, MPI::DOUBLE, from, 30);
	      MPI::COMM_WORLD.Recv(packagedWtaueta, sizeOfData, MPI::DOUBLE, from, 31);
	      MPI::COMM_WORLD.Recv(packagedWxx, sizeOfData, MPI::DOUBLE, from, 32);
	      MPI::COMM_WORLD.Recv(packagedWxy, sizeOfData, MPI::DOUBLE, from, 33);
	      MPI::COMM_WORLD.Recv(packagedWxeta, sizeOfData, MPI::DOUBLE, from, 34);
	      MPI::COMM_WORLD.Recv(packagedWyy, sizeOfData, MPI::DOUBLE, from, 35);
	      MPI::COMM_WORLD.Recv(packagedWyeta, sizeOfData, MPI::DOUBLE, from, 36);
	      MPI::COMM_WORLD.Recv(packagedWetaeta, sizeOfData, MPI::DOUBLE, from, 37);
	      //cerr << "Recv loop finished!" << endl;
	    }
	  }
	//cout << " receiving from the right on rank " << rank << endl;
	for(ix=0; ix<=nx; ix++)
	  {
	    for(iy=0; iy<=ny; iy++)
	      {
		position = ix + ((nx+1)*iy);
		//if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
		//cout << "Rneighbor[ix][iy][0].TJb[i][alpha][0]=" << package[position] << endl;
		Rneighbor_eps[ix][iy] = package[position];
		//      Rneighbor_utau[ix][iy] = packageutau[position];
		Rneighbor_ux[ix][iy] = packageux[position];
		Rneighbor_uy[ix][iy] = packageuy[position];
		Rneighbor_ueta[ix][iy] = packageueta[position];
		if(DATA->turn_on_rhob) Rneighbor_rhob[ix][iy] = packagerhob[position];
		else Rneighbor_rhob[ix][iy] = 0.;
		if(DATA->turn_on_rhob) Rneighbor_Pi_b[ix][iy] = packagePi_b[position];
		else Rneighbor_Pi_b[ix][iy] = 0.;
		if(DATA->turn_on_shear) 
		  {
		    //      Rneighbor_Wtautau[ix][iy] = packageWtautau[position];
		    //      Rneighbor_Wtaux[ix][iy] = packageWtaux[position];
		    //      Rneighbor_Wtauy[ix][iy] = packageWtauy[position];
		    //      Rneighbor_Wtaueta[ix][iy] = packageWtaueta[position];
		    Rneighbor_Wxx[ix][iy] = packageWxx[position];
		    Rneighbor_Wxy[ix][iy] = packageWxy[position];
		    Rneighbor_Wxeta[ix][iy] = packageWxeta[position];
		    Rneighbor_Wyy[ix][iy] = packageWyy[position];
		    Rneighbor_Wyeta[ix][iy] = packageWyeta[position];
		    //      Rneighbor_Wetaeta[ix][iy] = packageWetaeta[position];
		    if(DATA->fluctuatingHydroFlag == 2){
		      Rneighbor_dT[ix][iy] = packagedT[position];
		      Rneighbor_dutau[ix][iy] = packagedutau[position];
		      Rneighbor_dux[ix][iy] = packagedux[position];
		      Rneighbor_duy[ix][iy] = packageduy[position];
		      Rneighbor_dueta[ix][iy] = packagedueta[position];
		      Rneighbor_de[ix][iy] = packagede[position];
		      Rneighbor_dp[ix][iy] = packagedp[position];
		      Rneighbor_dWtautau[ix][iy] = packagedWtautau[position];
		      Rneighbor_dWtaux[ix][iy] = packagedWtaux[position];
		      Rneighbor_dWtauy[ix][iy] = packagedWtauy[position];
		      Rneighbor_dWtaueta[ix][iy] = packagedWtaueta[position];
		      Rneighbor_dWxx[ix][iy] = packagedWxx[position];
		      Rneighbor_dWxy[ix][iy] = packagedWxy[position];
		      Rneighbor_dWxeta[ix][iy] = packagedWxeta[position];
		      Rneighbor_dWyy[ix][iy] = packagedWyy[position];
		      Rneighbor_dWyeta[ix][iy] = packagedWyeta[position];
		      Rneighbor_dWetaeta[ix][iy] = packagedWetaeta[position];
		    }
		  }
		else
		  {
		    Rneighbor_Wxx[ix][iy] = 0.;
		    Rneighbor_Wxy[ix][iy] = 0.;
		    Rneighbor_Wxeta[ix][iy] = 0.;
		    Rneighbor_Wyy[ix][iy] = 0.;
		    Rneighbor_Wyeta[ix][iy] = 0.;
		  }
	      }
	  }
	//cerr << " done receiving from the right on rank " << rank << endl;
      }
    
    //   if (rank == size-1) maxEta = neta-fac;
    //   else maxEta = neta;
  
    //*******************************************
    //  end MPI code copied from FindFreezeOutSurface2
  }
  //   int maxEta;
  maxEta = neta-fac;
  
  FO = DATA->epsilonFreeze/hbarc;

  //cerr << "FO = " << FO << endl;

  //if (DATA->useEpsFO) FO = DATA->epsilonFreeze/hbarc;
  //else FO = DATA->TFO/hbarc;

  for(ix=0; ix<=nx; ix+=fac)
    {
      x = ix*(DATA->delta_x) - (DATA->x_size/2.0); 
      for(iy=0; iy<=ny; iy+=fac)
	{
	  y = iy*(DATA->delta_y) - (DATA->y_size/2.0);
	  for(ieta=0; ieta<=maxEta; ieta+=fac)
	    {
	      eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
	    
	      if (arena[ix][iy][ieta].epsilon >= FO) frozen = 0;
	      //if (DATA->useEpsFO && (arena[ix][iy][ieta].epsilon >= FO)) frozen = 0;
	      //else if (!DATA->useEpsFO && (eos->get_temperature(arena[ix][iy][ieta].epsilon,arena[ix][iy][ieta].rhob) >= FO)) frozen = 0;
	            
	      //         rhob = arena[ix][iy][ieta].rhob;
	      //         if (!DATA->useEpsFO && (1 == DATA->turn_on_rhob)) 
	      //   cerr << "Caution! Constant temperature freezeout with nonzero chemical potential doesn't yet work properly (need to implement rhob_prev)\n";

	      if(ix<nx)
		{
		  // check forward in x
		  xf = x + DX/2.;
		  yf = y;
		  etaf = eta;// no longer used
		  tauf = tau;

		  nix = ix+fac;
		  niy = iy;
		  nieta = ieta;
		  FULLSU[0]=0.;
		  FULLSU[1]=DTAU*DY*DETA/subsections;
		  FULLSU[2]=0.;
		  FULLSU[3]=0.;
		        
		  //if (DATA->useEpsFO)
		  //{
		  home = arena[ix][iy][ieta].epsilon;
		  neighbor = arena[nix][niy][nieta].epsilon;
		  //}
		  //else
		  //{
		  //  //   home = arena[ix][iy][ieta].T;
		  //  //   neighbor = arena[nix][niy][nieta].T;
		  //  home = eos->get_temperature(arena[ix][iy][ieta].epsilon,arena[ix][iy][ieta].rhob);
		  //  neighbor = eos->get_temperature(arena[nix][niy][nieta].epsilon,arena[nix][niy][nieta].rhob);
		  //}
		  if(((home > FO) && (neighbor <= FO)) || ((home <= FO) && (neighbor > FO)))
		    {
		      //cerr << "Is the 1st condition ever satisfied?" << endl;
		      if (home<=FO)
			SIG=-1.;
		      else SIG=1.;
		      for (int orient = 0; orient < 4; orient++) FULLSU[orient]*=SIG;
		      // simple linear interpolation of all independent variables
		      ux = 0.5*(arena[ix][iy][ieta].u[0][1] + arena[nix][niy][nieta].u[0][1]);
		      uy = 0.5*(arena[ix][iy][ieta].u[0][2] + arena[nix][niy][nieta].u[0][2]);
		      ueta = 0.5*(arena[ix][iy][ieta].u[0][3] + arena[nix][niy][nieta].u[0][3]);
		      utau = sqrt(1 + ux*ux + uy*uy + ueta*ueta);
		      iepsFO = 0.5*(arena[ix][iy][ieta].epsilon + arena[nix][niy][nieta].epsilon);
		      rhob = 0.5*(arena[ix][iy][ieta].rhob + arena[nix][niy][nieta].rhob);
		      double iP = eos->p_func(iepsFO, rhob);
		      double sFO = eos->s_func(iepsFO, iP, rhob, DATA);
		      //pi_b = 0.5*(arena[ix][iy][ieta].pi_b[0] + arena[nix][niy][nieta].pi_b[0]);
		      //TFO = eos->get_temperature(iepsFO, rhob);
		      //muB = eos->get_mu(iepsFO, rhob);
		      TFO = eos->interpolate2(iepsFO, 0., 1);
		      muB = 0.;
		      //eps_plus_p_over_T_FO=(iepsFO+eos->get_pressure(iepsFO, rhob))/TFO;
		      Wxx = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][1] + arena[nix][niy][nieta].Wmunu[0][1][1]);
		      Wxy = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][2] + arena[nix][niy][nieta].Wmunu[0][1][2]);
		      Wxeta = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][3] + arena[nix][niy][nieta].Wmunu[0][1][3]);
		      Wyy = 0.5*(arena[ix][iy][ieta].Wmunu[0][2][2] + arena[nix][niy][nieta].Wmunu[0][2][2]);
		      Wyeta = 0.5*(arena[ix][iy][ieta].Wmunu[0][2][3] + arena[nix][niy][nieta].Wmunu[0][2][3]);

		      Wetaeta = (2.*(ux*uy*Wxy 
                            + ux*ueta*Wxeta
				     + uy*ueta*Wyeta )
				 -( utau*utau - ux*ux )*Wxx
				 -( utau*utau - uy*uy )*Wyy
				 )/( utau*utau - ueta*ueta ) ;
		      Wtaux = (ux*Wxx + uy*Wxy + ueta*Wxeta)/utau;
		      Wtauy = (ux*Wxy + uy*Wyy + ueta*Wyeta)/utau;
		      Wtaueta = (ux*Wxeta + uy*Wyeta + ueta*Wetaeta)/utau;
		      Wtautau = (ux*Wtaux + uy*Wtauy + ueta*Wtaueta)/utau;
		        
		      if(DATA->fluctuatingHydroFlag == 2){
			//if(arena[ix][iy][ieta].deltaP[0] == 0.){dT = 0.;}
			//else{dT = 0.5*(arena[ix][iy][ieta].deltaP[0]
			//	       /eos->s_func(arena[ix][iy][ieta].epsilon, 
			//			    eos->p_func(arena[ix][iy][ieta].epsilon, 
			//					arena[ix][iy][ieta].rhob), 
			//			    arena[ix][iy][ieta].rhob, DATA));
			//}
			//if(arena[nix][niy][nieta].deltaP[0] != 0.){
			//dT += 0.5*(arena[nix][niy][nieta].deltaP[0]
			//	     /eos->s_func(arena[nix][niy][nieta].epsilon, 
			//			  eos->p_func(arena[nix][niy][nieta].epsilon,
			//				      arena[nix][niy][nieta].rhob), 
			//			  arena[nix][niy][nieta].rhob, DATA));}
			//dutau = 0.5*(arena[ix][iy][ieta].deltaU[0][0] + arena[nix][niy][nieta].deltaU[0][0]);
			dux = 0.5*(arena[ix][iy][ieta].deltaU[0][1] + arena[nix][niy][nieta].deltaU[0][1]);
			duy = 0.5*(arena[ix][iy][ieta].deltaU[0][2] + arena[nix][niy][nieta].deltaU[0][2]);
			dueta = 0.5*(arena[ix][iy][ieta].deltaU[0][3] + arena[nix][niy][nieta].deltaU[0][3]);
			dutau = (dux*ux + duy*uy + dueta*ueta)/utau;
			//de = 0.5*(arena[ix][iy][ieta].deltaP[0]/eos->get_dpOverde2(fmax(0.03, arena[ix][iy][ieta].epsilon), 
			//							   arena[ix][iy][ieta].rhob) 
			//	  + arena[nix][niy][nieta].deltaP[0]/eos->get_dpOverde2(fmax(0.03, arena[nix][niy][nieta].epsilon),
			//								arena[nix][niy][nieta].rhob) );
			dp = 0.5*(arena[ix][iy][ieta].deltaP[0] + arena[nix][niy][nieta].deltaP[0]);
			dT = dp/eos->s_func(iepsFO, eos->p_func(iepsFO, rhob), rhob, DATA);
			de = dp/eos->get_dpOverde2(iepsFO, rhob);
			//dWtautau = 0.5*(arena[ix][iy][ieta].deltaW[0][0][0] + arena[nix][niy][nieta].deltaW[0][0][0]);
			//dWtaux = 0.5*(arena[ix][iy][ieta].deltaW[0][0][1] + arena[nix][niy][nieta].deltaW[0][0][1]);
			//dWtauy = 0.5*(arena[ix][iy][ieta].deltaW[0][0][2] + arena[nix][niy][nieta].deltaW[0][0][2]);
			//dWtaueta = 0.5*(arena[ix][iy][ieta].deltaW[0][0][3] + arena[nix][niy][nieta].deltaW[0][0][3]);
			dWxx = 0.5*(arena[ix][iy][ieta].deltaW[0][1][1] + arena[nix][niy][nieta].deltaW[0][1][1]);
			dWxy = 0.5*(arena[ix][iy][ieta].deltaW[0][1][2] + arena[nix][niy][nieta].deltaW[0][1][2]);
			dWxeta = 0.5*(arena[ix][iy][ieta].deltaW[0][1][3] + arena[nix][niy][nieta].deltaW[0][1][3]);
			dWyy = 0.5*(arena[ix][iy][ieta].deltaW[0][2][2] + arena[nix][niy][nieta].deltaW[0][2][2]);
			dWyeta = 0.5*(arena[ix][iy][ieta].deltaW[0][2][3] + arena[nix][niy][nieta].deltaW[0][2][3]);
			dWetaeta = 0.5*(arena[ix][iy][ieta].deltaW[0][3][3] + arena[nix][niy][nieta].deltaW[0][3][3]);
			dWtaux = (ux*dWxx + uy*dWxy + ueta*dWxeta -dutau*Wtaux + dux*Wxx + duy*Wxy + dueta*Wxeta)/utau;
			dWtauy = (ux*dWxy + uy*dWyy + ueta*dWyeta -dutau*Wtauy + dux*Wxy + duy*Wyy + dueta*Wyeta)/utau;
			dWtaueta = (ux*dWxeta + uy*dWyeta + ueta*dWetaeta -dutau*Wtaueta + dux*Wxeta + duy*Wyeta + dueta*Wetaeta)/utau;
			dWtautau = (ux*dWtaux + uy*dWtauy + ueta*dWtaueta -dutau*Wtautau + dux*Wtaux + duy*Wtauy + dueta*Wtaueta)/utau;
		      }

		      for (int i = 0; i < subsections; i++)
			{
			  etaf = eta - DETA/2 + DETA/2/subsections + (i*DETA)/(subsections);
			  s_file << std::setprecision(10) << tauf << " " << xf << " " << yf << " " << etaf << " " 
				 << FULLSU[0] << " " <<FULLSU[1] << " " <<FULLSU[2] << " " <<FULLSU[3] 
				 << " " <<  utau << " " << ux << " " << uy << " " << ueta << " " 
				 << iepsFO << " " << TFO << " " << muB << " " << sFO << " " 
				 << Wtautau << " " << Wtaux << " " << Wtauy << " " << Wtaueta << " " 
				 << Wxx << " " << Wxy << " " << Wxeta << " " << Wyy << " " << Wyeta << " " << Wetaeta; 
			  if(DATA->turn_on_bulk) s_file << " " << pi_b;
			  s_file << endl;  
			  if(DATA->fluctuatingHydroFlag == 2){
			    ds_file << std::setprecision(10) << dT << " " << dutau << " " << dux << " " << duy << " " << dueta << " "
				    << de << " " << dp << " "
				    << dWtautau << " " << dWtaux << " " << dWtauy << " " << dWtaueta << " "
				    << dWxx << " " << dWxy << " " << dWxeta << " "
				    << dWyy << " " << dWyeta << " "
				    << dWetaeta << endl;
			  }
			}
		    } // if grid pair straddles freeze out density
		} //if(ix<nx)
	      
	      if(iy<ny)
		{
		  // check forward in y
		  xf = x;
		  yf = y + DY/2.;
		  etaf = eta; // no longer used
		  tauf = tau;

		  nix = ix;
		  niy = iy+fac;
		  nieta = ieta;
		  FULLSU[0]=0.;
		  FULLSU[1]=0.;
		  FULLSU[2]=DTAU*DX*DETA/subsections;
		  FULLSU[3]=0.;
		  
		  //if (DATA->useEpsFO)
		  //{
		  home = arena[ix][iy][ieta].epsilon;
		  neighbor = arena[nix][niy][nieta].epsilon;
		  //}
		  //else
		  //{
		  //  //   home = arena[ix][iy][ieta].T;
		  //  //   neighbor = arena[nix][niy][nieta].T;
		  //  home = eos->get_temperature(arena[ix][iy][ieta].epsilon,arena[ix][iy][ieta].rhob);
		  //  neighbor = eos->get_temperature(arena[nix][niy][nieta].epsilon,arena[nix][niy][nieta].rhob);
		  //}
		  if(((home > FO) && (neighbor <= FO)) || ((home <= FO) && (neighbor > FO)))
		    {
		      //cerr << "Is the 2nd condition ever satisfied?" << endl;
		      if (home<=FO)
			SIG=-1.;
		      else SIG=1.;
		      for (int orient = 0; orient < 4; orient++) FULLSU[orient]*=SIG;
		      // simple linear interpolation of all independent variables
		      ux = 0.5*(arena[ix][iy][ieta].u[0][1] + arena[nix][niy][nieta].u[0][1]);
		      uy = 0.5*(arena[ix][iy][ieta].u[0][2] + arena[nix][niy][nieta].u[0][2]);
		      ueta = 0.5*(arena[ix][iy][ieta].u[0][3] + arena[nix][niy][nieta].u[0][3]);
		      utau = sqrt(1 + ux*ux + uy*uy + ueta*ueta);
		      iepsFO = 0.5*(arena[ix][iy][ieta].epsilon + arena[nix][niy][nieta].epsilon);
		      rhob = 0.5*(arena[ix][iy][ieta].rhob + arena[nix][niy][nieta].rhob);
		      double iP = eos->p_func(iepsFO, rhob);
		      double sFO = eos->s_func(iepsFO, iP, rhob, DATA);
		      //pi_b = 0.5*(arena[ix][iy][ieta].pi_b[0] + arena[nix][niy][nieta].pi_b[0]);
		      //TFO = eos->get_temperature(iepsFO, rhob);
		      //muB = eos->get_mu(iepsFO, rhob);
		      TFO = eos->interpolate2(iepsFO, 0., 1);
		      muB = 0.;
		      //eps_plus_p_over_T_FO=(iepsFO+eos->get_pressure(iepsFO, rhob))/TFO;
		      Wxx = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][1] + arena[nix][niy][nieta].Wmunu[0][1][1]);
		      Wxy = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][2] + arena[nix][niy][nieta].Wmunu[0][1][2]);
		      Wxeta = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][3] + arena[nix][niy][nieta].Wmunu[0][1][3]);
		      Wyy = 0.5*(arena[ix][iy][ieta].Wmunu[0][2][2] + arena[nix][niy][nieta].Wmunu[0][2][2]);
		      Wyeta = 0.5*(arena[ix][iy][ieta].Wmunu[0][2][3] + arena[nix][niy][nieta].Wmunu[0][2][3]);

		      Wetaeta = (2.*(ux*uy*Wxy 
                            + ux*ueta*Wxeta
				     + uy*ueta*Wyeta )
				 -( utau*utau - ux*ux )*Wxx
				 -( utau*utau - uy*uy )*Wyy
				 )/( utau*utau - ueta*ueta ) ;
		      Wtaux = (ux*Wxx + uy*Wxy + ueta*Wxeta)/utau;
		      Wtauy = (ux*Wxy + uy*Wyy + ueta*Wyeta)/utau;
		      Wtaueta = (ux*Wxeta + uy*Wyeta + ueta*Wetaeta)/utau;
		      Wtautau = (ux*Wtaux + uy*Wtauy + ueta*Wtaueta)/utau;
		        
		      if(DATA->fluctuatingHydroFlag == 2){
			//if(arena[ix][iy][ieta].deltaP[0] == 0.){dT = 0.;}
			//else{dT = 0.5*(arena[ix][iy][ieta].deltaP[0]
			//	       /eos->s_func(arena[ix][iy][ieta].epsilon, 
			//			    eos->p_func(arena[ix][iy][ieta].epsilon, 
			//					arena[ix][iy][ieta].rhob), 
			//			    arena[ix][iy][ieta].rhob, DATA));}
			//if(arena[nix][niy][nieta].deltaP[0] != 0.){
			//dT += 0.5*(arena[nix][niy][nieta].deltaP[0]
			//	     /eos->s_func(arena[nix][niy][nieta].epsilon, 
			//			  eos->p_func(arena[nix][niy][nieta].epsilon, 
			//				      arena[nix][niy][nieta].rhob), 
			//			  arena[nix][niy][nieta].rhob, DATA));}
			//dutau = 0.5*(arena[ix][iy][ieta].deltaU[0][0] + arena[nix][niy][nieta].deltaU[0][0]);
			dux = 0.5*(arena[ix][iy][ieta].deltaU[0][1] + arena[nix][niy][nieta].deltaU[0][1]);
			duy = 0.5*(arena[ix][iy][ieta].deltaU[0][2] + arena[nix][niy][nieta].deltaU[0][2]);
			dueta = 0.5*(arena[ix][iy][ieta].deltaU[0][3] + arena[nix][niy][nieta].deltaU[0][3]);
			dutau = (dux*ux + duy*uy + dueta*ueta)/utau;
			//de = 0.5*(arena[ix][iy][ieta].deltaP[0]/eos->get_dpOverde2(fmax(0.03, arena[ix][iy][ieta].epsilon), 
			//							   arena[ix][iy][ieta].rhob) 
			//	  + arena[nix][niy][nieta].deltaP[0]/eos->get_dpOverde2(fmax(0.03, arena[nix][niy][nieta].epsilon),
			//								arena[nix][niy][nieta].rhob) );
			dp = 0.5*(arena[ix][iy][ieta].deltaP[0] + arena[nix][niy][nieta].deltaP[0]);
			dT = dp/eos->s_func(iepsFO, eos->p_func(iepsFO, rhob), rhob, DATA);
			de = dp/eos->get_dpOverde2(iepsFO, rhob);
			//dWtautau = 0.5*(arena[ix][iy][ieta].deltaW[0][0][0] + arena[nix][niy][nieta].deltaW[0][0][0]);
			//dWtaux = 0.5*(arena[ix][iy][ieta].deltaW[0][0][1] + arena[nix][niy][nieta].deltaW[0][0][1]);
			//dWtauy = 0.5*(arena[ix][iy][ieta].deltaW[0][0][2] + arena[nix][niy][nieta].deltaW[0][0][2]);
			//dWtaueta = 0.5*(arena[ix][iy][ieta].deltaW[0][0][3] + arena[nix][niy][nieta].deltaW[0][0][3]);
			dWxx = 0.5*(arena[ix][iy][ieta].deltaW[0][1][1] + arena[nix][niy][nieta].deltaW[0][1][1]);
			dWxy = 0.5*(arena[ix][iy][ieta].deltaW[0][1][2] + arena[nix][niy][nieta].deltaW[0][1][2]);
			dWxeta = 0.5*(arena[ix][iy][ieta].deltaW[0][1][3] + arena[nix][niy][nieta].deltaW[0][1][3]);
			dWyy = 0.5*(arena[ix][iy][ieta].deltaW[0][2][2] + arena[nix][niy][nieta].deltaW[0][2][2]);
			dWyeta = 0.5*(arena[ix][iy][ieta].deltaW[0][2][3] + arena[nix][niy][nieta].deltaW[0][2][3]);
			dWetaeta = 0.5*(arena[ix][iy][ieta].deltaW[0][3][3] + arena[nix][niy][nieta].deltaW[0][3][3]);
			dWtaux = (ux*dWxx + uy*dWxy + ueta*dWxeta -dutau*Wtaux + dux*Wxx + duy*Wxy + dueta*Wxeta)/utau;
			dWtauy = (ux*dWxy + uy*dWyy + ueta*dWyeta -dutau*Wtauy + dux*Wxy + duy*Wyy + dueta*Wyeta)/utau;
			dWtaueta = (ux*dWxeta + uy*dWyeta + ueta*dWetaeta -dutau*Wtaueta + dux*Wxeta + duy*Wyeta + dueta*Wetaeta)/utau;
			dWtautau = (ux*dWtaux + uy*dWtauy + ueta*dWtaueta -dutau*Wtautau + dux*Wtaux + duy*Wtauy + dueta*Wtaueta)/utau;
		      }

		      for (int i = 0; i < subsections; i++)
			{
			  etaf = eta - DETA/2 + DETA/2/subsections + (i*DETA)/(subsections);
			      
			  s_file << std::setprecision(10) << tauf << " " << xf << " " << yf << " " << etaf << " " 
				 << FULLSU[0] << " " <<FULLSU[1] << " " <<FULLSU[2] << " " <<FULLSU[3] 
				 << " " <<  utau << " " << ux << " " << uy << " " << ueta << " " 
				 << iepsFO << " " << TFO << " " << muB << " " << sFO << " " 
				 << Wtautau << " " << Wtaux << " " << Wtauy << " " << Wtaueta << " " 
				 << Wxx << " " << Wxy << " " << Wxeta << " " << Wyy << " " << Wyeta << " " << Wetaeta; 
			  if(DATA->turn_on_bulk) s_file << " " << pi_b;
			  s_file << endl; 
			  if(DATA->fluctuatingHydroFlag == 2){
			    ds_file << std::setprecision(10) << dT << " " << dutau << " " << dux << " " << duy << " " << dueta << " "
				    << de << " " << dp << " "
				    << dWtautau << " " << dWtaux << " " << dWtauy << " " << dWtaueta << " "
				    << dWxx << " " << dWxy << " " << dWxeta << " "
				    << dWyy << " " << dWyeta << " "
				    << dWetaeta << endl;
			  }
			}
		    } // if grid pair straddles freeze out density
		} //if(iy<ny)

	      if(ieta<(maxEta))
		{
		  // check forward in eta
		  xf = x;
		  yf = y;
		  etaf = eta + DETA/2;
		  tauf = tau;

		  nix = ix;
		  niy = iy;
		  nieta = ieta + fac;
		  FULLSU[0]=0.;
		  FULLSU[1]=0.;
		  FULLSU[2]=0.;
		  FULLSU[3]=DTAU*DX*DY;
		        
		  //if (DATA->useEpsFO)
		  //{
		  home = arena[ix][iy][ieta].epsilon;
		  neighbor = arena[nix][niy][nieta].epsilon;
		  //}
		  //else
		  //{
		  //  //   home = arena[ix][iy][ieta].T;
		  //  //   neighbor = arena[nix][niy][nieta].T;
		  //  home = eos->get_temperature(arena[ix][iy][ieta].epsilon,arena[ix][iy][ieta].rhob);
		  //  neighbor = eos->get_temperature(arena[nix][niy][nieta].epsilon,arena[nix][niy][nieta].rhob);
		  //}
		  if(((home > FO) && (neighbor <= FO)) || ((home <= FO) && (neighbor > FO)))
		    {
		      //cerr << "Is the 3rd condition ever satisfied?" << endl;
		      if (home<=FO)
			SIG=-1.;
		      else SIG=1.;
		      for (int orient = 0; orient < 4; orient++) FULLSU[orient]*=SIG;
		      // simple linear interpolation of all independent variables
		      ux = 0.5*(arena[ix][iy][ieta].u[0][1] + arena[nix][niy][nieta].u[0][1]);
		      uy = 0.5*(arena[ix][iy][ieta].u[0][2] + arena[nix][niy][nieta].u[0][2]);
		      ueta = 0.5*(arena[ix][iy][ieta].u[0][3] + arena[nix][niy][nieta].u[0][3]);
		      utau = sqrt(1 + ux*ux + uy*uy + ueta*ueta);
		      iepsFO = 0.5*(arena[ix][iy][ieta].epsilon + arena[nix][niy][nieta].epsilon);
		      rhob = 0.5*(arena[ix][iy][ieta].rhob + arena[nix][niy][nieta].rhob);
		      double iP = eos->p_func(iepsFO, rhob);
		      double sFO = eos->s_func(iepsFO, iP, rhob, DATA);
		      //pi_b = 0.5*(arena[ix][iy][ieta].pi_b[0] + arena[nix][niy][nieta].pi_b[0]);
		      //TFO = eos->get_temperature(iepsFO, rhob);
		      //muB = eos->get_mu(iepsFO, rhob);
		      TFO = eos->interpolate2(iepsFO, 0., 1);
		      muB = 0.;
		      //eps_plus_p_over_T_FO=(iepsFO+eos->get_pressure(iepsFO, rhob))/TFO;
		      Wxx = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][1] + arena[nix][niy][nieta].Wmunu[0][1][1]);
		      Wxy = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][2] + arena[nix][niy][nieta].Wmunu[0][1][2]);
		      Wxeta = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][3] + arena[nix][niy][nieta].Wmunu[0][1][3]);
		      Wyy = 0.5*(arena[ix][iy][ieta].Wmunu[0][2][2] + arena[nix][niy][nieta].Wmunu[0][2][2]);
		      Wyeta = 0.5*(arena[ix][iy][ieta].Wmunu[0][2][3] + arena[nix][niy][nieta].Wmunu[0][2][3]);

		      Wetaeta = (2.*(ux*uy*Wxy 
                            + ux*ueta*Wxeta
				     + uy*ueta*Wyeta )
				 -( utau*utau - ux*ux )*Wxx
				 -( utau*utau - uy*uy )*Wyy
				 )/( utau*utau - ueta*ueta ) ;
		      Wtaux = (ux*Wxx + uy*Wxy + ueta*Wxeta)/utau;
		      Wtauy = (ux*Wxy + uy*Wyy + ueta*Wyeta)/utau;
		      Wtaueta = (ux*Wxeta + uy*Wyeta + ueta*Wetaeta)/utau;
		      Wtautau = (ux*Wtaux + uy*Wtauy + ueta*Wtaueta)/utau;
		        
		      if(DATA->fluctuatingHydroFlag == 2){
			//if(arena[ix][iy][ieta].deltaP[0] == 0.){dT = 0.;}
			//else{dT = 0.5*(arena[ix][iy][ieta].deltaP[0]
			//	       /eos->s_func(arena[ix][iy][ieta].epsilon, 
			//			    eos->p_func(arena[ix][iy][ieta].epsilon, 
			//					arena[ix][iy][ieta].rhob), 
			//			    arena[ix][iy][ieta].rhob, DATA));}
			//if(arena[nix][niy][nieta].deltaP[0] != 0.){
			//dT += 0.5*(arena[nix][niy][nieta].deltaP[0]
			//	     /eos->s_func(arena[nix][niy][nieta].epsilon, 
			//			  eos->p_func(arena[nix][niy][nieta].epsilon, 
			//				      arena[nix][niy][nieta].rhob), 
			//			  arena[nix][niy][nieta].rhob, DATA));}
			//dutau = 0.5*(arena[ix][iy][ieta].deltaU[0][0] + arena[nix][niy][nieta].deltaU[0][0]);
			dux = 0.5*(arena[ix][iy][ieta].deltaU[0][1] + arena[nix][niy][nieta].deltaU[0][1]);
			duy = 0.5*(arena[ix][iy][ieta].deltaU[0][2] + arena[nix][niy][nieta].deltaU[0][2]);
			dueta = 0.5*(arena[ix][iy][ieta].deltaU[0][3] + arena[nix][niy][nieta].deltaU[0][3]);
			dutau = (dux*ux + duy*uy + dueta*ueta)/utau;
			de = 0.5*(arena[ix][iy][ieta].deltaP[0]/eos->get_dpOverde2(fmax(0.03, arena[ix][iy][ieta].epsilon), 
										   arena[ix][iy][ieta].rhob) 
				  + arena[nix][niy][nieta].deltaP[0]/eos->get_dpOverde2(fmax(0.03, arena[nix][niy][nieta].epsilon),
											arena[nix][niy][nieta].rhob) );
			dp = 0.5*(arena[ix][iy][ieta].deltaP[0] + arena[nix][niy][nieta].deltaP[0]);
			dT = dp/eos->s_func(iepsFO, eos->p_func(iepsFO, rhob), rhob, DATA);
			de = dp/eos->get_dpOverde2(iepsFO, rhob);
			//dWtautau = 0.5*(arena[ix][iy][ieta].deltaW[0][0][0] + arena[nix][niy][nieta].deltaW[0][0][0]);
			//dWtaux = 0.5*(arena[ix][iy][ieta].deltaW[0][0][1] + arena[nix][niy][nieta].deltaW[0][0][1]);
			//dWtauy = 0.5*(arena[ix][iy][ieta].deltaW[0][0][2] + arena[nix][niy][nieta].deltaW[0][0][2]);
			//dWtaueta = 0.5*(arena[ix][iy][ieta].deltaW[0][0][3] + arena[nix][niy][nieta].deltaW[0][0][3]);
			dWxx = 0.5*(arena[ix][iy][ieta].deltaW[0][1][1] + arena[nix][niy][nieta].deltaW[0][1][1]);
			dWxy = 0.5*(arena[ix][iy][ieta].deltaW[0][1][2] + arena[nix][niy][nieta].deltaW[0][1][2]);
			dWxeta = 0.5*(arena[ix][iy][ieta].deltaW[0][1][3] + arena[nix][niy][nieta].deltaW[0][1][3]);
			dWyy = 0.5*(arena[ix][iy][ieta].deltaW[0][2][2] + arena[nix][niy][nieta].deltaW[0][2][2]);
			dWyeta = 0.5*(arena[ix][iy][ieta].deltaW[0][2][3] + arena[nix][niy][nieta].deltaW[0][2][3]);
			dWetaeta = 0.5*(arena[ix][iy][ieta].deltaW[0][3][3] + arena[nix][niy][nieta].deltaW[0][3][3]);
			dWtaux = (ux*dWxx + uy*dWxy + ueta*dWxeta -dutau*Wtaux + dux*Wxx + duy*Wxy + dueta*Wxeta)/utau;
			dWtauy = (ux*dWxy + uy*dWyy + ueta*dWyeta -dutau*Wtauy + dux*Wxy + duy*Wyy + dueta*Wyeta)/utau;
			dWtaueta = (ux*dWxeta + uy*dWyeta + ueta*dWetaeta -dutau*Wtaueta + dux*Wxeta + duy*Wyeta + dueta*Wetaeta)/utau;
			dWtautau = (ux*dWtaux + uy*dWtauy + ueta*dWtaueta -dutau*Wtautau + dux*Wtaux + duy*Wtauy + dueta*Wtaueta)/utau;
		      }

		      s_file << std::setprecision(10) << tauf << " " << xf << " " << yf << " " << etaf << " " 
			     << FULLSU[0] << " " <<FULLSU[1] << " " <<FULLSU[2] << " " <<FULLSU[3] 
			     << " " <<  utau << " " << ux << " " << uy << " " << ueta << " " 
			     << iepsFO << " " << TFO << " " << muB << " " << sFO << " " 
			     << Wtautau << " " << Wtaux << " " << Wtauy << " " << Wtaueta << " " 
			     << Wxx << " " << Wxy << " " << Wxeta << " " << Wyy << " " << Wyeta << " " << Wetaeta; 
		      if(DATA->turn_on_bulk) s_file << " " << pi_b;
		      s_file << endl; 

		      if(DATA->fluctuatingHydroFlag == 2){
			ds_file << std::setprecision(10) << dT << " " << dutau << " " << dux << " " << duy << " " << dueta << " "
				<< de << " " << dp << " "
				<< dWtautau << " " << dWtaux << " " << dWtauy << " " << dWtaueta << " "
				<< dWxx << " " << dWxy << " " << dWxeta << " "
				<< dWyy << " " << dWyeta << " "
				<< dWetaeta << endl;
		      }
		    }// if grid pair straddles freeze out density
		}// if (ieta<maxEta)	      
	      
	      // check backward in time
	      xf = x;
	      yf = y;
	      etaf = eta;
	      tauf = tau-DTAU/2.;

	      nix = ix;
	      niy = iy;
	      nieta = ieta;
	      FULLSU[0]=DX*DY*DETA/subsections;
	      FULLSU[1]=0.;
	      FULLSU[2]=0.;
	      FULLSU[3]=0.;
	      
	      neighbor = home;
	      //if (DATA->useEpsFO)
	      //{
	      home = arena[ix][iy][ieta].epsilon_prev;
	      //}
	      //else
	      //{
	      //  //   home = arena[ix][iy][ieta].T;
	      //  //   home = eos->get_temperature(arena[ix][iy][ieta].epsilon_prev,arena[ix][iy][ieta].rhob);
	      //    
	      //  // Need temperature at previous time step. Use value of rhob at previous time step,
	      //  // or else use T_prev if that is available
	      //  home = eos->get_temperature(arena[ix][iy][ieta].epsilon_prev,arena[ix][iy][ieta].rhob_prev);
	      //}
	      if(((home > FO) && (neighbor <= FO)) || ((home <= FO) && (neighbor > FO)))
		{
		  //cerr << "Is the 4th condition ever satisfied?" << endl;
		  if (home<=FO)
		    SIG=-1.;
		  else SIG=1.;
		  for (int orient = 0; orient < 4; orient++) FULLSU[orient]*=SIG;
		  // simple linear interpolation of all independent variables
		  ux = 0.5*(arena[ix][iy][ieta].u[0][1] + arena[nix][niy][nieta].u_prev[1]);
		  uy = 0.5*(arena[ix][iy][ieta].u[0][2] + arena[nix][niy][nieta].u_prev[2]);
		  ueta = 0.5*(arena[ix][iy][ieta].u[0][3] + arena[nix][niy][nieta].u_prev[3]);
		  utau = sqrt(1 + ux*ux + uy*uy + ueta*ueta);
		  iepsFO = 0.5*(arena[ix][iy][ieta].epsilon + arena[ix][iy][ieta].epsilon_prev);
		  double iP = eos->p_func(iepsFO, rhob);
		  double sFO = eos->s_func(iepsFO, iP, rhob, DATA);
		  //I don't seem to have accesss to value of rhob from previous timestep.  Use value at current timestep.
		  //   rhob = 0.5*(arena[ix][iy][ieta].rhob + arena[nix][niy][nieta].rhob);
		  rhob = 0.5*(arena[ix][iy][ieta].rhob + arena[nix][niy][nieta].rhob_prev);
		  //pi_b = 0.5*(arena[ix][iy][ieta].pi_b[0] + arena[nix][niy][nieta].pi_b_prev);
		  //TFO = eos->get_temperature(iepsFO, rhob);
		  //muB = eos->get_mu(iepsFO, rhob);
		  TFO = eos->interpolate2(iepsFO, 0., 1);
		  muB = 0.;
		  //eps_plus_p_over_T_FO=(iepsFO+eos->get_pressure(iepsFO, rhob))/TFO;
		  Wxx = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][1] + arena[nix][niy][nieta].W_prev[1][1]);
		  Wxy = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][2] + arena[nix][niy][nieta].W_prev[1][2]);
		  Wxeta = 0.5*(arena[ix][iy][ieta].Wmunu[0][1][3] + arena[nix][niy][nieta].W_prev[1][3]);
		  Wyy = 0.5*(arena[ix][iy][ieta].Wmunu[0][2][2] + arena[nix][niy][nieta].W_prev[2][2]);
		  Wyeta = 0.5*(arena[ix][iy][ieta].Wmunu[0][2][3] + arena[nix][niy][nieta].W_prev[2][3]);

		  Wetaeta = (2.*(ux*uy*Wxy 
                            + ux*ueta*Wxeta
				 + uy*ueta*Wyeta )
			     -( utau*utau - ux*ux )*Wxx
			     -( utau*utau - uy*uy )*Wyy
			     )/( utau*utau - ueta*ueta ) ;
		  Wtaux = (ux*Wxx + uy*Wxy + ueta*Wxeta)/utau;
		  Wtauy = (ux*Wxy + uy*Wyy + ueta*Wyeta)/utau;
		  Wtaueta = (ux*Wxeta + uy*Wyeta + ueta*Wetaeta)/utau;
		  Wtautau = (ux*Wtaux + uy*Wtauy + ueta*Wtaueta)/utau;
		    
		      if(DATA->fluctuatingHydroFlag == 2){
			//if(arena[nix][niy][nieta].deltaP_prev == 0.){dT = 0.;}
			//else{dT = 0.5*(arena[nix][niy][nieta].deltaP_prev
			//	  /eos->s_func(arena[nix][niy][nieta].epsilon_prev, 
			//		       eos->p_func(arena[nix][niy][nieta].epsilon_prev, 
			//				   arena[nix][niy][nieta].rhob_prev), 
			//		       arena[nix][niy][nieta].rhob_prev, DATA));}
			//if(arena[ix][iy][ieta].deltaP[0] != 0.){
			//dT += 0.5*(arena[ix][iy][ieta].deltaP[0]
			//        /eos->s_func(arena[ix][iy][ieta].epsilon, 
			//		       eos->p_func(arena[ix][iy][ieta].epsilon, 
			//				   arena[ix][iy][ieta].rhob), 
			//		       arena[ix][iy][ieta].rhob, DATA));}
			//dutau = 0.5*(arena[ix][iy][ieta].deltaU[0][0] + arena[nix][niy][nieta].deltaU_prev[0]);
			dux = 0.5*(arena[ix][iy][ieta].deltaU[0][1] + arena[nix][niy][nieta].deltaU_prev[1]);
			duy = 0.5*(arena[ix][iy][ieta].deltaU[0][2] + arena[nix][niy][nieta].deltaU_prev[2]);
			dueta = 0.5*(arena[ix][iy][ieta].deltaU[0][3] + arena[nix][niy][nieta].deltaU_prev[3]);
			dutau = (dux*ux + duy*uy + dueta*ueta)/utau;
			//de = 0.5*(arena[ix][iy][ieta].deltaP[0]/eos->get_dpOverde2(fmax(0.03, arena[ix][iy][ieta].epsilon), 
			//							   arena[ix][iy][ieta].rhob) 
			//	  + arena[nix][niy][nieta].deltaP_prev/eos->get_dpOverde2(fmax(0.03, arena[nix][niy][nieta].epsilon_prev),
			//								arena[nix][niy][nieta].rhob_prev) );
			dp = 0.5*(arena[ix][iy][ieta].deltaP[0] + arena[nix][niy][nieta].deltaP_prev);
			dT = dp/eos->s_func(iepsFO, eos->p_func(iepsFO, rhob), rhob, DATA);
			de = dp/eos->get_dpOverde2(iepsFO, rhob);
			//dWtautau = 0.5*(arena[ix][iy][ieta].deltaW[0][0][0] + arena[nix][niy][nieta].deltaW_prev[0][0]);
			//dWtaux = 0.5*(arena[ix][iy][ieta].deltaW[0][0][1] + arena[nix][niy][nieta].deltaW_prev[0][1]);
			//dWtauy = 0.5*(arena[ix][iy][ieta].deltaW[0][0][2] + arena[nix][niy][nieta].deltaW_prev[0][2]);
			//dWtaueta = 0.5*(arena[ix][iy][ieta].deltaW[0][0][3] + arena[nix][niy][nieta].deltaW_prev[0][3]);
			dWxx = 0.5*(arena[ix][iy][ieta].deltaW[0][1][1] + arena[nix][niy][nieta].deltaW_prev[1][1]);
			dWxy = 0.5*(arena[ix][iy][ieta].deltaW[0][1][2] + arena[nix][niy][nieta].deltaW_prev[1][2]);
			dWxeta = 0.5*(arena[ix][iy][ieta].deltaW[0][1][3] + arena[nix][niy][nieta].deltaW_prev[1][3]);
			dWyy = 0.5*(arena[ix][iy][ieta].deltaW[0][2][2] + arena[nix][niy][nieta].deltaW_prev[2][2]);
			dWyeta = 0.5*(arena[ix][iy][ieta].deltaW[0][2][3] + arena[nix][niy][nieta].deltaW_prev[2][3]);
			dWetaeta = 0.5*(arena[ix][iy][ieta].deltaW[0][3][3] + arena[nix][niy][nieta].deltaW_prev[3][3]);
			dWtaux = (ux*dWxx + uy*dWxy + ueta*dWxeta -dutau*Wtaux + dux*Wxx + duy*Wxy + dueta*Wxeta)/utau;
			dWtauy = (ux*dWxy + uy*dWyy + ueta*dWyeta -dutau*Wtauy + dux*Wxy + duy*Wyy + dueta*Wyeta)/utau;
			dWtaueta = (ux*dWxeta + uy*dWyeta + ueta*dWetaeta -dutau*Wtaueta + dux*Wxeta + duy*Wyeta + dueta*Wetaeta)/utau;
			dWtautau = (ux*dWtaux + uy*dWtauy + ueta*dWtaueta -dutau*Wtautau + dux*Wtaux + duy*Wtauy + dueta*Wtaueta)/utau;
		      }

		  for (int i = 0; i < subsections; i++)
		    {
		      etaf = eta - DETA/2 + DETA/2/subsections + (i*DETA)/(subsections);
		          
		      s_file << std::setprecision(10) << tauf << " " << xf << " " << yf << " " << etaf << " " 
			     << FULLSU[0] << " " <<FULLSU[1] << " " <<FULLSU[2] << " " <<FULLSU[3] 
			     << " " <<  utau << " " << ux << " " << uy << " " << ueta << " " 
			     << iepsFO << " " << TFO << " " << muB << " " << sFO << " " 
			     << Wtautau << " " << Wtaux << " " << Wtauy << " " << Wtaueta << " " 
			     << Wxx << " " << Wxy << " " << Wxeta << " " << Wyy << " " << Wyeta << " " << Wetaeta; 
		      if(DATA->turn_on_bulk) s_file << " " << pi_b;
		      s_file << endl; 
		      if(DATA->fluctuatingHydroFlag == 2){
			ds_file << std::setprecision(10) << dT << " " << dutau << " " << dux << " " << duy << " " << dueta << " "
				<< de << " " << dp << " "
				<< dWtautau << " " << dWtaux << " " << dWtauy << " " << dWtaueta << " "
				<< dWxx << " " << dWxy << " " << dWxeta << " "
				<< dWyy << " " << dWyeta << " "
				<< dWetaeta << endl;
		      }
		    }
		  //cerr << "4th set of statements exited" << endl;
		} // if grid pair straddles density
	            
	      //       if(shown==0)
	      // {
	      //   fprintf(stderr,"DTAU=%f\n", DTAU);
	      //   fprintf(stderr,"DX=%f\n", DX);
	      //   fprintf(stderr,"DY=%f\n", DY);
	      //   fprintf(stderr,"DETA=%f\n", DETA);
	      //   shown=1;
	      // }
	      
	    }// eta loop

	  //  At the end of the eta grid in the +eta direction, check neighboring cell in eta, 
	  //  received from neighboring processor.  Only check if surface is in between cells in eta direction
	  //  (steps in x, y, and tau are checked on the neighboring processor)
	  if (rank != size-1 && size > 1) 
	    {
	      //cerr << "rank = " << rank << endl;

	      eta = (DATA->delta_eta)*(neta+DATA->neta*rank) - (DATA->eta_size)/2.0;
	      
	      // from Rneighbor, check backward in eta
	      xf = x;
	      yf = y;
	      etaf = eta - DETA/2;
	      tauf = tau;

	      nix = ix;
	      niy = iy;
	      nieta = maxEta;
	      FULLSU[0]=0.;
	      FULLSU[1]=0.;
	      FULLSU[2]=0.;
	      FULLSU[3]=DTAU*DX*DY;
	      
	      home = arena[ix][iy][nieta].epsilon;
	      neighbor = Rneighbor_eps[nix][niy];

	      if(((home > FO) && (neighbor <= FO)) || ((home <= FO) && (neighbor > FO)))
		{
		  //cerr << "Is the 5th condition ever satisfied?" << endl;
		  if (home<=FO)
		    SIG=-1.;
		  else SIG=1.;
		  for(int orient = 0; orient < 4; orient++){FULLSU[orient]*=SIG;}
		  // simple linear interpolation of all independent variables
		  ux = 0.5*(Rneighbor_ux[ix][iy] + arena[nix][niy][nieta].u[0][1]);
		  uy = 0.5*(Rneighbor_uy[ix][iy] + arena[nix][niy][nieta].u[0][2]);
		  ueta = 0.5*(Rneighbor_ueta[ix][iy] + arena[nix][niy][nieta].u[0][3]);
		  utau = sqrt(1 + ux*ux + uy*uy + ueta*ueta);
		  iepsFO = 0.5*(arena[ix][iy][nieta].epsilon + Rneighbor_eps[nix][niy]);
		  rhob = 0.5*(Rneighbor_rhob[ix][iy] + arena[nix][niy][nieta].rhob);
		  double iP = eos->p_func(iepsFO, rhob);
		  double sFO = eos->s_func(iepsFO, iP, rhob, DATA);
		  //pi_b = 0.5*(Rneighbor_Pi_b[ix][iy] + arena[nix][niy][nieta].pi_b[0]);
		  //TFO = eos->get_temperature(iepsFO, rhob);
		  //muB = eos->get_mu(iepsFO, rhob);
		  TFO = eos->interpolate2(iepsFO, 0., 1);
		  muB = 0.;
		  //eps_plus_p_over_T_FO=(iepsFO+eos->get_pressure(iepsFO, rhob))/TFO;
		  Wxx = 0.5*(Rneighbor_Wxx[ix][iy] + arena[nix][niy][nieta].Wmunu[0][1][1]);
		  Wxy = 0.5*(Rneighbor_Wxy[ix][iy] + arena[nix][niy][nieta].Wmunu[0][1][2]);
		  Wxeta = 0.5*(Rneighbor_Wxeta[ix][iy] + arena[nix][niy][nieta].Wmunu[0][1][3]);
		  Wyy = 0.5*(Rneighbor_Wyy[ix][iy] + arena[nix][niy][nieta].Wmunu[0][2][2]);
		  Wyeta = 0.5*(Rneighbor_Wyeta[ix][iy] + arena[nix][niy][nieta].Wmunu[0][2][3]);

		  Wetaeta = (2.*(ux*uy*Wxy 
                            + ux*ueta*Wxeta
				 + uy*ueta*Wyeta )
			     -( utau*utau - ux*ux )*Wxx
			     -( utau*utau - uy*uy )*Wyy
			     )/( utau*utau - ueta*ueta ) ;
		  Wtaux = (ux*Wxx + uy*Wxy + ueta*Wxeta)/utau;
		  Wtauy = (ux*Wxy + uy*Wyy + ueta*Wyeta)/utau;
		  Wtaueta = (ux*Wxeta + uy*Wyeta + ueta*Wetaeta)/utau;
		  Wtautau = (ux*Wtaux + uy*Wtauy + ueta*Wtaueta)/utau;

		      if(DATA->fluctuatingHydroFlag == 2){
			//dT = 0.5*Rneighbor_dT[ix][iy];
			//if(arena[nix][niy][nieta].deltaP[0] != 0){
			// dT += 0.5*arena[nix][niy][nieta].deltaP[0]
			//  /eos->s_func(arena[nix][niy][nieta].epsilon, 
			//		 eos->p_func(arena[nix][niy][nieta].epsilon, 
			//			     arena[nix][niy][nieta].rhob), 
			//		 arena[nix][niy][nieta].rhob, DATA);}
			//dutau = 0.5*(Rneighbor_dutau[ix][iy] + arena[nix][niy][nieta].deltaU[0][0]);
			dux = 0.5*(Rneighbor_dux[ix][iy] + arena[nix][niy][nieta].deltaU[0][1]);
			duy = 0.5*(Rneighbor_duy[ix][iy] + arena[nix][niy][nieta].deltaU[0][2]);
			dueta = 0.5*(Rneighbor_dueta[ix][iy] + arena[nix][niy][nieta].deltaU[0][3]);
			dutau = (dux*ux + duy*uy + dueta*ueta)/utau;
			//de = 0.5*(Rneighbor_dp[ix][iy]/eos->get_dpOverde2(fmax(0.03, Rneighbor_eps[ix][iy]), Rneighbor_rhob[ix][iy]) 
			//	  + arena[nix][niy][nieta].deltaP[0]/eos->get_dpOverde2(fmax(0.03, arena[nix][niy][nieta].epsilon),
			//								arena[nix][niy][nieta].rhob) );
			dp = 0.5*(arena[ix][iy][ieta].deltaP[0] + arena[nix][niy][nieta].deltaP[0]);
			dT = dp/eos->s_func(iepsFO, eos->p_func(iepsFO, rhob), rhob, DATA);
			de = dp/eos->get_dpOverde2(iepsFO, rhob);
			//dWtautau = 0.5*(Rneighbor_dWtautau[ix][iy] + arena[nix][niy][nieta].deltaW[0][0][0]);
			//dWtaux = 0.5*(Rneighbor_dWtaux[ix][iy] + arena[nix][niy][nieta].deltaW[0][0][1]);
			//dWtauy = 0.5*(Rneighbor_dWtauy[ix][iy] + arena[nix][niy][nieta].deltaW[0][0][2]);
			//dWtaueta = 0.5*(Rneighbor_dWtaueta[ix][iy] + arena[nix][niy][nieta].deltaW[0][0][3]);
			dWxx = 0.5*(Rneighbor_dWxx[ix][iy] + arena[nix][niy][nieta].deltaW[0][1][1]);
			dWxy = 0.5*(Rneighbor_dWxy[ix][iy] + arena[nix][niy][nieta].deltaW[0][1][2]);
			dWxeta = 0.5*(Rneighbor_dWxeta[ix][iy] + arena[nix][niy][nieta].deltaW[0][1][3]);
			dWyy = 0.5*(Rneighbor_dWyy[ix][iy] + arena[nix][niy][nieta].deltaW[0][2][2]);
			dWyeta = 0.5*(Rneighbor_dWyeta[ix][iy] + arena[nix][niy][nieta].deltaW[0][2][3]);
			dWetaeta = 0.5*(Rneighbor_dWetaeta[ix][iy] + arena[nix][niy][nieta].deltaW[0][3][3]);
			dWtaux = (ux*dWxx + uy*dWxy + ueta*dWxeta -dutau*Wtaux + dux*Wxx + duy*Wxy + dueta*Wxeta)/utau;
			dWtauy = (ux*dWxy + uy*dWyy + ueta*dWyeta -dutau*Wtauy + dux*Wxy + duy*Wyy + dueta*Wyeta)/utau;
			dWtaueta = (ux*dWxeta + uy*dWyeta + ueta*dWetaeta -dutau*Wtaueta + dux*Wxeta + duy*Wyeta + dueta*Wetaeta)/utau;
			dWtautau = (ux*dWtaux + uy*dWtauy + ueta*dWtaueta -dutau*Wtautau + dux*Wtaux + duy*Wtauy + dueta*Wtaueta)/utau;
		      }

		  s_file << std::setprecision(10) << tauf << " " << xf << " " << yf << " " << etaf << " " 
			 << FULLSU[0] << " " <<FULLSU[1] << " " <<FULLSU[2] << " " <<FULLSU[3] 
			 << " " <<  utau << " " << ux << " " << uy << " " << ueta << " " 
			 << iepsFO << " " << TFO << " " << muB << " " << sFO << " " 
			 << Wtautau << " " << Wtaux << " " << Wtauy << " " << Wtaueta << " " 
			 << Wxx << " " << Wxy << " " << Wxeta << " " << Wyy << " " << Wyeta << " " << Wetaeta; 
		  if(DATA->turn_on_bulk) s_file << " " << pi_b;
		  s_file << endl; 

		  if(DATA->fluctuatingHydroFlag == 2){
		    ds_file << std::setprecision(10) << dT << " " << dutau << " " << dux << " " << duy << " " << dueta << " "
			    << de << " " << dp << " "
			    << dWtautau << " " << dWtaux << " " << dWtauy << " " << dWtaueta << " "
			    << dWxx << " " << dWxy << " " << dWxeta << " "
			    << dWyy << " " << dWyeta << " "
			    << dWetaeta << endl;
		  }
		}// if grid pair straddles freeze out density
	      //cerr << "ix = " << ix << ", iy = " << iy << endl;

	    }// if (rank != size-1)
	}// y loop
    }// x loop
  //   delete []package;
  s_file.close();
  freeze_file.close();
  
  //cerr << "OK closing the files" << endl;

  //MPI_Barrier(MPI_COMM_WORLD);

  cerr << "rank = " << rank << endl;

  // check if all cells are frozen out 
  if(size > 1){
    if (rank!=0) //send to rank 0
      {
	//cerr << "About to send frozen" << endl;
	MPI::COMM_WORLD.Send(&frozen, 1, MPI::INT, 0, 38);
	//cerr << "frozen sent" << endl;
	//cerr << "frozen = " << frozen << endl;
      }
    if (rank==0) // receive from all ranks >0
      {  
	int frozenFrom = 1;
	//allfrozen = frozen;
	for (from = 1; from<size; from++)
	  {
	    //cerr << "frozen = " << frozen << endl;
	    //cerr << "About to receive frozenFrom" << endl;
	    //MPI::COMM_WORLD.Recv(&frozen, 1, MPI::INT, from, 38);
	    //allfrozen = (allfrozen && frozen);
	    MPI::COMM_WORLD.Recv(&frozenFrom, 1, MPI::INT, from, 38);
	    //cerr << "Received frozenFrom" << endl;
	    allfrozen = (allfrozen && frozenFrom);
	  }
	//for (from = 1; from<size; from ++)
	//{
	//  MPI::COMM_WORLD.Send(&allfrozen,1,MPI::INT,from, 39);
	//}
	if(allfrozen)
	  //{
	  //cout << "All cells frozen out. Exiting." << endl;
	  //exit(rank);
	  //}
	  {
	    cout << "All cells frozen out. Exiting hydro evolution." << endl;
	    exit(rank);
	    // write OSCAR header if wanted:  (seems to be unused at the moment, so I'll keep the code, but comment out)
	    //   FILE *oout_file;
	    //   char* oout_name = "OSCARheader.dat";
	    //   oout_file = fopen(oout_name, "w");
	    //   if(DATA->outputEvolutionData)
	    //     {
	    //       if(DATA->viscosity_flag==1)
	    // fprintf(oout_file,"OSCAR2008H viscous history\n");
	    //       else
	    // fprintf(oout_file,"OSCAR2008H ideal history\n");
	    // 
	    //       fprintf(oout_file,"INIT: MUSIC %s+%s at b=%f fm, Glauber\n",DATA->Target.c_str(), DATA->Projectile.c_str(),DATA->b);
	    //       fprintf(oout_file,"INIT: \n");
	    //       
	    //       if (DATA->whichEOS==0)
	    // {
	    //   fprintf(oout_file,"EOS: ideal gas EOS \n");
	    // }
	    //       else if (DATA->whichEOS==1)
	    // {
	    //     fprintf(oout_file,"EOS: EOS-Q from AZHYDRO \n");
	    // }
	    //       else if (DATA->whichEOS==2)
	    // {
	    //     fprintf(oout_file,"EOS: lattice EOS from Huovinen/Petreczky \n");
	    // }
	    //       else if (DATA->whichEOS==3)
	    // {
	    //     fprintf(oout_file,"EOS: lattice EOS from Huovinen/Petreczky with partial chemical equilibrium (PCE) chem. f.o. at 150 MeV\n");
	    // }
	    //       else if (DATA->whichEOS==4)
	    // {
	    //     fprintf(oout_file,"EOS: lattice EOS from Huovinen/Petreczky with partial chemical equilibrium (PCE) chem. f.o. at 155 MeV\n");
	    // }
	    //       else if (DATA->whichEOS==5)
	    // {
	    //     fprintf(oout_file,"EOS: lattice EOS from Huovinen/Petreczky with partial chemical equilibrium (PCE) chem. f.o. at 160 MeV\n");
	    // }
	    //       else if (DATA->whichEOS==6)
	    // {
	    //     fprintf(oout_file,"EOS: lattice EOS from Huovinen/Petreczky with partial chemical equilibrium (PCE) chem. f.o. at 165 MeV\n");
	    // }
	    // 
	    //       if(DATA->turn_on_rhob==0)
	    // fprintf(oout_file,"CHARGES: none\n");
	    //       else
	    // fprintf(oout_file,"CHARGES: baryon\n");
	    // 
	    //       fprintf(oout_file,"HYPER: full evolution\n");
	    //  
	    //       fprintf(oout_file,"GEOM: 3d\n");
	    //        
	    //       fprintf(oout_file,"GRID: Euler\n");
	    //  
	    //       fprintf(oout_file,"%d %d %d %d %d %d %d\n", static_cast<int>(((tau-DATA->tau0)/static_cast<double>(DATA->delta_tau)+1)/10.)+1, 
	    //       DATA->nx, DATA->ny, DATA->neta, 0, 0, 0);
	    //     
	    //       fprintf(oout_file,"%f %f %f %f %f %f %f %f\n", DATA->tau0, tau, -DATA->x_size/2., DATA->x_size/2., -DATA->y_size/2., DATA->y_size/2.,
	    //       -DATA->eta_size/2., DATA->eta_size/2.);
	    // 
	    //       if(DATA->viscosity_flag==1)
	    // {
	    //   if(DATA->turn_on_shear && DATA->turn_on_bulk)
	    //     {
	    //       fprintf(oout_file,"VISCOSITY: shear and bulk viscosity\n");
	    //       fprintf(oout_file,"VISCOSITY: eta/s = %f, zeta/s = %f\n", DATA->shear_to_s, DATA->bulk_to_s);
	    //     }
	    //   else if(DATA->turn_on_shear==1 && DATA->turn_on_bulk==0)
	    //     {
	    //       fprintf(oout_file,"VISCOSITY: shear viscosity only\n");
	    //       fprintf(oout_file,"VISCOSITY: eta/s = %f\n", DATA->shear_to_s);
	    //     }
	    //   else if(DATA->turn_on_shear==0 && DATA->turn_on_bulk==1)
	    //     {
	    //       fprintf(oout_file,"VISCOSITY: bulk viscosity only\n");
	    //       fprintf(oout_file,"VISCOSITY: zeta/s = %f\n", DATA->bulk_to_s);
	    //     }
	    //   else if(DATA->turn_on_shear==0 && DATA->turn_on_bulk==0)
	    //     {
	    //       fprintf(oout_file,"VISCOSITY: none\n");
	    //       fprintf(oout_file,"VISCOSITY: \n");
	    //     }
	    // }
	    //       else
	    // {
	    //   fprintf(oout_file,"VISCOSITY: none\n");
	    //   fprintf(oout_file,"VISCOSITY: \n");
	    // }
	    //       
	    //       fprintf(oout_file,"COMM:\n");
	    // 
	    //       fprintf(oout_file,"END_OF_HEADER\n");
	    //     }
	    //     
	    //   fclose(oout_file);
	    //   
	    //   int check=system ("cat OSCARheader.dat OSCAR.dat >OSCARoutput.dat");
	    //   if(check==0)
	    //     {
	    //       system ("rm OSCAR.dat");
	    //       system ("rm OSCARheader.dat");
	    //     }
	    //   
	    
	    //cerr << "OK before the system call" << endl;
	    
	    //system("cat surface?.dat surface??.dat > surface.dat");
	    //system("rm surface?.dat surface??.dat 2> /dev/null");
	    
	    //cerr << "OK after the system call" << endl;
	    
	    //   MPI::Finalize();
	    //    exit(1);
	    //   return 1;
	  }
	
      }
  }
  else{
    if(frozen){
      cout << "All cells frozen out. Exiting hydro evolution." << endl;
      exit(1);
    }
  }
  //if (rank!=0)
  //{
  //MPI::COMM_WORLD.Recv(&allfrozen, 1, MPI::INT, 0, 39);
  //if (allfrozen)
  //  {
  //	//   cout << "All cells frozen out. Exiting." << endl;
  //	//   MPI::Finalize();
  //	//    exit(1);
  //	//   return 1;
  //    }
  //}
  
  cerr << "OK up to here..." << endl;

  delete[] packageWxx;
  delete[] packageWxy;
  delete[] packageWxeta;
  delete[] packageWyy;
  delete[] packageWyeta; 
  
  delete[] package;
  delete[] packageux;
  delete[] packageuy;
  delete[] packageueta;
  delete[] packagerhob;
  
  util->mtx_free(Rneighbor_eps,nx+1,ny+1);
  util->mtx_free(Rneighbor_ux,nx+1,ny+1);
  util->mtx_free(Rneighbor_uy,nx+1,ny+1);
  util->mtx_free(Rneighbor_ueta,nx+1,ny+1);
  util->mtx_free(Rneighbor_rhob,nx+1,ny+1);
  
  util->mtx_free(Rneighbor_Wxx,nx+1,ny+1);
  util->mtx_free(Rneighbor_Wxy,nx+1,ny+1);
  util->mtx_free(Rneighbor_Wxeta,nx+1,ny+1);
  util->mtx_free(Rneighbor_Wyy,nx+1,ny+1);
  util->mtx_free(Rneighbor_Wyeta,nx+1,ny+1);
  
  //cerr << "Hypercubic freeze-out finished" << endl;
  
  return allfrozen;
}

