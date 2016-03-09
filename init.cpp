#include "util.h"
#include "grid.h"
#include "init.h"
#include "eos.h"
#include <string>
#include <sstream>

using namespace std;

Init::Init(EOS *eosIn, Glauber *glauberIn, int gsl_seed)
{
  eos = new EOS;
  eos = eosIn;
  util = new Util;
  glauber = new Glauber;
  glauber = glauberIn;
  random = new Random;

  gsl_ran = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (gsl_ran, gsl_seed);
}

// destructor
Init::~Init()
{
  delete random;
  delete eos;
  delete util;
  delete glauber;
}

void Init::InitArena(InitData *DATA, Grid ****arena, Grid ****Lneighbor, Grid ****Rneighbor, int size, int rank)
{
  Grid *helperGrid;
  helperGrid = new Grid;
  *arena = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, DATA->neta);
  *Lneighbor = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, 2);
  *Rneighbor = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, 2);
  
  while (InitTJb(DATA, arena, Lneighbor, Rneighbor, size, rank)==0);

  cout << "rank " << rank << " ok." << endl;
  LinkNeighbors(DATA, arena, size, rank);
  delete helperGrid;
}/* InitArena */


void Init::LinkNeighbors(InitData *DATA, Grid ****arena, int size, int rank)
{
 int ix, iy, ieta, mu, nu, nx, ny, neta;
 int x_nbr, y_nbr, eta_nbr;

 nx = DATA->nx;
 ny = DATA->ny;
 neta = DATA->neta;

/* allocate memory */
 for(ix=0; ix<=nx; ix++)
  {
   for(iy=0; iy<=ny; iy++)
    {
     for(ieta=0; ieta<neta; ieta++)
      {
       (*arena)[ix][iy][ieta].nbr_p_1 
               = (Grid **) malloc (sizeof(Grid *)*4); 
       (*arena)[ix][iy][ieta].nbr_m_1 
               = (Grid **) malloc (sizeof(Grid *)*4); 
       	 
       (*arena)[ix][iy][ieta].position[1] = ix;
       (*arena)[ix][iy][ieta].position[2] = iy;
       (*arena)[ix][iy][ieta].position[3] = ieta;
      }/* ieta */
    }/* iy */
  }/* ix */

 for(ix=0; ix<=nx; ix++)
  {
   for(iy=0; iy<=ny; iy++)
    {
     for(ieta=0; ieta<neta; ieta++)
      {
	if(ix != nx)
	  { (*arena)[ix][iy][ieta].nbr_p_1[1] = &(*arena)[ix+1][iy][ieta]; 
	  }
	else 
	  { //(*arena)[ix][iy][ieta].nbr_p_1[1] = NULL; 
	    (*arena)[ix][iy][ieta].nbr_p_1[1] = &(*arena)[0][iy][ieta]; 
	  }

	if(ix != 0)
          { (*arena)[ix][iy][ieta].nbr_m_1[1] = &(*arena)[ix-1][iy][ieta]; 
	  }
	else
	  { //(*arena)[ix][iy][ieta].nbr_m_1[1] = NULL;
	    (*arena)[ix][iy][ieta].nbr_m_1[1] = &(*arena)[nx][iy][ieta]; 
	  }
         
	if(iy != ny)
          { (*arena)[ix][iy][ieta].nbr_p_1[2] = &(*arena)[ix][iy+1][ieta]; 
	  }
	else
          { //(*arena)[ix][iy][ieta].nbr_p_1[2] = NULL;
	    (*arena)[ix][iy][ieta].nbr_p_1[2] = &(*arena)[ix][0][ieta];
	  }

	if(iy != 0)
          { (*arena)[ix][iy][ieta].nbr_m_1[2] = &(*arena)[ix][iy-1][ieta]; 
          }
	else
          { //(*arena)[ix][iy][ieta].nbr_m_1[2] = NULL;
	    (*arena)[ix][iy][ieta].nbr_m_1[2] = &(*arena)[ix][ny][ieta]; 
          }

	// do not care which rank it is - that is dealt with in evolve.cpp
	if(ieta != neta-1)
          { (*arena)[ix][iy][ieta].nbr_p_1[3] = &(*arena)[ix][iy][ieta+1]; 
	  }
	else 
	  {
	    if(size!=1){
	      (*arena)[ix][iy][ieta].nbr_p_1[3] = NULL;
	    }
	    else{
	      //Only to be used if size = 1:
	      (*arena)[ix][iy][ieta].nbr_p_1[3] = &(*arena)[ix][iy][0];
	    }
	  }

	if(ieta != 0)
	  { (*arena)[ix][iy][ieta].nbr_m_1[3] = &(*arena)[ix][iy][ieta-1]; 
	  }
	else
	  {
	    if(size!=1){
	      (*arena)[ix][iy][ieta].nbr_m_1[3] = NULL; 
	    }
	    else{
	      //Only to be used if size = 1:
	      (*arena)[ix][iy][ieta].nbr_m_1[3] = &(*arena)[ix][iy][neta-1];
	    }
	  }
      }/* ieta */
    }/* iy */
  }/* ix */
}/* LinkNeighbors */


void Init::sampleTA()
{
  int A = static_cast<int>(glauber->nucleusA());

  //  cout << "A=" << A << endl;
  ReturnValue rv, rv2;

  for (int i = 0; i < A; i++) // get all nucleon coordinates
    {
      rv = glauber->SampleTARejection(random);
      rv2 = glauber->SampleTARejection(random);
      nucleusA.push_back(rv);
      nucleusB.push_back(rv2);
    }
}  

int Init::InitTJb(InitData *DATA, Grid ****arena, Grid ****Lneighbor, Grid ****Rneighbor, int size, int rank)
{
 double epsilon0, testEpsilon, epsilon, R_A, a_A, p, T, h, u[4], x, y, eta, rho, volume;
 double cosheta, sinheta, rhob, exparg1, exparg;
 double eta0, eta_fall_off, eta_flat, a_short, a_long;
 int ix, iy, ieta, ietar, mu, nu, rk_order;
 int initializeEntropy = DATA->initializeEntropy;
 stringstream numIn;
 numIn << DATA->seed;
 FILE* e_file;
 //char* e_name;
 string e_name;
 if (rank ==0) e_name = "e_profile_"+numIn.str()+".dat"; 
 else e_name = "e_profile2_"+numIn.str()+".dat";
 e_file = fopen(e_name.c_str(), "w");
 //e_file = fopen(e_name, "w");
 FILE* e2_file;
 //char* e2_name = "e_x_profile.dat"; 
 string e2_name = "e_x_profile_"+numIn.str()+".dat"; 
 //e2_file = fopen(e2_name, "w");
 e2_file = fopen(e2_name.c_str(), "w");
 FILE* e3_file;
 //char* e3_name = "e_y_profile.dat"; 
 string e3_name = "e_y_profile_"+numIn.str()+".dat"; 
 //e3_file = fopen(e3_name, "w");
 e3_file = fopen(e3_name.c_str(), "w");
 FILE* e4_file;
 //char* e4_name = "e_x_y_profile.dat"; 
 string e4_name = "e_x_y_profile_"+numIn.str()+".dat"; 
 //e4_file = fopen(e4_name, "w");
 e4_file = fopen(e4_name.c_str(), "w");
 FILE* geometry_file;
 //char* geometry_name = "geometry.dat"; 
 string geometry_name = "geometry_"+numIn.str()+".dat"; 
 //geometry_file = fopen(geometry_name, "w");
 geometry_file = fopen(geometry_name.c_str(), "w");

 R_A = DATA->R_A;
 a_A = DATA->a_A;
 a_short = DATA->a_short;
 a_long = DATA->a_long;
 rk_order = DATA->rk_order;
 eta_fall_off = DATA->eta_fall_off;
 eta_flat = DATA->eta_flat;
 eta0 = 0.0; /* modify this to Hirano's transvere row-on-row shift */
 epsilon0 = DATA->epsilon0; /* this is in GeV/fm^3 */
 if (initializeEntropy==0)
 epsilon0 /= hbarc; /* everything is in fm now. eps is in fm^-4 */
 
/*  for(ix=0; ix<=10; ix++) */
/*   { */
/*     fprintf(stderr,"TA(%i,b=%f)=%f\n",ix,DATA->b,InterNuPInSP(ix)); */
/*   } */

 ////Determine sigmaV = sqrt(dV), needed for sampling the values of S:
 double sigmaV;
 if(DATA->Initial_profile != 7){
   sigmaV = sqrt( (DATA->tau0 + 0.5*DATA->delta_tau) * DATA->delta_x * DATA->delta_y * DATA->delta_eta * DATA->delta_tau * (double)(DATA->fluctuatingGridFactorTau) );
 }
 else{
   sigmaV = sqrt( DATA->delta_x * DATA->delta_x * DATA->delta_y * DATA->delta_tau );
 }
 cerr << "sigmaV = " << sigmaV << endl;
 if(DATA->fluctuatingHydroFlag == 1){
   for(ix=0; ix<=DATA->nx; ix++){
     for(iy=0; iy<=DATA->ny; iy++){
       (*Lneighbor)[ix][iy][0].S = util->mtx_malloc(3,3);
       (*Rneighbor)[ix][iy][0].S = util->mtx_malloc(3,3);
       (*Lneighbor)[ix][iy][1].S = util->mtx_malloc(3,3);
       (*Rneighbor)[ix][iy][1].S = util->mtx_malloc(3,3);

       //(*Lneighbor)[ix][iy][0].prev_S = util->mtx_malloc(3,3);
       //(*Rneighbor)[ix][iy][0].prev_S = util->mtx_malloc(3,3);
       //(*Lneighbor)[ix][iy][1].prev_S = util->mtx_malloc(3,3);
       //(*Rneighbor)[ix][iy][1].prev_S = util->mtx_malloc(3,3);

       (*Lneighbor)[ix][iy][0].Xi = util->cube_malloc(rk_order+1, 4,4);
       (*Rneighbor)[ix][iy][0].Xi = util->cube_malloc(rk_order+1, 4,4);
       (*Lneighbor)[ix][iy][1].Xi = util->cube_malloc(rk_order+1, 4,4);
       (*Rneighbor)[ix][iy][1].Xi = util->cube_malloc(rk_order+1, 4,4);

       (*Lneighbor)[ix][iy][0].prev_Xi0 = util->vector_malloc(4);
       (*Rneighbor)[ix][iy][0].prev_Xi0 = util->vector_malloc(4);
       (*Lneighbor)[ix][iy][1].prev_Xi0 = util->vector_malloc(4);
       (*Rneighbor)[ix][iy][1].prev_Xi0 = util->vector_malloc(4);

       for(ieta=0; ieta < DATA->neta; ieta++){
	 (*arena)[ix][iy][ieta].S = util->mtx_malloc(3,3);
	 //(*arena)[ix][iy][ieta].prev_S = util->mtx_malloc(3,3);
	 (*arena)[ix][iy][ieta].Xi = util->cube_malloc(rk_order+1, 4,4);
	 (*arena)[ix][iy][ieta].prev_Xi0 = util->vector_malloc(4);
	 //Initializing S:
	 (*arena)[ix][iy][ieta].initS(gsl_ran, sigmaV, DATA);
       }
     }
   }
 }	     
 
 if(DATA->fluctuatingHydroFlag == 2){
   for(ix=0; ix<=DATA->nx; ix++){
     for(iy=0; iy<=DATA->ny; iy++){
       //cerr << ix << " " << iy << endl;
       (*Lneighbor)[ix][iy][0].S = util->mtx_malloc(3,3);
       (*Rneighbor)[ix][iy][0].S = util->mtx_malloc(3,3);
       (*Lneighbor)[ix][iy][1].S = util->mtx_malloc(3,3);
       (*Rneighbor)[ix][iy][1].S = util->mtx_malloc(3,3);

       (*Lneighbor)[ix][iy][0].deltaU = util->mtx_malloc(rk_order+1, 4);
       (*Rneighbor)[ix][iy][0].deltaU = util->mtx_malloc(rk_order+1, 4);
       (*Lneighbor)[ix][iy][1].deltaU = util->mtx_malloc(rk_order+1, 4);
       (*Rneighbor)[ix][iy][1].deltaU = util->mtx_malloc(rk_order+1, 4);

       (*Lneighbor)[ix][iy][0].prev_deltaU = util->vector_malloc(4);
       (*Rneighbor)[ix][iy][0].prev_deltaU = util->vector_malloc(4);
       (*Lneighbor)[ix][iy][1].prev_deltaU = util->vector_malloc(4);
       (*Rneighbor)[ix][iy][1].prev_deltaU = util->vector_malloc(4);

       (*Lneighbor)[ix][iy][0].deltaP = util->vector_malloc(rk_order+1);
       (*Rneighbor)[ix][iy][0].deltaP = util->vector_malloc(rk_order+1);
       (*Lneighbor)[ix][iy][1].deltaP = util->vector_malloc(rk_order+1);
       (*Rneighbor)[ix][iy][1].deltaP = util->vector_malloc(rk_order+1);
 
       (*Lneighbor)[ix][iy][0].deltaT = util->cube_malloc(rk_order+1, 4, 4);
       (*Rneighbor)[ix][iy][0].deltaT = util->cube_malloc(rk_order+1, 4, 4);
       (*Lneighbor)[ix][iy][1].deltaT = util->cube_malloc(rk_order+1, 4, 4);
       (*Rneighbor)[ix][iy][1].deltaT = util->cube_malloc(rk_order+1, 4, 4);

       (*Lneighbor)[ix][iy][0].deltaW = util->cube_malloc(rk_order+1, 4, 4);
       (*Rneighbor)[ix][iy][0].deltaW = util->cube_malloc(rk_order+1, 4, 4);
       (*Lneighbor)[ix][iy][1].deltaW = util->cube_malloc(rk_order+1, 4, 4);
       (*Rneighbor)[ix][iy][1].deltaW = util->cube_malloc(rk_order+1, 4, 4);

       (*Lneighbor)[ix][iy][0].prev_deltaW = util->mtx_malloc(4, 4);
       (*Rneighbor)[ix][iy][0].prev_deltaW = util->mtx_malloc(4, 4);
       (*Lneighbor)[ix][iy][1].prev_deltaW = util->mtx_malloc(4, 4);
       (*Rneighbor)[ix][iy][1].prev_deltaW = util->mtx_malloc(4, 4);

       //(*Lneighbor)[ix][iy][0].Xi = util->cube_malloc(rk_order+1, 4, 4);
       //(*Rneighbor)[ix][iy][0].Xi = util->cube_malloc(rk_order+1, 4, 4);
       //(*Lneighbor)[ix][iy][1].Xi = util->cube_malloc(rk_order+1, 4, 4);
       //(*Rneighbor)[ix][iy][1].Xi = util->cube_malloc(rk_order+1, 4, 4);
       //
       //(*Lneighbor)[ix][iy][0].prev_Xi = util->mtx_malloc(4, 4);
       //(*Rneighbor)[ix][iy][0].prev_Xi = util->mtx_malloc(4, 4);
       //(*Lneighbor)[ix][iy][1].prev_Xi = util->mtx_malloc(4, 4);
       //(*Rneighbor)[ix][iy][1].prev_Xi = util->mtx_malloc(4, 4);

       for(ieta=0; ieta < DATA->neta; ieta++){
	 (*arena)[ix][iy][ieta].S = util->mtx_malloc(3,3);
	 (*arena)[ix][iy][ieta].deltaU = util->mtx_malloc(rk_order+1, 4);
	 (*arena)[ix][iy][ieta].prev_deltaU = util->vector_malloc(4);
	 (*arena)[ix][iy][ieta].deltaP = util->vector_malloc(rk_order+1);
	 (*arena)[ix][iy][ieta].deltaT = util->cube_malloc(rk_order+1, 4, 4);
	 (*arena)[ix][iy][ieta].deltaW = util->cube_malloc(rk_order+1, 4, 4);
	 (*arena)[ix][iy][ieta].prev_deltaW = util->mtx_malloc(4, 4);
	 //(*arena)[ix][iy][ieta].Xi = util->cube_malloc(rk_order+1, 4, 4);
	 //(*arena)[ix][iy][ieta].prev_Xi = util->mtx_malloc(4, 4);
	 //Initializing S:
	 (*arena)[ix][iy][ieta].initS(gsl_ran, sigmaV, DATA);
       }

       //Finally, these have trivial initial values:
       (*Lneighbor)[ix][iy][0].deltaP[0] = 0.;
       (*Lneighbor)[ix][iy][1].deltaP[0] = 0.;
       (*Rneighbor)[ix][iy][0].deltaP[0] = 0.;
       (*Rneighbor)[ix][iy][1].deltaP[0] = 0.;

       (*Lneighbor)[ix][iy][0].prev_deltaP = 0.;
       (*Lneighbor)[ix][iy][1].prev_deltaP = 0.;
       (*Rneighbor)[ix][iy][0].prev_deltaP = 0.;
       (*Rneighbor)[ix][iy][1].prev_deltaP = 0.;

       for(ieta=0; ieta<DATA->neta; ieta++){

	 (*arena)[ix][iy][ieta].deltaP[0] = 0.;
	 (*arena)[ix][iy][ieta].prev_deltaP = 0.;

       }
       for(int ii=0; ii<4; ii++){

	 (*Lneighbor)[ix][iy][0].deltaU[0][ii] = 0.;
	 (*Lneighbor)[ix][iy][1].deltaU[0][ii] = 0.;
	 (*Rneighbor)[ix][iy][0].deltaU[0][ii] = 0.;
	 (*Rneighbor)[ix][iy][1].deltaU[0][ii] = 0.;

	 (*Lneighbor)[ix][iy][0].prev_deltaU[ii] = 0.;
	 (*Lneighbor)[ix][iy][1].prev_deltaU[ii] = 0.;
	 (*Rneighbor)[ix][iy][0].prev_deltaU[ii] = 0.;
	 (*Rneighbor)[ix][iy][1].prev_deltaU[ii] = 0.;

	 for(ieta=0; ieta < DATA->neta; ieta++){
	   (*arena)[ix][iy][ieta].deltaU[0][ii] = 0.;
	   (*arena)[ix][iy][ieta].prev_deltaU[ii] = 0.;
	 }
	 for(int ij=0; ij<4; ij++){

	   (*Lneighbor)[ix][iy][0].deltaT[0][ii][ij] = 0.;
	   (*Lneighbor)[ix][iy][1].deltaT[0][ii][ij] = 0.;
	   (*Rneighbor)[ix][iy][0].deltaT[0][ii][ij] = 0.;
	   (*Rneighbor)[ix][iy][1].deltaT[0][ii][ij] = 0.;

	   (*Lneighbor)[ix][iy][0].deltaW[0][ii][ij] = 0.;
	   (*Lneighbor)[ix][iy][1].deltaW[0][ii][ij] = 0.;
	   (*Rneighbor)[ix][iy][0].deltaW[0][ii][ij] = 0.;
	   (*Rneighbor)[ix][iy][1].deltaW[0][ii][ij] = 0.;

	   (*Lneighbor)[ix][iy][0].prev_deltaW[ii][ij] = 0.;
	   (*Lneighbor)[ix][iy][1].prev_deltaW[ii][ij] = 0.;
	   (*Rneighbor)[ix][iy][0].prev_deltaW[ii][ij] = 0.;
	   (*Rneighbor)[ix][iy][1].prev_deltaW[ii][ij] = 0.;

	   //(*Lneighbor)[ix][iy][0].Xi[0][ii][ij] = 0.;
	   //(*Lneighbor)[ix][iy][1].Xi[0][ii][ij] = 0.;
	   //(*Rneighbor)[ix][iy][0].Xi[0][ii][ij] = 0.;
	   //(*Rneighbor)[ix][iy][1].Xi[0][ii][ij] = 0.;
	   
	   //(*Lneighbor)[ix][iy][0].prev_Xi[ii][ij] = 0.;
	   //(*Lneighbor)[ix][iy][1].prev_Xi[ii][ij] = 0.;
	   //(*Rneighbor)[ix][iy][0].prev_Xi[ii][ij] = 0.;
	   //(*Rneighbor)[ix][iy][1].prev_Xi[ii][ij] = 0.;

	   for(ieta=0; ieta < DATA->neta; ieta++){

	     (*arena)[ix][iy][ieta].deltaW[0][ii][ij] = 0.;

	     (*arena)[ix][iy][ieta].prev_deltaW[ii][ij] = 0.;

	     //(*arena)[ix][iy][ieta].Xi[0][ii][ij] = 0.;
	     
	     //(*arena)[ix][iy][ieta].prev_Xi[ii][ij] = 0.;
	   }
	 }
       }

     }
   }
   cerr << "OK initializing the fluctuating variables..." << endl;
 }

 if(DATA->Initial_profile==0) // sangyong's initial version
   {
     for(ix=0; ix<=DATA->nx; ix++)
       {
	 x = (DATA->delta_x)*ix - (DATA->x_size)/2.0;
	 for(iy=0; iy<=DATA->ny; iy++)
	   {
	     y = (DATA->delta_y)*iy - (DATA->y_size)/2.0;
	     rho = sqrt(x*x/(a_short*a_short) + y*y/(a_long*a_long));
	     rho *= sqrt(a_short*a_short + a_long*a_long);
	     (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
	     
	     (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     
	     //(*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     
	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     
	     //(*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     
	     for(ieta=0; ieta<DATA->neta; ieta++)
	       {
		 eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
		 
		 //cout << " on rank " << rank << " do eta = " << eta << endl;
		 
		 exparg1 = (fabs(eta-eta0) - eta_flat/2.0)/eta_fall_off;
		 exparg = exparg1*exparg1/2.0;
		 
		 epsilon = epsilon0*exp(-exparg*theta(exparg1));
		 epsilon *= 1.0/(1.0 + exp((rho-R_A)/a_A));
		 
		 rhob = (1e-6)*epsilon; /* in fm^-3 */
		 
		 p = eos->p_func(epsilon, rhob);
		 
		 if(DATA->whichEOS == 1){
		   T = eos->interpolate(epsilon, rhob, 0);
		 }
		 else {
		   T = eos->interpolate2(epsilon, rhob, 1);
		 }
		 
		 (*arena)[ix][iy][ieta].epsilon = epsilon;
		 (*arena)[ix][iy][ieta].rhob = rhob;
		 (*arena)[ix][iy][ieta].p = p;
		 (*arena)[ix][iy][ieta].trouble = 0;

		 //The temporary values must be instantiated because getXi uses these. -CFY
		 (*arena)[ix][iy][ieta].epsilon_t = epsilon;
		 (*arena)[ix][iy][ieta].rhob_t = rhob;
		 (*arena)[ix][iy][ieta].p_t = p;
		 (*arena)[ix][iy][ieta].T_t = T;

		 //Finally, the previous values are set to be the current ones at the initial step. -CFY
		 (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
		 (*arena)[ix][iy][ieta].prev_p = p;
		
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 //(*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 //(*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
		 
		 /* 0 is actual, any others are RK intermediate values */
		 /* for a shock test */
		 //       u[0] = (*arena)[ix][iy][ieta].u[0][0] = cosh(10.0);
		 //       u[3] = (*arena)[ix][iy][ieta].u[0][3] = sinh(10.0);
		 
		 /* for HIC */
		 u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
		 u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
		 u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
		 u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][2] = 0.0;
		 
		 //(*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;

		 for(mu=0; mu<4; mu++)
		   {
		     /* baryon density */
		     (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
		     
		     for(nu=0; nu<4; nu++)
		       {
			 (*arena)[ix][iy][ieta].TJb[0][nu][mu] 
			   = (epsilon + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];
		 
			 (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevWmunu[0][nu][mu] = (double) 0.0;
			 
			 //(*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
			 
		       }/* nu */
		   }/* mu */
		 
		 if(DATA->fluctuatingHydroFlag == 1){
		   //With T, epsilon, p, and rho determined, Xi can be determined as well:
		   (*arena)[ix][iy][ieta].initXi(util, eos, DATA);
		   //Initially, prev_Xi is equal to Xi[0]:
		   for(mu=0; mu<4; mu++){
		     (*arena)[ix][iy][ieta].prev_Xi0[mu] = (*arena)[ix][iy][ieta].Xi[0][0][mu];
		   }
		 }
		 
	       }}}/* ix, iy, ieta */
   }
 else if (DATA->Initial_profile==1) //full average initial conditions using Glauber
   {
     //cerr << "OK starting when DATA->Initial_profile == 1..." << endl;
     double va, test;
     //va = findRoot(&ssolve, 0., 110., 1.15*0.+0.001,1500.,0.001);
     //fprintf(stderr,"test=%f\n",ssolve(va, 0., 110.));
     //fprintf(stderr,"T(e)=%f\n",interpolate(va, 0., 0));
     //fprintf(stderr,"e=%f\n",va);
     //fprintf(stderr,"P=%f\n",interpolate_pressure(va,0.));
     
     double s, r1, r2, W, nBinary, nWounded, W0;
     // impact parameter:
     double b=DATA->b;
     int i;
     //normalization of TATarget and TAProjectile is given by Target.A and Projectie.A, respectively
     double normT = glauber->LexusData.Target.A;
     double normP = glauber->LexusData.Projectile.A;
     // the fraction of the hard (binary collisions) contribution (will be made a parameter)
     double hard;
     hard=DATA->hard;
     
//      //secret program bonus:
//      //*************************************************************************************************
//      //compute average b: (kinda rough but good enough)

//      double sigmaAAtot, sigmaAA;
//      double integral, integralb;
//      double bvar, s1, s2;
//      double TAB[100], b1, b2, ratio, deltab, deltas;
//      int ib1, ib2;
//      integral = 0.;
//      deltab = 0.025;
//      deltas = 2.;
//      cout << " computing sigmaAAtot" << endl;
//      for (int ib=0; ib<800.; ib++)
//        {
// 	 bvar=ib*deltab;
// 	 TAB[ib] = 0.;
// 	 for (int ix=0; ix<20; ix++)
// 	   for (int iy=0; iy<20; iy++)
// 	     {
// 	       x = -20.+ix*deltas;
// 	       y = -20.+iy*deltas;
// 	       s1=sqrt(pow(x+bvar/2.,2.)+y*y);
// 	       s2=sqrt(pow(x-bvar/2.,2.)+y*y);
// 	       TAB[ib] += deltas*deltas*TATarget(DATA,s1)*TAProjectile(DATA,s2);
// 	     }
// 	 integral += deltab*2*PI*bvar*(1.-exp(-TAB[ib]*(DATA->SigmaNN/10.)));
// 	 cout << "b=" << bvar << ", TAB=" << TAB[ib] << ", " << deltab*2*PI*bvar*bvar*(1.-exp(-TAB[ib]*(DATA->SigmaNN/10.))) << " " << integral << endl;
//        }
//      sigmaAAtot = integral;

// //      for (ib2=1; ib2<50; ib2++)
// //        {
// // 	 integral = 0.;
// // 	 integralb = 0.;
// // 	 b2 = ib2*deltab;
// // 	 cout << "b_2=" << b2 << " fm" << endl;
// // 	 for (int ib=0; ib<ib2; ib++)
// // 	   {
// // 	     bvar=ib*deltab;
// // 	     integral += deltab*2*PI*bvar*(1.-exp(-TAB[ib]*(DATA->SigmaNN/10.)));
// // 	     integralb += deltab*2*PI*bvar*bvar*(1.-exp(-TAB[ib]*(DATA->SigmaNN/10.)));
// // 	   }
// // 	 sigmaAA = integral;
// // 	 ratio = sigmaAA/sigmaAAtot;
// // 	 cout << "percentage=" << ratio*100 << " % at ib2= " << ib2 << endl;
// // 	 cout << "<b>=" << integralb/sigmaAA << " fm" << endl;
// //        }

//      // fix lower limit b1 and vary upper limit b2 until the desired percentage is reached.

//      for (ib1=0; ib1<1; ib1++) // same value here for both , e.g.     for (ib1=191; ib1<192; ib1++)
//        {
// 	 b1=ib1*deltab;
// 	 for (ib2=0; ib2<500; ib2++)
// 	   {
// 	     integral = 0.;
// 	     integralb = 0.;
// 	     b2 = ib2*deltab;
// 	     if (ib1>=ib2) continue;
// 	     for (int ib=ib1; ib<ib2; ib++)
// 	       {
// 		 bvar=ib*deltab;
// 		 integral += deltab*2*PI*bvar*(1.-exp(-TAB[ib]*(DATA->SigmaNN/10.)));
// 		 integralb += deltab*2*PI*bvar*bvar*(1.-exp(-TAB[ib]*(DATA->SigmaNN/10.)));
// 	       }
// 	     sigmaAA = integral;
// 	     ratio = sigmaAA/sigmaAAtot;
// 	     cout << "percentage=" << ratio*100 << " % at ib2=" << ib2 << endl;
// 	     cout << "b1=" << b1 << " fm, b2=" << b2 << " fm, <b>=" << integralb/sigmaAA << " fm" << endl;
// 	   }
//        }



//      exit(1);


//      // *************************************************************************************************


     // The value of W at x=y=0 for normalization of W in central collisions:
     W0 = hard*(TATarget(DATA, 0.+0./2.)*TAProjectile(DATA, 0.+0./2.)*DATA->SigmaNN/10.)
       + (1.-hard)*(TATarget(DATA, 0.+0./2.)*(1.-pow((1.-(DATA->SigmaNN/10.)*TAProjectile(DATA, 0.+0./2.)/normP),normP))
		    + TAProjectile(DATA, 0.+0./2.)*(1.-pow((1.-(DATA->SigmaNN/10.)*TATarget(DATA, 0.+0./2.)/normT),normT))); 
     
     // loop over the whole lattice and initialize values:
     cerr << "OK starting the loop in init..." << endl;
     for(ix=0; ix<=DATA->nx; ix++)
       {
	 x = (DATA->delta_x)*ix - (DATA->x_size)/2.0;
	 for(iy=0; iy<=DATA->ny; iy++)
	   {
	     y = (DATA->delta_y)*iy - (DATA->y_size)/2.0;
	     if(DATA->rotateBy45degrees ==0)
	       {
		 r1 = sqrt((x+b/2.)*(x+b/2.)+(y*y));
		 r2 = sqrt((x-b/2.)*(x-b/2.)+(y*y));
	       }
	     else
	       {
		 r1 = sqrt((1/sqrt(2.)*(x+y)+b/2.)*(1/sqrt(2.)*(x+y)+b/2.)+0.5*(-x+y)*(-x+y));
		 r2 = sqrt((1/sqrt(2.)*(x+y)-b/2.)*(1/sqrt(2.)*(x+y)-b/2.)+0.5*(-x+y)*(-x+y));
	       }
	     //number of binary collisions:
	     nBinary = TAProjectile(DATA, r2)*TATarget(DATA, r1)*DATA->SigmaNN/10.;
	     //number of wounded nucleons:
	     nWounded = TATarget(DATA, r1)*(1.-pow((1.-(DATA->SigmaNN/10.)*TAProjectile(DATA, r2)/normP),normP))
	       + TAProjectile(DATA, r2)*(1.-pow((1.-(DATA->SigmaNN/10.)*TATarget(DATA, r1)/normT),normT));
	     //distribution in the transverse plane, normalized so that maximum value is 1:
	     W = ((1.-hard)*nWounded + hard*nBinary)/W0;

	     (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);

	     //(*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);

	     //(*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     int fakeEnvelope=0;
	     if(fakeEnvelope==1)
	       {
		 if (ix==0 && iy==0) cout << "WARNING: USING THE FAKE TESTING ENVELOPE!!!!!!" << endl;
		 //		 W*=(1.3+1.2*tanh((sqrt(0.001*pow(y,4.))))-1.2*tanh((sqrt(0.1*pow(x,4.))))); //(0.85+1.2*tanh((sqrt(0.001*pow(y,4.)))));
		 //	 W*=(1.5+1.4*tanh((sqrt(0.002*pow(y,4.))))-1.5*tanh((sqrt(1.*pow(x,4.))))); //(0.85+1.2*tanh((sqrt(0.001*pow(y,4.)))));
		 W*=(1.6+1.6*tanh((sqrt(0.001*pow(y,4.))))-1.6*tanh(sqrt(0.2*pow(x,4))));
	       }
	      for(ieta=0; ieta<DATA->neta; ieta++)
	       {
		 eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;

		 if ( initializeEntropy==0 )
		   {
		     exparg1 = (fabs(eta-eta0) - eta_flat/2.0)/eta_fall_off;
		     exparg = exparg1*exparg1/2.0;
		     
		     // distribution of the initial energy density in eta:
		     epsilon = epsilon0*exp(-exparg*theta(exparg1));
		     // and x,y:
		     epsilon *= W;

		     // distribution of initial baryon density:
		     rhob = DATA->rhoB0/epsilon0*epsilon; /* in fm^-3 */
		     if (x==0 && y==0) 
		       {
			 fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
			 fprintf(e_file,"%e %e\n",eta,epsilon*hbarc);
		       }
		     if (y==0 && eta==0) 
		       {
			 fprintf(e2_file,"%e %e\n",x,epsilon*hbarc);
		       }
		     if (x==0 && eta==0) 
		       {
			 fprintf(e3_file,"%e %e\n",y,epsilon*hbarc);
		       }
		   }
		 else if ( initializeEntropy==1 )
		   {
		     s = epsilon0*W;
		     exparg1 = (fabs(eta-eta0) - eta_flat/2.0)/eta_fall_off;
		     exparg = exparg1*exparg1/2.0;
		     s *= exp(-exparg*theta(exparg1));
		     rhob = DATA->rhoB0/epsilon0*s;
		     if ( DATA->whichEOS==1 && s>96.288 )
		       {
			 epsilon=(pow((4./169./pow(PI,2.)),1./3.)/30.*pow(22.5,(4./3.))*
				  pow(s,(4./3.))*0.197+0.0)/hbarc;
		       }
		     else if (s>0.2)
		       {
			 epsilon = eos->findRoot(&EOS::ssolve, rhob, s, 1.15*rhob+0.001, 300.,0.001);
		       }
		     else 
		       {
			 epsilon=0;
		       }
		     if (x==0 && y==0) 
		       {
			 fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
			 fprintf(e_file,"%e %e %e %e\n",eta,s,rhob, epsilon*hbarc);
		       }
		   }
		 // intial pressure distribution:
		 p = eos->p_func(epsilon, rhob);

		 // set all values in the grid element:
		 (*arena)[ix][iy][ieta].epsilon = epsilon;
		 (*arena)[ix][iy][ieta].rhob = rhob;
		 (*arena)[ix][iy][ieta].p = p;
		 (*arena)[ix][iy][ieta].trouble = 0;

		 if (DATA->whichEOS==1)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate(epsilon, rho, 0);
		     (*arena)[ix][iy][ieta].mu = eos->interpolate(epsilon, rho, 1);
		   }
		 else if (DATA->whichEOS==2)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate2(epsilon, rho, 1);
		     (*arena)[ix][iy][ieta].mu = 0.0;
		   }
		 
		 //The temporary values must be instantiated because getXi uses these. -CFY
		 (*arena)[ix][iy][ieta].epsilon_t = epsilon;
		 (*arena)[ix][iy][ieta].rhob_t = rhob;
		 (*arena)[ix][iy][ieta].p_t = p;
		 (*arena)[ix][iy][ieta].T_t = (*arena)[ix][iy][ieta].T;
		
		 //Finally, the previous values are set to be the current ones at the initial step. -CFY
		 (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
		 (*arena)[ix][iy][ieta].prev_p = p;
		
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 //(*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 //(*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
		 int taui;
		 //for (taui=0; taui<3; taui++)
		 //{
		 //  (*arena)[ix][iy][ieta].tauF[taui] = 0.;
		 //}
		 /* 0 is actual, any others are RK intermediate values */
		 /* for a shock test */
		 //       u[0] = (*arena)[ix][iy][ieta].u[0][0] = cosh(10.0);
		 //       u[3] = (*arena)[ix][iy][ieta].u[0][3] = sinh(10.0);
		 
		 /* for HIC */
		 u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
		 u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
		 u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
		 u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;

		 (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][2] = 0.0;
		 
		 //(*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
		 
		 for(mu=0; mu<4; mu++)
		   {
		     /* baryon density */
		     (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
		     
		     for(nu=0; nu<4; nu++)
		       {
			 (*arena)[ix][iy][ieta].TJb[0][nu][mu] 
			   = (epsilon + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];

			 (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevWmunu[0][nu][mu] = (double) 0.0;
			 
			 //(*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 
		 if(DATA->fluctuatingHydroFlag == 1){
		   //With T, epsilon, p, and rho determined, Xi can be determined as well:
		   //cerr << "OK before initXi..." << endl;
		   (*arena)[ix][iy][ieta].initXi(util, eos, DATA);
		   //cerr << "OK after initXi..." << endl;
		   //Initially, prev_Xi is equal to Xi[0]:
		   for(mu=0; mu<4; mu++){
		     (*arena)[ix][iy][ieta].prev_Xi0[mu] = (*arena)[ix][iy][ieta].Xi[0][0][mu];
		   }
		 }
	       }}
       }/* ix, iy, ieta */
     cerr << "OK initializing..." << endl;
   }   
 else if (DATA->Initial_profile==2) // this is a test scenario for testing the freeze-out surface algorithms
   {
     double va, test;
     va = eos->findRoot(&EOS::ssolve, 0., 110., 1.15*0.+0.001,1500.,0.001);
     //fprintf(stderr,"test=%f\n",ssolve(va, 0., 110.));
     //fprintf(stderr,"T(e)=%f\n",interpolate(va, 0., 0));
     //fprintf(stderr,"e=%f\n",va);
     //fprintf(stderr,"P=%f\n",interpolate_pressure(va,0.));
     
     double s, r1, r2, W, nBinary, nWounded, W0;
     // impact parameter:
     double b=DATA->b;
     int i;
         //normalization of TATarget and TAProjectile is given by Target.A and Projectie.A, respectively
     
     // loop over the whole lattice and initialize values:
     for(ix=0; ix<=DATA->nx; ix++)
       {
	 x = (DATA->delta_x)*ix - (DATA->x_size)/2.0;
	 for(iy=0; iy<=DATA->ny; iy++)
	   {
	     y = (DATA->delta_y)*iy - (DATA->y_size)/2.0;
	     (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);

	     //(*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);

	     //(*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);

	     for(ieta=0; ieta<DATA->neta; ieta++)
	       {
		 eta = (DATA->delta_eta)*ieta - (DATA->eta_size)/2.0;
		 if ((iy)<DATA->ny/2+1)
		   epsilon=DATA->epsilonFreeze/hbarc+1.;
		 else
		   epsilon = DATA->epsilonFreeze/hbarc-1.;
		 
		 // intial pressure distribution:
		 p = eos->p_func(epsilon, rhob);
		 
		 // set all values in the grid element:
		 (*arena)[ix][iy][ieta].epsilon = epsilon;
		 (*arena)[ix][iy][ieta].rhob = rhob;
		 (*arena)[ix][iy][ieta].p = p;
		 (*arena)[ix][iy][ieta].trouble = 0;

		 if (DATA->whichEOS==1)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate(epsilon, rho, 0);
		     (*arena)[ix][iy][ieta].mu = eos->interpolate(epsilon, rho, 1);
		   }
		 else if (DATA->whichEOS==2)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate2(epsilon, rho, 1);
		     (*arena)[ix][iy][ieta].mu = 0.0;
		   }
	
		 //The temporary values must be instantiated because getXi uses these. -CFY
		 (*arena)[ix][iy][ieta].epsilon_t = epsilon;
		 (*arena)[ix][iy][ieta].rhob_t = rhob;
		 (*arena)[ix][iy][ieta].p_t = p;
		 (*arena)[ix][iy][ieta].T_t = (*arena)[ix][iy][ieta].T;
		
		 //Finally, the previous values are set to be the current ones at the initial step. -CFY
		 (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
		 (*arena)[ix][iy][ieta].prev_p = p;
		
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 //(*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 //(*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
		 int taui;
		 //for (taui=0; taui<3; taui++)
		 //{
		 //  (*arena)[ix][iy][ieta].tauF[taui] = 0.;
		 //}
		 /* 0 is actual, any others are RK intermediate values */
		 /* for a shock test */
		 //       u[0] = (*arena)[ix][iy][ieta].u[0][0] = cosh(10.0);
		 //       u[3] = (*arena)[ix][iy][ieta].u[0][3] = sinh(10.0);
		 
		 /* for HIC */
		 u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
		 u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
		 u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
		 u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][2] = 0.0;
		 
		 //(*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
      		 
		 for(mu=0; mu<4; mu++)
		   {
		     /* baryon density */
		     (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
		     
		     for(nu=0; nu<4; nu++)
		       {
			 (*arena)[ix][iy][ieta].TJb[0][nu][mu] 
			   = (epsilon + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];
			 (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevWmunu[0][nu][mu] = (double) 0.0;
			 
			 //(*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 
		 if(DATA->fluctuatingHydroFlag == 1){
		   //With T, epsilon, p, and rho determined, Xi can be determined as well:
		   (*arena)[ix][iy][ieta].initXi(util, eos, DATA);
		   //Initially, prev_Xi is equal to Xi[0]:
		   for(mu=0; mu<4; mu++){
		     (*arena)[ix][iy][ieta].prev_Xi0[mu] = (*arena)[ix][iy][ieta].Xi[0][0][mu];
		   }
		 }
		 
	       }}}/* ix, iy, ieta */
   }
 else if (DATA->Initial_profile==3) // sample initial distribution with a Glauber Monte-Carlo
   {

     //keep this in mind: from Nara and Hirano: 0904.2966
//      In our MC-Glauber model, we find the defaultWoods-Saxon
//      distribution is reproduced by a larger radius parameter
//      r0 = 6.42 (4.28) fm and a smaller diffuseness parameter
//      d = 0.44 (0.50) fm (i.e., sharper boundary of a nucleus)
//      for a gold (copper) nucleus than the default parameters.
     
     // sample Nu (which gives the number of nucleons at radial position r) for both nuclei
     // for now use the same nucleus for both
     //double rnum = 5;
     double rnum=time(0)+DATA->seed; // initialize random seed using the system time
     cout << "random seed=" << rnum-DATA->seed << " + " << DATA->seed << " = " << rnum << endl;
     random->init_genrand64(rnum); 
     double width = DATA->sigma0; // width of Gaussian in fm around the center of the nucleon-nucleon collision;
     double norm = 0.2357518037*(0.4/width)*(0.4/width); // norm to get the right average energy density in the center. can be determined bu using option Initial_Profile = 4
     double A;
     double Z;
     double bArray[1];
     double AxArray[300], AyArray[300]; // arrays for sending when using MPI
     double BxArray[300], ByArray[300];
     ReturnValue rv, rv2; //auxiliary nucleons
     double eccentricity2, Psi2;
     double eccentricity3, Psi3;
     double x, y, xm, ym, dx, dy, dij;
     double xA, yA, rA;
     double phiA;
     double avrSq;
     int i;
     double b;

     // minimal and maximal impact parameter in fm
     double bmin = DATA->bmin;
     double bmax = DATA->bmax;
   

     // set up number of nucleons, free arrays with nucleons in nucleus A and B
     A=glauber->nucleusA();
     if (A>300) 
       {
	 cerr << "A>300, change the code to make it work" << endl;
	 exit(1);
       }
     nucleusA.clear();
     nucleusB.clear();


     // sample the distribution of nucleons when on rank 0, send the result to the other ranks
     if (rank==0) //let only one rank do the sampling and then pass the same distribution to all other ranks
       {
	 // first sample a b from the distribution
	 // b db = 2b/(b_max^2-b_min^2)
	 double xb = random->genrand64_real1(); // uniformly distributed random variable
	 b = sqrt((bmax*bmax-bmin*bmin)*xb+bmin*bmin);
	 // that's all, we have b.
	 // find a testing routine at the end of this file
	 sampleTA();                                  // populate the lists nucleusA and nucleusB with position data of the nucleons
	 for (i = 0; i<A; i++) 
	   {
	     bArray[0] = b;
	     AxArray[i] = nucleusA.at(i).x;
	     AyArray[i] = nucleusA.at(i).y;
	     BxArray[i] = nucleusB.at(i).x;
	     ByArray[i] = nucleusB.at(i).y;
	   }

	 for(i=1; i<size; i++) // send to all other ranks
	   {
	     MPI::COMM_WORLD.Send(bArray,1,MPI::DOUBLE,i,5);
	     MPI::COMM_WORLD.Send(AxArray,300,MPI::DOUBLE,i,1);
	     MPI::COMM_WORLD.Send(AyArray,300,MPI::DOUBLE,i,2);
	     MPI::COMM_WORLD.Send(BxArray,300,MPI::DOUBLE,i,3);
	     MPI::COMM_WORLD.Send(ByArray,300,MPI::DOUBLE,i,4);
	   }
       }

     if (rank!=0) //let only one rank do the sampling and then pass the same distribution to all other ranks
       {
	 MPI::COMM_WORLD.Recv(bArray,1,MPI::DOUBLE,0,5);
	 MPI::COMM_WORLD.Recv(AxArray,300,MPI::DOUBLE,0,1);
	 MPI::COMM_WORLD.Recv(AyArray,300,MPI::DOUBLE,0,2);
	 MPI::COMM_WORLD.Recv(BxArray,300,MPI::DOUBLE,0,3);
	 MPI::COMM_WORLD.Recv(ByArray,300,MPI::DOUBLE,0,4);
	 
       	 for (int i = 0; i<A; i++) 
	   {
	     b = bArray[0];
	     rv.x= AxArray[i];
	     rv.y = AyArray[i];
	     rv2.x = BxArray[i];
	     rv2.y = ByArray[i];

	     nucleusA.push_back(rv);
	     nucleusB.push_back(rv2);
     	   }
       }

     cout << "rank " << rank << ": Using impact parameter b=" << b << endl;

     // define the radius^2 in fm^2:
     double d2 = DATA->SigmaNN/(PI*10.);          // in fm^2
     
     // determine whether the nucleons collide, depending on their position/overlap:
     for (int i = 0; i<A; i++) // shift the nuclei's position by =/- b/2 or +b/2 respectively
       {
	 nucleusB.at(i).x=nucleusB.at(i).x+b/2.;
	 nucleusA.at(i).x=nucleusA.at(i).x-b/2.;
       }   
     
     int ibinColl = 0;      // counter for binary collisions
     int NbinColl = 0.;     // total number of binary collisions - will be counted below
     double xbinColl[500];  // x coordinate of binary collisions
     double ybinColl[500];  // y coordinate
     double WbinColl; // W binary collision scaled
     
     for (int i = 0; i<A; i++) 
       {
	 for (int j = 0 ; j<A ;j++) 
	   {
	     dx = (nucleusB.at(j).x)-(nucleusA.at(i).x); // one guy shifted to the right by b/2 the other to the left
	     dy = nucleusB.at(j).y-nucleusA.at(i).y;
	     dij = dx*dx+dy*dy;
	     if (dij < d2) 
	       {
		 nucleusB.at(j).collided=1;
		 nucleusA.at(i).collided=1;
		 xbinColl[ibinColl] = (nucleusA.at(i).x+nucleusB.at(j).x)/2.;
		 ybinColl[ibinColl] = (nucleusA.at(i).y+nucleusB.at(j).y)/2.;
		 //cout << xbinColl[ibinColl] << endl;
		 //cout << ybinColl[ibinColl] << endl;
		 ibinColl+=1;
	       }
	   }
       }

     NbinColl = ibinColl;
     cout << "number of binary collisions = " << NbinColl << endl;
     
     if(rank==0)
       {
	 // determine the eccentricity, participant reaction plane, triangularity, triangularity axis Psi_3:
	 
	 // first shift the position of all participants such that <x>=<y>=0:
	 // -----------------------------------------------------------------
	 double avx = 0.;
	 double avy = 0.;
	 int participants = 0;
	 
	 for (int i = 0; i<A; i++) 
	   {
	     if ( nucleusA.at(i).collided == 1)
	       {
		 xA = nucleusA.at(i).x;
		 yA = nucleusA.at(i).y;
		 
		 avx += xA;
		 avy += yA;
		 participants ++;
	       }
	     if ( nucleusB.at(i).collided == 1)
	       {
		 xA = nucleusB.at(i).x;
		 yA = nucleusB.at(i).y;
		 
		 avx += xA;
		 avy += yA;
		 participants ++;
	       }
	   }
	 
	 avx/=participants;
	 avy/=participants;
	 
	 cout << " average x before shift=" << avx << endl;
	 cout << " average y before shift=" << avy << endl;
	 
	 for (int i = 0; i<A; i++) 
	   {
	     nucleusA.at(i).x-=avx;
	     nucleusB.at(i).x-=avx;
	     nucleusA.at(i).y-=avy;
	     nucleusB.at(i).y-=avy;
	   }
	 // initial positions are shifted such that <x>=<y>=0 for the participants
	 // -----------------------------------------------------------------
	 
	 // test start -----------------------------------
	 avx = 0.;
	 avy = 0.;
	 
	 for (int i = 0; i<A; i++) 
	   {
	     if ( nucleusA.at(i).collided == 1)
	       {
		 xA = nucleusA.at(i).x;
		 yA = nucleusA.at(i).y;
		 
		 avx += xA;
		 avy += yA;
	       }
	     if ( nucleusB.at(i).collided == 1)
	       {
		 xA = nucleusB.at(i).x;
		 yA = nucleusB.at(i).y;
		 
		 avx += xA;
		 avy += yA;
	       }
	   }
	 
	 avx/=participants;
	 avy/=participants;
	 cout << " average x after shift=" << avx << endl;
	 cout << " average y after shift=" << avy << endl;
	 // test end -------------------------------------
	 
	 
	 // next compute the eccentricities, triangularity and axes
	 // I use Alver and Roland's definition of ecc_3 and Psi_3, which is questioned by Qin who wants r^2 replaced by r^3 in there
	 // -----------------------------------------------------------------
	 avrSq = 0.;
	 double avcos = 0.;
	 double avsin = 0.;
	 double avcos3 = 0.;
	 double avsin3 = 0.;
	 
	 for (int i = 0; i<A; i++) 
	   {
	     if ( nucleusA.at(i).collided == 1)
	       {
		 xA = nucleusA.at(i).x;
		 yA = nucleusA.at(i).y;
		 rA = sqrt( xA*xA + yA*yA );
		 phiA = atan(yA/xA);
		 
		 avrSq += rA*rA; // compute average r^2
		 avcos += rA*rA*cos(2.*phiA);
		 avsin += rA*rA*sin(2.*phiA);
		 avcos3 += rA*rA*cos(3.*phiA);
		 avsin3 += rA*rA*sin(3.*phiA);
	       }
	     
	     if ( nucleusB.at(i).collided == 1)
	       {
		 xA = nucleusB.at(i).x;
		 yA = nucleusB.at(i).y;
		 rA = sqrt( xA*xA + yA*yA );
		 phiA = atan(yA/xA);
		 
		 avrSq += rA*rA; // compute average r^2
		 avcos += rA*rA*cos(2.*phiA);
		 avsin += rA*rA*sin(2.*phiA);
		 avcos3 += rA*rA*cos(3.*phiA);
		 avsin3 += rA*rA*sin(3.*phiA);
	       }
	   }
	 
	 // compute average r^2 and average terms in the numerator:
	 avrSq/=participants;
	 avcos/=participants;
	 avsin/=participants;
	 avcos3/=participants;
	 avsin3/=participants;
	 
	 Psi2 = (atan2(avsin, avcos)+PI)/2.;
	 eccentricity2 = sqrt(avcos*avcos+avsin*avsin)/avrSq;
	 Psi3 = (atan2(avsin3, avcos3)+PI)/3.;
	 eccentricity3 = sqrt(avcos3*avcos3+avsin3*avsin3)/avrSq;
	 
	 cout << " N-part = " << participants << endl;
	 cout << " the participant eccentricity is eps_2=" << eccentricity2 << endl;
	 cout << " the minor axis of the ellipse defined by the participants is Psi_2=" << Psi2 << endl;
	 cout << " the participant triangularity is eps_3=" << eccentricity3 << endl;
	 cout << " the minor axis of the participant triangularity is Psi_3=" << Psi3 << endl;
	 
	 // now I need to save Psi_2 and Psi_3 (and also the eccentricities) to compute v_2 and v_3 later
	 // v_2 = <cos(2(phi-Psi_2))>
	 // v_3 = <cos(3(phi-Psi_3))>
	 
	 fprintf(geometry_file,"%e %e %e %e\n",eccentricity2, Psi2, eccentricity3, Psi3);
       } 
       
     
     //make a distribution function from the initial positions by putting Gaussians on top of each other
     //use positions of wounded nucleons only:

     double phipart, rpart; // angular coordinates of the participants
     double va, test;
     double s, r1, r2, W, nBinary, nWounded, W0;
     // impact parameter:
     //normalization of TATarget and TAProjectile is given by Target.A and Projectie.A, respectively
     double normT = glauber->LexusData.Target.A;
     double normP = glauber->LexusData.Projectile.A;
     // the fraction of the hard (binary collisions) contribution (will be made a parameter)
     double hard;
     hard=DATA->hard;

     // loop over the whole lattice and initialize values:
     for(ix=0; ix<=DATA->nx; ix++)
       {
	 x = (DATA->delta_x)*ix - (DATA->x_size)/2.0;
	 for(iy=0; iy<=DATA->ny; iy++)
	   {
	     y = (DATA->delta_y)*iy - (DATA->y_size)/2.0;
	     
	     W = 0.;
	     WbinColl = 0.;
	     
	     int wns=0;
	     for (int i = 0; i<A; i++) 
	       {
		 if (nucleusB.at(i).collided==1)
		   {
		     xm = nucleusB.at(i).x;
		     ym = nucleusB.at(i).y;
		     W += norm*exp((-(x-xm)*(x-xm)-(y-ym)*(y-ym))/(2.*width*width)); // normalize later to energy density I need
		   }
		 if (nucleusA.at(i).collided==1)
		   {
		     xm = nucleusA.at(i).x;
		     ym = nucleusA.at(i).y;
		     W += norm*exp((-(x-xm)*(x-xm)-(y-ym)*(y-ym))/(2.*width*width)); // normalize later to energy density I need
		   }
	       }
		
	     for (int i = 0; i<NbinColl; i++) 
	       {
		   {
		     xm = xbinColl[i];
		     ym = ybinColl[i];
		     WbinColl += norm*exp((-(x-xm)*(x-xm)-(y-ym)*(y-ym))/(2.*width*width)); // normalize later to energy density I need
		   }
	       }
	     
	     double Wfull = hard * WbinColl + (1.-hard) * W;
	   
	     W= Wfull;
	     if (hard==0.25) W = Wfull*0.645;

	     (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);

	     //(*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);

	     //(*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	 	 
	     epsilon = epsilon0*W;
	     
	  //    int check[size];
// 	     check[rank] = 0;
	     
// 	     if (epsilon*hbarc>311.)
// 	       {
// 		 if (rank==0) 
// 		   {		
// 		     //communicate so that rank 0 does not run away by itself to get everything stuck...
// 		     for(int i = 1; i<size; i++)
// 		       {
// 			 MPI::COMM_WORLD.Recv(check,size,MPI::INT,i,1);
// 			 if ( check[i]>0 )
// 			   {
// 			     cout << "check[" << i << "]=" << check[i] << endl;
// 			   }
// 			 cout << "i=" << i << endl;
// 		       }
// 		     for(int i = 1; i<size; i++)
// 		       {
// 			 if (check[i]>0)
// 			   return 0;
// 		       }
// 		   }
		 
// 		 if (rank!=0) 
// 		   {
// 		     cout << "energy density out of range for equation of state..." << endl;
// 		     cout << x << " " << y << " eps=" << epsilon*hbarc << endl;
// 		     cout << "redoing initialization on rank " << rank << endl;
		     
// 		     check[rank] = rank;
// 		     MPI::COMM_WORLD.Send(check,size,MPI::INT,0,1);
// 		     return 0;
// 		   }
// 	       }
	   
		     
	     for(ieta=0; ieta<DATA->neta; ieta++)
	       {
		 eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
		 if ( initializeEntropy==0 )
		   {
		     exparg1 = (fabs(eta-eta0) - eta_flat/2.0)/eta_fall_off;
		     exparg = exparg1*exparg1/2.0;
		     
		     epsilon = epsilon0*W;
				     
		     epsilon *= exp(-exparg*theta(exparg1));


		     // distribution of initial baryon density:
		     rhob = DATA->rhoB0/epsilon0*epsilon; /* in fm^-3 */
		    
		     if (x==0 && y==0) 
		       {
			 fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
			 fprintf(e_file,"%e %e\n",eta,epsilon*hbarc);
		       }
		     if (y==0 && eta==0) 
		       {
			 fprintf(e2_file,"%e %e\n",x,epsilon*hbarc);
		       }
		     if (x==0 && eta==0) 
		       {
			 fprintf(e3_file,"%e %e\n",y,epsilon*hbarc);
		       }
		   }
		 else if ( initializeEntropy==1 )
		   {
		     s = epsilon0*W;
		     exparg1 = (fabs(eta-eta0) - eta_flat/2.0)/eta_fall_off;
		     exparg = exparg1*exparg1/2.0;
		     s *= exp(-exparg*theta(exparg1));
		     rhob = DATA->rhoB0/epsilon0*s;
		     if ( DATA->whichEOS==1 && s>96.288 )
		       {
			 epsilon=(pow((4./169./pow(PI,2.)),1./3.)/30.*pow(22.5,(4./3.))*
				  pow(s,(4./3.))*0.197+0.0)/hbarc;
		       }
		     else if (s>0.2)
		       {
			 epsilon = eos->findRoot(&EOS::ssolve, rhob, s, 1.15*rhob+0.001, 300.,0.001);
		       }
		     else 
		       {
			 epsilon=0;
		       }
		     if (x==0 && y==0) 
		       {
			 fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
			 fprintf(e_file,"%e %e %e %e\n",eta,s,rhob, epsilon*hbarc);
		       }
		   }
		 // intial pressure distribution:
		 p = eos->p_func(epsilon, rhob);

		 // set all values in the grid element:
		 (*arena)[ix][iy][ieta].epsilon = epsilon;
		 (*arena)[ix][iy][ieta].rhob = rhob;
		 (*arena)[ix][iy][ieta].p = p;
		 (*arena)[ix][iy][ieta].trouble = 0;

		 if (DATA->whichEOS==1)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate(epsilon, rho, 0);
		     (*arena)[ix][iy][ieta].mu = eos->interpolate(epsilon, rho, 1);
		   }
		 else if (DATA->whichEOS==2)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate2(epsilon, rho, 1);
		     (*arena)[ix][iy][ieta].mu = 0.0;
		   }
		 
		 //The temporary values must be instantiated because getXi uses these. -CFY
		 (*arena)[ix][iy][ieta].epsilon_t = epsilon;
		 (*arena)[ix][iy][ieta].rhob_t = rhob;
		 (*arena)[ix][iy][ieta].p_t = p;
		 (*arena)[ix][iy][ieta].T_t = (*arena)[ix][iy][ieta].T;
		
		 //Finally, the previous values are set to be the current ones at the initial step. -CFY
		 (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
		 (*arena)[ix][iy][ieta].prev_p = p;
		
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 //(*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 //(*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
				
		 int taui;
		 //for (taui=0; taui<3; taui++)
		 //{
		 //  (*arena)[ix][iy][ieta].tauF[taui] = 0.;
		 //}
		 /* 0 is actual, any others are RK intermediate values */
		 /* for a shock test */
		 //       u[0] = (*arena)[ix][iy][ieta].u[0][0] = cosh(10.0);
		 //       u[3] = (*arena)[ix][iy][ieta].u[0][3] = sinh(10.0);
		 
		 /* for HIC */
		 u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
		 u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
		 u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
		 u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;

		 (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][2] = 0.0;
		 
		 //(*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
		 
		 for(mu=0; mu<4; mu++)
		   {
		     /* baryon density */
		     (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
		     
		     for(nu=0; nu<4; nu++)
		       {
			 (*arena)[ix][iy][ieta].TJb[0][nu][mu] 
			   = (epsilon + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];

			 (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevWmunu[0][nu][mu] = (double) 0.0;
			 
			 //(*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 if(DATA->fluctuatingHydroFlag == 1){
		   //With T, epsilon, p, and rho determined, Xi can be determined as well:
		   (*arena)[ix][iy][ieta].initXi(util, eos, DATA);
		   //Initially, prev_Xi is equal to Xi[0]:
		   for(mu=0; mu<4; mu++){
		     (*arena)[ix][iy][ieta].prev_Xi0[mu] = (*arena)[ix][iy][ieta].Xi[0][0][mu];
		   }
		 }
		 
	       }}}/* ix, iy, ieta */
   }
 else if (DATA->Initial_profile==4) // sample initial distribution with a Glauber Monte-Carlo and average over them to compare to average as in 1
   {
     // sample Nu (which gives the number of nucleons at radial position r) for both nuclei
     // for now use the same nucleus for both
     int nevents = 300;
     
     double rnum=time(0);
     //double rnum = 1;
     random->init_genrand64(rnum);

     double width = 0.4; // width of Gaussian in fm around the center;
     double norm = 0.2357518037;
    
     double A;
     double Z;
     
     double AxArray[300], AyArray[300];
     double BxArray[300], ByArray[300];
     
     ReturnValue rv, rv2;
     
     double sigma, epsSigma;
     double a, x, y, xm, ym, dx, dy, dij;
     double b = DATA->b;

     int posx;
     int posy;
     int n1=0;
     int n2=0;
     int i;
     double W;     
     double r = 1.2*pow(glauber->nucleusA(),1./3.);
     double s;
     double WSq=0.;
     A=glauber->nucleusA();
     if (A>300) 
       {
	 cerr << "A>300, change the code to make it work" << endl;
	 exit(1);
       }
     double WFull = 0.;
     for(int runs=0; runs<nevents; runs++)
       {		
	 W = 0.;
	 nucleusA.clear();
	 nucleusB.clear();
	 
	 sampleTA();                                  // populate the lists nucleusA and nucleusB with position data of the nucleons
	 
	 double d2 = DATA->SigmaNN/(PI*10.);          // in fm^2
	 
	 // for each of the A nucleons in nucleus B
	 
	 for (int i = 0; i<A; i++) // shift the nuclei's position by =/- b/2
	   {
	     nucleusB.at(i).x=nucleusB.at(i).x+b/2.;
	     nucleusA.at(i).x=nucleusA.at(i).x-b/2.;
	   }   
	 
	 for (int i = 0; i<A; i++) 
	   {
	     for (int j = 0 ; j<A ;j++) 
	       {
		 dx = (nucleusB.at(j).x)-(nucleusA.at(i).x); // one guy shifted to the right by b/2 the other to the left
		 dy = nucleusB.at(j).y-nucleusA.at(i).y;
		 dij = dx*dx+dy*dy;
		 if (dij < d2) 
		   {
		     nucleusB.at(j).collided=1;
		     nucleusA.at(i).collided=1;
		   }
	       }
	   }
	 
	 //make a distribution function from the initial positions by putting Gaussians on top of each other
	 // use positions of wounded nucleons only:
	 
	 
	 double va, test;
	 double r1, r2, nBinary, nWounded, W0;
	 // impact parameter:
	 //normalization of TATarget and TAProjectile is given by Target.A and Projectie.A, respectively
	 double normT = glauber->LexusData.Target.A;
	 double normP = glauber->LexusData.Projectile.A;
	 // the fraction of the hard (binary collisions) contribution (will be made a parameter)
	 double hard;
	 hard=DATA->hard;
	 
	 // The value of W at x=y=0 for normalization of W in central collisions:
	 W0 = hard*(TATarget(DATA, 0.+0./2.)*TAProjectile(DATA, 0.+0./2.)*DATA->SigmaNN/10.)
	   + (1.-hard)*(TATarget(DATA, 0.+0./2.)*(1.-pow((1.-(DATA->SigmaNN/10.)*TAProjectile(DATA, 0.+0./2.)/normP),normP))
			+ TAProjectile(DATA, 0.+0./2.)*(1.-pow((1.-(DATA->SigmaNN/10.)*TATarget(DATA, 0.+0./2.)/normT),normT))); 
	 
	 // loop over the whole lattice and initialize values:
	 
	 
	 ix = floor(static_cast<double>(DATA->nx)/2.+0.000001);
	 iy = floor(static_cast<double>(DATA->ny)/2.+0.000001);
	 
	 x = 0.;
	 y = 0.;
	 eta = 0.;
	 
	 for (int i = 0; i<A; i++) 
	   {
	     if (nucleusB.at(i).collided==1)
	       {
		 xm = nucleusB.at(i).x;
		 ym = nucleusB.at(i).y;
		 W += norm*exp((-(x-xm)*(x-xm)-(y-ym)*(y-ym))/(2.*width*width)); // normalize later to energy density I need
	       }
	     if (nucleusA.at(i).collided==1)
	       {
		 xm = nucleusA.at(i).x;
		 ym = nucleusA.at(i).y;
		 W += norm*exp((-(x-xm)*(x-xm)-(y-ym)*(y-ym))/(2.*width*width)); // normalize later to energy density I need
	       }
	   }
	 WFull += W;
	 WSq += W*W;
       }

     WFull/=static_cast<double>(nevents);
     WSq/=static_cast<double>(nevents);

     sigma = sqrt(WSq-WFull*WFull)/sqrt(static_cast<double>(nevents));

     cout << "WFull=" << WFull << endl;
     cout << "sigma=" << sigma << endl;

     if ( initializeEntropy==0 )
       {
	 exparg1 = (fabs(eta-eta0) - eta_flat/2.0)/eta_fall_off;
	 exparg = exparg1*exparg1/2.0;
	 
	 // distribution of the initial energy density in eta:
	 epsilon = epsilon0*exp(-exparg*theta(exparg1));
	 // and x,y:
	 epsilon *= WFull;
	 epsSigma = epsilon / WFull *sigma;
       }
     else if ( initializeEntropy==1 )
       {
	 s = epsilon0*WFull;
	 exparg1 = (fabs(eta-eta0) - eta_flat/2.0)/eta_fall_off;
	 exparg = exparg1*exparg1/2.0;
	 s *= exp(-exparg*theta(exparg1));
	 rhob = DATA->rhoB0/epsilon0*s;
	 if ( DATA->whichEOS==1 && s>96.288 )
	   {
	     epsilon=(pow((4./169./pow(PI,2.)),1./3.)/30.*pow(22.5,(4./3.))*
		      pow(s,(4./3.))*0.197+0.0)/hbarc;
	   }
	 else if (s>0.2)
	   {
	     epsilon = eos->findRoot(&EOS::ssolve, rhob, s, 1.15*rhob+0.001, 300.,0.001);
	   }
	 else 
	   {
	     epsilon=0;
	   }
       }
 
     cout << "the desired initial epsilon0 = " << DATA->epsilon0 << " GeV/fm^3" << endl;
     cout << "the average central epsilon = " << epsilon*hbarc << " GeV/fm^3" << endl;
     cout << "the standard deviation is = " << epsSigma*hbarc << " GeV/fm^3" << endl;
     cout << "that is a " <<sigma/WFull*100. << " % error" << endl; 
     
     cout << "Multiply 'norm' by " << DATA->epsilon0/(epsilon*hbarc) << " to get the correct average energy density." << endl;
   
     cout << "This routine only computes the scaling factor for the energy distribution and does not run hydro. Change Initial_Distibution to 1 or 3 to run the hydro evolution"
	  << endl;

     
     exit(1);
   }
 else if (DATA->Initial_profile==5) //something like pp
   {
     double va, test;
     double s, r1, r2, W, nBinary, nWounded, W0;
     // impact parameter:
     double b=DATA->b;
     int i;
     //normalization of TATarget and TAProjectile is given by Target.A and Projectie.A, respectively
     double normT = glauber->LexusData.Target.A;
     double normP = glauber->LexusData.Projectile.A;
     // the fraction of the hard (binary collisions) contribution (will be made a parameter)
     double hard;
     hard=DATA->hard;

     // loop over the whole lattice and initialize values:
     for(ix=0; ix<=DATA->nx; ix++)
       {
	 x = (DATA->delta_x)*ix - (DATA->x_size)/2.0;
	 for(iy=0; iy<=DATA->ny; iy++)
	   {
	     y = (DATA->delta_y)*iy - (DATA->y_size)/2.0;
	  
	     //number of binary collisions:
	     nBinary = TAProjectile(DATA, r2)*TATarget(DATA, r1)*DATA->SigmaNN/10.;
	     //number of wounded nucleons:
	     nWounded = TATarget(DATA, r1)*(1.-pow((1.-(DATA->SigmaNN/10.)*TAProjectile(DATA, r2)/normP),normP))
	       + TAProjectile(DATA, r2)*(1.-pow((1.-(DATA->SigmaNN/10.)*TATarget(DATA, r1)/normT),normT));
	     //distribution in the transverse plane, normalized so that maximum value is 1:
	     
	     W = exp(-x*x/(2.*0.25*0.25)-y*y/(2.*0.25*0.25))/(2.*PI*0.25*0.25);

	     (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);

	     //(*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);

	     //(*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	   
	     for(ieta=0; ieta<DATA->neta; ieta++)
	       {
		 eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
		 if ( initializeEntropy==0 )
		   {
		     exparg1 = (fabs(eta-eta0) - eta_flat/2.0)/eta_fall_off;
		     exparg = exparg1*exparg1/2.0;
		     
		     // distribution of the initial energy density in eta:
		     epsilon = epsilon0*exp(-exparg*theta(exparg1));
		     // and x,y:
		     epsilon *= W;

		     // distribution of initial baryon density:
		     rhob = DATA->rhoB0/epsilon0*epsilon; /* in fm^-3 */
		     if (x==0 && y==0) 
		       {
			 fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
			 fprintf(e_file,"%e %e\n",eta,epsilon*hbarc);
		       }
		     if (y==0 && eta==0) 
		       {
			 fprintf(e2_file,"%e %e\n",x,epsilon*hbarc);
		       }
		     if (x==0 && eta==0) 
		       {
			 fprintf(e3_file,"%e %e\n",y,epsilon*hbarc);
		       }
		   }
		 else if ( initializeEntropy==1 )
		   {
		     cout << "not implemented" << endl;
		     exit(1);
		   }
		 // intial pressure distribution:
		 p = eos->p_func(epsilon, rhob);

		 // set all values in the grid element:
		 (*arena)[ix][iy][ieta].epsilon = epsilon;
		 (*arena)[ix][iy][ieta].rhob = rhob;
		 (*arena)[ix][iy][ieta].p = p;
		 (*arena)[ix][iy][ieta].trouble = 0;

		 if (DATA->whichEOS==1)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate(epsilon, rho, 0);
		     (*arena)[ix][iy][ieta].mu = eos->interpolate(epsilon, rho, 1);
		   }
		 else if (DATA->whichEOS==2)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate2(epsilon, rho, 1);
		     (*arena)[ix][iy][ieta].mu = 0.0;
		   }
		 
		 //The temporary values must be instantiated because getXi uses these. -CFY
		 (*arena)[ix][iy][ieta].epsilon_t = epsilon;
		 (*arena)[ix][iy][ieta].rhob_t = rhob;
		 (*arena)[ix][iy][ieta].p_t = p;
		 (*arena)[ix][iy][ieta].T_t = (*arena)[ix][iy][ieta].T;
		
		 //Finally, the previous values are set to be the current ones at the initial step. -CFY
		 (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
		 (*arena)[ix][iy][ieta].prev_p = p;
		
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 //(*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 //(*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
		 int taui;
		 //for (taui=0; taui<3; taui++)
		 //{
		 //  (*arena)[ix][iy][ieta].tauF[taui] = 0.;
		 //}
		 /* 0 is actual, any others are RK intermediate values */
		 /* for a shock test */
		 //       u[0] = (*arena)[ix][iy][ieta].u[0][0] = cosh(10.0);
		 //       u[3] = (*arena)[ix][iy][ieta].u[0][3] = sinh(10.0);
		 
		 /* for HIC */
		 u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
		 u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
		 u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
		 u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;

		 (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][2] = 0.0;
		 
		 //(*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
		 
		 for(mu=0; mu<4; mu++)
		   {
		     /* baryon density */
		     (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
		     
		     for(nu=0; nu<4; nu++)
		       {
			 (*arena)[ix][iy][ieta].TJb[0][nu][mu] 
			   = (epsilon + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];

			 (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevWmunu[0][nu][mu] = (double) 0.0;
			 
			 //(*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 
		 if(DATA->fluctuatingHydroFlag == 1){
		   //With T, epsilon, p, and rho determined, Xi can be determined as well:
		   (*arena)[ix][iy][ieta].initXi(util, eos, DATA);
		   //Initially, prev_Xi is equal to Xi[0]:
		   for(mu=0; mu<4; mu++){
		     (*arena)[ix][iy][ieta].prev_Xi0[mu] = (*arena)[ix][iy][ieta].Xi[0][0][mu];
		   }
		 }
		 
	       }}}/* ix, iy, ieta */
   }   
 else if (DATA->Initial_profile==6) //An energy density uniform in each cell
   {
     // Loop over the whole lattice and initialize values:
     for(ix=0; ix<=DATA->nx; ix++)
       {
	 //x = (DATA->delta_x)*ix - (DATA->x_size)/2.0;
	 for(iy=0; iy<=DATA->ny; iy++)
	   {
	     //y = (DATA->delta_y)*iy - (DATA->y_size)/2.0;
	  
	     ////number of binary collisions:
	     //nBinary = TAProjectile(DATA, r2)*TATarget(DATA, r1)*DATA->SigmaNN/10.;
	     ////number of wounded nucleons:
	     //nWounded = TATarget(DATA, r1)*(1.-pow((1.-(DATA->SigmaNN/10.)*TAProjectile(DATA, r2)/normP),normP))
	     //+ TAProjectile(DATA, r2)*(1.-pow((1.-(DATA->SigmaNN/10.)*TATarget(DATA, r1)/normT),normT));
	     ////distribution in the transverse plane, normalized so that maximum value is 1:
	     //
	     //W = exp(-x*x/(2.*0.25*0.25)-y*y/(2.*0.25*0.25))/(2.*PI*0.25*0.25);

	     (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);

	     //(*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);

	     //(*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	   
	     for(ieta=0; ieta<DATA->neta; ieta++)
	       {
		 //eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
		 //exparg1 = (fabs(eta-eta0) - eta_flat/2.0)/eta_fall_off;
		 //exparg = exparg1*exparg1/2.0;
		 
		 // distribution of the initial energy density in eta:
		 epsilon = epsilon0;
		 //// and x,y:
		 //epsilon *= W;
		 
		 // distribution of initial baryon density:
		 rhob = DATA->rhoB0/epsilon0*epsilon0; /* in fm^-3 */
		 if (x==0 && y==0) 
		   {
		     fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
		     fprintf(e_file,"%e %e\n",eta,epsilon*hbarc);
		   }
		 if (y==0 && eta==0) 
		   {
		     fprintf(e2_file,"%e %e\n",x,epsilon*hbarc);
		   }
		 if (x==0 && eta==0) 
		   {
		     fprintf(e3_file,"%e %e\n",y,epsilon*hbarc);
		   }

		 // intial pressure distribution:
		 p = eos->p_func(epsilon0, rhob);

		 // set all values in the grid element:
		 (*arena)[ix][iy][ieta].epsilon = epsilon0;
		 (*arena)[ix][iy][ieta].rhob = rhob;
		 (*arena)[ix][iy][ieta].p = p;
		 (*arena)[ix][iy][ieta].trouble = 0;

		 if (DATA->whichEOS==1)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate(epsilon0, rhob, 0);
		     (*arena)[ix][iy][ieta].mu = eos->interpolate(epsilon0, rhob, 1);
		   }
		 else if (DATA->whichEOS==2)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate2(epsilon0, rhob, 1);
		     (*arena)[ix][iy][ieta].mu = 0.0;
		   }
		 
		 //The temporary values must be instantiated because getXi uses these. -CFY
		 (*arena)[ix][iy][ieta].epsilon_t = epsilon;
		 (*arena)[ix][iy][ieta].rhob_t = rhob;
		 (*arena)[ix][iy][ieta].p_t = p;
		 (*arena)[ix][iy][ieta].T_t = (*arena)[ix][iy][ieta].T;
		
		 //Finally, the previous values are set to be the current ones at the initial step. -CFY
		 (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
		 (*arena)[ix][iy][ieta].prev_p = p;
		
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 //(*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 //(*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
		 int taui;
		 //for (taui=0; taui<3; taui++)
		 //{
		 //  (*arena)[ix][iy][ieta].tauF[taui] = 0.;
		 //}
		 /* 0 is actual, any others are RK intermediate values */
		 /* for a shock test */
		 //       u[0] = (*arena)[ix][iy][ieta].u[0][0] = cosh(10.0);
		 //       u[3] = (*arena)[ix][iy][ieta].u[0][3] = sinh(10.0);
		 
		 /* for HIC */
		 u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
		 u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
		 u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
		 u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;

		 (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][2] = 0.0;
		 
		 //(*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
		 
		 for(mu=0; mu<4; mu++)
		   {
		     /* baryon density */
		     (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
		     
		     for(nu=0; nu<4; nu++)
		       {
			 (*arena)[ix][iy][ieta].TJb[0][nu][mu] 
			   = (epsilon0 + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];

			 (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevWmunu[0][nu][mu] = (double) 0.0;
			 
			 //(*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 
		 if(DATA->fluctuatingHydroFlag == 1){
		   //With T, epsilon, p, and rho determined, Xi can be determined as well:
		   (*arena)[ix][iy][ieta].initXi(util, eos, DATA);
		   //Initially, prev_Xi is equal to Xi[0]:
		   for(mu=0; mu<4; mu++){
		     (*arena)[ix][iy][ieta].prev_Xi0[mu] = (*arena)[ix][iy][ieta].Xi[0][0][mu];
		   }
		 }
		 
	       }}}/* ix, iy, ieta */
   }
 else if (DATA->Initial_profile==7) //An energy density uniform in each cell with no fluid flow in the rest frame:
   {
     // Loop over the whole lattice and initialize values:
     for(ix=0; ix<=DATA->nx; ix++)
       {
	 //x = (DATA->delta_x)*ix - (DATA->x_size)/2.0;
	 for(iy=0; iy<=DATA->ny; iy++)
	   {
	     //y = (DATA->delta_y)*iy - (DATA->y_size)/2.0;
	  
	     ////number of binary collisions:
	     //nBinary = TAProjectile(DATA, r2)*TATarget(DATA, r1)*DATA->SigmaNN/10.;
	     ////number of wounded nucleons:
	     //nWounded = TATarget(DATA, r1)*(1.-pow((1.-(DATA->SigmaNN/10.)*TAProjectile(DATA, r2)/normP),normP))
	     //+ TAProjectile(DATA, r2)*(1.-pow((1.-(DATA->SigmaNN/10.)*TATarget(DATA, r1)/normT),normT));
	     ////distribution in the transverse plane, normalized so that maximum value is 1:
	     //
	     //W = exp(-x*x/(2.*0.25*0.25)-y*y/(2.*0.25*0.25))/(2.*PI*0.25*0.25);

	     (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);

	     //(*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     //(*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);

	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);

	     //(*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     //(*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     //(*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	   
	     for(ieta=0; ieta<DATA->neta; ieta++)
	       {
		 eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
		 //cerr << "eta = " << eta << endl;
		 //exparg1 = (fabs(eta-eta0) - eta_flat/2.0)/eta_fall_off;
		 //exparg = exparg1*exparg1/2.0;
		 
		 // distribution of the initial energy density in eta:
		 epsilon = epsilon0;
		 //// and x,y:
		 //epsilon *= W;
		 
		 // distribution of initial baryon density:
		 rhob = DATA->rhoB0/epsilon0*epsilon0; /* in fm^-3 */
		 if (x==0 && y==0) 
		   {
		     //fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
		     fprintf(e_file,"%e %e\n",eta,epsilon*hbarc);
		   }
		 if (y==0 && eta==0) 
		   {
		     fprintf(e2_file,"%e %e\n",x,epsilon*hbarc);
		   }
		 if (x==0 && eta==0) 
		   {
		     fprintf(e3_file,"%e %e\n",y,epsilon*hbarc);
		   }

		 // intial pressure distribution:
		 p = eos->p_func(epsilon0, rhob);

		 // set all values in the grid element:
		 (*arena)[ix][iy][ieta].epsilon = epsilon0;
		 (*arena)[ix][iy][ieta].rhob = rhob;
		 (*arena)[ix][iy][ieta].p = p;
		 (*arena)[ix][iy][ieta].trouble = 0;

		 if (DATA->whichEOS==1)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate(epsilon0, rhob, 0);
		     (*arena)[ix][iy][ieta].mu = eos->interpolate(epsilon0, rhob, 1);
		   }
		 else if (DATA->whichEOS==2)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate2(epsilon0, rhob, 1);
		     (*arena)[ix][iy][ieta].mu = 0.0;
		   }
		 
		 if(x==0 && y==0 && eta==0){
		   cerr << "T = " << eos->interpolate2(epsilon0, rhob, 1) << ", e = " << epsilon0 << ", p = " << p << ", dedp = " << 1./eos->get_dpOverde2(epsilon0, rhob) << endl;
		 }

		 //The temporary values must be instantiated because getXi uses these. -CFY
		 (*arena)[ix][iy][ieta].epsilon_t = epsilon;
		 (*arena)[ix][iy][ieta].rhob_t = rhob;
		 (*arena)[ix][iy][ieta].p_t = p;
		 (*arena)[ix][iy][ieta].T_t = (*arena)[ix][iy][ieta].T;
		
		 //Finally, the previous values are set to be the current ones at the initial step. -CFY
		 (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
		 (*arena)[ix][iy][ieta].prev_p = p;
		
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 //(*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 //(*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 //(*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 //(*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
		 int taui;
		 //for (taui=0; taui<3; taui++)
		 //{
		 //  (*arena)[ix][iy][ieta].tauF[taui] = 0.;
		 //}
		 /* 0 is actual, any others are RK intermediate values */
		 /* for a shock test */
		 //       u[0] = (*arena)[ix][iy][ieta].u[0][0] = cosh(10.0);
		 //       u[3] = (*arena)[ix][iy][ieta].u[0][3] = sinh(10.0);
		 
		 /* for HIC */
		 for(int i=0; i<=rk_order; i++){
		   //u[0] = (*arena)[ix][iy][ieta].u[i][0] = cosh(eta);
		   //u[3] = (*arena)[ix][iy][ieta].u[i][3] = -sinh(eta);
		   //u[1] = (*arena)[ix][iy][ieta].u[i][1] = 0.0;
		   //u[2] = (*arena)[ix][iy][ieta].u[i][2] = 0.0;
		   u[0] = (*arena)[ix][iy][ieta].u[i][0] = 1.;
		   u[3] = (*arena)[ix][iy][ieta].u[i][3] = 0.0;
		   u[1] = (*arena)[ix][iy][ieta].u[i][1] = 0.0;
		   u[2] = (*arena)[ix][iy][ieta].u[i][2] = 0.0;
		 }
		 
		 //(*arena)[ix][iy][ieta].prev_u[0][0] = cosh(eta);
		 //(*arena)[ix][iy][ieta].prev_u[0][3] = -sinh(eta);
		 //(*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
		 //(*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
		 //(*arena)[ix][iy][ieta].pprev_u[0][0] = cosh(eta);
		 //(*arena)[ix][iy][ieta].pprev_u[0][3] = -sinh(eta);
		 //(*arena)[ix][iy][ieta].pprev_u[0][1] = 0.0;
		 //(*arena)[ix][iy][ieta].pprev_u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][0] = 1.;
		 (*arena)[ix][iy][ieta].prev_u[0][3] = 0.;
		 (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][0] = 1.;
		 (*arena)[ix][iy][ieta].pprev_u[0][3] = 0.;
		 (*arena)[ix][iy][ieta].pprev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_u[0][2] = 0.0;

		 //(*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 //(*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
		 
		 for(mu=0; mu<4; mu++)
		   {
		     /* baryon density */
		     (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
		     
		     for(nu=0; nu<4; nu++)
		       {

			 for(int i=0; i<=rk_order; i++){
			   (*arena)[ix][iy][ieta].TJb[i][nu][mu] 
			     = (epsilon0 + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];
			   
			   (*arena)[ix][iy][ieta].Wmunu[i][nu][mu] = (double) 0.0;
			 }
			 (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevWmunu[0][nu][mu] = (double) 0.0;

			 //(*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 //(*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 
		 if(DATA->fluctuatingHydroFlag == 1){
		   //With T, epsilon, p, and rho determined, Xi can be determined as well:
		   (*arena)[ix][iy][ieta].initXi(util, eos, DATA);
		   //Initially, prev_Xi is equal to Xi[0]:
		   for(mu=0; mu<4; mu++){
		     (*arena)[ix][iy][ieta].prev_Xi0[mu] = (*arena)[ix][iy][ieta].Xi[0][0][mu];
		   }
		 }
		 
	       }}}/* ix, iy, ieta */
   }
 
 fclose(e_file);
 fclose(e2_file);
 fclose(e3_file);
 fclose(e4_file);
 fclose(geometry_file);
 cerr << "InitTJb ended successfully..." << endl;
 return 1;
}/* InitTJb */
 
double Init::TATarget(InitData *DATA, double r)
{
  return glauber->InterNuTInST(r)/(DATA->SigmaNN/10.);
}

double Init::TAProjectile(InitData *DATA, double r)
{
  return glauber->InterNuPInSP(r)/(DATA->SigmaNN/10.);
}

	 // test the b sampling:
// 	      int bpos;
// 	      int bArray[20];
// 	      for (int i=0; i<20; i++)
// 	        {
// 	 	 bArray[i]=0;
// 	        }
// 	      for (int i=0; i<1000000; i++)
// 	        {
// 	 	 xb = random->genrand64_real1(); // uniformly distributed random variable
// 	 	 b = sqrt((bmax*bmax-bmin*bmin)*xb);
// 	 	 bpos = floor(b/(bmax-bmin)*20);
// 	 	 bArray[bpos]+=1;
// 	        }
// 	      for (int i=0; i<20; i++)
// 	        {
// 	 	 b = (bmax-bmin)*static_cast<double>(i)/20.;
// 	 	 cout << b << " " << bArray[i] << endl;
// 	        }
	 // done testing b sampling - works fine

 
