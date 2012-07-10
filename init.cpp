#include "util.h"
#include "grid.h"
#include "init.h"
#include "eos.h"

using namespace std;

Init::Init(EOS *eosIn, Glauber *glauberIn)
{
  eos = new EOS;
  eos = eosIn;
  util = new Util;
  glauber = new Glauber;
  glauber = glauberIn;
  random = new Random;
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
  cout << "initArena" << endl;

  if (DATA->Initial_profile==6 || DATA->Initial_profile==7)
    {
      string dummy;
      int nx, ny, neta;
      cout << DATA->initName <<endl;
      ifstream profile(DATA->initName.c_str());
      double deta, dx, dy, dummy2;
      //     profile = fopen(DATA->initName, "r");
      // read the first line with general info
      
      profile >> dummy >> dummy >> dummy2 >> dummy >> neta >> dummy >> nx >> dummy >> ny >> dummy >> deta >> dummy >> dx >> dummy >> dy ;
      cout << "neta = " << neta << endl;
   
      //   bytes_read=fscanf(profile,"%s %s %lf %s %d %s %d %s %d %s %lf %s %lf %s %lf",
      //		       &dummy,&dummy,&dummy2,&dummy,&neta,&dummy,&nx,&dummy,&ny,&dummy,&deta,&dummy,&dx,&dummy,&dy);
      cout << "Using Initial_profile=" << DATA->Initial_profile << ". Overwriting lattice dimensions:" << endl;
      DATA->neta = neta/size;
      DATA->nx = nx-1;
      DATA->ny = ny-1;
      DATA->delta_eta = deta;
      DATA->delta_x = dx;
      DATA->delta_y = dy;
      cout << "neta=" << neta << ", nx=" << nx << ", ny=" << ny << endl;
      cout << "deta=" << DATA->delta_eta << ", dx=" << DATA->delta_x << ", dy=" << DATA->delta_y << endl;
      //      cout << DATA->neta << " " << DATA->nx << " " << DATA->ny << endl;
      profile.close();
   }    
  
  // rewrite with new

  *arena = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, DATA->neta);
  *Lneighbor = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, 2);
  *Rneighbor = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, 2);
  
  cout << "Grid allocated." << endl;
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
	//       (*arena)[ix][iy][ieta].nbr_p_1 =  (Grid **) malloc (sizeof(Grid *)*4); 
	//       (*arena)[ix][iy][ieta].nbr_m_1 = (Grid **) malloc (sizeof(Grid *)*4); 
	(*arena)[ix][iy][ieta].nbr_p_1 = new Grid *[4];
	(*arena)[ix][iy][ieta].nbr_m_1 = new Grid *[4];
	
	for(int i=1; i<4; i++)
	  {
	    (*arena)[ix][iy][ieta].nbr_p_1[i] = new Grid;
	    (*arena)[ix][iy][ieta].nbr_m_1[i] = new Grid;
	  }

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
          { (*arena)[ix][iy][ieta].nbr_p_1[1] = NULL; 
	  }
         
	 if(ix != 0)
          { (*arena)[ix][iy][ieta].nbr_m_1[1] = &(*arena)[ix-1][iy][ieta]; 
	  }
         else
          { (*arena)[ix][iy][ieta].nbr_m_1[1] = NULL; 
	  }
         
	 if(iy != ny)
          { (*arena)[ix][iy][ieta].nbr_p_1[2] = &(*arena)[ix][iy+1][ieta]; 
	  }
         else
          { (*arena)[ix][iy][ieta].nbr_p_1[2] = NULL; 
	  }
         
	 if(iy != 0)
          { (*arena)[ix][iy][ieta].nbr_m_1[2] = &(*arena)[ix][iy-1][ieta]; 
          }
         else
          { (*arena)[ix][iy][ieta].nbr_m_1[2] = NULL; 
          }

	 // do not care which rank it is - that is dealt with in evolve.cpp
	 if(ieta != neta-1)
          { (*arena)[ix][iy][ieta].nbr_p_1[3] = &(*arena)[ix][iy][ieta+1]; 
	  }
         else 
	   {
	     (*arena)[ix][iy][ieta].nbr_p_1[3] = NULL; 
	   }
         
	 if(ieta != 0)
	   { (*arena)[ix][iy][ieta].nbr_m_1[3] = &(*arena)[ix][iy][ieta-1]; 
	   }
         else
	   {
	     (*arena)[ix][iy][ieta].nbr_m_1[3] = NULL; 
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
 double epsilon0, testEpsilon, epsilon, R_A, a_A, p, h, u[4], x, y, eta, rho, volume;
 double cosheta, sinheta, rhob, exparg1, exparg;
 double eta0, eta_fall_off, eta_flat, a_short, a_long;
 int ix, iy, ieta, ietar, mu, nu, rk_order;
 int initializeEntropy = DATA->initializeEntropy;


//  FILE* e_file;
//  char* e_name;
//  if (rank ==0) e_name = "e_profile.dat"; 
//  else e_name = "e_profile2.dat";
//  e_file = fopen(e_name, "w");
//  FILE* e2_file;
//  char* e2_name = "e_x_profile.dat"; 
//  e2_file = fopen(e2_name, "w");
//  FILE* e3_file;
//  char* e3_name = "e_y_profile.dat"; 
//  e3_file = fopen(e3_name, "w");
//  FILE* e4_file;
//  char* e4_name = "e_x_y_profile.dat"; 
//  e4_file = fopen(e4_name, "w");
//  FILE* e5_file;
//  char* e5_name = "e_x_y_profile_05.dat"; 
//  e5_file = fopen(e5_name, "w");

 ofstream geometry_file("geometry.dat");
 ofstream geometry_file_rn("geometry_rn.dat");

 ofstream twodplot_file("2d_positions.dat");
 ofstream twodplotA_file("2dA_positions.dat");
 ofstream twodplotB_file("2dB_positions.dat");

 ofstream psi_file("psi.gp");
 
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

 cout << "rk_order=" << rk_order << endl;
 
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
	     (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	  		     
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
		 
		 (*arena)[ix][iy][ieta].epsilon = epsilon;
		 (*arena)[ix][iy][ieta].rhob = rhob;
		 (*arena)[ix][iy][ieta].p = p;
		 (*arena)[ix][iy][ieta].trouble = 0;
		
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

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
		 
		 (*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
		 

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
			 
			 (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
			 
		       }/* nu */
		   }/* mu */
		 
	       }}}/* ix, iy, ieta */
   }
 else if (DATA->Initial_profile==1) //full average initial conditions using Glauber
   {
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
     
// //      //secret program bonus:
// //      //*************************************************************************************************
// //      //compute average b: (kinda rough but good enough)

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

//      for (ib1=490; ib1<491; ib1++) // same value here for both , e.g.     for (ib1=191; ib1<192; ib1++)
//        {
// 	 b1=ib1*deltab;
// 	 for (ib2=0; ib2<5000; ib2++)
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
	     (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
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

		     if (epsilon<0.00000000001)
		       epsilon = 0.00000000001;
		     
		     // distribution of initial baryon density:
		     rhob = DATA->rhoB0/epsilon0*epsilon; /* in fm^-3 */
		  //    if (x==0 && y==0) 
// 		       {
// 			 fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
// 			 fprintf(e_file,"%e %e\n",eta,epsilon*hbarc);
// 		       }
// 		     if (y==0 && eta==0) 
// 		       {
// 			 fprintf(e2_file,"%e %e\n",x,epsilon*hbarc);
// 		       }
// 		     if (x==0 && eta==0) 
// 		       {
// 			 fprintf(e3_file,"%e %e\n",y,epsilon*hbarc);
// 		       }
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
	// 	     if (x==0 && y==0) 
// 		       {
// 			 fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
// 			 fprintf(e_file,"%e %e %e %e\n",eta,s,rhob, epsilon*hbarc);
// 		       }
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
		     (*arena)[ix][iy][ieta].T = eos->interpolate(epsilon, rhob, 0);
		     (*arena)[ix][iy][ieta].mu = eos->interpolate(epsilon, rhob, 1);
		   }
		 else if (DATA->whichEOS==2)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate2(epsilon, rhob, 1);
		     (*arena)[ix][iy][ieta].mu = 0.0;
		   }
		 
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
		//  int taui;
// 		 for (taui=0; taui<3; taui++)
// 		   {
// 		     (*arena)[ix][iy][ieta].tauF[taui] = 0.;
// 		   }
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
		 
		 (*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;

		 //(*arena)[ix][iy][ieta].correction = 1.0;

		 
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
			 
			 (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 
	       }}}/* ix, iy, ieta */
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
	     (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
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
		     (*arena)[ix][iy][ieta].T = eos->interpolate(epsilon, rhob, 0);
		     (*arena)[ix][iy][ieta].mu = eos->interpolate(epsilon, rhob, 1);
		   }
		 else if (DATA->whichEOS==2)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate2(epsilon, rhob, 1);
		     (*arena)[ix][iy][ieta].mu = 0.0;
		   }
	
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
// 		 int taui;
// 		 for (taui=0; taui<3; taui++)
// 		   {
// 		     (*arena)[ix][iy][ieta].tauF[taui] = 0.;
// 		   }
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
		 
		 (*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
      		 
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
			 
			 (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 
	       }}}/* ix, iy, ieta */
   }
 else if (DATA->Initial_profile==3) // sample initial distribution with a Glauber Monte-Carlo (event-by-event initial condition)
   {

     //keep this in mind: from Nara and Hirano: 0904.2966
//      In our MC-Glauber model, we find the defaultWoods-Saxon
//      distribution is reproduced by a larger radius parameter
//      r0 = 6.42 (4.28) fm and a smaller diffuseness parameter
//      d = 0.44 (0.50) fm (i.e., sharper boundary of a nucleus)
//      for a gold (copper) nucleus than the default parameters.
     
     // sample Nu (which gives the number of nucleons at radial position r) for both nuclei
     // for now use the same nucleus for both
     //double rnum = 6;
     double rnum=time(0)+DATA->seed; // initialize random seed using the system time
     cout << "random seed=" << rnum-DATA->seed << " + " << DATA->seed << " = " << rnum << endl;
     random->init_genrand64(rnum); 
     double width = DATA->sigma0; // width of Gaussian in fm around the center of the nucleon-nucleon collision;
     double norm = 0.2357518037*(0.4/width)*(0.4/width); // norm to get the right average energy density in the center. can be determined bu using option Initial_Profile = 4
     double A;
     double avxBinary, avyBinary;
     double Z;
     double bArray[1];
     double AxArray[300], AyArray[300]; // arrays for sending when using MPI
     double BxArray[300], ByArray[300];
     ReturnValue rv, rv2; //auxiliary nucleons
     double eccentricity2, Psi2, Psi2a, Psi2b;
     double eccentricity3, Psi3;
     double eccentricity4, Psi4;
     double eccentricity5, Psi5;
     double eccentricity6, Psi6;
     double eccentricity3_rn, Psi3_rn;
     double eccentricity4_rn, Psi4_rn;
     double eccentricity5_rn, Psi5_rn;
     double eccentricity6_rn, Psi6_rn;
     double x, y, xm, ym, dx, dy, dij;
     double xA, yA, rA;
     double phiA;
     double avrSq, avr3, avr4, avr5, avr6;
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
     double xbinColl[5000];  // x coordinate of binary collisions
     double ybinColl[5000];  // y coordinate
     double WbinColl; // W binary collision scaled
     
     for (int i = 0; i<A; i++) 
       {
	 for (int j = 0 ; j<A ;j++) 
	   {
	     dx = nucleusB.at(j).x-nucleusA.at(i).x;
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

	 avxBinary=avx;
	 avyBinary=avy;
	 
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
	 // I use both Alver and Roland's definition of ecc_3 and Psi_3, and that by Qin et al who want r^2 replaced by r^3 in there
	 // -----------------------------------------------------------------
	 avrSq = 0.;
	 double avcos = 0.;
	 double avcosa = 0.;
	 double avcosb = 0.;
	 double avsin = 0.;
	 double avcos3 = 0.;
	 double avsin3 = 0.;
	 double avcos4 = 0.;
	 double avsin4 = 0.;
	 double avcos5 = 0.;
	 double avsin5 = 0.;
	 double avcos6 = 0.;
	 double avsin6 = 0.;
	
	 avr3 = 0.;
	 avr4 = 0.;
	 avr5 = 0.;
	 avr6 = 0.;
	 double avcos3_rn = 0.;
	 double avsin3_rn = 0.;
	 double avcos4_rn = 0.;
	 double avsin4_rn = 0.;
	 double avcos5_rn = 0.;
	 double avsin5_rn = 0.;
	 double avcos6_rn = 0.;
	 double avsin6_rn = 0.;
	 double newx, newxn3, newxn3rn, newxn4, newxn4rn, newxn5, newxn5rn, newxn6, newxn6rn;
	 double newy, newyn3, newyn3rn, newyn4, newyn4rn, newyn5, newyn5rn, newyn6, newyn6rn;
	 double avgx=0;
	 double avgy=0;
	 double avgx2=0; // average x squared
	 double avgy2=0;
	 double avgxAbs=0; 
	 double avgyAbs=0;
	 double avgxn3=0;
	 double avgyn3=0;
	 double avgx2n3=0; // average x squared
	 double avgy2n3=0;
	 double avgxAbsn3=0; 
	 double avgyAbsn3=0;
	 double avgxn3rn=0;
	 double avgyn3rn=0;
	 double avgx2n3rn=0; // average x squared
	 double avgy2n3rn=0;
	 double avgxAbsn3rn=0; 
	 double avgyAbsn3rn=0;
	 double avgxn4=0;
	 double avgyn4=0;
	 double avgx2n4=0; // average x squared
	 double avgy2n4=0;
	 double avgxAbsn4=0; 
	 double avgyAbsn4=0;
	 double avgxn4rn=0;
	 double avgyn4rn=0;
	 double avgx2n4rn=0; // average x squared
	 double avgy2n4rn=0;
	 double avgxAbsn4rn=0; 
	 double avgyAbsn4rn=0;
	 double avgxn5=0;
	 double avgyn5=0;
	 double avgx2n5=0; // average x squared
	 double avgy2n5=0;
	 double avgxAbsn5=0; 
	 double avgyAbsn5=0;
	 double avgxn5rn=0;
	 double avgyn5rn=0;
	 double avgx2n5rn=0; // average x squared
	 double avgy2n5rn=0;
	 double avgxAbsn5rn=0; 
	 double avgyAbsn5rn=0;
	 double avgxn6=0;
	 double avgyn6=0;
	 double avgx2n6=0; // average x squared
	 double avgy2n6=0;
	 double avgxAbsn6=0; 
	 double avgyAbsn6=0;
	 double avgxn6rn=0;
	 double avgyn6rn=0;
	 double avgx2n6rn=0; // average x squared
	 double avgy2n6rn=0;
	 double avgxAbsn6rn=0; 
	 double avgyAbsn6rn=0;

	 double sigmax2;
	 double sigmay2;
	 double sigmax2n3;
	 double sigmay2n3;
	 double sigmax2n3rn;
	 double sigmay2n3rn;
	 double sigmax2n4;
	 double sigmay2n4;
	 double sigmax2n4rn;
	 double sigmay2n4rn;
	 double sigmax2n5;
	 double sigmay2n5;
	 double sigmax2n5rn;
	 double sigmay2n5rn;
	 double sigmax2n6;
	 double sigmay2n6;
	 double sigmax2n6rn;
	 double sigmay2n6rn;

	 double Rbar2;
	 double Rbar3;
	 double Rbar3rn;
	 double Rbar4;
	 double Rbar4rn;
	 double Rbar5;
	 double Rbar5rn;
	 double Rbar6;
	 double Rbar6rn;

	 double DeltaL2;
	 double DeltaL3;
	 double DeltaL3rn;
	 double DeltaL4;
	 double DeltaL4rn;
	 double DeltaL5;
	 double DeltaL5rn;
	 double DeltaL6;
	 double DeltaL6rn;

	 double randalpha;
	 int counter=0;
	 double alpha =0.14;
	 
       	 for (int i = 0; i<A; i++) 
	   {
	     if ( nucleusA.at(i).collided == 1)
	       {
		 // randalpha=random->genrand64_real1();	
		 
		 //if(randalpha>alpha)
		 //  {
		 //xA = xbinColl[i];
		 xA=nucleusA.at(i).x;
		 //yA = ybinColl[i] ;
		 yA=nucleusA.at(i).y;
		 rA = sqrt( xA*xA + yA*yA );
		 if (xA>=0)
		   {
		     phiA = atan(yA/xA);
		   }
		 else
		   {
		     phiA = atan(yA/xA)+PI;
		   }
		 
		 avrSq += rA*rA; // compute average r^2
		 avr3  += rA*rA*rA;
		 avr4  += rA*rA*rA*rA;
		 avr5  += rA*rA*rA*rA*rA;
		 avr6  += rA*rA*rA*rA*rA*rA;
		 
		 avcos  += rA*rA*cos(2.*phiA);
		 avsin  += rA*rA*sin(2.*phiA);
		 avcos3 += rA*rA*cos(3.*phiA);
		 avsin3 += rA*rA*sin(3.*phiA);
		 avcos4 += rA*rA*cos(4.*phiA);
		 avsin4 += rA*rA*sin(4.*phiA);
		 avcos5 += rA*rA*cos(5.*phiA);
		 avsin5 += rA*rA*sin(5.*phiA);
		 avcos6 += rA*rA*cos(6.*phiA);
		 avsin6 += rA*rA*sin(6.*phiA);
		 
		 avcos3_rn += rA*rA*rA*cos(3.*phiA);
		 avsin3_rn += rA*rA*rA*sin(3.*phiA);
		 avcos4_rn += rA*rA*rA*rA*cos(4.*phiA);
		 avsin4_rn += rA*rA*rA*rA*sin(4.*phiA);
		 avcos5_rn += rA*rA*rA*rA*rA*cos(5.*phiA);
		 avsin5_rn += rA*rA*rA*rA*rA*sin(5.*phiA);
		 avcos6_rn += rA*rA*rA*rA*rA*rA*cos(6.*phiA);
		 avsin6_rn += rA*rA*rA*rA*rA*rA*sin(6.*phiA);
		 counter++;
	       }
	     // 		 else
	     // 		   {
	     // 		     xA = (nucleusA.at(i).x+nucleusB.at(i).x)/2.;
	     // 		     yA = (nucleusA.at(i).y+nucleusA.at(i).y)/2.;
	     // 		     rA = sqrt( xA*xA + yA*yA );
	     // 		     phiA = atan(yA/xA);
	     
	     // 		     avrSq += rA*rA; // compute average r^2
	     // 		     avr3 += rA*rA*rA;
	     // 		     avr4 += rA*rA*rA*rA;
	     // 		     avr5 += rA*rA*rA*rA*rA;
	     // 		     avr6 += rA*rA*rA*rA*rA*rA;
	     
	     // 		     avcos += rA*rA*cos(2.*phiA);
	     // 		     avsin += rA*rA*sin(2.*phiA);
	     // 		     avcos3 += rA*rA*cos(3.*phiA);
	     // 		     avsin3 += rA*rA*sin(3.*phiA);
	     // 		     avcos4 += rA*rA*cos(4.*phiA);
	     // 		     avsin4 += rA*rA*sin(4.*phiA);
	     // 		     avcos5 += rA*rA*cos(5.*phiA);
	     // 		     avsin5 += rA*rA*sin(5.*phiA);
	     // 		     avcos6 += rA*rA*cos(6.*phiA);
	     // 		     avsin6 += rA*rA*sin(6.*phiA);
	     
	     // 		     avcos3_rn += rA*rA*rA*cos(3.*phiA);
	     // 		     avsin3_rn += rA*rA*rA*sin(3.*phiA);
	     // 		     avcos4_rn += rA*rA*rA*rA*cos(4.*phiA);
	     // 		     avsin4_rn += rA*rA*rA*rA*sin(4.*phiA);
	     // 		     avcos5_rn += rA*rA*rA*rA*rA*cos(5.*phiA);
	     // 		     avsin5_rn += rA*rA*rA*rA*rA*sin(5.*phiA);
	     // 		     avcos6_rn += rA*rA*rA*rA*rA*rA*cos(6.*phiA);
	     // 		     avsin6_rn += rA*rA*rA*rA*rA*rA*sin(6.*phiA);
	     // 		     counter++;
	     // 		   }
	     
	     if ( nucleusB.at(i).collided == 1)
	       {
// 		 //		 randalpha=random->genrand64_real1();	
		 
// 		 //if(randalpha>alpha)
		 //{
		 xA = nucleusB.at(i).x;
		 yA = nucleusB.at(i).y;
		 rA = sqrt( xA*xA + yA*yA );
		 if (xA>=0)
		   {
		     phiA = atan(yA/xA);
		   }
		 else
		   {
		     phiA = atan(yA/xA)+PI;
		   }
		 
		 avrSq += rA*rA; // compute average r^2
		 avr3 += rA*rA*rA;
		 avr4 += rA*rA*rA*rA;
		 avr5 += rA*rA*rA*rA*rA;
		 avr6 += rA*rA*rA*rA*rA*rA;
		 
		 avcos += rA*rA*cos(2.*phiA);
		 avsin += rA*rA*sin(2.*phiA);
		 avcos3 += rA*rA*cos(3.*phiA);
		 avsin3 += rA*rA*sin(3.*phiA);
		 avcos4 += rA*rA*cos(4.*phiA);
		 avsin4 += rA*rA*sin(4.*phiA);
		 avcos5 += rA*rA*cos(5.*phiA);
		 avsin5 += rA*rA*sin(5.*phiA);
		 avcos6 += rA*rA*cos(6.*phiA);
		 avsin6 += rA*rA*sin(6.*phiA);
		 
		 avcos3_rn += rA*rA*rA*cos(3.*phiA);
		 avsin3_rn += rA*rA*rA*sin(3.*phiA);
		 avcos4_rn += rA*rA*rA*rA*cos(4.*phiA);
		 avsin4_rn += rA*rA*rA*rA*sin(4.*phiA);
		 avcos5_rn += rA*rA*rA*rA*rA*cos(5.*phiA);
		 avsin5_rn += rA*rA*rA*rA*rA*sin(5.*phiA);
		 avcos6_rn += rA*rA*rA*rA*rA*rA*cos(6.*phiA);
		 avsin6_rn += rA*rA*rA*rA*rA*rA*sin(6.*phiA);
		 counter++;
		 //		   }
	       }
	     
	   }	     
       
	 cout << " N-part = " << participants << endl;
	 int Npart = participants;
	 
	 participants=counter;
	 
	 // compute average r^2 and average terms in the numerator:
	 avrSq/=participants;
	 avr3/=participants;
	 avr4/=participants;
	 avr5/=participants;
	 avr6/=participants;
	 
	 avcos/=participants;
	 avsin/=participants;
	 avcos3/=participants;
	 avsin3/=participants;
	 avcos4/=participants;
	 avsin4/=participants;
	 avcos5/=participants;
	 avsin5/=participants;
	 avcos6/=participants;
	 avsin6/=participants;
	 
	 avcos3_rn/=participants;
	 avsin3_rn/=participants;
	 avcos4_rn/=participants;
	 avsin4_rn/=participants;
	 avcos5_rn/=participants;
	 avsin5_rn/=participants;
	 avcos6_rn/=participants;
	 avsin6_rn/=participants;
	 
	 //Psi2 = (atan2(avsin, avcos))/2.;
	 //Psi2 = (atan2(avsin, avcos)+PI)/2.;
	 Psi2 = (atan(avsin/avcos))/2.;
	 eccentricity2 = sqrt(avcos*avcos+avsin*avsin)/avrSq;
	 //Psi3 = (atan2(avsin3, avcos3))/3.;
	 Psi3 = (atan2(avsin3, avcos3))/3.;
	 // Psi3 = (atan(avsin3/avcos3))/3.;
	 eccentricity3 = sqrt(avcos3*avcos3+avsin3*avsin3)/avrSq;
	 Psi4 = (atan2(avsin4, avcos4))/4.;
	 eccentricity4 = sqrt(avcos4*avcos4+avsin4*avsin4)/avrSq;
	 Psi5 = (atan2(avsin5, avcos5))/5.;
	 eccentricity5 = sqrt(avcos5*avcos5+avsin5*avsin5)/avrSq;
	 Psi6 = (atan2(avsin6, avcos6))/6.;
	 eccentricity6 = sqrt(avcos6*avcos6+avsin6*avsin6)/avrSq;
	
	 //Psi3_rn = (atan2(avsin3_rn, avcos3_rn))/3.;
	 //Psi3_rn = (atan2(avsin3_rn, avcos3_rn)+PI)/3.;
	 Psi3_rn = (atan(avsin3_rn/avcos3_rn))/3.;
	 eccentricity3_rn = sqrt(avcos3_rn*avcos3_rn+avsin3_rn*avsin3_rn)/avr3;
 	 Psi4_rn = (atan2(avsin4_rn, avcos4_rn))/4.;
	 eccentricity4_rn = sqrt(avcos4_rn*avcos4_rn+avsin4_rn*avsin4_rn)/avr4;
	 Psi5_rn = (atan2(avsin5_rn, avcos5_rn))/5.;
	 eccentricity5_rn = sqrt(avcos5_rn*avcos5_rn+avsin5_rn*avsin5_rn)/avr5;
	 Psi6_rn = (atan2(avsin6_rn, avcos6_rn))/6.;
	 eccentricity6_rn = sqrt(avcos6_rn*avcos6_rn+avsin6_rn*avsin6_rn)/avr6;
	 
	 cout << " the participant eccentricity is eps_2=" << eccentricity2 << endl;
	 cout << " the minor axis of the ellipse defined by the participants is Psi_2=" << Psi2 << endl;
	 cout << " the participant triangularity is eps_3=" << eccentricity3_rn << endl;
	 cout << " the minor axis of the participant triangularity is Psi_3=" << Psi3_rn << endl;


	 // compute Roy's quantities
	 for (int i = 0; i<A; i++) 
	   {
	     if ( nucleusA.at(i).collided == 1)
	       {
		 xA = nucleusA.at(i).x;
		 yA = nucleusA.at(i).y;
		 
		 newx = xA*cos(Psi2)+yA*sin(Psi2);
		 newy = -xA*sin(Psi2)+yA*cos(Psi2);

		 avgx += newx;
		 avgy += newy;
		 avgxAbs += abs(newx);
		 avgyAbs += abs(newy);
		 avgx2 += newx*newx;
		 avgy2 += newy*newy;

		 newxn3 = xA*cos(Psi3)+yA*sin(Psi3);
		 //what is called y is really x shifted by another 60 degrees (for n=3)
		 newyn3 = xA*cos(Psi3+1.047197551)+yA*sin(Psi3+1.047197551);

		 avgxn3 += newxn3;
		 avgyn3 += newyn3;
		 avgxAbsn3 += abs(newxn3);
		 avgyAbsn3 += abs(newyn3);
		 avgx2n3 += newxn3*newxn3;
		 avgy2n3 += newyn3*newyn3;

		 newxn3rn = xA*cos(Psi3_rn)+yA*sin(Psi3_rn);
		 newyn3rn = xA*cos(Psi3_rn+1.047197551)+yA*sin(Psi3_rn+1.047197551);

		 avgxn3rn += newxn3rn;
		 avgyn3rn += newyn3rn;
		 avgxAbsn3rn += abs(newxn3rn);
		 avgyAbsn3rn += abs(newyn3rn);
		 avgx2n3rn += newxn3rn*newxn3rn;
		 avgy2n3rn += newyn3rn*newyn3rn;

		 newxn4 = xA*cos(Psi4)+yA*sin(Psi4);
		 newyn4 = xA*cos(Psi4+0.785398163)+yA*sin(Psi4+0.785398163);

		 avgxn4 += newxn4;
		 avgyn4 += newyn4;
		 avgxAbsn4 += abs(newxn4);
		 avgyAbsn4 += abs(newyn4);
		 avgx2n4 += newxn4*newxn4;
		 avgy2n4 += newyn4*newyn4;

		 newxn4rn = xA*cos(Psi4_rn)+yA*sin(Psi4_rn);
		 newyn4rn = xA*cos(Psi4_rn+0.785398163)+yA*sin(Psi4_rn+0.785398163);

		 avgxn4rn += newxn4rn;
		 avgyn4rn += newyn4rn;
		 avgxAbsn4rn += abs(newxn4rn);
		 avgyAbsn4rn += abs(newyn4rn);
		 avgx2n4rn += newxn4rn*newxn4rn;
		 avgy2n4rn += newyn4rn*newyn4rn;

		 newxn5 = xA*cos(Psi5)+yA*sin(Psi5);
		 newyn5 = xA*cos(Psi5+0.62831853)+yA*sin(Psi5+0.62831853);

		 avgxn5 += newxn5;
		 avgyn5 += newyn5;
		 avgxAbsn5 += abs(newxn5);
		 avgyAbsn5 += abs(newyn5);
		 avgx2n5 += newxn5*newxn5;
		 avgy2n5 += newyn5*newyn5;

		 newxn5rn = xA*cos(Psi5_rn)+yA*sin(Psi5_rn);
		 newyn5rn = xA*cos(Psi5_rn+0.62831853)+yA*sin(Psi5_rn+0.62831853);

		 avgxn5rn += newxn5rn;
		 avgyn5rn += newyn5rn;
		 avgxAbsn5rn += abs(newxn5rn);
		 avgyAbsn5rn += abs(newyn5rn);
		 avgx2n5rn += newxn5rn*newxn5rn;
		 avgy2n5rn += newyn5rn*newyn5rn;

		 newxn6 = xA*cos(Psi6)+yA*sin(Psi6);
		 newyn6 = xA*cos(Psi6+0.523598775)+yA*sin(Psi6+0.523598775);

		 avgxn6 += newxn6;
		 avgyn6 += newyn6;
		 avgxAbsn6 += abs(newxn6);
		 avgyAbsn6 += abs(newyn6);
		 avgx2n6 += newxn6*newxn6;
		 avgy2n6 += newyn6*newyn6;

		 newxn6rn = xA*cos(Psi6_rn)+yA*sin(Psi6_rn);
		 newyn6rn = xA*cos(Psi6_rn+0.523598775)+yA*sin(Psi6_rn+0.523598775);

		 avgxn6rn += newxn6rn;
		 avgyn6rn += newyn6rn;
		 avgxAbsn6rn += abs(newxn6rn);
		 avgyAbsn6rn += abs(newyn6rn);
		 avgx2n6rn += newxn6rn*newxn6rn;
		 avgy2n6rn += newyn6rn*newyn6rn;

	       }

	     if ( nucleusB.at(i).collided == 1)
	       {
		 xA = nucleusB.at(i).x;
		 yA = nucleusB.at(i).y;

		 newx = xA*cos(Psi2)+yA*sin(Psi2);
		 newy = -xA*sin(Psi2)+yA*cos(Psi2);

		 //	 cout << "psi_2=" << Psi2 << endl;
		 // cout << "x=" << xA << ", y=" << yA << ", new x=" << newx << ", new y=" << newy << endl;

		 avgx += newx;
		 avgy += newy;
		 avgxAbs += abs(newx);
		 avgyAbs += abs(newy);
		 avgx2 += newx*newx;
		 avgy2 += newy*newy;

		 newxn3 = xA*cos(Psi3)+yA*sin(Psi3);
		 //what is called y is really x shifted by another 60 degrees (for n=3)
		 newyn3 = xA*cos(Psi3+1.047197551)+yA*sin(Psi3+1.047197551);

		 avgxn3 += newxn3;
		 avgyn3 += newyn3;
		 avgxAbsn3 += abs(newxn3);
		 avgyAbsn3 += abs(newyn3);
		 avgx2n3 += newxn3*newxn3;
		 avgy2n3 += newyn3*newyn3;

		 newxn3rn = xA*cos(Psi3_rn)+yA*sin(Psi3_rn);
		 newyn3rn = xA*cos(Psi3_rn+1.047197551)+yA*sin(Psi3_rn+1.047197551);

		 avgxn3rn += newxn3rn;
		 avgyn3rn += newyn3rn;
		 avgxAbsn3rn += abs(newxn3rn);
		 avgyAbsn3rn += abs(newyn3rn);
		 avgx2n3rn += newxn3rn*newxn3rn;
		 avgy2n3rn += newyn3rn*newyn3rn;

		 newxn4 = xA*cos(Psi4)+yA*sin(Psi4);
		 newyn4 = xA*cos(Psi4+0.785398163)+yA*sin(Psi4+0.785398163);

		 avgxn4 += newxn4;
		 avgyn4 += newyn4;
		 avgxAbsn4 += abs(newxn4);
		 avgyAbsn4 += abs(newyn4);
		 avgx2n4 += newxn4*newxn4;
		 avgy2n4 += newyn4*newyn4;

		 newxn4rn = xA*cos(Psi4_rn)+yA*sin(Psi4_rn);
		 newyn4rn = xA*cos(Psi4_rn+0.785398163)+yA*sin(Psi4_rn+0.785398163);

		 avgxn4rn += newxn4rn;
		 avgyn4rn += newyn4rn;
		 avgxAbsn4rn += abs(newxn4rn);
		 avgyAbsn4rn += abs(newyn4rn);
		 avgx2n4rn += newxn4rn*newxn4rn;
		 avgy2n4rn += newyn4rn*newyn4rn;

		 newxn5 = xA*cos(Psi5)+yA*sin(Psi5);
		 newyn5 = xA*cos(Psi5+0.62831853)+yA*sin(Psi5+0.62831853);

		 avgxn5 += newxn5;
		 avgyn5 += newyn5;
		 avgxAbsn5 += abs(newxn5);
		 avgyAbsn5 += abs(newyn5);
		 avgx2n5 += newxn5*newxn5;
		 avgy2n5 += newyn5*newyn5;

		 newxn5rn = xA*cos(Psi5_rn)+yA*sin(Psi5_rn);
		 newyn5rn = xA*cos(Psi5_rn+0.62831853)+yA*sin(Psi5_rn+0.62831853);

		 avgxn5rn += newxn5rn;
		 avgyn5rn += newyn5rn;
		 avgxAbsn5rn += abs(newxn5rn);
		 avgyAbsn5rn += abs(newyn5rn);
		 avgx2n5rn += newxn5rn*newxn5rn;
		 avgy2n5rn += newyn5rn*newyn5rn;

		 newxn6 = xA*cos(Psi6)+yA*sin(Psi6);
		 newyn6 = xA*cos(Psi6+0.523598775)+yA*sin(Psi6+0.523598775);

		 avgxn6 += newxn6;
		 avgyn6 += newyn6;
		 avgxAbsn6 += abs(newxn6);
		 avgyAbsn6 += abs(newyn6);
		 avgx2n6 += newxn6*newxn6;
		 avgy2n6 += newyn6*newyn6;

		 newxn6rn = xA*cos(Psi6_rn)+yA*sin(Psi6_rn);
		 newyn6rn = xA*cos(Psi6_rn+0.523598775)+yA*sin(Psi6_rn+0.523598775);

		 avgxn6rn += newxn6rn;
		 avgyn6rn += newyn6rn;
		 avgxAbsn6rn += abs(newxn6rn);
		 avgyAbsn6rn += abs(newyn6rn);
		 avgx2n6rn += newxn6rn*newxn6rn;
		 avgy2n6rn += newyn6rn*newyn6rn;

	       }
	   }
	 
	 avgx/=Npart;
	 avgy/=Npart;
	 avgx2/=Npart;
	 avgy2/=Npart;
	 avgxn3/=Npart;
	 avgyn3/=Npart;
	 avgx2n3/=Npart;
	 avgy2n3/=Npart;
	 avgxn3rn/=Npart;
	 avgyn3rn/=Npart;
	 avgx2n3rn/=Npart;
	 avgy2n3rn/=Npart;
	 avgxn4/=Npart;
	 avgyn4/=Npart;
	 avgx2n4/=Npart;
	 avgy2n4/=Npart;
	 avgxn4rn/=Npart;
	 avgyn4rn/=Npart;
	 avgx2n4rn/=Npart;
	 avgy2n4rn/=Npart;
	 avgxn5/=Npart;
	 avgyn5/=Npart;
	 avgx2n5/=Npart;
	 avgy2n5/=Npart;
	 avgxn5rn/=Npart;
	 avgyn5rn/=Npart;
	 avgx2n5rn/=Npart;
	 avgy2n5rn/=Npart;
	 avgxn6/=Npart;
	 avgyn6/=Npart;
	 avgx2n6/=Npart;
	 avgy2n6/=Npart;
	 avgxn6rn/=Npart;
	 avgyn6rn/=Npart;
	 avgx2n6rn/=Npart;
	 avgy2n6rn/=Npart;
	 avgxAbs/=Npart;
	 avgyAbs/=Npart;
	 avgxAbsn3/=Npart;
	 avgyAbsn3/=Npart;
	 avgxAbsn3rn/=Npart;
	 avgyAbsn3rn/=Npart;
	 avgxAbsn4/=Npart;
	 avgyAbsn4/=Npart;
	 avgxAbsn4rn/=Npart;
	 avgyAbsn4rn/=Npart;
	 avgxAbsn5/=Npart;
	 avgyAbsn5/=Npart;
	 avgxAbsn5rn/=Npart;
	 avgyAbsn5rn/=Npart;
	 avgxAbsn6/=Npart;
	 avgyAbsn6/=Npart;
	 avgxAbsn6rn/=Npart;
	 avgyAbsn6rn/=Npart;



	 sigmax2 = avgx2-avgx*avgx;
	 sigmay2 = avgy2-avgy*avgy;
	 sigmax2n3 = avgx2n3-avgxn3*avgxn3;
	 sigmay2n3 = avgy2n3-avgyn3*avgyn3;
	 sigmax2n3rn = avgx2n3rn-avgxn3rn*avgxn3rn;
	 sigmay2n3rn = avgy2n3rn-avgyn3rn*avgyn3rn;
	 sigmax2n4 = avgx2n4-avgxn4*avgxn4;
	 sigmay2n4 = avgy2n4-avgyn4*avgyn4;
	 sigmax2n4rn = avgx2n4rn-avgxn4rn*avgxn4rn;
	 sigmay2n4rn = avgy2n4rn-avgyn4rn*avgyn4rn;
	 sigmax2n5 = avgx2n5-avgxn5*avgxn5;
	 sigmay2n5 = avgy2n5-avgyn5*avgyn5;
	 sigmax2n5rn = avgx2n5rn-avgxn5rn*avgxn5rn;
	 sigmay2n5rn = avgy2n5rn-avgyn5rn*avgyn5rn;
	 sigmax2n6 = avgx2n6-avgxn6*avgxn6;
	 sigmay2n6 = avgy2n6-avgyn6*avgyn6;
	 sigmax2n6rn = avgx2n6rn-avgxn6rn*avgxn6rn;
	 sigmay2n6rn = avgy2n6rn-avgyn6rn*avgyn6rn;

	 Rbar2 = 1./sqrt(1./sigmax2+1./sigmay2);
	 Rbar3 = 1./sqrt(1./sigmax2n3+1./sigmay2n3);
	 Rbar3rn = 1./sqrt(1./sigmax2n3rn+1./sigmay2n3rn);
	 Rbar4 = 1./sqrt(1./sigmax2n4+1./sigmay2n4);
	 Rbar4rn = 1./sqrt(1./sigmax2n4rn+1./sigmay2n4rn);
	 Rbar5 = 1./sqrt(1./sigmax2n5+1./sigmay2n5);
	 Rbar5rn = 1./sqrt(1./sigmax2n5rn+1./sigmay2n5rn);
	 Rbar6 = 1./sqrt(1./sigmax2n6+1./sigmay2n6);
	 Rbar6rn = 1./sqrt(1./sigmax2n6rn+1./sigmay2n6rn);
	 
	 DeltaL2 = avgxAbs-avgyAbs;
	 DeltaL3 = avgxAbsn3-avgyAbsn3;
	 DeltaL3rn = avgxAbsn3rn-avgyAbsn3rn;
	 DeltaL4 = avgxAbsn4-avgyAbsn4;
	 DeltaL4rn = avgxAbsn4rn-avgyAbsn4rn;
	 DeltaL5 = avgxAbsn5-avgyAbsn5;
	 DeltaL5rn = avgxAbsn5rn-avgyAbsn5rn;
	 DeltaL6 = avgxAbsn6-avgyAbsn6;
	 DeltaL6rn = avgxAbsn6rn-avgyAbsn6rn;
 


	 cout << "sigma_x^2=" << sigmax2 << ", sigma_y^2=" << sigmay2 << endl;
	 cout << "sigma_x^2(3)=" << sigmax2n3 << ", sigma_y^2(3)=" << sigmay2n3 << endl;



// 	 // alternative eps
// 	 counter = 0;
// 	 avcos=0.;
// 	 avcos3=0.;
// 	 avcos4=0.;
// 	 avcos5=0.;
// 	 avcos6=0.;

// 	 for (int i = 0; i<A; i++) 
// 	   {
// 	     if ( nucleusA.at(i).collided == 1)
// 	       {
// 		 //		 randalpha=random->genrand64_real1();	
// 		 //if(randalpha>alpha)
// 		 //{
// 		     xA = nucleusA.at(i).x;
// 		     yA = nucleusA.at(i).y;
// 		     rA = sqrt( xA*xA + yA*yA );
// 		     if (xA>=0)
// 		       {
// 			 phiA = atan(yA/xA);
// 		       }
// 		     else
// 		       {
// 			 phiA = atan(yA/xA)+PI;
// 		       }
		       
// 		     avcos += rA*rA*cos(2.*(phiA-Psi2));
// 		     avcos3 += rA*rA*rA*cos(3.*(phiA-Psi3_rn));
// 		     avcos4 += rA*rA*rA*rA*cos(4.*(phiA-Psi4_rn));
// 		     avcos5 += rA*rA*rA*rA*rA*cos(5.*(phiA-Psi5_rn));
// 		     avcos6 += rA*rA*rA*rA*rA*rA*cos(6.*(phiA-Psi6_rn));
		     
// 		     counter++;
// 		   }
// // 		 else
// // 		   {
// // 		     xA = (nucleusA.at(i).x+nucleusB.at(i).x)/2.;
// // 		     yA = (nucleusA.at(i).y+nucleusA.at(i).y)/2.;
// // 		     rA = sqrt( xA*xA + yA*yA );
// // 		     phiA = atan(yA/xA);
		     
// // 		     avcos += cos(2.*(phiA-Psi2));
// // 		     avcos3 += cos(3.*(phiA-Psi3_rn));
// // 		     avcos4 += cos(4.*(phiA-Psi4_rn));
// // 		     avcos5 += cos(5.*(phiA-Psi5_rn));
// // 		     avcos6 += cos(6.*(phiA-Psi6_rn));
		     
// // 		     counter++;
// // 		   }
// //	       }
// 	 if ( nucleusB.at(i).collided == 1)
// 	   {
// 	     //	     randalpha=random->genrand64_real1();	
	     
// 	     // if(randalpha>alpha)
// 	     //   {
// 	     xA = nucleusB.at(i).x;
// 	     yA = nucleusB.at(i).y;
// 	     rA = sqrt( xA*xA + yA*yA );
// 	     if (xA>=0)
// 	       {
// 		 phiA = atan(yA/xA);
// 	       }
// 	     else
// 	       {
// 		 phiA = atan(yA/xA)+PI;
// 	       }
	     
// 	     avcos += rA*rA*cos(2.*(phiA-Psi2));
// 	     avcos3 += rA*rA*rA*cos(3.*(phiA-Psi3_rn));
// 	     avcos4 += rA*rA*rA*rA*cos(4.*(phiA-Psi4_rn));
// 	     avcos5 += rA*rA*rA*rA*rA*cos(5.*(phiA-Psi5_rn));
// 	     avcos6 += rA*rA*rA*rA*rA*rA*cos(6.*(phiA-Psi6_rn));
	     
// 	     counter++;
// 	   }
//        }
 
// 	 double eccentricity2a,eccentricity2b;
	 
// 	 eccentricity2 = avcos/counter/avrSq;
// 	 eccentricity3_rn = avcos3/counter/avr3;
// 	 eccentricity4_rn = avcos4/counter/avr4;

// 	 cout << " the ALTERNATIVE participant eccentricity is eps_2=" << eccentricity2 << endl;
// 	 cout << " the ALTERNATIVE participant triangularity is eps_3=" << eccentricity3_rn << endl;
// 	 cout << " the ALTERNATIVE participant eps_4=" << eccentricity4_rn << endl;

 
// 	 // now I need to save Psi_2 and Psi_3 (and also the eccentricities) to compute v_2 and v_3 later
// 	 // v_2 = <cos(2(phi-Psi_2))>
// 	 // v_3 = <cos(3(phi-Psi_3))>
	 
	 
	 geometry_file << "2 " << eccentricity2 << " " << Psi2 << " " << Rbar2 << " " << -DeltaL2 << " " << sigmax2 << " " << sigmay2 << endl;
	 geometry_file << "3 " << eccentricity3 << " " << Psi3 << " " << Rbar3 << " " << -DeltaL3 << " " << sigmax2n3 << " " << sigmay2n3 << endl;
	 geometry_file << "4 " << eccentricity4 << " " << Psi4 << " " << Rbar4 << " " << -DeltaL4 << " " << sigmax2n4 << " " << sigmay2n4 << endl;
	 geometry_file << "5 " << eccentricity5 << " " << Psi5 << " " << Rbar5 << " " << -DeltaL5 << " " << sigmax2n5 << " " << sigmay2n5 << endl;
	 geometry_file << "6 " << eccentricity6 << " " << Psi6 << " " << Rbar6 << " " << -DeltaL6 << " " << sigmax2n6 << " " << sigmay2n6 << endl;

	 geometry_file_rn << "2 " << eccentricity2 << " " << Psi2 << " " << Rbar2 << " " << -DeltaL2 << " " << sigmax2 << " " << sigmay2 << endl;
	 geometry_file_rn << "3 " << eccentricity3_rn << " " << Psi3_rn << " " << Rbar3rn << " " 
			  << -DeltaL3rn << " " << sigmax2n3rn << " " << sigmay2n3rn << endl;
	 geometry_file_rn << "4 " << eccentricity4_rn << " " << Psi4_rn << " " << Rbar4rn << " " 
			  << -DeltaL4rn << " " << sigmax2n4rn << " " << sigmay2n4rn << endl;
	 geometry_file_rn << "5 " << eccentricity5_rn << " " << Psi5_rn << " " << Rbar5rn << " " 
			  << -DeltaL5rn << " " << sigmax2n5rn << " " << sigmay2n5rn << endl;
	 geometry_file_rn << "6 " << eccentricity6_rn << " " << Psi6_rn << " " << Rbar6rn << " " 
			  << -DeltaL6rn << " " << sigmax2n6rn << " " << sigmay2n6rn << endl;


	 psi_file << "a = " << Psi2 << endl;
	 psi_file << "b = " << Psi3 << endl;
	 
	 for(int i=0; i<A; i++)
	   {
	     if ( nucleusA.at(i).collided == 1)
	       twodplot_file << nucleusA.at(i).x << " " << nucleusA.at(i).y << endl; 
	     if ( nucleusB.at(i).collided == 1)
	       twodplot_file << nucleusB.at(i).x << " " << nucleusB.at(i).y << endl; 
	   }

	 for(int i=0; i<A; i++)
	   {
	     twodplotA_file << nucleusA.at(i).x << " " << nucleusA.at(i).y << endl; 
	   }

	 for(int i=0; i<A; i++)
	   {
	     twodplotB_file << nucleusB.at(i).x << " " << nucleusB.at(i).y << endl; 
	   }

       } 

     //     printf("hallo\n");
     

     //exit(1);
     
     //make a distribution function from the initial positions by putting Gaussians on top of each other
     double phipart, rpart; // angular coordinates of the participants
     double va, test;
     double s, r1, r2, W, nBinary, nWounded, W0;
     // impact parameter:
     //normalization of TATarget and TAProjectile is given by Target.A and Projectie.A, respectively
     double normT = glauber->LexusData.Target.A;
     double normP = glauber->LexusData.Projectile.A;
     // the fraction of the hard (binary collisions) contribution
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
		     xm = xbinColl[i]-avxBinary;
		     ym = ybinColl[i]-avyBinary;
		     WbinColl += norm*exp((-(x-xm)*(x-xm)-(y-ym)*(y-ym))/(2.*width*width)); // normalize later to energy density I need
		   }
	       }
	     
	     double Wfull = hard * WbinColl + (1.-hard) * W;
	   
	     W= Wfull;
	     //	     if (hard==0.25) W = Wfull*0.645;

	     (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	 	 
	     epsilon = epsilon0*W;
	     if (epsilon<0.00000000001)
	       epsilon = 0.00000000001;
		     
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
		    
	// 	     if (x==0 && y==0) 
// 		       {
// 			 fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
// 			 fprintf(e_file,"%e %e\n",eta,epsilon*hbarc);
// 		       }
// 		     if (y==0 && eta==0) 
// 		       {
// 			 fprintf(e2_file,"%e %e\n",x,epsilon*hbarc);
// 		       }
// 		     if (x==0 && eta==0) 
// 		       {
// 			 fprintf(e3_file,"%e %e\n",y,epsilon*hbarc);
// 		       }
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
		//      if (x==0 && y==0) 
// 		       {
// 			 fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
// 			 fprintf(e_file,"%e %e %e %e\n",eta,s,rhob, epsilon*hbarc);
// 		       }
		   }
		 // intial pressure distribution:
		 p = eos->p_func(epsilon, rhob);

		 // set all values in the grid element:
		 (*arena)[ix][iy][ieta].epsilon = epsilon;
		 // cout <<  &(*arena)[ix][iy][ieta].epsilon << endl;
		 (*arena)[ix][iy][ieta].rhob = rhob;
		 (*arena)[ix][iy][ieta].p = p;
		 (*arena)[ix][iy][ieta].trouble = 0;

		 if (DATA->whichEOS==1)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate(epsilon, rhob, 0);
		     (*arena)[ix][iy][ieta].mu = eos->interpolate(epsilon, rhob, 1);
		   }
		 else if (DATA->whichEOS==2)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate2(epsilon, rhob, 1);
		     (*arena)[ix][iy][ieta].mu = 0.0;
		   }
		 
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
				
// 		 int taui;
// 		 for (taui=0; taui<3; taui++)
// 		   {
// 		     (*arena)[ix][iy][ieta].tauF[taui] = 0.;
// 		   }
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
		 
		 (*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
		 
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
			 
			 (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
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
	     (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
	     (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
	     (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	     (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	   
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
		 //     if (x==0 && y==0) 
// 		       {
// 			 fprintf(stderr,"e(%f)=%f GeV/fm^3, rhob=%f\n",eta,epsilon*hbarc, rhob);
// 			 fprintf(e_file,"%e %e\n",eta,epsilon*hbarc);
// 		       }
// 		     if (y==0 && eta==0) 
// 		       {
// 			 fprintf(e2_file,"%e %e\n",x,epsilon*hbarc);
// 		       }
// 		     if (x==0 && eta==0) 
// 		       {
// 			 fprintf(e3_file,"%e %e\n",y,epsilon*hbarc);
// 		       }
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
		     (*arena)[ix][iy][ieta].T = eos->interpolate(epsilon, rhob, 0);
		     (*arena)[ix][iy][ieta].mu = eos->interpolate(epsilon, rhob, 1);
		   }
		 else if (DATA->whichEOS==2)
		   {
		     (*arena)[ix][iy][ieta].T = eos->interpolate2(epsilon, rhob, 1);
		     (*arena)[ix][iy][ieta].mu = 0.0;
		   }
		 
		 (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta].prev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta].pprevPimunu = util->cube_malloc(1, 5,4);

		 (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(5,4);
// 		 int taui;
// 		 for (taui=0; taui<3; taui++)
// 		   {
// 		     (*arena)[ix][iy][ieta].tauF[taui] = 0.;
// 		   }
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
		 
		 (*arena)[ix][iy][ieta].pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta].prev_pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta].pprev_pi_b[0] = 0.0;
		 
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
			 
			 (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 
	       }}}/* ix, iy, ieta */
   }   
 else if (DATA->Initial_profile==6) //read in the profile from file - work with Ohio group
   {
    double va, test;
     double s, r1, r2, W, nBinary, nWounded, W0;
     // impact parameter:
     double b=DATA->b;
     int i;
     int bytes_read;
     string dummy;
     int nx, ny, neta;
     int y_size, x_size;
     double dx, dy, deta;
     int ieta2;
     size = DATA->size;
     cout << "size=" << size << endl;
     
     cout << " ----- information on initial distribution -----" << endl;
     cout << "file name used: " << DATA->initName << endl;
  
     ifstream profile(DATA->initName.c_str());

     //     profile = fopen(DATA->initName, "r");
     // read the first line with general info
    
     //          exit(1);

     for(ix=0; ix<=DATA->nx; ix++)
       {
	 // x = (DATA->delta_x)*ix - (DATA->x_size)/2.0;
	 for(iy=0; iy<=DATA->ny; iy++)
	   {
	     //  cout << "ix=" << ix << " iy=" << iy << endl;
 	     (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
 	     (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
 	     (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
 	     (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
  	     (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
  	     (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
  	     (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
  	     (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
 	     (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
 	     (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
 	     (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
 	     (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
 	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
 	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
 	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
 	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
 	     (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
 	     (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
 	     (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
 	     (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	   }
	 //	 sleep(1);
       }
  
     double density, dummy1, dummy2, dummy3, eta0, x0, y0;
     // loop over the whole lattice and initialize values:
     //     bytes_read=fscanf(profile,"%s %s %s %s %d %s %d %s %d %s %lf %s %lf %s %lf",
     //		       &dummy,&dummy,&dummy,&dummy,&neta,&dummy,&nx,&dummy,&ny,&dummy,
     //		       &deta,&dummy,&dx,&dummy,&dy);
   
     profile >> dummy >>dummy >>dummy >>dummy >>neta >> dummy >> nx >> dummy >> ny >> dummy >> deta >> dummy >> dx >> dummy >> dy;

     cout << "DATA->nx=" << DATA->nx << endl;

     cout << "neta=" << DATA->neta << ", nx=" << nx << ", ny=" << ny << ", deta=" << deta << ", dx=" << dx << ", dy=" << dy << endl;
     for(ieta=0; ieta<DATA->neta*size; ieta++)
       {
	 //	 cout << "ieta=" << ieta << endl;
	 ieta2 = ieta - rank*DATA->neta;
	 for(ix=0; ix<=DATA->nx; ix++)
	   {
	     for(iy=0; iy<=DATA->ny; iy++)
	       {

		 profile >> dummy1 >> dummy2 >> dummy3 >> density;
		 //bytes_read=fscanf(profile,"%lf %lf %lf %lf",&dummy1,&dummy2,&dummy3,&density);
		 
		 if (ieta==0 && ix == 0 && iy == 0)
		   {
		     DATA->eta_size = -dummy1*2;
		     DATA->x_size = -dummy2*2;
		     DATA->y_size = -dummy3*2;
		     cout << "eta_size=" << -dummy1 << ", x_size=" << -dummy2 << ", y_size=" << -dummy3 << endl;
		   }

// 		 if (ix==0 && iy==0)
// 		   {
// 		     cout << "eta=" << dummy1 << ", x=" << dummy2 << ", y=" << dummy3 << endl; 
// 		     cout << "density=" << density << endl;
// 		   }

		 if ( ieta<DATA->neta*rank || ieta>DATA->neta*(rank+1)-1 )
		   {
		     //if (ix==0 && iy==0) cout << "skipping ieta=" << ieta << " on rank " << rank << ". size=" << size << endl;
		     continue;
		   }

		 
		//  if (ieta==1 && ix==1 && iy==1)
// 		   {
// 		     DATA->delta_eta = dummy1-eta0;
// 		     DATA->delta_x = dummy2-x0;
// 		     DATA->delta_y = dummy3-y0;
// 		     cout << "x_max=" << DATA->x_size << ", y_max=" << DATA->y_size 
// 			  << ", eta_max=" << DATA->eta_size << endl;
// 		     cout << "delta_x_=" << DATA->delta_x << ", delta_y=" << DATA->delta_y << ", delta_eta=" 
// 			  << DATA->delta_eta << endl;
// 		     //exit(1);
// 		   }

		 s = density*DATA->sFactor;
		
		 //factors:
		 //AuAu@200, N = 5.16347  
		 //Auau@63, N = 5.75164
		 //CuCu@63 N = 5.700148477
		 //CuCu@200, N = 5.11928
		 //PbPb@2760, N = 5.03831
		 //PbPb@5500, N = 5.15075
		 
		 //rhob = DATA->rhoB0/epsilon0*s;
		 rhob=0.;
		 if ( DATA->whichEOS==1 && s>96.288 )
		   {
		     epsilon=(pow((4./169./pow(PI,2.)),1./3.)/30.*pow(22.5,(4./3.))*
			      pow(s,(4./3.))*0.197+0.0)/hbarc;
		   }
		 else if (s>0.2)
		   {
		     epsilon = eos->findRoot(&EOS::ssolve, rhob, s, 1.15*rhob+0.001, 300.,0.001);
	   //cout << epsilon << " " <<  eos->findRoot(&EOS::Tsolve, rhob, eos->interpolate2(epsilon, rho, 1), 1.15*rhob+0.001, 300.,0.001) << endl;
		   }
		 else 
		   {
		     epsilon=0.00000000001;
		   }
 	 
		 if (epsilon<0.00000000001)
		   epsilon = 0.00000000001;

		 // intial pressure distribution:
		 p = eos->p_func(epsilon, rhob);
		 
		 // set all values in the grid element:
		 (*arena)[ix][iy][ieta2].epsilon = epsilon;
		 (*arena)[ix][iy][ieta2].rhob = rhob;
		 (*arena)[ix][iy][ieta2].p = p;
		 (*arena)[ix][iy][ieta2].trouble = 0;
		 
		 if (DATA->whichEOS==1)
		   {
		     (*arena)[ix][iy][ieta2].T = eos->interpolate(epsilon, rhob, 0);
		     (*arena)[ix][iy][ieta2].mu = eos->interpolate(epsilon, rhob, 1);
		   }
		 else if (DATA->whichEOS>1)
		   {
		     (*arena)[ix][iy][ieta2].T = eos->interpolate2(epsilon, rhob, 1);
		     (*arena)[ix][iy][ieta2].mu = 0.0;
		   }
		 
		//  if (ix==DATA->nx/2+1 && iy==DATA->ny/2+1 && ieta2==DATA->neta/2)
// 		   {
// 		     cout << "s=" << s << ", eps=" << epsilon*0.1973 << ", P=" << p*0.1973 << ", T=" << eos->interpolate2(epsilon, rho, 1)*0.1973 << endl; 
// 		     cout << eos->interpolate2(epsilon, rho, 1)*0.1973*s-p*0.1973-epsilon*0.1973 << endl;
// 		     // exit(1);
// 		   }

		 (*arena)[ix][iy][ieta2].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta2].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta2].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta2].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta2].theta_u = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta2].pi_b = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta2].prev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta2].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta2].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta2].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta2].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta2].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta2].pprevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta2].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta2].prevPimunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta2].pprevPimunu = util->cube_malloc(1, 5,4);
		 
		 (*arena)[ix][iy][ieta2].W_prev = util->mtx_malloc(5,4);
// 		 int taui;
// 		 for (taui=0; taui<3; taui++)
// 		   {
// 		     (*arena)[ix][iy][ieta2].tauF[taui] = 0.;
// 		   }
		 /* 0 is actual, any others are RK intermediate values */
		 /* for a shock test */
		 //u[0] = (*arena)[ix][iy][ieta2].u[0][0] = cosh(10.0);
		 //u[3] = (*arena)[ix][iy][ieta2].u[0][3] = sinh(10.0);
		 
		 /* for HIC */
		 u[0] = (*arena)[ix][iy][ieta2].u[0][0] = 1.0;
		 u[3] = (*arena)[ix][iy][ieta2].u[0][3] = 0.0;
		 u[1] = (*arena)[ix][iy][ieta2].u[0][1] = 0.0;
		 u[2] = (*arena)[ix][iy][ieta2].u[0][2] = 0.0;
		 
		 (*arena)[ix][iy][ieta2].prev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta2].prev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta2].prev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta2].prev_u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta2].pprev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta2].pprev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta2].pprev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta2].pprev_u[0][2] = 0.0;
		 
		 (*arena)[ix][iy][ieta2].pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta2].prev_pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta2].pprev_pi_b[0] = 0.0;
		 
		 for(mu=0; mu<4; mu++)
		   {
		     /* baryon density */
		     (*arena)[ix][iy][ieta2].TJb[0][4][mu] = rhob*u[mu];
		     
		     for(nu=0; nu<4; nu++)
		       {
			 (*arena)[ix][iy][ieta2].TJb[0][nu][mu] 
			   = (epsilon + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];
			 
			 (*arena)[ix][iy][ieta2].Wmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta2].prevWmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta2].pprevWmunu[0][nu][mu] = (double) 0.0;
			 
			 (*arena)[ix][iy][ieta2].Pimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta2].prevPimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta2].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 
	       }}
	 //cout << "ok " << ieta << " is less than " << DATA->neta*size << endl;
	 //cout << "DATA->neta=" << DATA->neta << " , size=" << size << endl; 
	 //exit(1);
       }/* ix, iy, ieta */
   }

else if (DATA->Initial_profile==7) //read in the profile from file - IPSat initial conditions
   {
     double va, test;
     double s, r1, r2, W, nBinary, nWounded, W0;
     // impact parameter:
     double b=DATA->b;
     int i;
     int bytes_read;
     string dummy;
     int nx, ny, neta;
     int y_size, x_size;
     double dx, dy, deta;
     int ieta2;
     size = DATA->size;
     cout << "size=" << size << endl;
     
     cout << " ----- information on initial distribution -----" << endl;
     cout << "file name used: " << DATA->initName << endl;
  
     ifstream profile(DATA->initName.c_str());

     //     profile = fopen(DATA->initName, "r");
     // read the first line with general info
    
     //          exit(1);

     for(ix=0; ix<=DATA->nx; ix++)
       {
	 // x = (DATA->delta_x)*ix - (DATA->x_size)/2.0;
	 for(iy=0; iy<=DATA->ny; iy++)
	   {
	     //  cout << "ix=" << ix << " iy=" << iy << endl;
 	     (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
 	     (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
 	     (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
 	     (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
  	     (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
  	     (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 5,4);
  	     (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
  	     (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 5,4);
 	     (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
 	     (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 5,4);
 	     (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
 	     (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 5,4);
 	     (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
 	     (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
 	     (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
 	     (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
 	     (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
 	     (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
 	     (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
 	     (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
	   }
	 //	 sleep(1);
       }
  
     double density, dummy1, dummy2, dummy3, eta0, x0, y0;
     // loop over the whole lattice and initialize values:
     //     bytes_read=fscanf(profile,"%s %s %s %s %d %s %d %s %d %s %lf %s %lf %s %lf",
     //		       &dummy,&dummy,&dummy,&dummy,&neta,&dummy,&nx,&dummy,&ny,&dummy,
     //		       &deta,&dummy,&dx,&dummy,&dy);
   
     profile >> dummy >>dummy >>dummy >>dummy >>neta >> dummy >> nx >> dummy >> ny >> dummy >> deta >> dummy >> dx >> dummy >> dy;

     cout << "DATA->nx=" << DATA->nx << endl;

     cout << "neta=" << DATA->neta << ", nx=" << nx << ", ny=" << ny << ", deta=" << deta << ", dx=" << dx << ", dy=" << dy << endl;
     for(ieta=0; ieta<DATA->neta*size; ieta++)
       {
	 //	 cout << "ieta=" << ieta << endl;
	 ieta2 = ieta - rank*DATA->neta;
	 for(ix=0; ix<=DATA->nx; ix++)
	   {
	     for(iy=0; iy<=DATA->ny; iy++)
	       {

		 //		 bytes_read=fscanf(profile,"%lf %lf %lf %lf",&dummy1,&dummy2,&dummy3,&density);
		 profile >> dummy1 >> dummy2 >> dummy3 >> density;
		 //cout << density << endl;
	
		 if (ieta==0 && ix == 0 && iy == 0)
		   {
		     DATA->eta_size = -dummy1*2;
		     DATA->x_size = -dummy2*2;
		     DATA->y_size = -dummy3*2;
		     cout << "eta_size=" << -dummy1 << ", x_size=" << -dummy2 << ", y_size=" << -dummy3 << endl;
		   }

// 		 if (ix==0 && iy==0)
// 		   {
// 		     cout << "eta=" << dummy1 << ", x=" << dummy2 << ", y=" << dummy3 << endl; 
// 		     cout << "density=" << density << endl;
// 		   }

		 if ( ieta<DATA->neta*rank || ieta>DATA->neta*(rank+1)-1 )
		   {
		     //if (ix==0 && iy==0) cout << "skipping ieta=" << ieta << " on rank " << rank << ". size=" << size << endl;
		     continue;
		   }

		 
		//  if (ieta==1 && ix==1 && iy==1)
// 		   {
// 		     DATA->delta_eta = dummy1-eta0;
// 		     DATA->delta_x = dummy2-x0;
// 		     DATA->delta_y = dummy3-y0;
// 		     cout << "x_max=" << DATA->x_size << ", y_max=" << DATA->y_size 
// 			  << ", eta_max=" << DATA->eta_size << endl;
// 		     cout << "delta_x_=" << DATA->delta_x << ", delta_y=" << DATA->delta_y << ", delta_eta=" 
// 			  << DATA->delta_eta << endl;
// 		     //exit(1);
// 		   }

		 epsilon= density*DATA->sFactor/hbarc;
		
		 //factors:
		 // here the factors should be tuned by hand
		 
		 rhob=0.;
		 if (epsilon<0.00000000001)
		   epsilon = 0.00000000001;
		 
		 // intial pressure distribution:
		 p = eos->p_func(epsilon, rhob);
		 
		 // set all values in the grid element:
		 (*arena)[ix][iy][ieta2].epsilon = epsilon;
		 (*arena)[ix][iy][ieta2].rhob = rhob;
		 (*arena)[ix][iy][ieta2].p = p;
		 (*arena)[ix][iy][ieta2].trouble = 0;
		 
		 if (DATA->whichEOS==1)
		   {
		     (*arena)[ix][iy][ieta2].T = eos->interpolate(epsilon, rhob, 0);
		     (*arena)[ix][iy][ieta2].mu = eos->interpolate(epsilon, rhob, 1);
		   }
		 else if (DATA->whichEOS>1)
		   {
		     (*arena)[ix][iy][ieta2].T = eos->interpolate2(epsilon, rhob, 1);
		     (*arena)[ix][iy][ieta2].mu = 0.0;
		   }
		 
		//  if (ix==DATA->nx/2+1 && iy==DATA->ny/2+1 && ieta2==DATA->neta/2)
// 		   {
// 		     cout << "s=" << s << ", eps=" << epsilon*0.1973 << ", P=" << p*0.1973 << ", T=" << eos->interpolate2(epsilon, rho, 1)*0.1973 << endl; 
// 		     cout << eos->interpolate2(epsilon, rho, 1)*0.1973*s-p*0.1973-epsilon*0.1973 << endl;
// 		     // exit(1);
// 		   }

		 (*arena)[ix][iy][ieta2].TJb = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta2].dUsup = util->cube_malloc(rk_order+1, 4,4);
		 (*arena)[ix][iy][ieta2].u = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta2].a = util->mtx_malloc(rk_order+1, 4);
		 (*arena)[ix][iy][ieta2].theta_u = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta2].pi_b = util->vector_malloc(rk_order+1);
		 (*arena)[ix][iy][ieta2].prev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta2].pprev_pi_b = util->vector_malloc(1);
		 (*arena)[ix][iy][ieta2].prev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta2].pprev_u = util->mtx_malloc(1, 4);
		 (*arena)[ix][iy][ieta2].Wmunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta2].prevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta2].pprevWmunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta2].Pimunu = util->cube_malloc(rk_order+1, 5,4);
		 (*arena)[ix][iy][ieta2].prevPimunu = util->cube_malloc(1, 5,4);
		 (*arena)[ix][iy][ieta2].pprevPimunu = util->cube_malloc(1, 5,4);
		 
		 (*arena)[ix][iy][ieta2].W_prev = util->mtx_malloc(5,4);
// 		 int taui;
// 		 for (taui=0; taui<3; taui++)
// 		   {
// 		     (*arena)[ix][iy][ieta2].tauF[taui] = 0.;
// 		   }
		 /* 0 is actual, any others are RK intermediate values */
		 /* for a shock test */
		 //u[0] = (*arena)[ix][iy][ieta2].u[0][0] = cosh(10.0);
		 //u[3] = (*arena)[ix][iy][ieta2].u[0][3] = sinh(10.0);
		 
		 /* for HIC */
		 u[0] = (*arena)[ix][iy][ieta2].u[0][0] = 1.0;
		 u[3] = (*arena)[ix][iy][ieta2].u[0][3] = 0.0;
		 u[1] = (*arena)[ix][iy][ieta2].u[0][1] = 0.0;
		 u[2] = (*arena)[ix][iy][ieta2].u[0][2] = 0.0;
		 
		 (*arena)[ix][iy][ieta2].prev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta2].prev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta2].prev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta2].prev_u[0][2] = 0.0;
		 (*arena)[ix][iy][ieta2].pprev_u[0][0] = 1.0;
		 (*arena)[ix][iy][ieta2].pprev_u[0][3] = 0.0;
		 (*arena)[ix][iy][ieta2].pprev_u[0][1] = 0.0;
		 (*arena)[ix][iy][ieta2].pprev_u[0][2] = 0.0;
		 
		 (*arena)[ix][iy][ieta2].pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta2].prev_pi_b[0] = 0.0;
		 (*arena)[ix][iy][ieta2].pprev_pi_b[0] = 0.0;
		 
		 for(mu=0; mu<4; mu++)
		   {
		     /* baryon density */
		     (*arena)[ix][iy][ieta2].TJb[0][4][mu] = rhob*u[mu];
		     
		     for(nu=0; nu<4; nu++)
		       {
			 (*arena)[ix][iy][ieta2].TJb[0][nu][mu] 
			   = (epsilon + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];
			 
			 (*arena)[ix][iy][ieta2].Wmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta2].prevWmunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta2].pprevWmunu[0][nu][mu] = (double) 0.0;
			 
			 (*arena)[ix][iy][ieta2].Pimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta2].prevPimunu[0][nu][mu] = (double) 0.0;
			 (*arena)[ix][iy][ieta2].pprevPimunu[0][nu][mu] = (double) 0.0;
		       }/* nu */
		   }/* mu */
		 
	       }}
	 //cout << "ok " << ieta << " is less than " << DATA->neta*size << endl;
	 //cout << "DATA->neta=" << DATA->neta << " , size=" << size << endl; 
	 //exit(1);
       }/* ix, iy, ieta */
   }


//  fclose(e_file);
//  fclose(e2_file);
//  fclose(e3_file);
//  fclose(e4_file);
//  fclose(e5_file);
 psi_file.close();
 twodplot_file.close();
 twodplotA_file.close();
 twodplotB_file.close();
 geometry_file.close();
 geometry_file_rn.close();

 //system("gnuplot 2dplot.plt");
 cout << "initial distribution done." << endl;
 return 1;
}/* InitTJb*/
 

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

 
