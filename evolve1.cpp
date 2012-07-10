#include "evolve.h"

using namespace std;

Evolve::Evolve(EOS *eosIn)
{
  eos = new EOS;
  eos = eosIn;
  grid = new Grid;
  reconst = new Reconst(eosIn, grid);
  util = new Util;
}

// destructor
Evolve::~Evolve()
{
  delete eos;
  delete reconst;
  delete grid;
  delete util;
}

int Evolve::EvolveIt(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank)
{
/* implement Kurganov-Tadmor */
 int ix, iy, ieta, nx, ny, neta, it, itmax, rk_flag, flag, cent_eta;
 double dt, tau0, tau, x;

 FILE *v2_file;
 char* v2_name = "aniso.dat";
 v2_file = fopen(v2_name, "w");
 fprintf(v2_file,"tau [fm], momentum anisotropy\n");
 fclose(v2_file);
 FILE *ecc_file;
 char* ecc_name = "eccentricity.dat";
 ecc_file = fopen(ecc_name, "w");
 fprintf(ecc_file,"tau [fm], spatial anisotropy\n");
 fclose(ecc_file);
 FILE *t_file;
 char* t_name = "tauf.dat";
 t_file = fopen(t_name, "w");
 fprintf(t_file,"x_f [fm], tau_f [fm], Sx[fm^3], Stau [fm^3] \n");
 fclose(t_file);
 FILE *t2_file;
 char* t2_name = "taufx.dat";
 t2_file = fopen(t2_name, "w");
 fprintf(t2_file,"x_f [fm], tau_f [fm], Sx[fm^3], Stau [fm^3] \n");
 fclose(t2_file);
 FILE *t3_file;
 char* t3_name = "taufy.dat";
 t3_file = fopen(t3_name, "w");
 fprintf(t3_file,"y_f [fm], tau_f [fm], Sy[fm^3], Stau [fm^3] \n");
 fclose(t3_file);

 char *buf;
 buf = util->char_malloc(40);
 FILE *s_file;
 char* s_name;
 s_name = util->char_malloc(100);
 
 sprintf (buf, "%d", rank);
 strcat(s_name, "surface");
 strcat(s_name,buf);
 strcat(s_name, ".dat");
 
 s_file = fopen(s_name, "w");
 fclose(s_file);

 FILE *out_file;
 char* out_name = "evolution.dat";
 out_file = fopen(out_name, "w");
 fprintf(out_file,"");
 fclose(out_file);

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
   tau = tau0 + dt*it;
   //fprintf(stderr, "Starting time step %d/%d on rank %d.\n", it, itmax, rank);
   if(it==0) 
     {
       // storePreviousT(tau, DATA, arena);
       storePreviousEpsilon2(tau, DATA, arena);
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
       grid->PrintxEpsilon(arena, DATA, tau, size);
       grid->ComputeEccentricity(DATA, arena, tau);
       grid->ComputeAnisotropy(DATA, arena, tau);
       grid->OutputXY(arena, DATA, eos, tau, size, rank);
     }

   if(it%10==0 && DATA->outputEvolutionData) 
     {
       grid->OutputEvolutionDataXYZ(arena, DATA, eos, tau, size, rank); 
     }
   
   /* execute rk steps */
   flag = AdvanceRK(tau, DATA, arena, Lneighbor, Rneighbor, size, rank);
   
   UpdateArena(tau, DATA, arena);
   
   //check energy conservation
   //grid->ComputeEnergyConservation(DATA, arena, tau);
   //storePreviousT(tau, DATA, arena);

   //determine freeze-out surface
    if (it%facTau==0 && it>0) 
      {
	if (DATA->freezeOutMethod == 1)
	  FindFreezeOutSurface(tau, DATA, arena, size, rank);
	else if (DATA->freezeOutMethod == 2)
	  FindFreezeOutSurface2(tau, DATA, arena, size, rank);
	storePreviousEpsilon2(tau, DATA, arena);
      } 
    
    if (rank == 0) fprintf(stderr, "Done time step %d/%d.\n", it, itmax);
 
    if (it==itmax) grid->PrintAxy(DATA, arena, tau);

  }/* it */ 
 fprintf(stderr,"SUM=%f\n", SUM);
 return 1; /* successful */

}/* Evolve */

void Evolve::MPISendReceive(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag)
{
  // this sends and receives information from neighboring cells in the next processor in the eta direction 
  // and stores it in Lneighbor and Rneighbor (unless the processor is really at the edge of the total grid

  int ix, iy, ieta, nx, ny, neta, i, alpha, iflag;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;
  int sizeOfData = 5*2*(nx+1)*(ny+1);
  int position;
  double *package;
  double *package2;
  
  package = new double[sizeOfData];
  package2 = new double[sizeOfData];

  // receive from the right / send to the left
  int from = rank+1;
  int to = rank-1;
  // packing the package to send
  if ( rank != 0 )
    {
      //      cout << " sending to the left on rank " << rank << endl;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(alpha=0; alpha<5; alpha++)
		{
		  for(i=0; i<2; i++)
		    {
		      //if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
		      //cout << "arena[ix][iy][0].TJb[i][alpha][0]=" << arena[ix][iy][0].TJb[i][alpha][0] << endl;
		      position = (i+2*(alpha+(5*(ix + (nx*iy)))));
		      package[position] = arena[ix][iy][0].TJb[i][alpha][0];
		      package2[position] = arena[ix][iy][1].TJb[i][alpha][0];
		    }
		}
	    }
	}
      MPI::COMM_WORLD.Send(package,sizeOfData,MPI::DOUBLE,to,10+rk_flag);
      MPI::COMM_WORLD.Send(package2,sizeOfData,MPI::DOUBLE,to,20+rk_flag);
      //cout << " done sending to the left on rank " << rank << endl;
    }
  // receiving and unwrapping the package
  if ( rank != size-1 )
    {  
      MPI::COMM_WORLD.Recv(package,sizeOfData,MPI::DOUBLE,from,10+rk_flag);
      MPI::COMM_WORLD.Recv(package2,sizeOfData,MPI::DOUBLE,from,20+rk_flag);
      //cout << " receiving from the right on rank " << rank << endl;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(alpha=0; alpha<5; alpha++)
		{
		  for(i=0; i<2; i++)
		    {
		      position = (i+2*(alpha+(5*(ix + (nx*iy)))));
		      //	      if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
		      //cout << "Rneighbor[ix][iy][0].TJb[i][alpha][0]=" << package[position] << endl;
		      Rneighbor[ix][iy][0].TJb[i][alpha][0] = package[position];
		      Rneighbor[ix][iy][1].TJb[i][alpha][0] = package2[position];
		    }
		}
	    }
	}
      //cout << " done receiving from the right on rank " << rank << endl;
    }


  // receive from the left / send to the right
  from = rank-1;
  to = rank+1;
  // packing the package to send
  if ( rank != size-1 )
    {
      //cout << " **"<<rank<<"** "<<  " sending to the right " << endl;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(alpha=0; alpha<5; alpha++)
		{
		  for(i=0; i<2; i++)
		    {
		      position = (i+2*(alpha+(5*(ix + (nx*iy)))));
		      package[position] = arena[ix][iy][neta-1].TJb[i][alpha][0];
		      package2[position] = arena[ix][iy][neta-2].TJb[i][alpha][0];
		    }
		}
	    }
	}
      MPI::COMM_WORLD.Send(package,sizeOfData,MPI::DOUBLE,to,30+rk_flag);
      MPI::COMM_WORLD.Send(package2,sizeOfData,MPI::DOUBLE,to,40+rk_flag);
      //cout << " **"<<rank<<"** "<<  " done sending to the right " << endl;
    }
  // receiving and unwrapping the package
  if ( rank != 0 )
    {  
      //cout << " ++"<<rank<<"++ "<<  " receiving from the left " << endl;
      MPI::COMM_WORLD.Recv(package,sizeOfData,MPI::DOUBLE,from,30+rk_flag);
      MPI::COMM_WORLD.Recv(package2,sizeOfData,MPI::DOUBLE,from,40+rk_flag);
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(alpha=0; alpha<5; alpha++)
		{
		  for(i=0; i<2; i++)
		    {
		      position = (i+2*(alpha+(5*(ix + (nx*iy)))));
		      Lneighbor[ix][iy][0].TJb[i][alpha][0] = package[position];
		      Lneighbor[ix][iy][1].TJb[i][alpha][0] = package2[position];
		    }
		}
	    }
	}
      //cout << " ++"<<rank<<"++ "<<  " receiving from the left " << endl;
    }
  delete[] package;
  delete[] package2;
}//end MPISendReceive


void Evolve::storePreviousEpsilon(double tau, InitData *DATA, Grid ***arena)
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
	      if (tau>tau0+8*DATA->delta_tau)
		{
		  arena[ix][iy][ieta].epsilon_prev[9]=arena[ix][iy][ieta].epsilon_prev[8];
		  arena[ix][iy][ieta].u_prev[9][0]=arena[ix][iy][ieta].u_prev[8][0];
		  arena[ix][iy][ieta].u_prev[9][1]=arena[ix][iy][ieta].u_prev[8][1];
		  arena[ix][iy][ieta].u_prev[9][2]=arena[ix][iy][ieta].u_prev[8][2];
		  arena[ix][iy][ieta].u_prev[9][3]=arena[ix][iy][ieta].u_prev[8][3];
		}
	      if (tau>tau0+7*DATA->delta_tau)
		{
		  arena[ix][iy][ieta].epsilon_prev[8]=arena[ix][iy][ieta].epsilon_prev[7];
		  arena[ix][iy][ieta].u_prev[8][0]=arena[ix][iy][ieta].u_prev[7][0];
		  arena[ix][iy][ieta].u_prev[8][1]=arena[ix][iy][ieta].u_prev[7][1];
		  arena[ix][iy][ieta].u_prev[8][2]=arena[ix][iy][ieta].u_prev[7][2];
		  arena[ix][iy][ieta].u_prev[8][3]=arena[ix][iy][ieta].u_prev[7][3];
		}
	      if (tau>tau0+6*DATA->delta_tau)
		{
		  arena[ix][iy][ieta].epsilon_prev[7]=arena[ix][iy][ieta].epsilon_prev[6];
		  arena[ix][iy][ieta].u_prev[7][0]=arena[ix][iy][ieta].u_prev[6][0];
		  arena[ix][iy][ieta].u_prev[7][1]=arena[ix][iy][ieta].u_prev[6][1];
		  arena[ix][iy][ieta].u_prev[7][2]=arena[ix][iy][ieta].u_prev[6][2];
		  arena[ix][iy][ieta].u_prev[7][3]=arena[ix][iy][ieta].u_prev[6][3];
		}
	      if (tau>tau0+5*DATA->delta_tau)
		{
		  arena[ix][iy][ieta].epsilon_prev[6]=arena[ix][iy][ieta].epsilon_prev[5];
		  arena[ix][iy][ieta].u_prev[6][0]=arena[ix][iy][ieta].u_prev[5][0];
		  arena[ix][iy][ieta].u_prev[6][1]=arena[ix][iy][ieta].u_prev[5][1];
		  arena[ix][iy][ieta].u_prev[6][2]=arena[ix][iy][ieta].u_prev[5][2];
		  arena[ix][iy][ieta].u_prev[6][3]=arena[ix][iy][ieta].u_prev[5][3];
		}
	      if (tau>tau0+4*DATA->delta_tau)
		{
		  arena[ix][iy][ieta].epsilon_prev[5]=arena[ix][iy][ieta].epsilon_prev[4];
	     	  arena[ix][iy][ieta].u_prev[5][0]=arena[ix][iy][ieta].u_prev[4][0];
		  arena[ix][iy][ieta].u_prev[5][1]=arena[ix][iy][ieta].u_prev[4][1];
		  arena[ix][iy][ieta].u_prev[5][2]=arena[ix][iy][ieta].u_prev[4][2];
		  arena[ix][iy][ieta].u_prev[5][3]=arena[ix][iy][ieta].u_prev[4][3];
		}
	      if (tau>tau0+3*DATA->delta_tau)
		{
		  arena[ix][iy][ieta].epsilon_prev[4]=arena[ix][iy][ieta].epsilon_prev[3];
	     	  arena[ix][iy][ieta].u_prev[4][0]=arena[ix][iy][ieta].u_prev[3][0];
		  arena[ix][iy][ieta].u_prev[4][1]=arena[ix][iy][ieta].u_prev[3][1];
		  arena[ix][iy][ieta].u_prev[4][2]=arena[ix][iy][ieta].u_prev[3][2];
		  arena[ix][iy][ieta].u_prev[4][3]=arena[ix][iy][ieta].u_prev[3][3];
		}
	      if (tau>tau0+2*DATA->delta_tau)
		{
		  arena[ix][iy][ieta].epsilon_prev[3]=arena[ix][iy][ieta].epsilon_prev[2];
	     	  arena[ix][iy][ieta].u_prev[3][0]=arena[ix][iy][ieta].u_prev[2][0];
		  arena[ix][iy][ieta].u_prev[3][1]=arena[ix][iy][ieta].u_prev[2][1];
		  arena[ix][iy][ieta].u_prev[3][2]=arena[ix][iy][ieta].u_prev[2][2];
		  arena[ix][iy][ieta].u_prev[3][3]=arena[ix][iy][ieta].u_prev[2][3];
		}
	      if (tau>tau0+DATA->delta_tau)
		{
		  arena[ix][iy][ieta].epsilon_prev[2]=arena[ix][iy][ieta].epsilon_prev[1];
	     	  arena[ix][iy][ieta].u_prev[2][0]=arena[ix][iy][ieta].u_prev[1][0];
		  arena[ix][iy][ieta].u_prev[2][1]=arena[ix][iy][ieta].u_prev[1][1];
		  arena[ix][iy][ieta].u_prev[2][2]=arena[ix][iy][ieta].u_prev[1][2];
		  arena[ix][iy][ieta].u_prev[2][3]=arena[ix][iy][ieta].u_prev[1][3];
		}
	      if (tau>tau0)
		{
		  arena[ix][iy][ieta].epsilon_prev[1]=arena[ix][iy][ieta].epsilon_prev[0];
	      	  arena[ix][iy][ieta].u_prev[1][0]=arena[ix][iy][ieta].u_prev[0][0];
		  arena[ix][iy][ieta].u_prev[1][1]=arena[ix][iy][ieta].u_prev[0][1];
		  arena[ix][iy][ieta].u_prev[1][2]=arena[ix][iy][ieta].u_prev[0][2];
		  arena[ix][iy][ieta].u_prev[1][3]=arena[ix][iy][ieta].u_prev[0][3];
		}
	      arena[ix][iy][ieta].epsilon_prev[0]=arena[ix][iy][ieta].epsilon;
	      arena[ix][iy][ieta].u_prev[0][0]=arena[ix][iy][ieta].u[0][0];
	      arena[ix][iy][ieta].u_prev[0][1]=arena[ix][iy][ieta].u[0][1];
	      arena[ix][iy][ieta].u_prev[0][2]=arena[ix][iy][ieta].u[0][2];
	      arena[ix][iy][ieta].u_prev[0][3]=arena[ix][iy][ieta].u[0][3];
	    }
	}
    }
}

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
	      arena[ix][iy][ieta].epsilon_prev[facTau-1]=arena[ix][iy][ieta].epsilon;
	      arena[ix][iy][ieta].u_prev[facTau-1][0]=arena[ix][iy][ieta].u[0][0];
	      arena[ix][iy][ieta].u_prev[facTau-1][1]=arena[ix][iy][ieta].u[0][1];
	      arena[ix][iy][ieta].u_prev[facTau-1][2]=arena[ix][iy][ieta].u[0][2];
	      arena[ix][iy][ieta].u_prev[facTau-1][3]=arena[ix][iy][ieta].u[0][3];
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


int Evolve::UpdateArena(double tau, InitData *DATA, Grid ***arena)
{
 int rk_flag, ix, iy, ieta, nx, ny, neta, flag, alpha, mu, rk_order;
 
 nx = DATA->nx;
 ny = DATA->ny;
 neta = DATA->neta;
 rk_order = DATA->rk_order;

 for(ix=0; ix<=nx; ix++)
  {
   for(iy=0; iy<=ny; iy++)
    {
     for(ieta=0; ieta<neta; ieta++)
      {
       arena[ix][iy][ieta].p = arena[ix][iy][ieta].p_t;
       arena[ix][iy][ieta].epsilon = arena[ix][iy][ieta].epsilon_t;
       arena[ix][iy][ieta].rhob = arena[ix][iy][ieta].rhob_t;
       for(mu=0; mu<4; mu++)
        {
	 arena[ix][iy][ieta].u[0][mu] = 
	               arena[ix][iy][ieta].u[rk_order][mu]; 
	 for(alpha=0; alpha<5; alpha++)
	 {
	  arena[ix][iy][ieta].TJb[0][alpha][mu] = 
	               arena[ix][iy][ieta].TJb[rk_order][alpha][mu]; 
	 }}/* mu, alpha */
      }}}/* ix, iy, ieta */

}/* update arena */

int Evolve::AdvanceRK(double tau, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank)
{
 int rk_flag, ix, iy, ieta, nx, ny, neta, flag;
 
 nx = DATA->nx;
 ny = DATA->ny;
 neta = DATA->neta;

   /* this executes RK steps */
   for(rk_flag=0; rk_flag<DATA->rk_order; rk_flag++) 
    {
      //fprintf(stderr, "Evolve rk_flag = %d on rank %d\n", rk_flag,rank);
      // send and receive the two neighboring layers of cells
      MPISendReceive(DATA, arena, Lneighbor, Rneighbor, size, rank,rk_flag);
      
      //      cout << "going through Advance rk_flag=" << rk_flag << " rank " << rank << endl;
      
      for(ix=0; ix<=nx; ix++)
       {
        for(iy=0; iy<=ny; iy++)
         {
          for(ieta=0; ieta<neta; ieta++)
           {
	     flag = Advance(tau, DATA, &(arena[ix][iy][ieta]), &(Lneighbor[ix][iy][0]), &(Rneighbor[ix][iy][0])
			    , &(Lneighbor[ix][iy][1]), &(Rneighbor[ix][iy][1]), rk_flag, size, rank);
	    if(flag==0) 
	     { 
              //PrintArena(DATA, arena, ieta, tau); 
	      //PrintdEdEta(DATA, arena);
	      return 0; 
	     }
           }/* ieta */
         }/*iy */
       }/* ix */
      //      fprintf(stderr, "Done %d-th rk step on rank %d.\n", rk_flag,rank);
    }/* rk_flag */
 return 1; /* successful */
}/* AdvanceRK */


int Evolve::Advance(double tau_it, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, 
		    int rk_flag, int size, int rank)
{
 static Grid grid_rk;
 static double **qirk, *qi, *rhs;
 double tau_now, tau_next;
 int trk_flag, rk_order_m1, flag;
 int i, alpha, mu;
 static int ind=0;

 ind++;
 if(ind == 1)
 {
  qirk = util->mtx_malloc(5,4);
  qi = util->vector_malloc(5);
  rhs = util->vector_malloc(5);
  grid_rk.TJb = util->cube_malloc(DATA->rk_order,5,4);
  grid_rk.u = util->mtx_malloc(DATA->rk_order,4);
 } 
 
 if(rk_flag == 0) /* first rk step */
  {
   tau_now = tau_it;
   tau_next = tau_it + (DATA->delta_tau);
   /* this calculates k1*dt = f(t, u)*dt */
   MakeDeltaQI(tau_now, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, qi, rhs, DATA, rk_flag, size, rank);
  
   for(alpha=0; alpha<5; alpha++)
    {
     qirk[alpha][0] = qi[alpha] + rhs[alpha];
     if(!finite(qirk[alpha][0]))
      {
       fprintf(stderr, "qirk[%d][0] = %e is a nan.\n", alpha, qirk[alpha][0]);
       fprintf(stderr, "qi[%d] = %e\n", alpha, qi[alpha]);
       fprintf(stderr, "rhs[%d] = %e\n", alpha, rhs[alpha]);
      }
    }
   
   //if (rank==0) cout << "rhob =" <<  grid_pt->rhob << endl;
   flag=
     reconst->ReconstIt(&grid_rk, 0, tau_next, qirk, grid_pt,
			grid_pt->epsilon, grid_pt->rhob, DATA, rk_flag); 
   
   if(flag==0)
     {
       reconst->ReconstError("grid_rk", 0, rk_flag, qi, qirk, grid_pt);
    }/* flag == 0 */
   else if(flag != 0) 
     {
       if (flag == -1 && grid_rk.epsilon>1 )
	 {
	   cout << "on rank " << rank << " got eps=" << grid_rk.TJb[0][0][0] << endl;
	   cout << "on rank " << rank << " had eps=" << grid_pt->TJb[0][0][0] << endl;
	 }	   
       UpdateTJbRK(&grid_rk, grid_pt, DATA, rk_flag); 
     }


  }/* rk_flag == 0 */
 else /* rk step is finished since we only implement 2nd order rk */
  {
   tau_now = tau_it + (DATA->delta_tau);
   /* this calculates k2*dt = f(t+dt, u0+k1*dt)*dt */
   MakeDeltaQI(tau_now, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, qi, rhs, DATA, rk_flag, size, rank);
 
  /* at this point qi = u0 + k1*h */
  /* I need u0 + (k1*h + k2*h)/2 */
  /* qi + u0 + k2*h = 2u0 + (k1+k2)*h */

   for(alpha=0; alpha<5; alpha++)
   {
    qirk[alpha][0] = qi[alpha];
    qirk[alpha][0] += grid_pt->TJb[0][alpha][0]*tau_it + rhs[alpha];
    qirk[alpha][0] *= 0.5;
   }
   
   flag=reconst->ReconstIt(&grid_rk, 0, tau_now, qirk, grid_pt,
			   grid_pt->epsilon, grid_pt->rhob, DATA, rk_flag); 
   
   if(flag==0)
    {
      reconst->ReconstError("grid_rk", 0, rk_flag, qi, qirk, grid_pt);
    }/* flag == 0 */
   else if(flag != 0) UpdateTJbRK(&grid_rk, grid_pt, DATA, rk_flag); 

  }/* rk_flag == 1 */

 
 return 1; /* if successful */
}/* Advance */


/* Calculate the right-hand-side */
/*
 du/dt = (1/Delta x)(F_imh - F_iph) + D
 F_iph(a) = (1/2)(f_a(qiphR) + f_a(qiphL))- (1/2)aiph(qiphR - qiphL)
 F_imh(a) = (1/2)(f_a(qimhR) + f_a(qimhL))- (1/2)aimh(qimhR - qimhL)
 RK 2nd order Heun's rules:

 k1 = f(t, u);
 u1 = u + k1*h
 k2 = f(t+h, u1);
 u2 = u1 + k2*h = u + (k1+k2)*h;
 u_next = u + (1/2)(k1+k2)*h
        = (1/2)(u + u2);
*/
void Evolve::MakeDeltaQI(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			 Grid *Lneighbor2, Grid *Rneighbor2, double *qi, double *rhs, 
			 InitData *DATA, int rk_flag, int size, int rank) 
{
  double xjet, yjet, xtrigger, ytrigger;
  static double delta[4], sumf;
  static double tau_fac[4];
  static int alpha, i, rk_order_m1, nmax[4], flag;
  static double **DFmmp;
  static NbrQs NbrCells;
  static BdryCells HalfwayCells;
  static int ind=0;
   double x=grid_pt->position[1]*DATA->delta_x-DATA->x_size/2;
 double y=grid_pt->position[2]*DATA->delta_y-DATA->y_size/2;
 double eta;
 if(size>1)
   eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2+static_cast<double>(rank)/static_cast<double>(size)*DATA->eta_size;
 else
   eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2;
   
 ind++;
 if(ind==1)
  {
   InitTempGrids(&HalfwayCells, DATA->rk_order); 
   InitNbrQs(&NbrCells);
   DFmmp = util->mtx_malloc(5,4);
  }

 delta[1] = DATA->delta_x;
 delta[2] = DATA->delta_y;
 delta[3] = DATA->delta_eta;

 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta-1;

/* \partial_tau (tau Ttautau) + \partial_eta Tetatau 
           + \partial_x (tau Txtau) + \partial_y (tau Tytau) + Tetaeta = 0 */
/* \partial_tau (tau Ttaueta) + \partial_eta Teteta 
           + \partial_x (tau Txeta) + \partial_y (tau Txeta) + Tetatau = 0 */
/* \partial_tau (tau Txtau) + \partial_eta Tetax + \partial_x tau T_xx
  + \partial_y tau Tyx = 0 */

     tau_fac[1] = tau;
     tau_fac[2] = tau;
     tau_fac[3] = 1.0;
 
 for(alpha=0; alpha<5; alpha++) 
  {
   qi[alpha] = grid_pt->TJb[rk_flag][alpha][0]*tau;
  }/* get qi first */

 /* implement Kurganov-Tadmor scheme */
 GetQIs(tau, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, qi, &NbrCells, rk_flag, DATA, size, rank);
 
 flag = 
   MakeQIHalfs(qi, &NbrCells, &HalfwayCells, grid_pt, DATA);
 
 flag = 
   ConstHalfwayCells(tau, &HalfwayCells, qi, grid_pt, DATA, rk_flag, size, rank);
 
 MakeKTCurrents(tau, DFmmp, grid_pt, &HalfwayCells, rk_flag);

 for(alpha=0; alpha<5; alpha++) 
  {
   sumf = 0.0; 
   for(i=1; i<=3; i++)
    {
     sumf += DFmmp[alpha][i]/delta[i];
    }/* i */
    
   if(alpha==0)
    {
      sumf -= grid_pt->TJb[rk_flag][3][3];
    }
   else if(alpha==3)
    {
     sumf -= grid_pt->TJb[rk_flag][3][0];
    }

   //   cout << tau << " " << x << " " << " " << y << " " << eta << endl;

   int jetPosition=DATA->includeJet;
   int includeTrigger = DATA->includeTrigger;

   if ( DATA->includeJet>=1 )
     {
       double sTau;
       double etaWidth;
       etaWidth = asinh(0.2/tau);
       if (rk_flag==0)
	 sTau=tau;
       else if(rk_flag==1)
	 sTau=tau+DATA->delta_tau;
       
       if( jetPosition==1)
	 {
	   // adding for both RK steps
	   if((alpha==0 || alpha==1) && eta>-1 && eta < 1)
	     {
	       // initial position x_0=-6.5 fm, y_0=0 fm, EOS-L
	       xjet=sTau-0.4-6.5;
	       sumf += exp(-0.002900958784179574*pow (sTau,4))*(4.400274607409827*sTau - 8.208280106812962*pow (sTau,0.5)
								- 0.017616443714506047*pow (sTau,2) + 4.787674665037308*
								sin (1.6438127597829693 - 0.43552254762923354*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+y*y/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	 }
       else if (jetPosition==2)
	 {
	   // initial position x_0=-6 fm, y_0=+1 fm, EOS-L
	   xjet=sTau-0.4-6.25;
	   yjet=1.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0029509200139146478*pow (sTau,4))*(5.278151242285277*sTau - 8.444559533760918
								 *pow (sTau,0.5) - 0.13145550184759378*pow (sTau,2) 
								 + 5.132513932354498*sin (1.3176771285339544 - 0.4119506018253078*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0029491413291687354*pow (sTau,4))*(5.229129413110761*sTau - 8.565901304561994*pow (sTau,0.5) - 0.11421321941761275*pow (sTau,2) + 5.247913321697691*sin (1.7971957760975532 + 0.4081129909441432*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0018074622135051813*pow (sTau,4))*(-0.23598811781595871*sTau + 0.5541626880497722*pow (sTau,0.5) + 0.00032718366971014106*pow (sTau,2) - 0.27463077002612946*sin (2.1291457914893175 - 0.38798771916344515*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	 }
       else if (jetPosition==3)
	 {
	   xjet=sTau-0.4-6.;
	   yjet=2.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0035906166780202307*pow (sTau,4))*(5.0640915858911235*sTau - 8.837228675391408*pow (sTau,0.5) - 0.05234351009920327*pow (sTau,2) + 5.159927517626032*sin (1.6053207443255721 + 0.44520302760905406*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.003509062975494955*pow (sTau,4))*(5.828299397705474*sTau - 9.292956048608538*pow (sTau,0.5) - 0.1363354443218772*pow (sTau,2) + 5.564340280765995*sin (1.8041345571841076 + 0.428573490402628*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.003996847829075283*pow (sTau,4))*(-0.09149886641543215*sTau + 0.861042590511802*pow (sTau,0.5) - 1.3673488106527625*sin (0.2771896626053987 + 0.239904224781499*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	 }
       else if (jetPosition==4)
	 {
	   xjet=sTau-0.4-5.5;
	   yjet=3.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0022070817909056214*pow (sTau,4))*(1.2225746479646122*sTau - 3.3388320834061833*pow (sTau,0.5) - 0.026762107563405337*pow (sTau,2) + 2.977467633816383*sin (1.6216984712182159 - 0.1535019209578384*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.005091846471240839*pow (sTau,4))*(3.75523839883194*sTau - 8.95538490791062*pow (sTau,0.5) + 0.22394351008974536*pow (sTau,2) + 5.837401017993541*sin (1.1927699719121576 + 0.4891469058626414*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.003864430680285909*pow (sTau,4))*(2.8260675143096687 + 2.797578650088276*sTau - 0.4274871788208842*pow (sTau,2) - 7.311363420959585*sin (0.41333182808593466 + 0.3513714115809187*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	 }
       else if (jetPosition==5)
	 {
	   xjet=sTau-0.4-4.9;
	   yjet=4.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.002952204681235363*pow (sTau,4))*(0.980538761527065*sTau - 2.879485450981131*pow (sTau,0.5) - 0.043264051912027625*pow (sTau,2) + 2.7505739320219744*sin (1.567916771125228 + 0.00846187190357157*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0014176497960131133*pow (sTau,4))*(1.8293748566705894*sTau - 3.6900275499412873*pow (sTau,0.5) - 0.023590887484978628*pow (sTau,2) + 3.4430240378773633*sin (1.978813217520836 + 0.2471300228847286*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.005935461648444595*pow (sTau,4))*(0.5361892466090209 + 3.1232593963980597*sTau - 0.518309841204045*pow (sTau,2) + 5.096209631102419*sin (3.289770626761662 + 0.46875605498704376*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	 }
       else if (jetPosition==6)
	 {
	   xjet=sTau-0.4-4.;
	   yjet=5.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.003062415535182853*pow (sTau,4))*(1.7666475366776815*sTau - 4.026269551910104*pow (sTau,0.5) + 0.004166676592246689*pow (sTau,2) + 3.339782964137437*sin (1.4136239734383476 - 0.2781007470319588*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.002942876969252957*pow (sTau,4))*(2.0626741503281187*sTau - 4.231599576452959*pow (sTau,0.5) - 0.014612972396351212*pow (sTau,2) + 3.6417916923616764*sin (1.2448269477167029 - 0.2751141289045627*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.011754747085468394*pow (sTau,4))*(-0.7623185279608651 + 3.6514675381277626*sTau - 0.6980718072605644*pow (sTau,2) + 3.798715759890013*sin (3.0362952821637506 + 0.6380506141557616*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	 }
       else if (jetPosition==7)
	 {
	   xjet=sTau-0.4-3.;
	   yjet=6.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.06864751588043007*pow (sTau,4))*(-13.743120718017126*sTau + 11.116936151187065*pow (sTau,0.5) + 2.8747006909362116*pow (sTau,2) - 0.9329035581311574*sin (7.6723700300717415 + 2.070332200524868*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.06864751588043007*pow (sTau,4))*(-13.743120718017126*sTau + 11.116936151187065*pow (sTau,0.5) + 2.8747006909362116*pow (sTau,2) - 0.9329035581311574*sin (7.6723700300717415 + 2.070332200524868*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.04475541692624296*pow (sTau,4))*(-1.5384153966930805 + 5.001776521463333*sTau - 1.2641085658650624*pow (sTau,2) + 3.1365932521524305*sin (2.8379398966235247 + 0.9665279157420975*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	 }
       else if (jetPosition==8)
	 {
	   xjet=sTau-0.4-1.;
	   yjet=7.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 0.
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 0.
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	 }
       else if (jetPosition==9)
	 {
	   xjet=sTau-0.4-5.5;
	   yjet=0.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.00294276431756867*pow (sTau,4))*(0.40543626547798606*sTau - 6.066521450106897*pow (sTau,0.5) + 0.42696219332946833*pow (sTau,2) + 6.756244273644252*sin (2.3741241238376887 - 0.4248126062148173*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0029266103802329095*pow (sTau,4))*(0.5261692898305094*sTau - 6.529510085959178*pow (sTau,0.5) + 0.43572136244182547*pow (sTau,2) + 7.141239170695929*sin (2.3642116208529296 - 0.4192923683313183*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }

	   if (includeTrigger==1)
	     {
	       xtrigger=-5.5-sTau+0.4;
	       ytrigger=0.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.985272760379311*pow (sTau,4))*(-12.714544798547672*sTau - 11.352078820732894*pow (sTau,0.5)
								 + 21.185236148526357*pow (sTau,2) 
								 + 10.5548749883366*sin (0.45026013602488324 + 2.3167468118844567*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-1.007505093427274*pow (sTau,4))*(15.429504345313838*sTau + 17.384984572674774*pow (sTau,0.5) - 27.229770851187904*pow (sTau,2) - 14.70991393388861*sin (2.670871176613221 - 2.2053827398438592*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }

	 }
       else if (jetPosition==10)
	 {
	   xjet=sTau-0.4-5.25;
	   yjet=1.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.001913558624126035*pow (sTau,4))*(1.5459813961207587*sTau - 4.284252971753834*pow (sTau,0.5) - 0.06759638067063292*pow (sTau,2) + 3.9142035233877785*sin (1.564395039254696 + 0.004525763896476653*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0028304168723991925*pow (sTau,4))*(4.730252160845711*sTau - 9.311695746532726*pow (sTau,0.5) - 0.00928347433856184*pow (sTau,2) + 6.187616312129986*sin (1.5444175977562364 - 0.370469655675186*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0022955058518539873*pow (sTau,4))*(0.062454235564433774*sTau + 0.12278914599327842*pow (sTau,0.5) - 0.016205004054166215*pow (sTau,2) - 0.2324533195815475*sin (0.24278947417294008 + 0.41848781382546957*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   
	   if (includeTrigger==1)
	     {
	       xtrigger=-5.25-sTau+0.4;
	       ytrigger=1.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.8127084353935011*pow (sTau,4))*(-13.450020289755795*sTau - 15.380058425584462*pow (sTau,0.5) + 22.89270127691148*pow (sTau,2) + 13.53888833054219*sin (0.4831687569725234 + 2.081527998284279*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.8099755500100739*pow (sTau,4))*(10.916895894492189*sTau + 13.910082762393523*pow (sTau,0.5) - 20.440129840054848*pow (sTau,2) - 11.904424326141177*sin (2.6349294709655253 - 2.1202791303724475*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==11)
	 {
	   xjet=sTau-0.4-5.;
	   yjet=2.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0020788605369378627*pow (sTau,4))*(1.391140336454321*sTau - 3.980553791418186*pow (sTau,0.5) - 0.03695043702387419*pow (sTau,2) + 3.7396961892944596*sin (1.5124985949970577 + 0.11995380807294627*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.00342532394303531*pow (sTau,4))*(3.6890671064376903*sTau - 8.312825381563028*pow (sTau,0.5) + 0.09426714892616761*pow (sTau,2) + 5.791657213833773*sin (1.7432626624031993 - 0.39510716928255707*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.003776029815679842*pow (sTau,4))*(0.7545652652847962 + 0.014261787683823312*sTau + 0.46911022764044474*pow (sTau,0.5) - 1.9319696697437254*sin (2.6195396538874434 - 0.13814672660505822*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-5.-sTau+0.4;
	       ytrigger=2.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.9650537872598234*pow (sTau,4))*(-9.996226964616918*sTau - 26.520846838100706*pow (sTau,0.5) + 24.748559974849023*pow (sTau,2) + 19.238636942065675*sin (0.5450349494376384 + 1.8337241544971552*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.8033881190653333*pow (sTau,4))*(2.0872660286885396*sTau + 5.071408958413009*pow (sTau,0.5) - 6.2874909807136525*pow (sTau,2) - 4.7168666848261935*sin (2.445102487613561 - 1.9590308275378376*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==12)
	 {
	   xjet=sTau-0.4-4.5;
	   yjet=3.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0027672015224750727*pow (sTau,4))*(1.579868744085268*sTau - 4.355255691593592*pow (sTau,0.5) + 0.001787891873849861*pow (sTau,2) + 3.909184380107434*sin (1.4799447657693332 + 0.2187624105734425*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.004507345590933594*pow (sTau,4))*(4.070680827772322*sTau - 8.880979460528396*pow (sTau,0.5) + 0.1028781070350226*pow (sTau,2) + 6.014300524889755*sin (1.7347449181862455 - 0.4143091544265476*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0038928286349715923*pow (sTau,4))*(2.607872837877081 - 0.46792726277662056*sTau + 0.7992949117864787*pow (sTau,0.5) - 3.0350924071982948*sin (1.35077705227024 - 0.12656815720992065*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-4.5-sTau+0.4;
	       ytrigger=3.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }

	 }
       else if (jetPosition==13)
	 {
	   xjet=sTau-0.4-3.9;
	   yjet=4.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0038928286349715923*pow (sTau,4))*(2.607872837877081 - 0.46792726277662056*sTau + 0.7992949117864787*pow (sTau,0.5) - 3.0350924071982948*sin (1.35077705227024 - 0.12656815720992065*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.007345852179482*pow (sTau,4))*(2.7521306017777416*sTau - 7.850387786019089*pow (sTau,0.5) + 0.328993277901861*pow (sTau,2) + 5.877957222689557*sin (1.991102663926878 - 0.4855694909212403*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.005584314579207876*pow (sTau,4))*(-0.2582522444760471 - 0.21620252567191448*sTau + 0.6001041476633*pow (sTau,0.5) - 0.0015819227819163317*sin (7.535818088682731 - 10.200834737073285*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-3.9-sTau+0.4;
	       ytrigger=4.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }

	 }
       else if (jetPosition==14)
	 {
	   xjet=sTau-0.4-3.;
	   yjet=5.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.005584314579207876*pow (sTau,4))*(-0.2582522444760471 - 0.21620252567191448*sTau + 0.6001041476633*pow (sTau,0.5) - 0.0015819227819163317*sin (7.535818088682731 - 10.200834737073285*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.009462099666702754*pow (sTau,4))*(1.4531441132857945*sTau - 4.0008577042702935*pow (sTau,0.5) - 0.07605156480868927*pow (sTau,2) + 3.6962908472582447*sin (1.613902367489429 - 0.012662120256958903*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.013418778818857923*pow (sTau,4))*(-1.3422574926343736 - 0.06394454375124038*sTau + 1.5101778289182743*pow (sTau,0.5) - 1.1663080069033012*sin (18.30936835494525 + 0.45775491282333725*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-3.-sTau+0.4;
	       ytrigger=5.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }

	 }
       else if (jetPosition==15)
	 {
	   xjet=sTau-0.4-2.;
	   yjet=6.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.04786478740658437*pow (sTau,4))*(2.791349433324794*sTau - 4.443661128955812*pow (sTau,0.5) + 0.2456605815128899*pow (sTau,2) + 3.48410205301595*sin (1.7858586420329274 + 0.6549330329206062*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.04643883508340708*pow (sTau,4))*(3.5583912833558986*sTau - 5.202268415395203*pow (sTau,0.5) + 0.19151939143460414*pow (sTau,2) + 3.854227221251983*sin (1.8237285767976308 + 0.6634774383799322*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.03551099360640726*pow (sTau,4))*(-0.6924002285423925 - 0.7219058428452331*sTau + 1.6919924124128298*pow (sTau,0.5) + 0.016257159228193784*sin (2.7507517283724874 - 8.368088572061957*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-2.-sTau+0.4;
	       ytrigger=6.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }

	 }
       else if (jetPosition==16)
	 {
	   xjet=sTau-0.4-0.;
	   yjet=7.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 0.
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 0.
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-sTau+0.4;
	       ytrigger=7.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }

	 }
       else if (jetPosition==17)
	 {
	   xjet=sTau-0.4-4.5;
	   yjet=0.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.002902751000018112*pow (sTau,4))*(0.08868492437195961*sTau - 5.505729415090479*pow (sTau,0.5) + 0.36906973443894664*pow (sTau,2) + 6.905884249459013*sin (2.3757215745552664 - 0.37253537536081105*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0028929736606367815*pow (sTau,4))*(0.6524226527514703*sTau - 6.393786013776885*pow (sTau,0.5) + 0.340858895489331*pow (sTau,2) + 7.140028288646667*sin (2.305993365464475 - 0.3662848457962157*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-4.5-sTau+0.4;
	       ytrigger=0.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.2240631925743222*pow (sTau,4))*(32.725522589815*sTau - 40.05622815581742*pow (sTau,0.5) + 1.2417172332276807*pow (sTau,2) + 14.397429198583213*sin (1.3099149912664074 + 1.268359751562469*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.05109616429075409*pow (sTau,4))*(5.823939768790299*sTau - 6.4909687061759715*pow (sTau,0.5) - 0.6893155082300996*pow (sTau,2) - 0.02805457827623055*sin (5.894800163675925 - 10.636182004977723*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==18)
	 {
	   xjet=sTau-0.4-4.25;
	   yjet=1.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0029118546832658876*pow (sTau,4))*(2.6855806896808847*sTau - 6.861154426971757*pow (sTau,0.5) + 0.08405958135016414*pow (sTau,2) + 5.387851899037637*sin (1.7457016539918702 - 0.32974388711486646*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.002922733342658492*pow (sTau,4))*(2.919466409749543*sTau - 7.569925591486737*pow (sTau,0.5) + 0.09250952032283395*pow (sTau,2) + 5.875396983479187*sin (1.773213186690503 - 0.32949575533559483*sTau)) 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0010690135804432253*pow (sTau,4))*(0.05976895307991918*sTau + 0.0995046403053554*pow (sTau,0.5) - 1.1596991664037817*sin (0.03887628790564185 + 0.08446366182522316*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-4.25-sTau+0.4;
	       ytrigger=1.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.05273062507375901*pow (sTau,4))*(6.77767238307406*sTau - 9.517047870703166*pow (sTau,0.5) - 0.02096480044799833*pow (sTau,2) + 5.55700644921248*sin (1.551559513789632 + 0.7581316274651319*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.08785735636953232*pow (sTau,4))*(-3.0104886376548046*sTau + 6.27923533901553*pow (sTau,0.5) + 0.06538658243511675*pow (sTau,2) - 4.904216978985655*sin (1.3929572101040348 + 0.42166506504974754*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==19)
	 {
	   xjet=sTau-0.4-4.;
	   yjet=2.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.00400655654968671*pow (sTau,4))*(-4.906063135530068*sTau - 2.2019366668616236*pow (sTau,0.5) + 0.9836841907332459*pow (sTau,2) + 10.385323057700257*sin (2.7880933074297993 - 0.43531518485089726*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.00393910308289692*pow (sTau,4))*(-4.154239356495493*sTau - 3.246405214225931*pow (sTau,0.5) + 0.9185666471852633*pow (sTau,2) + 10.541130810594293*sin (2.736897611128812 - 0.4232156141451517*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0013604316443845536*pow (sTau,4))*(-0.0010329765635707859*sTau + 0.23885772490959412*pow (sTau,0.5) - 0.7003764804837042*sin (0.1583389373462131 + 0.13555642700760426*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-4.-sTau+0.4;
	       ytrigger=2.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.1876616080422508*pow (sTau,4))*(18.61730925134063*sTau - 28.04318403229763*pow (sTau,0.5) + 2.7486503612172415*pow (sTau,2) + 11.874176780440402*sin (1.169826323359127 + 1.205494809764435*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.1883157781121874*pow (sTau,4))*(-17.283203491827674*sTau + 27.823088416467527*pow (sTau,0.5) - 3.2349288746324127*pow (sTau,2) - 12.276260113895768*sin (2.004876070519177 - 1.1940954001579402*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.15252249721386935*pow (sTau,4))*(0.2146957962268962*sTau - 0.10760874773298813*pow (sTau,0.5) - 0.03344893342113466*pow (sTau,2) + 0.013586267407018753*sin (7.208778260565238 - 6.137262734340389*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==20)
	 {
	   xjet=sTau-0.4-3.5;
	   yjet=3.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.004673145785568747*pow (sTau,4))*(3.6784707067485356*sTau - 7.895126250763775*pow (sTau,0.5) + 0.05016777504754546*pow (sTau,2) + 5.705037238159215*sin (1.5940444839073495 - 0.3572970879562185*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.004731996603181068*pow (sTau,4))*(3.9226603879421793*sTau - 8.661365759368922*pow (sTau,0.5) + 0.06825388452021515*pow (sTau,2) + 6.186336410058382*sin (1.6439912124181604 - 0.36122457104965733*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0036602175118534905*pow (sTau,4))*(-0.2139723463717529*sTau + 0.4746277454360757*pow (sTau,0.5) - 0.2090877756544745*sin (1.4612211145440397 + 0.374002023649221*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-3.5-sTau+0.4;
	       ytrigger=3.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.049834071347335175*pow (sTau,4))*(4.828672659834893*sTau - 7.689332316755995*pow (sTau,0.5) + 0.15264554453521245*pow (sTau,2) + 4.971572509297799*sin (1.4667722881528054 + 0.7319268030706501*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.08198009973792314*pow (sTau,4))*(-2.8037619914674536*sTau + 5.580235247480415*pow (sTau,0.5) + 0.38340542988426224*pow (sTau,2) - 4.593076274918189*sin (1.4343069470688408 + 0.0012072290792002525*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==21)
	 {
	   xjet=sTau-0.4-2.9;
	   yjet=4.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.007921799212089945*pow (sTau,4))*(2.7292121164887684*sTau - 7.993500885608372*pow (sTau,0.5) + 0.3017596173450703*pow (sTau,2) + 6.156966237789564*sin (1.9907541344769812 - 0.45356155465520703*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.007879568636504682*pow (sTau,4))*(3.5330385109713913*sTau - 8.923393099138014*pow (sTau,0.5) + 0.2406832385947608*pow (sTau,2) + 6.434372174658485*sin (1.881849532930337 - 0.4422321111038994*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.004550573490057557*pow (sTau,4))*(0.6178959718775526*pow (sTau,0.5) - 0.03015637903208191*pow (sTau,2) - 0.7072375427651219*sin (0.4409304949258095 + 0.32861032114370187*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-2.9-sTau+0.4;
	       ytrigger=4.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.17083858793803347*pow (sTau,4))*(17.219858108328395*sTau - 23.47201380703663*pow (sTau,0.5) + 1.442478497912089*pow (sTau,2) + 9.819275750576576*sin (1.3052251069532794 + 1.1454017216909782*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.08550181070139284*pow (sTau,4))*(-1.4541383939492034*sTau + 3.6306419443333082*pow (sTau,0.5) + 0.20411944994292452*pow (sTau,2) - 4.135567985728974*sin (1.1202455018133615 + 0.00247302133844981*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==22)
	 {
	   xjet=sTau-0.4-2;
	   yjet=5.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.017497100808295727*pow (sTau,4))*(3.0912358279713703*sTau - 8.77280120046681*pow (sTau,0.5) + 0.5232802633822561*pow (sTau,2) + 6.33119992547775*sin (2.038240121809075 - 0.5736999285362875*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.017420420281954944*pow (sTau,4))*(3.5516973860531698*sTau - 9.711100792819272*pow (sTau,0.5) + 0.5433737524585578*pow (sTau,2) + 6.840639990393554*sin (2.0264657912886395 - 0.572502079872621*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.012276453674695894*pow (sTau,4))*(-0.5192594143801705 - 0.395031520035804*sTau + 1.076290982337544*pow (sTau,0.5) + 0.002165689641501208*sin (1.0574592916669154 + 7.298591826508869*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-2.-sTau+0.4;
	       ytrigger=5.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==23)
	 {
	   xjet=sTau-0.4-1;
	   yjet=6.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.06263886155697326*pow (sTau,4))*(4.611994365422411*sTau - 6.271654548366204*pow (sTau,0.5) + 0.1408712039293767*pow (sTau,2) + 4.185626290907246*sin (1.7828512866776576 + 0.7237402285710154*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.062494557431210865*pow (sTau,4))*(4.55837328776617*sTau - 7.1563033550787205*pow (sTau,0.5) + 0.33815138962124197*pow (sTau,2) + 4.505772553177724*sin (1.5691868145762282 + 0.757178441995617*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.05681006509193217*pow (sTau,4))*(-1.5467761379793241 - 1.1093462724287124*sTau + 2.9361588834287384*pow (sTau,0.5) - 0.10678835105181621*sin (3.7923142282106297 + 2.257113949788277*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-1.-sTau+0.4;
	       ytrigger=6.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 0.
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==24)
	 {
	   xjet=sTau-0.4-3.5;
	   yjet=0.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0030334597905961494*pow (sTau,4))*(1.141229864555664*sTau - 5.548694848770199*pow (sTau,0.5) + 0.1518638211655223*pow (sTau,2) + 5.754453728646889*sin (2.1169002270193227 - 0.2874298572279565*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0030362735867955387*pow (sTau,4))*(1.828095763580888*sTau - 6.493079280335862*pow (sTau,0.5) + 0.11692552477321608*pow (sTau,2) + 5.935703393351542*sin (1.9722219997244084 - 0.28165042620595837*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-3.5-sTau+0.4;
	       ytrigger=0.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.027044309795761348*pow (sTau,4))*(6.060981314408092*sTau - 12.86451984574572*pow (sTau,0.5) + 0.401839123343581*pow (sTau,2) + 7.999683711824608*sin (2.0099615494257983 - 0.6472469157633255*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.010529021569601825*pow (sTau,4))*(4.632006918391199*sTau - 5.6915918133690235*pow (sTau,0.5) - 0.4497065146208425*pow (sTau,2) + 0.024770007163109575*sin (0.8207811012801189 + 8.37398160769294*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==25)
	 {
	   xjet=sTau-0.4-3.25;
	   yjet=1.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.003523579361815024*pow (sTau,4))*(2.804798623985102*sTau - 7.305184763038861*pow (sTau,0.5) + 0.09572921227397486*pow (sTau,2) + 5.789675971394264*sin (1.7390164516970634 - 0.3213110080970541*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.003562112377789245*pow (sTau,4))*(2.8505586902981395*sTau - 8.034970489360276*pow (sTau,0.5) + 0.1333067862725772*pow (sTau,2) + 6.41712763334741*sin (1.8433965932314873 - 0.3294210486562722*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.004221357376357216*pow (sTau,4))*(0.9336484071499895 - 0.06549939156428349*sTau + 0.19045019524957732*pow (sTau,0.5) - 1.0541308785691728*sin (1.8021146136346435 - 0.11148885907413109*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-3.25-sTau+0.4;
	       ytrigger=1.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.004221357376357216*pow (sTau,4))*(0.9336484071499895 - 0.06549939156428349*sTau + 0.19045019524957732*pow (sTau,0.5) - 1.0541308785691728*sin (1.8021146136346435 - 0.11148885907413109*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.008522281602190951*pow (sTau,4))*(4.383991451916347*sTau - 5.481880483932806*pow (sTau,0.5) - 0.41159527695159664*pow (sTau,2) - 0.022340611229409453*sin (4.109702681443229 + 7.91026410275091*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.04695119403578456*pow (sTau,4))*(0.875489280789133*sTau - 0.4197794730835628*pow (sTau,0.5) - 0.19202168849756018*pow (sTau,2) + 0.35041718892176466*sin (0.23229472728333989 - 0.9737471876463205*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==26)
	 {
	   xjet=sTau-0.4-3.;
	   yjet=2.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0046373820942023835*pow (sTau,4))*(3.763707513942669*sTau - 8.210948628318384*pow (sTau,0.5) + 0.07462526072638705*pow (sTau,2) + 6.0425824487611255*sin (1.5712439970535894 - 0.35180756139140573*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.004725328524525551*pow (sTau,4))*(4.56777889272573*sTau - 9.351634557577139*pow (sTau,0.5) + 0.04666848047077209*pow (sTau,2) + 6.677795690450636*sin (1.4974798638567597 - 0.3548088668253056*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.006588373127720927*pow (sTau,4))*(0.8822147789428378 - 0.11650377349372705*sTau + 0.31689527258101363*pow (sTau,0.5) - 1.0578287419854382*sin (1.852648366975483 - 0.17037776897577556*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-3.5-sTau+0.4;
	       ytrigger=0.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.0645952324959411*pow (sTau,4))*(-1.3522816798080468*sTau - 9.354515750379464*pow (sTau,0.5) + 2.670054332600475*pow (sTau,2) + 9.499831189312836*sin (2.4135097532702976 - 0.8069549729157073*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.06654503821107935*pow (sTau,4))*(0.5731607544892123*sTau + 10.970669623768421*pow (sTau,0.5) - 2.815125058109428*pow (sTau,2) + 10.225031712437389*sin (5.518587751739419 - 0.8099952024996363*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.04240835146565604*pow (sTau,4))*(0.12873857878247968*sTau + 0.07768957715163734*pow (sTau,0.5) - 0.05343497498409488*pow (sTau,2) - 0.10151852387395735*sin (2.2463440085206634 - 0.7310307006616044*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==27)
	 {
	   xjet=sTau-0.4-2.5;
	   yjet=3.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.005705053326956571*pow (sTau,4))*(0.6457710769995819*sTau - 6.13157944049133*pow (sTau,0.5) + 0.34850008487549017*pow (sTau,2) + 6.825430815056587*sin (2.2871218348164777 - 0.37297584561949737*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.005647998709257811*pow (sTau,4))*(1.9144298956047414*sTau - 7.342992542984457*pow (sTau,0.5) + 0.23316477941343*pow (sTau,2) + 6.5978438952235985*sin (2.075880111274593 - 0.3529584053896968*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.008057416135930551*pow (sTau,4))*(0.6262854277970386 + 0.05542756512548674*sTau + 0.4366798991732166*pow (sTau,0.5) - 1.6960535163423422*sin (2.6340042519326294 - 0.17549647167879534*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-2.5-sTau+0.4;
	       ytrigger=3.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.06286866807654191*pow (sTau,4))*(0.22646965205389902*sTau - 10.990357463983583*pow (sTau,0.5) + 2.564218357905855*pow (sTau,2) + 9.639204163050122*sin (2.3536864158270205 - 0.8183588042968982*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.008412857974835639*pow (sTau,4))*(4.320109255126841*sTau - 5.3810049414883485*pow (sTau,0.5) - 0.41033926459486547*pow (sTau,2) + 0.021354652798506873*sin (0.9125217219845126 + 8.220790586478646*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.04494644689037953*pow (sTau,4))*(0.8071024878731886*sTau - 0.5590572294527016*pow (sTau,0.5) - 0.1548819976090776*pow (sTau,2) + 0.07137433494760724*sin (1.5189176778906328 - 1.6911089889084943*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==28)
	 {
	   xjet=sTau-0.4-1.9;
	   yjet=4.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.009690249115348491*pow (sTau,4))*(2.9559591860255456*sTau - 7.454245142447624*pow (sTau,0.5) + 0.19288846058908318*pow (sTau,2) + 5.655635132634967*sin (1.8056759651839085 - 0.4091867926969617*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.009922472605293591*pow (sTau,4))*(3.348846083808456*sTau - 8.295835363222213*pow (sTau,0.5) + 0.21421092633125322*pow (sTau,2) + 6.145263494838894*sin (1.7982863539074145 - 0.4164095899128632*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.011720025995420436*pow (sTau,4))*(1.531328055013795 + 0.16344108671437438*sTau + 0.5758362307945345*pow (sTau,0.5) - 3.6418401572804107*sin (2.625880181760986 - 0.12957311587097853*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-1.9-sTau+0.4;
	       ytrigger=4.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.05479969925453919*pow (sTau,4))*(-1.3376466926968762*sTau - 5.98268138885275*pow (sTau,0.5) + 1.4004113617834026*pow (sTau,2) + 7.457174737887126*sin (2.446276883692137 - 0.6619429848468878*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.07906830693546153*pow (sTau,4))*(22.1113650506217*sTau - 17.635242069596917*pow (sTau,0.5) - 5.026578782728813*pow (sTau,2) - 1.6829694371474142*sin (4.578779128801658 - 1.9877966657131847*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.03920089658879859*pow (sTau,4))*(0.08332517348697588*sTau + 0.3562972309083527*pow (sTau,0.5) - 0.06196716327385682*pow (sTau,2) - 0.37669270825330614*sin (0.7444133671095676 - 3.6830778569240643e-6*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==29)
	 {
	   xjet=sTau-0.4-1.;
	   yjet=5.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.02542393696791293*pow (sTau,4))*(5.1689067258781325*sTau - 11.523560814205702*pow (sTau,0.5) + 0.5653808916130942*pow (sTau,2) + 7.238736539125021*sin (1.9390249039248098 - 0.6251584891157057*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.025291387067865857*pow (sTau,4))*(7.185944457653998*sTau - 13.338037443306224*pow (sTau,0.5) + 0.35076168370066946*pow (sTau,2) + 7.688092454532429*sin (1.7726519509510852 - 0.6110321379505738*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.01755911948733197*pow (sTau,4))*(-0.5376358025996556 - 0.3959882084397944*sTau + 1.0888983865739965*pow (sTau,0.5) - 0.0025873284505280878*sin (9.003917638715105 - 5.3138873408890905*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-1.-sTau+0.4;
	       ytrigger=5.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.05764733091059714*pow (sTau,4))*(3.982436163920859*sTau - 9.818620979568749*pow (sTau,0.5) + 0.7266980893347043*pow (sTau,2) + 6.551059583955836*sin (2.02670032147639 - 0.6890711606366067*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.011869018252653414*pow (sTau,4))*(4.25045044784334*sTau - 5.113967015034581*pow (sTau,0.5) - 0.43463476634627957*pow (sTau,2) - 0.008324718061378591*sin (19.231848841548288 - 27.686594847257144*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.06583533755357535*pow (sTau,4))*(1.67967774475678*sTau - 1.2160883736539134*pow (sTau,0.5) - 0.29790321907218653*pow (sTau,2) - 0.14411106209090518*sin (4.121656918491851 + 1.9972909997202604*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==30)
	 {
	   xjet=sTau-0.4;
	   yjet=6.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.10936406253059605*pow (sTau,4))*(9.098502608563937*sTau - 10.091693909897613*pow (sTau,0.5) - 0.18896692136362572*pow (sTau,2) + 5.326423998586413*sin (1.8188990006226227 + 0.8957773951061764*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.10797725738696093*pow (sTau,4))*(10.779174157654008*sTau - 11.70866107205591*pow (sTau,0.5) - 0.3341728400787949*pow (sTau,2) + 5.972053688567354*sin (1.8409179650483574 + 0.8976894675336098*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.11599216654447202*pow (sTau,4))*(-0.6950062674530637 - 0.3960373931272013*sTau + 1.3554648691541817*pow (sTau,0.5) - 0.06511998054802402*sin (1.037357119488681 + 4.342844934971587*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-sTau+0.4;
	       ytrigger=6.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.11451622266341746*pow (sTau,4))*(7.938248274201743*sTau - 11.981814123341952*pow (sTau,0.5) + 0.6784772081696436*pow (sTau,2) + 6.077388237390741*sin (1.391729504362421 + 0.9417568052649801*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.08444829221235442*pow (sTau,4))*(-1.0332519018784816*sTau + 7.049007051128628*pow (sTau,0.5) + 0.3042552525478686*pow (sTau,2) - 25.15892811262912*sin (25.331560223667758 + 0.0973818306635045*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.12301290561720846*pow (sTau,4))*(1.2278131464553501*sTau - 0.7651689306510236*pow (sTau,0.5) - 0.2151029728374755*pow (sTau,2) - 0.10881847242795059*sin (7.737476764056306 - 3.822119729438926*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==31)
	 {
	   xjet=-2.5+sTau-0.4;
	   yjet=0.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0014372262495466733*pow (sTau,4))*(3.0265648099984834*sTau - 7.815972522500506*pow (sTau,0.5) + 0.016383263731009248*pow (sTau,2) + 6.135577382190255*sin (1.7773203279759904 - 0.28487530952335266*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.004752850698514316*pow (sTau,4))*(3.8743527595349936*sTau - 8.98198370238839*pow (sTau,0.5) + 0.12307750522980652*pow (sTau,2) + 6.654606105144967*sin (1.6332791892251581 - 0.35484456573771334*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-2.5-sTau+0.4;
	       ytrigger=0.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.038823404968060395*pow (sTau,4))*(1.7512195635939352*sTau - 10.810336773138099*pow (sTau,0.5) + 1.7290220208164637*pow (sTau,2) + 8.792476908966865*sin (2.2117924558807247 - 0.7070286156276723*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.04062701226721497*pow (sTau,4))*(-0.5649363491274665*sTau + 11.6136961793932*pow (sTau,0.5) - 2.291255701102521*pow (sTau,2) - 10.239921801668467*sin (2.3039808200945746 - 0.7338799586689239*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==32)
	 {
	   xjet=-2.25+sTau-0.4;
	   yjet=1.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0014753032787306376*pow (sTau,4))*(2.4583348516154233*sTau - 7.593801937030575*pow (sTau,0.5) + 0.07295796791942533*pow (sTau,2) + 6.453195669632415*sin (1.9728110621696333 - 0.2969800521589791*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0015213939609284035*pow (sTau,4))*(3.198383973297099*sTau - 8.582030355904596*pow (sTau,0.5) + 0.032438578348265246*pow (sTau,2) + 6.747944492053158*sin (1.835869844336653 - 0.28959040822174614*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0038656962132657996*pow (sTau,4))*(0.8680226730793907 - 0.3193977063637287*sTau + 0.10039997923413957*pow (sTau,0.5) - 1.835184731415154*sin (0.5139900270873518 - 0.15959047639865695*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-2.25-sTau+0.4;
	       ytrigger=1.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.03679602231176933*pow (sTau,4))*(5.183905362655249*sTau - 14.56796670465583*pow (sTau,0.5) + 1.490289141750878*pow (sTau,2) + 9.570112259071575*sin (2.0544935619520714 - 0.7101746322613057*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.038244543127490475*pow (sTau,4))*(-5.265624451030655*sTau + 16.34203582231867*pow (sTau,0.5) - 1.8449348337498899*pow (sTau,2) + 10.897227712898507*sin (4.17442066044022 + 0.7271019571426925*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.022385385373970598*pow (sTau,4))*(0.18870405053300354*sTau - 0.08118636692434567*pow (sTau,0.5) - 0.041780859095448326*pow (sTau,2) + 0.06404497028759396*sin (0.0981521591402687 - 0.9632203685187715*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==33)
	 {
	   xjet=-2.+sTau-0.4;
	   yjet=2.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.004599979348719903*pow (sTau,4))*(2.106427237432534*sTau - 6.2976574174958415*pow (sTau,0.5) - 0.07490251930450051*pow (sTau,2) + 6.0969745761442615*sin (2.0773668234249874 - 0.10593129248053645*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.005383459331929371*pow (sTau,4))*(2.271872927017533*sTau - 7.2024356401501946*pow (sTau,0.5) + 0.1239806710592679*pow (sTau,2) + 6.1969849993640285*sin (1.8929739593139379 - 0.28702009773727755*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.007069210792196933*pow (sTau,4))*(0.8323934741594162 - 0.10418221680232768*sTau + 0.17844848417599818*pow (sTau,0.5) - 0.9276903354714578*sin (1.748480280327207 + 0.13394254462809974*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-2.-sTau+0.4;
	       ytrigger=2.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.018395357319584537*pow (sTau,4))*(3.7610109584248725*sTau - 9.976238247121444*pow (sTau,0.5) + 0.05023906133790709*pow (sTau,2) + 7.635499363894264*sin (2.086696979598251 - 0.3790115140825889*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.035115312332170634*pow (sTau,4))*(-4.1436657635017555*sTau + 12.393809946407249*pow (sTau,0.5) - 1.2431518214453643*pow (sTau,2) - 8.64257101962757*sin (2.034353960536658 - 0.6636804798530802*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.02442009909060968*pow (sTau,4))*(-0.35022987087474355*sTau + 0.4734996481441619*pow (sTau,0.5) + 0.01579120732371071*pow (sTau,2) - 0.18582441639179606*sin (1.3838852011694214 - 0.7975773008120983*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==34)
	 {
	   xjet=-1.5+sTau-0.4;
	   yjet=3.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.004919650312757161*pow (sTau,4))*(2.3222002954740444*sTau - 6.731367138826253*pow (sTau,0.5) - 0.05161636283142127*pow (sTau,2) + 5.92509162167337*sin (1.982944538589701 - 0.20140956086936307*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.004614763732048756*pow (sTau,4))*(2.760245053426462*sTau - 7.7358987697115*pow (sTau,0.5) - 0.011615871673029849*pow (sTau,2) + 6.387428880318173*sin (1.9201689607719694 - 0.25161303133769736*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.010650506608673596*pow (sTau,4))*(0.612262351331548 + 0.03897299647736938*sTau + 0.3796926600548727*pow (sTau,0.5) - 1.3838164736126017*sin (0.5977560185681787 + 0.19545812082528913*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-1.5-sTau+0.4;
	       ytrigger=3.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.033120044725285756*pow (sTau,4))*(5.681273611909777*sTau - 13.090486306382346*pow (sTau,0.5) + 0.8577772278144475*pow (sTau,2) + 8.281790561082987*sin (1.9225908135807444 - 0.6389267604568034*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.035252353532555465*pow (sTau,4))*(-2.6034020264554067*sTau + 13.056239802358686*pow (sTau,0.5) - 1.7633794192838086*pow (sTau,2) + 10.13081558824343*sin (4.065266141882272 + 0.6908354931020331*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.024589988421586956*pow (sTau,4))*(-0.45276166424473424*sTau + 0.7035268183497906*pow (sTau,0.5) + 0.006700188017706982*pow (sTau,2) - 0.28069915761905523*sin (1.5712513805234207 - 0.7850398369269365*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==35)
	 {
	   xjet=-0.9+sTau-0.4;
	   yjet=4.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.015016858979594813*pow (sTau,4))*(5.492018484701715*sTau - 10.936185660284307*pow (sTau,0.5) + 0.22886827829669834*pow (sTau,2) + 6.978975442091971*sin (1.698584231424903 - 0.4965627609701093*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.015320149753947758*pow (sTau,4))*(6.441280463975216*sTau - 12.128789183892398*pow (sTau,0.5) + 0.19029161718374765*pow (sTau,2) + 7.5472495934027215*sin (1.6265123044717564 - 0.4951787818805935*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.018035767758243088*pow (sTau,4))*(1.8593713843154713 + 0.13338971555013868*sTau + 0.7374907435567433*pow (sTau,0.5) - 3.542273545908911*sin (2.4789162167386185 - 0.1818787628974768*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-0.9-sTau+0.4;
	       ytrigger=4.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.03191145478377903*pow (sTau,4))*(3.316424638290932*sTau - 9.309755950403918*pow (sTau,0.5) + 0.6644091469867403*pow (sTau,2) + 6.7740870991532685*sin (1.9899262451384112 - 0.5788505404035744*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.03335371954250087*pow (sTau,4))*(-1.7964746174130608*sTau + 9.713315142481486*pow (sTau,0.5) - 1.156197784796767*pow (sTau,2) + 8.159772415353649*sin (5.344646787072439 - 0.619756053548155*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.033498758369693585*pow (sTau,4))*(-0.11764467049417265*sTau + 0.4200058464392235*pow (sTau,0.5) + 0.014841573022471284*pow (sTau,2) - 1.0634403424311756*sin (2.9325030349507113 - 0.003313293016223239*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==36)
	 {
	   xjet=sTau-0.4;
	   yjet=5.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.03641658664501125*pow (sTau,4))*(3.41818475905804*sTau - 8.876351834896278*pow (sTau,0.5) + 0.6670809620165384*pow (sTau,2) + 6.207328471780915*sin (1.9881081179530646 - 0.642510701976173*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.0371836040829966*pow (sTau,4))*(4.183045209775258*sTau - 10.285643800052082*pow (sTau,0.5) + 0.7232635316800676*pow (sTau,2) + 6.890442983213694*sin (1.9730879589574135 - 0.649326110980083*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.03258901630865071*pow (sTau,4))*(-0.519171005578741 - 0.34298223401082506*sTau + 1.02619821505636*pow (sTau,0.5) + 0.0023151172816239874*sin (1.1896944880247515 + 6.876127258786177*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-sTau+0.4;
	       ytrigger=5.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.03835761083236513*pow (sTau,4))*(5.926585204572508*sTau - 12.675412835255464*pow (sTau,0.5) + 0.6827845741787248*pow (sTau,2) + 7.658106651535916*sin (1.9460672425163175 - 0.6714955732285345*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.03937817986964633*pow (sTau,4))*(-5.060948750864896*sTau + 13.396929610011545*pow (sTau,0.5) - 1.0665015919954162*pow (sTau,2) + 8.710492608443767*sin (10.484744707167826 + 0.6904255969619694*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.034540566261759796*pow (sTau,4))*(1.0304728924194353*sTau - 0.7345505589119642*pow (sTau,0.5) - 0.1721812527784118*pow (sTau,2) + 0.08538967393121112*sin (0.8832670096749428 + 1.794129069252943*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==37)
	 {
	   xjet=-1.5+sTau-0.4;
	   yjet=0.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.006036580083576199*pow (sTau,4))*(1.9820365144411103*sTau - 6.736024588537131*pow (sTau,0.5) + 0.23595000117534565*pow (sTau,2) + 5.902406166466287*sin (1.886642304938082 - 0.3446947621724734*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.006556571465291123*pow (sTau,4))*(3.0772918948613377*sTau - 8.495374378473672*pow (sTau,0.5) + 0.2742532356526463*pow (sTau,2) + 6.653524055842914*sin (1.7953319898409197 - 0.3899551750121932*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-1.5-sTau+0.4;
	       ytrigger=0.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.017653276142393027*pow (sTau,4))*(2.5540027686973112*sTau - 6.804678840134507*pow (sTau,0.5) - 0.05850686422196289*pow (sTau,2) + 6.598948880862519*sin (2.126334836040546 - 0.03634140467721772*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.01885081880519051*pow (sTau,4))*(-2.369721952178173*sTau + 7.959244026048064*pow (sTau,0.5) - 0.3658021646080463*pow (sTau,2) - 6.682821563458861*sin (1.9290712616537062 - 0.37557576874362003*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==38)
	 {
	   xjet=-1.25+sTau-0.4;
	   yjet=1.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.006188981458820095*pow (sTau,4))*(1.7019033390457576*sTau - 5.581024371209775*pow (sTau,0.5) - 0.025059099422016697*pow (sTau,2) + 6.3710291278851185*sin (2.2118588893014266 - 0.05109220662751584*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.006110207364120221*pow (sTau,4))*(2.075119173619182*sTau - 6.4127748091563594*pow (sTau,0.5) - 0.04802624846553211*pow (sTau,2) + 6.662375492990042*sin (2.137733348632212 - 0.057304957977116655*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.009619758054052692*pow (sTau,4))*(0.8077872030401471 + 0.027906616663023142*sTau + 0.10280912938085073*pow (sTau,0.5) - 1.082032757831442*sin (0.9055764952547456 + 0.1256651967298902*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-1.25-sTau+0.4;
	       ytrigger=1.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.022334780651378852*pow (sTau,4))*(7.069759764469539*sTau - 15.0954833476681*pow (sTau,0.5) + 0.8313340979644008*pow (sTau,2) + 9.137684828614452*sin (1.8650574049395605 - 0.6387755634800332*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.022726857540027616*pow (sTau,4))*(-9.637113234380319*sTau + 17.444819721537176*pow (sTau,0.5) - 0.5733543413397832*pow (sTau,2) + 9.847948961145866*sin (4.8484756580576835 - 0.6277460452838047*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.019937635801130827*pow (sTau,4))*(0.9474050663241974 + 0.6139936012864207*sTau - 0.15742585714062612*pow (sTau,2) - 1.5210582789703448*sin (13.25354046447905 + 0.468645393152694*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==39)
	 {
	   xjet=-1.+sTau-0.4;
	   yjet=2.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.007111619491665004*pow (sTau,4))*(2.393373182349039*sTau - 6.598759847013704*pow (sTau,0.5) - 0.08911281313585782*pow (sTau,2) + 6.271059081740174*sin (2.064222353549638 - 0.05122759335699637*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.008758225232150727*pow (sTau,4))*(2.444713691394675*sTau - 7.9932345777685825*pow (sTau,0.5) + 0.34758067415547256*pow (sTau,2) + 6.626531696844264*sin (1.9046438201295472 - 0.3950676115665397*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.009526463037671531*pow (sTau,4))*(0.7934830760491649 - 0.04355949878137271*sTau + 0.30014631538177833*pow (sTau,0.5) - 1.1288383498649686*sin (2.169687743205054 - 0.16760248559899987*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-1.-sTau+0.4;
	       ytrigger=2.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.021379161911016582*pow (sTau,4))*(6.353411921631624*sTau - 13.376421651252473*pow (sTau,0.5) + 0.6380038054036753*pow (sTau,2) + 8.293808399018655*sin (1.8123837033089392 - 0.5968753448674866*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.02182704733263932*pow (sTau,4))*(-8.081729855907522*sTau + 15.328151974450419*pow (sTau,0.5) - 0.5305957908917311*pow (sTau,2) + 9.041644304467763*sin (4.857163122782111 - 0.5937395079304582*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.019014403257030684*pow (sTau,4))*(3.454780117656778 + 1.377823867992245*sTau - 0.3395807481748257*pow (sTau,2) + 4.796303473635836*sin (3.957976513933571 + 0.3816948452272641*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==40)
	 {
	   xjet=-0.5+sTau-0.4;
	   yjet=3.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.01074477698634356*pow (sTau,4))*(1.983857634687937*sTau - 6.126458786564379*pow (sTau,0.5) - 0.045494359772938885*pow (sTau,2) + 6.195624072084153*sin (2.1450743490031163 - 0.11310287338466851*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.010696095046045535*pow (sTau,4))*(2.261449897237186*sTau - 6.884097996305059*pow (sTau,0.5) - 0.05172861125045582*pow (sTau,2) + 6.648480809013435*sin (2.1181461834233377 - 0.12858715208463384*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.015409314439574625*pow (sTau,4))*(0.6133435657215833 + 0.05233158462515747*sTau + 0.44544257253332353*pow (sTau,0.5) - 1.4668447345942563*sin (2.5566453985225275 - 0.22952211627526*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-0.5-sTau+0.4;
	       ytrigger=3.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.0068445048761359645*pow (sTau,4))*(4.988293676322439*sTau - 11.47390914257041*pow (sTau,0.5) + 0.11401485448322002*pow (sTau,2) + 7.7616326791162225*sin (1.910115856775044 - 0.43587390379321217*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.020813656206150933*pow (sTau,4))*(-6.701016227478931*sTau + 13.674575255251407*pow (sTau,0.5) - 0.5101191758716368*pow (sTau,2) + 8.45615144571175*sin (4.90259686336129 - 0.5629468192737487*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.017626547008869304*pow (sTau,4))*(2.9215059847781872 + 1.3908126630375552*sTau - 0.36388833500888434*pow (sTau,2) - 4.162211563855551*sin (0.795464279618991 + 0.41922435344299164*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==41)
	 {
	   xjet=-0.5+sTau-0.4;
	   yjet=0.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.006992960998210004*pow (sTau,4))*(2.6444370563625723*sTau - 6.985309412800202*pow (sTau,0.5) - 0.10283752933685916*pow (sTau,2) + 6.442045951329669*sin (2.0327080664472987 - 0.02628289004240723*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.009903510847640313*pow (sTau,4))*(3.096880923810613*sTau - 9.735187556825972*pow (sTau,0.5) + 0.5705380742236034*pow (sTau,2) + 7.548054581805482*sin (1.9446429218930512 - 0.4849519973480446*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-0.5-sTau+0.4;
	       ytrigger=0.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.012200093260722404*pow (sTau,4))*(2.3596010535763714*sTau - 7.968027303653305*pow (sTau,0.5) + 0.37555998857692646*pow (sTau,2) + 6.639804265797814*sin (1.9339842291875744 - 0.39434896159690985*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.014282653681934066*pow (sTau,4))*(-3.8726956167474524*sTau + 10.752519194043137*pow (sTau,0.5) - 0.6346248579180311*pow (sTau,2) + 7.836974401807125*sin (4.395395160548572 + 0.5118053896153262*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==42)
	 {
	   xjet=-0.25+sTau-0.4;
	   yjet=1.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.009062508914462615*pow (sTau,4))*(2.468105173586902*sTau - 7.325038214001014*pow (sTau,0.5) + 0.20868988637033975*pow (sTau,2) + 6.071283684896168*sin (1.812150859778465 - 0.3332991978021835*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.011373287592985611*pow (sTau,4))*(3.595006102914763*sTau - 10.232855728852456*pow (sTau,0.5) + 0.5797515849458237*pow (sTau,2) + 7.606208686178364*sin (1.8997531586641176 - 0.50139973170975*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.010670953026393863*pow (sTau,4))*(0.6979730624900741 + 0.03222738287279124*sTau + 0.12430513001524042*pow (sTau,0.5) - 1.0599862631421735*sin (0.7989985028739367 + 0.12398117497412166*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-0.25-sTau+0.4;
	       ytrigger=1.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.010328116984713472*pow (sTau,4))*(2.573915048497882*sTau - 6.868355139048105*pow (sTau,0.5) - 0.07607068349813273*pow (sTau,2) + 6.452777290207822*sin (2.053748195938428 - 0.020374006562163573*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.010264773815524986*pow (sTau,4))*(-2.8849451095316705*sTau + 7.776535575927119*pow (sTau,0.5) + 0.10333337112046262*pow (sTau,2) + 7.7185350835441096*sin (4.087567922091075 + 0.04063633627756037*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.009450801615453898*pow (sTau,4))*(4.075966505053479 + 0.7088373432702433*sTau - 0.1807975206590012*pow (sTau,2) - 4.752706979311239*sin (13.600556308484684 + 0.2771978575371827*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
       else if (jetPosition==43)
	 {
	   xjet=sTau-0.4;
	   yjet=2.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.010048562155659663*pow (sTau,4))*(2.8058128252753654*sTau - 7.158710002056507*pow (sTau,0.5) - 0.11794755268638339*pow (sTau,2) + 6.403784379721178*sin (2.0221229760212256 - 0.033835243041232164*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.013667693188424822*pow (sTau,4))*(2.8171906167308354*sTau - 9.29783622598996*pow (sTau,0.5) + 0.658835491695092*pow (sTau,2) + 7.302517604657306*sin (1.9592246441234424 - 0.5103485617600766*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += exp(-0.014587818306744963*pow (sTau,4))*(0.7231343782702441 - 0.04796348028514608*sTau + 0.3774944667453521*pow (sTau,0.5) - 1.1475578797313286*sin (2.2686825398287396 - 0.2349461671076996*sTau))
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(0.3,2.)+pow((y-yjet),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		  *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=-sTau+0.4;
	       ytrigger=2.;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.010458791141194664*pow (sTau,4))*(2.5129369689643304*sTau - 6.66301941540482*pow (sTau,0.5) - 0.08222258174894903*pow (sTau,2) + 6.554656084444954*sin (2.1310318780384248 - 0.019778310693981067*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.010323737751631229*pow (sTau,4))*(-2.861056668178366*sTau + 7.561899084791621*pow (sTau,0.5) + 0.10965013106917074*pow (sTau,2) + 7.088029032163617*sin (4.176119371302653 + 0.04149238022513548*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += exp(-0.014970918869615761*pow (sTau,4))*(-1.581116500092839 + 0.14773636522567435*sTau + 0.04860894463762158*pow (sTau,2) - 1.6339826841775276*sin (4.977908178884162 + 0.26716340454923154*sTau))
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(0.3,2.)+pow((y-ytrigger),2.)/2/pow(0.3,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(0.3*0.3*asinh(0.2/tau));
		 }
	     }
	 }
     }
   rhs[alpha] = sumf*(DATA->delta_tau);
   
  }/* alpha */
   
 return;
}/* MakeDeltaQI */



/* %%%%%%%%%%%%%%%%% Kurganov-Tadmor %%%%%%%%%%%%%%%%% */


void Evolve::GetQIs
(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2,
 double *qi, NbrQs *NbrCells, int rk_flag, InitData *DATA, int size, int rank)
{
 int alpha, i;
 double tempg, T00, K00, tempf;
 int nmax[4];
 // auxiliary grids to store the received grids from the neighbor processor
 int from;
 int to;

 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta-1; // reduced all eta grids by one for parallel

 tempg = tau;

 for(alpha=0; alpha<5; alpha++)
  {
    /* qs from the neighbors */
    /* implement outflow boundary condition - simply set the two outside
     * values the same as at the boundary */
        
    for(i=1;i<=2;i++)
      {
	if(grid_pt->position[i] == nmax[i])/* plus edge */
	  {
	    NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] 
	      = grid_pt->nbr_m_1[i]->nbr_m_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qim2[alpha][i] *= tempg;  
	 
	    NbrCells->qip1[alpha][i] = qi[alpha];
	    
	    NbrCells->qip2[alpha][i] = qi[alpha];
	  }
	else if(grid_pt->position[i] == (nmax[i]-1))/* plus edge - 1 */
	  {
	    NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] 
	      = grid_pt->nbr_m_1[i]->nbr_m_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qim2[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] = qi[alpha];
	  }
	else if(grid_pt->position[i] == 0)/* minus edge */
	  {
	    NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] 
	      = grid_pt->nbr_p_1[i]->nbr_p_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qip2[alpha][i] *= tempg;  
	    
	    NbrCells->qim1[alpha][i] = qi[alpha];
	    
	    NbrCells->qim2[alpha][i] = qi[alpha];
	  }
	else if(grid_pt->position[i] == 1)/* minus edge + 1 */
	  {
	    NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] 
	      = grid_pt->nbr_p_1[i]->nbr_p_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qip2[alpha][i] *= tempg;  
	    
	    NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] = qi[alpha];
	  }
	else
	  {
	    NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] 
	      = grid_pt->nbr_p_1[i]->nbr_p_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qip2[alpha][i] *= tempg;  
	    
	    NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] 
	      = grid_pt->nbr_m_1[i]->nbr_m_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qim2[alpha][i] *= tempg;  
	  }
      }/* i */    
    
    i=3; //in eta direction where the lattice is cut
    if(grid_pt->position[i] == nmax[i])/* plus edge */
      {
	if(rank == size-1) // for the right most rank do boundary condition on the right
	  {
	    NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] = grid_pt->nbr_m_1[i]->nbr_m_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qim2[alpha][i] *= tempg;  
	    
	    NbrCells->qip1[alpha][i] = qi[alpha];
	    
	    NbrCells->qip2[alpha][i] = qi[alpha];
	    }
	else
	  {
	    NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] = grid_pt->nbr_m_1[i]->nbr_m_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qim2[alpha][i] *= tempg;  
	    
	    NbrCells->qip1[alpha][i] = Rneighbor->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] = Rneighbor2->TJb[rk_flag][alpha][0]; 
	    NbrCells->qip2[alpha][i] *= tempg;  
	  }
      }
    else if(grid_pt->position[i] == (nmax[i]-1))/* plus edge - 1 */
      {
	if(rank == size-1) // for the right most rank do boundary condition on the right
	  {
	    NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] = grid_pt->nbr_m_1[i]->nbr_m_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qim2[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] = qi[alpha];
	  }
	else
	  {
	    NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] = grid_pt->nbr_m_1[i]->nbr_m_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qim2[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] = Rneighbor->TJb[rk_flag][alpha][0]; 
	    NbrCells->qip2[alpha][i] *= tempg;  
	  }
      }
    
    else if(grid_pt->position[i] == 0)/* minus edge */
      {
	if(rank == 0) // for the left most rank do boundary condition on the left
	  {
	    NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] = grid_pt->nbr_p_1[i]->nbr_p_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qip2[alpha][i] *= tempg;  
	    
	    NbrCells->qim1[alpha][i] = qi[alpha];
	      
	    NbrCells->qim2[alpha][i] = qi[alpha];
	  }
	else
	  {
	    NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] = grid_pt->nbr_p_1[i]->nbr_p_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qip2[alpha][i] *= tempg;  
	    
	    NbrCells->qim1[alpha][i] = Lneighbor->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] = Lneighbor2->TJb[rk_flag][alpha][0]; 
	    NbrCells->qim2[alpha][i] *= tempg;  
	 //    if (grid_pt->position[1]==20 && grid_pt->position[2] == 20)
// 	      {
// 		cout << "x=" << grid_pt->position[1] << " y=" << grid_pt->position[2] 
// 		     << " Lneighbor->TJb[" << rk_flag << "][" << alpha << "][0]=" << Lneighbor->TJb[rk_flag][alpha][0] << endl;
// 		cout << "x=" << grid_pt->position[1] << " y=" << grid_pt->position[2] 
// 		     << " Lneighbor2->TJb[" << rk_flag << "][" << alpha << "][0]=" << Lneighbor2->TJb[rk_flag][alpha][0] << endl;
// 	      }
	  }
      }
    else if(grid_pt->position[i] == 1)/* minus edge + 1 */
      {
	if(rank == 0)
	  {
	    NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] = grid_pt->nbr_p_1[i]->nbr_p_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qip2[alpha][i] *= tempg;  
	    
	    NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] = qi[alpha];
	  }
	else
	  {
	    NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip1[alpha][i] *= tempg;  
	    
	    NbrCells->qip2[alpha][i] = grid_pt->nbr_p_1[i]->nbr_p_1[i]->TJb[rk_flag][alpha][0]; 
	    NbrCells->qip2[alpha][i] *= tempg;  
	    
	    NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim1[alpha][i] *= tempg;
	    
	    NbrCells->qim2[alpha][i] = Lneighbor->TJb[rk_flag][alpha][0]; 
	    NbrCells->qim2[alpha][i] *= tempg;  
	  }
      }
    else // all normal (not edge) cells
      {
	NbrCells->qip1[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	NbrCells->qip1[alpha][i] *= tempg;  
	
	NbrCells->qip2[alpha][i] = grid_pt->nbr_p_1[i]->nbr_p_1[i]->TJb[rk_flag][alpha][0]; 
	NbrCells->qip2[alpha][i] *= tempg;  
	
	NbrCells->qim1[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	NbrCells->qim1[alpha][i] *= tempg;
	
	NbrCells->qim2[alpha][i] = grid_pt->nbr_m_1[i]->nbr_m_1[i]->TJb[rk_flag][alpha][0]; 
	NbrCells->qim2[alpha][i] *= tempg;  
      }
  }/* alpha */
 
 return; 
}/* GetQIs */     


int Evolve::MakeQIHalfs
(double *qi, NbrQs *NbrCells, BdryCells *HalfwayCells, 
 Grid *grid_pt, InitData *DATA)
{
 int alpha, direc, nmax[4], flag;
 double fphL, fphR, fmhL, fmhR;
 double gphL, gphR, gmhL, gmhR;
 double x, y, eta, tempf;

  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;
  
  for(alpha=0; alpha<5; alpha++)
  {
   for(direc=1; direc<=3; direc++)
    {
      gphL = qi[alpha];
      fphL = 0.5*minmod_dx(NbrCells->qip1[alpha][direc], qi[alpha],
                           NbrCells->qim1[alpha][direc], DATA);
      
      gphR = NbrCells->qip1[alpha][direc];
      fphR = -0.5*minmod_dx(NbrCells->qip2[alpha][direc],
                            NbrCells->qip1[alpha][direc], qi[alpha], DATA);
      
      gmhL = NbrCells->qim1[alpha][direc];
      fmhL = 0.5*minmod_dx(qi[alpha], NbrCells->qim1[alpha][direc], 
                           NbrCells->qim2[alpha][direc], DATA);
      
      gmhR = qi[alpha];
      fmhR = -0.5*minmod_dx(NbrCells->qip1[alpha][direc], qi[alpha],
                            NbrCells->qim1[alpha][direc], DATA);

     tempf = HalfwayCells->qiphL[alpha][direc] = gphL + fphL;
     if(!finite(tempf))
      {
       fprintf(stderr, "qiphL is not finite with g = %e, f = %e\n", gphL, fphL);
       fprintf(stderr, "alpha = %d\n", alpha);
       fprintf(stderr, "direc = %d\n", direc);
       exit(0);
      }
     tempf = HalfwayCells->qiphR[alpha][direc] = gphR + fphR;
     if(!finite(tempf))
      {
       fprintf(stderr, "qiphR is not finite with g = %e, f = %e\n", gphR, fphR);
       fprintf(stderr, "alpha = %d\n", alpha);
       fprintf(stderr, "direc = %d\n", direc);
       exit(0);
      }
     tempf = HalfwayCells->qimhL[alpha][direc] = gmhL + fmhL;
     if(!finite(tempf))
      {
       fprintf(stderr, "qimhL is not finite with g = %e, f = %e\n", gmhL, fmhL);
       fprintf(stderr, "alpha = %d\n", alpha);
       fprintf(stderr, "direc = %d\n", direc);
       exit(0);
      }
     tempf = HalfwayCells->qimhR[alpha][direc] = gmhR + fmhR;
     if(!finite(tempf))
      {
       fprintf(stderr, "qimhR is not finite with g = %e, f = %e\n", gmhR, fmhR);
       fprintf(stderr, "alpha = %d\n", alpha);
       fprintf(stderr, "direc = %d\n", direc);
       exit(0);
      }
    
   }/* direc */
  }/* alpha */

 return 1; /* if successful */
}/* MakeQIHalfs */


int Evolve::ConstHalfwayCells(double tau, BdryCells *HalfwayCells, double *qi, Grid *grid_pt,
			      InitData *DATA, int rk_flag, int size, int rank)
{
 int direc, flag;
 double epsilon_init, rhob_init;

 epsilon_init = grid_pt->epsilon;
 rhob_init = grid_pt->rhob;

 for(direc=1; direc<=3; direc++)
 {
 /* for each direction, reconstruct half-way cells */

  flag=
    reconst->ReconstIt(&(HalfwayCells->grid_p_h_L[direc]), 
		       direc, tau, HalfwayCells->qiphL, grid_pt,
		       epsilon_init, rhob_init, DATA, rk_flag); 
  if(flag==0) {
    reconst->ReconstError("grid_p_h_L", direc, rk_flag, qi, HalfwayCells->qiphL, grid_pt);
    return 0;
   }/* if Reconst returns error */

  flag=
    reconst->ReconstIt(&(HalfwayCells->grid_p_h_R[direc]), 
		       direc, tau, HalfwayCells->qiphR, grid_pt,
		       epsilon_init, rhob_init, DATA, rk_flag); 
  if(flag==0) {
    reconst->ReconstError("grid_p_h_R", direc, rk_flag, qi, HalfwayCells->qiphR, grid_pt);
    return 0;
   }
  
  flag=
    reconst->ReconstIt(&(HalfwayCells->grid_m_h_L[direc]), 
		       direc, tau, HalfwayCells->qimhL, grid_pt,
		       epsilon_init, rhob_init, DATA, rk_flag); 
  if(flag==0) {
    reconst->ReconstError("grid_m_h_L", direc, rk_flag, qi, HalfwayCells->qimhL, grid_pt);
    return 0;
  }
  
  flag=
    reconst->ReconstIt(&(HalfwayCells->grid_m_h_R[direc]), 
		       direc, tau, HalfwayCells->qimhR, grid_pt,
		       epsilon_init, rhob_init, DATA, rk_flag); 
  if(flag==0) {
    reconst->ReconstError("grid_m_h_R", direc, rk_flag, qi, HalfwayCells->qimhR, grid_pt);
    return 0;
   }

  if (HalfwayCells->grid_m_h_R[direc].rhob>1.&& HalfwayCells->grid_m_h_R[direc].epsilon<0.0000001) 
    HalfwayCells->grid_m_h_R[direc].rhob=0.;
 }/* direc */

//  if (HalfwayCells->grid_p_h_R[direc].rhob>2) cout << "phR rhob =" <<  HalfwayCells->grid_p_h_R[direc].rhob << endl;
//  if (HalfwayCells->grid_m_h_R[direc].rhob>1.) cout << "mhR rhob =" <<  HalfwayCells->grid_m_h_R[direc].rhob 
// 						   << "mhR rhob_init=" << rhob_init << " eps=" << HalfwayCells->grid_m_h_R[direc].epsilon 
// 						   <<" on rank " << rank << endl;
//  if (HalfwayCells->grid_p_h_L[direc].rhob>2) cout << "phL rhob =" <<  HalfwayCells->grid_p_h_L[direc].rhob << endl;
//  if (HalfwayCells->grid_m_h_L[direc].rhob>2) cout << "mhL rhob =" <<  HalfwayCells->grid_m_h_L[direc].rhob << endl;

 
 return 1;/* upon successful execution */
}/* ConstHalfwayCells */


void Evolve::MakeKTCurrents
 (double tau, double **DFmmp, Grid *grid_pt, 
  BdryCells *HalfwayCells, int rk_flag)
{
 int i, alpha;
 double FiphL[5][4], FiphR[5][4], FimhL[5][4], FimhR[5][4];
 double Fiph[5][4], Fimh[5][4], delta[4];
 double aiph[4], aimh[4], tau_fac[4], tempf;

 MakeMaxSpeedAs(tau, HalfwayCells, aiph, aimh, rk_flag);

/* Current J^i is calculated from the halfway cells in the same
 * direction because what we need is the divergence d_i J^i */
 /* d_tau tau Jtau + d_x (tau Jx) + d_y (tau Jy) + d_eta Jeta */

   tau_fac[1] = tau;
   tau_fac[2] = tau;
   tau_fac[3] = 1.0;

   for(alpha=0; alpha<5; alpha++)
    {
     for(i=1; i<=3; i++)
     {
      /* x_i current for TJb[0][alpha][0] from reconstructed halfway cells.
         Reconst only uses TJb[0] */
     FiphL[alpha][i] = HalfwayCells->grid_p_h_L[i].TJb[0][alpha][i]*tau_fac[i];
     FiphR[alpha][i] = HalfwayCells->grid_p_h_R[i].TJb[0][alpha][i]*tau_fac[i];
     FimhL[alpha][i] = HalfwayCells->grid_m_h_L[i].TJb[0][alpha][i]*tau_fac[i]; 
     FimhR[alpha][i] = HalfwayCells->grid_m_h_R[i].TJb[0][alpha][i]*tau_fac[i];

/* KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2})/2 
       - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2 */

     Fiph[alpha][i] = 0.5*(FiphL[alpha][i] + FiphR[alpha][i]);
     Fiph[alpha][i] -= 0.5*aiph[i]*(HalfwayCells->qiphR[alpha][i] 
                                  - HalfwayCells->qiphL[alpha][i]);
     Fimh[alpha][i] = 0.5*(FimhL[alpha][i] + FimhR[alpha][i]);
     Fimh[alpha][i] -= 0.5*aimh[i]*(HalfwayCells->qimhR[alpha][i] 
                                  - HalfwayCells->qimhL[alpha][i]);
    
tempf=DFmmp[alpha][i] = (Fimh[alpha][i] - Fiph[alpha][i]);
if(!finite(tempf))
{
 fprintf(stderr, "DFmmp[%d][%d] is not finite.\n", alpha, i);
 fprintf(stderr, "FimhL[%d][%d] is %e.\n", alpha, i, FimhL[alpha][i]);
 fprintf(stderr, "FiphL[%d][%d] is %e.\n", alpha, i, FiphL[alpha][i]);
 fprintf(stderr, "FimhR[%d][%d] is %e.\n", alpha, i, FimhR[alpha][i]);
 fprintf(stderr, "FiphR[%d][%d] is %e.\n", alpha, i, FiphR[alpha][i]);
 fprintf(stderr, "Fimh[%d][%d] is %e.\n", alpha, i, Fimh[alpha][i]);
 fprintf(stderr, "Fiph[%d][%d] is %e.\n", alpha, i, Fiph[alpha][i]);
 fprintf(stderr, "aimh[%d] = %e\n", i, aimh[i]);
 fprintf(stderr, "aiph[%d] = %e\n", i, aiph[i]);
 exit(0);
}
     }/* i - loop over directions */
    }/* alpha */
}/* MakeKTCurrents */


/* Calculate the right-hand-side */
/*
 du/dt = (1/Delta x)(F_imh - F_iph) + D
 F_iph(a) = (1/2)(f_a(qiphR) + f_a(qiphL))- (1/2)aiph(qiphR - qiphL)
 F_imh(a) = (1/2)(f_a(qimhR) + f_a(qimhL))- (1/2)aimh(qimhR - qimhL)
 RK 2nd order Heun's rules:

 k1 = f(t, u);
 u1 = u + k1*h
 k2 = f(t+h, u1);
 u2 = u1 + k2*h = u + (k1+k2)*h;
 u_next = u + (1/2)(k1+k2)*h
        = (1/2)(u + u2);
*/


void Evolve::MakeMaxSpeedAs
(double tau, 
 BdryCells *HalfwayCells,
 double aiph[], double aimh[], int rk_flag)
{
/* Implement Kurganov-Tadmor */

 int i;
 double aiphL[4], aiphR[4], aimhL[4], aimhR[4];

       for(i=1; i<=3; i++)
        {
         aiphL[i] = MaxSpeed(tau, i, &(HalfwayCells->grid_p_h_L[i]), rk_flag);
         aiphR[i] = MaxSpeed(tau, i, &(HalfwayCells->grid_p_h_R[i]), rk_flag);
         aimhL[i] = MaxSpeed(tau, i, &(HalfwayCells->grid_m_h_L[i]), rk_flag);
         aimhR[i] = MaxSpeed(tau, i, &(HalfwayCells->grid_m_h_R[i]), rk_flag);

         aiph[i] = maxi(aiphL[i], aiphR[i]);
         aimh[i] = maxi(aimhL[i], aimhR[i]);

        }/* i */
 return;
}/* MakeMaxSpeedAs */




double Evolve::MaxSpeed (double tau, int direc, Grid *grid_p, int rk_flag)
{
 double f, den, num;
 double rhob, utau, utau2, ux2, ux, uy, ueta, p, eps, h; 
 double deriv_p_rho, deriv_p_eps, deriv_p_h, rho_p_rho, h_p_h;
 double ut2mux2, ut, pe, rho, rpr;
 
 rhob = grid_p->rhob;

/* grid_p = grid_p_h_L, grid_p_h_R, grid_m_h_L, grid_m_h_R
   these are reconstructed by Reconst which only uses u[0] and TJb[0] */
 
 utau = (grid_p->u[0][0]);
 utau2 = utau*utau;
 ut = utau;

 /* to get the maximum value */
 ux = fabs((grid_p->u[0][direc]));
 ux2 = ux*ux;

 ut2mux2 = utau2-ux2;
 
 p = grid_p->p;
 eps = grid_p->epsilon;
 h = p+eps;
 
 pe = eos->p_e_func(eps, rhob);
 rpr = rhob*eos->p_rho_func(eps, rhob);

 den = -((1.0 + pe)*rpr*(-1.0 + util->Power(ut,2))) 
       + h*(pe + util->Power(ut,2) - pe*util->Power(ut,2));
 
 num = 
 sqrt(-((h*pe + rpr + pe*rpr)*(h*(pe*(-1.0 + ut2mux2) - ut2mux2) 
          + (1.0 + pe)*rpr*(-1.0 + ut2mux2)))) 
	  - h*(-1.0 + pe)*ut*ux - rpr*ut*ux - pe*rpr*ut*ux;

 if(-((h*pe + rpr + pe*rpr)*(h*(pe*(-1.0 + ut2mux2) - ut2mux2) 
			     + (1.0 + pe)*rpr*(-1.0 + ut2mux2)))<0) 
   {
     if(pe<0.001) 
       {
	 rpr=0.;
	 num = 
	   sqrt(-((h*pe + rpr + pe*rpr)*(h*(pe*(-1.0 + ut2mux2) - ut2mux2) 
					 + (1.0 + pe)*rpr*(-1.0 + ut2mux2)))) 
	   - h*(-1.0 + pe)*ut*ux - rpr*ut*ux - pe*rpr*ut*ux;
       }
     else 
       {
	 //fprintf(stderr,"num=%lf for pr=%lf\n",num,p_rho_func(eps, rhob));
	 fprintf(stderr,"WARNING: in MaxSpeed. Expression under sqrt in num=%lf. \n",
		 -((h*pe + rpr + pe*rpr)*(h*(pe*(-1.0 + ut2mux2) - ut2mux2) + (1.0 + pe)*rpr*(-1.0 + ut2mux2))));
	 fprintf(stderr,"at value e=%lf. \n",eps);
	 fprintf(stderr,"at value p=%lf. \n",p);
	 //fprintf(stderr,"at value rho=%lf. \n",rhob);
	 //fprintf(stderr,"at value prho=%lf. \n",p_rho_func(eps, rhob));
	 //fprintf(stderr,"at value h=%lf. \n",h);
	 fprintf(stderr,"at value rpr=%lf. \n",rpr);
	 fprintf(stderr,"at value pe=%lf. \n",pe);
	 fprintf(stderr,"at value (h*pe + rpr + pe*rpr)=%lf. \n",(h*pe + rpr + pe*rpr));
	 exit(0);
       }
   
   }
// if(den == 0.0) den += SMALL; 

 if( (num == 0.0) && (den != 0.0) )
  {
   f = 0.0;
  }
 else if( (num != 0.0) && (den != 0.0) )
  {
   f = num/den;
  }
 else if( (num == 0.0) && (den == 0.0) )
  {
   /* nothing there */
   f = 0.0;
  }
 else /* (num != 0.0) && (den == 0.0) */
  {
   fprintf(stderr, "SpeedMax = is infinite.\n");
   fprintf(stderr, "Can't happen.\n");
   exit(0);
  }


 if(f < 0.0) 
  {
   fprintf(stderr, "SpeedMax = %e\n is negative.\n", f);
   fprintf(stderr, "Can't happen.\n");
   exit(0);
  }
 else if(f <  ux/ut) 
  {
    if(num != 0.0)
    {
     if(fabs(f-ux/ut)<0.0001) f = ux/ut;
     else 
       {
	 fprintf(stderr, "SpeedMax-v = %lf\n", f-ux/ut);
	 fprintf(stderr, "SpeedMax = %e\n is smaller than v = %e.\n", f, ux/ut);
	 fprintf(stderr, "Can't happen.\n");
	 exit(0);
       }
    }
  }
 else if(f >  1.0) 
  {
   fprintf(stderr, "SpeedMax = %e\n is bigger than 1.\n", f);
   fprintf(stderr, "Can't happen.\n");
   exit(0);
  }

 if(direc == 3) f /= tau;
 
 return f;
}/* MaxSpeed */


	
double Evolve::minmod_dx(double up1, double u, double um1, InitData *DATA)
{
 double theta, diffup, diffdown, diffmid;
 double tempf;
 
 theta = DATA->minmod_theta;
 
 diffup = up1 - u;
 diffup *= theta;

 diffdown = u - um1;
 diffdown *= theta;

 diffmid = up1 - um1;
 diffmid *= 0.5;

 if( (diffup > 0.0) && (diffdown > 0.0) && (diffmid > 0.0) )
  {
   tempf = mini(diffdown, diffmid);
   return mini(diffup, tempf);
  }
 else if( (diffup < 0.0) && (diffdown < 0.0) && (diffmid < 0.0) )
  {
   tempf = maxi(diffdown, diffmid);
   return maxi(diffup, tempf);
  }
 else return 0.0;

}/* minmod_dx */


void Evolve::InitNbrQs(NbrQs *NbrCells)
{
 (NbrCells->qip1) = util->mtx_malloc(5,4);
 (NbrCells->qip2) = util->mtx_malloc(5,4);
 (NbrCells->qim1) = util->mtx_malloc(5,4);
 (NbrCells->qim2) = util->mtx_malloc(5,4);
}/* InitNbrQs */


void Evolve::InitTempGrids(BdryCells *HalfwayCells, int rk_order)
{
 int direc;

 (HalfwayCells->grid_p_h_L) = grid->grid_v_malloc(4);
 (HalfwayCells->grid_p_h_R) = grid->grid_v_malloc(4);
 (HalfwayCells->grid_m_h_L) = grid->grid_v_malloc(4);
 (HalfwayCells->grid_m_h_R) = grid->grid_v_malloc(4);

 HalfwayCells->qiphL = util->mtx_malloc(5,4);
 HalfwayCells->qiphR = util->mtx_malloc(5,4);
 HalfwayCells->qimhL = util->mtx_malloc(5,4);
 HalfwayCells->qimhR = util->mtx_malloc(5,4);
 
 for(direc=0; direc<4; direc++)
  {
    (HalfwayCells->grid_p_h_L)[direc].TJb = util->cube_malloc(rk_order,5,4);
    (HalfwayCells->grid_p_h_R)[direc].TJb = util->cube_malloc(rk_order,5,4);
    (HalfwayCells->grid_m_h_L)[direc].TJb = util->cube_malloc(rk_order,5,4);
    (HalfwayCells->grid_m_h_R)[direc].TJb = util->cube_malloc(rk_order,5,4);
    
    (HalfwayCells->grid_p_h_L)[direc].u = util->mtx_malloc(rk_order,4);
    (HalfwayCells->grid_p_h_R)[direc].u = util->mtx_malloc(rk_order,4);
    (HalfwayCells->grid_m_h_L)[direc].u = util->mtx_malloc(rk_order,4);
    (HalfwayCells->grid_m_h_R)[direc].u = util->mtx_malloc(rk_order,4);
   }
 return;
}/* InitTempGrids */



void Evolve::UpdateTJbRK(Grid *grid_rk, Grid *grid_pt, InitData *DATA, int rk_flag)
{
 int trk_flag, mu, alpha;

 trk_flag = rk_flag+1;
 grid_pt->p_t = grid_rk->p;
 grid_pt->epsilon_t = grid_rk->epsilon;
 grid_pt->rhob_t = grid_rk->rhob;

/* reconstructed grid_rk uses rk_flag 0 only */
 for(mu=0; mu<4; mu++)
  {
   grid_pt->u[trk_flag][mu] = grid_rk->u[0][mu];
   for(alpha=0; alpha<5; alpha++)
    {
      grid_pt->TJb[trk_flag][alpha][mu] = grid_rk->TJb[0][alpha][mu];
    }/* alpha */
  }/* mu */
}/* UpdateTJbRK */




void Evolve::FindFreezeOutSurface(double tau, InitData *DATA, Grid ***arena, int size, int rank)
{	
  FILE *t_file;
  char* t_name = "tauf.dat";
  t_file = fopen(t_name, "a");
  char *buf;
  buf = util->char_malloc(40);
  FILE *s_file;
  char* s_name;
  s_name = util->char_malloc(100);
  sprintf (buf, "%d", rank);
  strcat(s_name, "surface");
  strcat(s_name,buf);
  strcat(s_name, ".dat");

  s_file = fopen(s_name, "a");
  //fprintf(t_file,"x [fm], tau_f(x) [fm] \n");
  int ix, iy, ieta, nx, ny, neta;
  double x, y, eta;
  double epsFO=DATA->epsilonFreeze/hbarc;
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
	       if((arena[ix][iy][ieta].epsilon-epsFO)*(arena[ix][iy][ieta].epsilon_prev[facTau-1]-epsFO)<0.)
		{
		  if (arena[ix][iy][ieta].epsilon>epsFO)
		    SIG=-1.;
		  else SIG=1.;
		  intersect = 1;
		  intersecttau = 1;
		  tauf = tau - DTAU + DTAU * (arena[ix][iy][ieta].epsilon_prev[facTau-1]-epsFO)
		    /(arena[ix][iy][ieta].epsilon_prev[facTau-1]-arena[ix][iy][ieta].epsilon);
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
	      
	      utauX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][0]);
	      utauX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][0]);
	      utauX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][0]);
	      utauX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][0]);
	      
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
	      
	      uxX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][1]);
	      uxX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][1]);
	      uxX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][1]);
	      uxX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][1]);
	      
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
	      
	      uyX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][2]);
	      uyX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][2]);
	      uyX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][2]);
	      uyX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][2]);
	      
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
	      
	      uetaX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][3]);
	      uetaX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][3]);
	      uetaX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][3]);
	      uetaX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][3]);
	      uetaY1 = ((1.-yfrac)*uetaX1+yfrac*uetaX2);
	      uetaY2 = ((1.-yfrac)*uetaX3+yfrac*uetaX4);
	      
	      ueta2 = ((1.-etafrac)*uetaY1+etafrac*uetaY2);
	      
	      ueta = (1.-taufrac)*ueta2+taufrac*ueta1;
	      
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
		  fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
			  tauf, x, y, eta, FULLSU[0], 0., 0., 0.,
			  utau, ux, uy, ueta, epsFO, TFO, muB);
		  if(fabs(x)<0.05 && (fabs(y)<0.05))
		    {
		      fprintf(t_file,"%f %f %f %f %f %f %f %f\n",tauf, x, y, eta, FULLSU[0],0.,0.,0.);
		    }
		}
	      if (intersectx)
		{
		  fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
			  tau, xf, y, eta, 0., FULLSU[1], 0., 0.,
			  utau, ux, uy, ueta, epsFO, TFO, muB);
		  if(fabs(x)<0.05 && (fabs(y)<0.05))
		    {
		      fprintf(t_file,"%f %f %f %f %f %f %f %f\n",tau, xf, y, eta, 0., FULLSU[1],0.,0.);
		    }
		}
	      if (intersecty)
		{
		  fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
			  tau, x, yf, eta, 0., 0., FULLSU[2], 0.,
			  utau, ux, uy, ueta, epsFO, TFO, muB);
		  if(fabs(x)<0.05 && (fabs(y)<0.05))
		    {
		      fprintf(t_file,"%f %f %f %f %f %f %f %f\n",tau, x, yf, eta, 0.,0.,FULLSU[2],0.);
		    }
		}
	      if (intersecteta)
		{
		  fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
			  tau, x, y, etaf, 0., 0., 0., FULLSU[3],
			  utau, ux, uy, ueta, epsFO, TFO, muB);
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


void Evolve::FindFreezeOutSurface2(double tau, InitData *DATA, Grid ***arena, int size, int rank)
{	
  FILE *t_file;
  char* t_name = "tauf.dat";
  t_file = fopen(t_name, "a");
  FILE *t2_file;
  char* t2_name = "taufx.dat";
  t2_file = fopen(t2_name, "a");
  FILE *t3_file;
  char* t3_name = "taufy.dat";
  t3_file = fopen(t3_name, "a");
  char *buf;
  buf = util->char_malloc(40);
  FILE *s_file;
  char* s_name;
  s_name = util->char_malloc(100);
  
  sprintf (buf, "%d", rank);
  strcat(s_name, "surface");
  strcat(s_name,buf);
  strcat(s_name, ".dat");

  s_file = fopen(s_name, "a");

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
  double rhob, utau, ux, uy, ueta, TFO, muB;
  double utauX1, utauX2, utauX3, utauX4, utauY1, utauY2, utau1, utau2;
  double rhobX1, rhobX2, rhobX3, rhobX4, rhobY1, rhobY2, rhob1, rhob2;
  double uxX1, uxX2, uxX3, uxX4, uxY1, uxY2, ux1, ux2;
  double uyX1, uyX2, uyX3, uyX4, uyY1, uyY2, uy1, uy2;
  double uetaX1, uetaX2, uetaX3, uetaX4, uetaY1, uetaY2, ueta1, ueta2;
  double xfrac, yfrac, etafrac, taufrac;
  // who is direct neighbor (1) or neighbor across the plane (2), more distant neighbor (3), (4)
  int intersections;
  int maxEta;
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

  double *packageutau_prev;
  double *packageux_prev;
  double *packageuy_prev;
  double *packageueta_prev;
  double *packagerhob_prev;
  
  double **Rneighbor_eps;
  double **Rneighbor_eps_prev;
  double **Rneighbor_utau;
  double **Rneighbor_ux;
  double **Rneighbor_uy;
  double **Rneighbor_ueta;
  double **Rneighbor_rhob;
 
  double **Rneighbor_utau_prev;
  double **Rneighbor_ux_prev;
  double **Rneighbor_uy_prev;
  double **Rneighbor_ueta_prev;
  double **Rneighbor_rhob_prev;

  package = new double[sizeOfData];
  package_prev = new double[sizeOfData];
  packageutau = new double[sizeOfData];
  packageux = new double[sizeOfData];
  packageuy = new double[sizeOfData];
  packageueta = new double[sizeOfData];
  packagerhob = new double[sizeOfData];
  
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
  
  Rneighbor_utau_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_ux_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_uy_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_ueta_prev = util->mtx_malloc(nx+1,ny+1);
  Rneighbor_rhob_prev = util->mtx_malloc(nx+1,ny+1);

  
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
	      package_prev[position] = arena[ix][iy][0].epsilon_prev[facTau-1];
	      packageutau[position] = arena[ix][iy][0].u[0][0];
	      packageux[position] = arena[ix][iy][0].u[0][1];
	      packageuy[position] = arena[ix][iy][0].u[0][2];
	      packageueta[position] = arena[ix][iy][0].u[0][3];
	      packagerhob[position] = arena[ix][iy][0].rhob;
	      packageutau_prev[position] = arena[ix][iy][0].u_prev[facTau-1][0];
	      packageux_prev[position] = arena[ix][iy][0].u_prev[facTau-1][1];
	      packageuy_prev[position] = arena[ix][iy][0].u_prev[facTau-1][2];
	      packageueta_prev[position] = arena[ix][iy][0].u_prev[facTau-1][3];
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
      //cout << " done sending to the left on rank " << rank << endl;
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
      //cout << " receiving from the right on rank " << rank << endl;
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
		  if((arena[ix+fac][iy+fac][ieta+fac].epsilon-epsFO)*(arena[ix][iy][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
		    if((arena[ix+fac][iy][ieta].epsilon-epsFO)*(arena[ix][iy+fac][ieta+fac].epsilon_prev[facTau-1]-epsFO)>0.)
		      if((arena[ix][iy+fac][ieta].epsilon-epsFO)*(arena[ix+fac][iy][ieta+fac].epsilon_prev[facTau-1]-epsFO)>0.)
			if((arena[ix][iy][ieta+fac].epsilon-epsFO)*(arena[ix+fac][iy+fac][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
			  if((arena[ix+fac][iy+fac][ieta].epsilon-epsFO)*(arena[ix][iy][ieta+fac].epsilon_prev[facTau-1]-epsFO)>0.)
			    if((arena[ix+fac][iy][ieta+fac].epsilon-epsFO)*(arena[ix][iy+fac][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
			      if((arena[ix][iy+fac][ieta+fac].epsilon-epsFO)*(arena[ix+fac][iy][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
				if((arena[ix][iy][ieta].epsilon-epsFO)*(arena[ix+fac][iy+fac][ieta+fac].epsilon_prev[facTau-1]-epsFO)>0.)
				  intersect=0;
		}
	      else // if this is the right most edge
		{
		  intersect=1;
		  if((Rneighbor_eps[ix+fac][iy+fac] - epsFO)*(arena[ix][iy][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
		    if((arena[ix+fac][iy][ieta].epsilon-epsFO)*(Rneighbor_eps_prev[ix][iy+fac]-epsFO)>0.)
		      if((arena[ix][iy+fac][ieta].epsilon-epsFO)*(Rneighbor_eps_prev[ix+fac][iy]-epsFO)>0.)
			if((Rneighbor_eps[ix][iy]-epsFO)*(arena[ix+fac][iy+fac][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
			  if((arena[ix+fac][iy+fac][ieta].epsilon-epsFO)*(Rneighbor_eps_prev[ix][iy]-epsFO)>0.)
			    if((Rneighbor_eps[ix+fac][iy]-epsFO)*(arena[ix][iy+fac][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
			      if((Rneighbor_eps[ix][iy+fac]-epsFO)*(arena[ix+fac][iy][ieta].epsilon_prev[facTau-1]-epsFO)>0.)
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
		      cube[0]=arena[ix][iy][ieta].epsilon_prev[facTau-1];
		      cube[1]=arena[ix+fac][iy][ieta].epsilon_prev[facTau-1];
		      cube[2]=arena[ix+fac][iy+fac][ieta].epsilon_prev[facTau-1];
		      cube[3]=arena[ix][iy+fac][ieta].epsilon_prev[facTau-1];
		      cube[4]=arena[ix][iy][ieta].epsilon;
		      cube[5]=arena[ix+fac][iy][ieta].epsilon;
		      cube[6]=arena[ix+fac][iy+fac][ieta].epsilon;
		      cube[7]=arena[ix][iy+fac][ieta].epsilon;
		      cube[8]=arena[ix][iy][ieta+fac].epsilon_prev[facTau-1];
		      cube[9]=arena[ix+fac][iy][ieta+fac].epsilon_prev[facTau-1];
		      cube[10]=arena[ix+fac][iy+fac][ieta+fac].epsilon_prev[facTau-1];
		      cube[11]=arena[ix][iy+fac][ieta+fac].epsilon_prev[facTau-1];
		      cube[12]=arena[ix][iy][ieta+fac].epsilon;
		      cube[13]=arena[ix+fac][iy][ieta+fac].epsilon;
		      cube[14]=arena[ix+fac][iy+fac][ieta+fac].epsilon;
		      cube[15]=arena[ix][iy+fac][ieta+fac].epsilon;
		    }
		  else
		    {
		      cube[0]=arena[ix][iy][ieta].epsilon_prev[facTau-1];
		      cube[1]=arena[ix+fac][iy][ieta].epsilon_prev[facTau-1];
		      cube[2]=arena[ix+fac][iy+fac][ieta].epsilon_prev[facTau-1];
		      cube[3]=arena[ix][iy+fac][ieta].epsilon_prev[facTau-1];
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
		      fprintf(stderr,"NUMBER OF SINGLE EDGES IS NOT A MULTIPLE OF 3, mmhhh... number=%d\n",countSingleEdges);
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
		      //continue; // don't add the flawed cells at all
		    }
	  
		  //fprintf(stderr,"Volume=%f\n", SUM);
		  //fprintf(stderr,"Cells done=%d\n", cells);
		  //fprintf(stderr,"Warnings=%d\n", warnings);
		  fprintf(stderr,"percent error=%f\n", 100*(warnings*1.0)/(cells*1.0));
		  
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
		
		  utauX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][0]);
		  utauX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][0] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][0]);
		  
		  if (ieta<neta-fac)
		    {
		      utauX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][0] 
				+ xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][0]);
		      utauX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][0] 
				+ xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][0]);
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
			  
		  uxX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][1]);
		  uxX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][1]);
		  

		  if (ieta<neta-fac)
		    {
		      uxX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][1]);
		      uxX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][1] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][1]);
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
			  
		  uyX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][2]);
		  uyX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][2]);

		  if (ieta<neta-fac)
		    {
		      uyX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][2]);
		      uyX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][2] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][2]);
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
			  
		  uetaX1 = ((1.-xfrac)*arena[ix][iy][ieta].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy][ieta].u_prev[facTau-1][3]);
		  uetaX2 = ((1.-xfrac)*arena[ix][iy+fac][ieta].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy+fac][ieta].u_prev[facTau-1][3]);

		  if (ieta<neta-fac)
		    {
		      uetaX3 = ((1.-xfrac)*arena[ix][iy][ieta+fac].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy][ieta+fac].u_prev[facTau-1][3]);
		      uetaX4 = ((1.-xfrac)*arena[ix][iy+fac][ieta+fac].u_prev[facTau-1][3] + xfrac*arena[ix+fac][iy+fac][ieta+fac].u_prev[facTau-1][3]);
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
		  fprintf(s_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",
			  tauf, xf, yf, etaf, FULLSU[0],FULLSU[1],FULLSU[2],FULLSU[3],
			  utau, ux, uy, ueta, epsFO, TFO, muB);
		  
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
      if (allIntersections==0)
	{
	  cout << "All cells frozen out. Exiting." << endl;
	  MPI::Finalize();
 	  exit(1);
	}
    }
  if (rank!=0)
    {
      MPI::COMM_WORLD.Recv(allIntersectionsArray,1,MPI::INT,0,2);
      if (allIntersectionsArray[0]==0)
	{
	  cout << "All cells frozen out. Exiting." << endl;
	  MPI::Finalize();
 	  exit(1);
	}
    }
  
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
