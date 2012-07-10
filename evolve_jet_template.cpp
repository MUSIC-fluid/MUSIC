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
   
   //plug in rates here:
   
   if ( DATA->includeJet>=1 )
     {
       double etaWidth = asinh(0.3/tau); // keep a fixed width in z
       double sTau;
       if (rk_flag==0)
	 sTau=tau;
       else if(rk_flag==1)
	 sTau=tau+DATA->delta_tau;
       
       double width = 0.3+0.00065*pow(sTau,2.5);
       double etaWidth = asinh(width/tau); // keep a fixed width in z

       //first jet (out of 2, actually four with trigger side) +px +py, trigger: -px +py
       
       if( jetPosition==1)
	 {
	   // adding for both RK steps
	   xjet=;
	   yjet=;
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==2)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==3)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf +=
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==4)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==5)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf +=
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==6)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==7)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==8)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
		 }
       else if (jetPosition==9)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==10)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==11)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==12)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }

	 }
       else if (jetPosition==13)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }

	 }
       else if (jetPosition==14)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }

	 }
       else if (jetPosition==15)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }

	 }
       else if (jetPosition==16)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }

	 }
       else if (jetPosition==17)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==18)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==19)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==20)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf +=
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==21)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==22)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==23)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==24)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==25)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==26)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==27)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==28)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==29)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==30)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==31)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==32)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==33)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==34)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==35)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==36)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==37)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;
	       
	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==38)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==39)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==40)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==41)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==42)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==43)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
     
       //second jet (out of 2, actually four with trigger side)  +px -py, trigger: -px -py
      
       if( jetPosition==1)
	 {
	   // adding for both RK steps
	   xjet=;
	   yjet=;
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==2)
	 {
	   xjet=sTau-0.4-6.25;
	   yjet=1.;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==3)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf +=
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==4)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==5)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf +=
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==6)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==7)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==8)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	 }
       else if (jetPosition==9)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==10)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==11)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==12)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }

	 }
       else if (jetPosition==13)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }

	 }
       else if (jetPosition==14)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }

	 }
       else if (jetPosition==15)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }

	 }
       else if (jetPosition==16)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }

	 }
       else if (jetPosition==17)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==18)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==19)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==20)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf +=
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==21)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==22)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==23)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==24)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==25)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==26)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==27)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==28)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==29)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==30)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==31)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==32)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==33)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==34)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==35)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==36)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==37)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==38)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==39)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==40)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==41)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==42)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	     }
	 }
       else if (jetPosition==43)
	 {
	   xjet=;
	   yjet=;
	   
	   if((alpha==0) && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==1 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   else if (alpha==2 && eta>-1 && eta < 1)
	     {
	       sumf += 
		 /hbarc*exp(-(pow((x-xjet),2.)/2/pow(width,2.)+pow((y-yjet),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		 *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
	     }
	   if (includeTrigger==1)
	     {
	       xtrigger=;
	       ytrigger=;

	       if((alpha==0) && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==1 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
		 }
	       else if (alpha==2 && eta>-1 && eta < 1)
		 {
		   sumf += 
		     /hbarc*exp(-(pow((x-xtrigger),2.)/2/pow(width,2.)+pow((y-ytrigger),2.)/2/pow(width,2.)+eta*eta/2/pow(etaWidth,2.)))
		     *1/pow((sqrt(2*PI)),3.)/(width*width*etaWidth);
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
 rpr = rhob*eos->p_rho_func(eps, rhob)/(1.+pe); // fixed this, June 15th 2010 times 1/(1.+pe)

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
