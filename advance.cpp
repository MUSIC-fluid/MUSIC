#include "util.h"
#include "data.h"
#include "grid.h"
#include "eos.h"
#include "evolve.h"
#include "advance.h"

using namespace std;

Advance::Advance(EOS *eosIn, Grid *gridIn)
{
  eos = new EOS;
  eos = eosIn;
  grid = new Grid;
  grid = gridIn;
  reconst = new Reconst(eosIn, grid);
  util = new Util;
  diss = new Diss(eosIn);
  minmod = new Minmod;
  u_derivative = new U_derivative;
}

// destructor
Advance::~Advance()
{
  delete reconst;
  delete util;
  delete diss;
  delete grid;
  delete eos;
  delete minmod;
  delete u_derivative;
}

int Advance::AdvanceIt(double tau, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int rk_flag, int size, int rank)
{
  int ix, iy, ieta;
//   int flag;
  //cout << "0 AdvanceIt Wmunu=" << (Lneighbor[1][1][0]).Wmunu[rk_flag][1][1] << endl;
  
  for(ix=0; ix<=DATA->nx; ix++)
    {
      for(iy=0; iy<=DATA->ny; iy++)
	{
          for(ieta=0; ieta<DATA->neta; ieta++)
	    {
	      AdvanceLocalT(tau, DATA, &(arena[ix][iy][ieta]), &(Lneighbor[ix][iy][0]), &(Rneighbor[ix][iy][0]), 
				   &(Lneighbor[ix][iy][1]), &(Rneighbor[ix][iy][1]), rk_flag, size, rank);
	    }/* ieta */
	}/*iy */
    }/* ix */
  
  MPISendReceiveT(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);

  if(DATA->viscosity_flag == 1)
    {
      for(ix=0; ix<=DATA->nx; ix++)
	{
	  for(iy=0; iy<=DATA->ny; iy++)
	    {
	      for(ieta=0; ieta<DATA->neta; ieta++)
		{
		  AdvanceLocalW(tau, DATA, &(arena[ix][iy][ieta]), &(Lneighbor[ix][iy][0]), &(Rneighbor[ix][iy][0]), 
				       &(Lneighbor[ix][iy][1]), &(Rneighbor[ix][iy][1]), rk_flag, size, rank);
		}/* ieta */
	    }/*iy */
	}/* ix */
    }/* if viscosity flag is set */
  
  MPISendReceiveW(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);
  
  return 1;
}/* AdvanceIt */

void Advance::MPISendReceiveT(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag)
{
  // this sends and receives information from neighboring cells in the next processor in the eta direction 
  // and stores it in Lneighbor and Rneighbor (unless the processor is really at the edge of the total grid

  int ix, iy, nx, ny, neta, i, alpha;
//   int ieta, iflag;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;
  int sizeOfData = 5*3*(nx+1)*(ny+1);
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
		  for(i=0; i<=DATA->rk_order; i++) // rk_order 2 means rk_flag = 0, 1, 2
		    {
		      //if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
		      //cout << "arena[ix][iy][0].TJb[i][alpha][0]=" << arena[ix][iy][0].TJb[i][alpha][0] << endl;
		      position = (i+3*(alpha+(5*(ix + ((nx+1)*iy)))));
		      package[position]  = arena[ix][iy][0].TJb[i][alpha][0];
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
		  for(i=0; i<=DATA->rk_order; i++)
		    {
		      position = (i+3*(alpha+(5*(ix + ((nx+1)*iy)))));
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
		  for(i=0; i<=DATA->rk_order; i++)
		    {
		      position = (i+3*(alpha+(5*(ix + ((nx+1)*iy)))));
		      package[position]  = arena[ix][iy][neta-1].TJb[i][alpha][0];
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
		  for(i=0; i<=DATA->rk_order; i++)
		    {
		      position = (i+3*(alpha+(5*(ix + ((nx+1)*iy)))));
		      Lneighbor[ix][iy][0].TJb[i][alpha][0] = package[position];
		      Lneighbor[ix][iy][1].TJb[i][alpha][0] = package2[position];
		    }
		}
	    }
	}
      //      cout << " ++"<<rank<<"++ "<<  " receiving from the left " << endl;
    }
  //  if(rank==0)
    //   cout << "on rank " <<rank<< " " <<  " MPISendReceiveT done " << endl;
  delete[] package;
  delete[] package2;
}//end MPISendReceive



void Advance::MPISendReceiveW(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag)
{
  // this sends and receives information from neighboring cells in the next processor in the eta direction 
  // and stores it in Lneighbor and Rneighbor (unless the processor is really at the edge of the total grid
  
  int ix, iy, nx, ny, neta, i, alpha, beta;
//   int ieta, iflag;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;

  int sizeOfDataDis = 5*3*(nx+1)*(ny+1)*4; // size of data package for W's and Pi's
  int sizeOfDataU = 4*3*(nx+1)*(ny+1); // size of data package for u's
  int sizeOfDataPi_b = 3*(nx+1)*(ny+1); // size of data package for pi_b's

//   int position;
  int positionDis, positionU, positionPi_b;

  double *packageW; // dissipative parts
  double *packagePi;
  double *packageW2; 
  double *packagePi2;
  double *packageU; 
  double *packageU2;
  double *packagePi_b;
  double *packagePi_b2;
  
  packageW = new double[sizeOfDataDis];
  packageW2 = new double[sizeOfDataDis];
  packagePi = new double[sizeOfDataDis];
  packagePi2 = new double[sizeOfDataDis];
  packageU = new double[sizeOfDataU];
  packageU2 = new double[sizeOfDataU];
  packagePi_b = new double[sizeOfDataPi_b];
  packagePi_b2 = new double[sizeOfDataPi_b];
  
  // receive from the right / send to the left
  int from = rank+1;
  int to = rank-1;
  // packing the package to send
  if ( rank != 0 )
    {
      //cout << " sending to the left on rank " << rank << endl;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(i=0; i<=DATA->rk_order; i++) // rk_flag iterator
		{
		  //if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
		  //cout << "arena[ix][iy][0].TJb[i][alpha][0]=" << arena[ix][iy][0].TJb[i][alpha][0] << endl;
		  positionPi_b = i+3*(ix + ((nx+1)*iy));
		  packagePi_b[positionPi_b] = arena[ix][iy][0].pi_b[i];
		  packagePi_b2[positionPi_b] = arena[ix][iy][1].pi_b[i];
		  for(beta=0; beta<4; beta++) // dissipative part
		    {
		      for(alpha=0; alpha<5; alpha++)
			{
			  positionDis = beta+4*(i+3*(alpha+(5*(ix + ((nx+1)*iy)))));
			  packageW[positionDis] = arena[ix][iy][0].Wmunu[i][alpha][beta];
			  packageW2[positionDis] = arena[ix][iy][1].Wmunu[i][alpha][beta];
			  packagePi[positionDis] = arena[ix][iy][0].Pimunu[i][alpha][beta];
			  packagePi2[positionDis] = arena[ix][iy][1].Pimunu[i][alpha][beta];
			}	  
		      positionU = (i+3*(beta+(4*(ix + ((nx+1)*iy)))));
		      packageU[positionU] = arena[ix][iy][0].u[i][beta];
		      packageU2[positionU] = arena[ix][iy][1].u[i][beta];
		    }
		}
	    }
	}
      
      MPI::COMM_WORLD.Send(packageW,sizeOfDataDis,MPI::DOUBLE,to,30);
      MPI::COMM_WORLD.Send(packageW2,sizeOfDataDis,MPI::DOUBLE,to,40);
      MPI::COMM_WORLD.Send(packagePi,sizeOfDataDis,MPI::DOUBLE,to,50);
      MPI::COMM_WORLD.Send(packagePi2,sizeOfDataDis,MPI::DOUBLE,to,60);
      MPI::COMM_WORLD.Send(packageU,sizeOfDataU,MPI::DOUBLE,to,70);
      MPI::COMM_WORLD.Send(packageU2,sizeOfDataU,MPI::DOUBLE,to,80);
      MPI::COMM_WORLD.Send(packagePi_b,sizeOfDataPi_b,MPI::DOUBLE,to,90);
      MPI::COMM_WORLD.Send(packagePi_b2,sizeOfDataPi_b,MPI::DOUBLE,to,100);
      //cout << " done sending to the left on rank " << rank << endl;
    }
  // receiving and unwrapping the package
  if ( rank != size-1 )
    {  
      MPI::COMM_WORLD.Recv(packageW,sizeOfDataDis,MPI::DOUBLE,from,30);
      MPI::COMM_WORLD.Recv(packageW2,sizeOfDataDis,MPI::DOUBLE,from,40);
      MPI::COMM_WORLD.Recv(packagePi,sizeOfDataDis,MPI::DOUBLE,from,50);
      MPI::COMM_WORLD.Recv(packagePi2,sizeOfDataDis,MPI::DOUBLE,from,60);
      MPI::COMM_WORLD.Recv(packageU,sizeOfDataU,MPI::DOUBLE,from,70);
      MPI::COMM_WORLD.Recv(packageU2,sizeOfDataU,MPI::DOUBLE,from,80);
      MPI::COMM_WORLD.Recv(packagePi_b,sizeOfDataPi_b,MPI::DOUBLE,from,90);
      MPI::COMM_WORLD.Recv(packagePi_b2,sizeOfDataPi_b,MPI::DOUBLE,from,100);
      //cout << " receiving from the right on rank " << rank << endl;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(i=0; i<=DATA->rk_order; i++)
		{
		  positionPi_b = i+3*(ix + ((nx+1)*iy));
		  Rneighbor[ix][iy][0].pi_b[i] = packagePi_b[positionPi_b];
		  Rneighbor[ix][iy][1].pi_b[i] = packagePi_b2[positionPi_b];
		  for(beta=0; beta<4; beta++) // dissipative part
		    {
		      for(alpha=0; alpha<5; alpha++)
			{
			  positionDis = beta+4*(i+3*(alpha+(5*(ix + ((nx+1)*iy)))));
			  Rneighbor[ix][iy][0].Wmunu[i][alpha][beta] = packageW[positionDis];
			  Rneighbor[ix][iy][1].Wmunu[i][alpha][beta] = packageW2[positionDis];
			  Rneighbor[ix][iy][0].Pimunu[i][alpha][beta] = packagePi[positionDis];
			  Rneighbor[ix][iy][1].Pimunu[i][alpha][beta] = packagePi2[positionDis];
			  //cout << "received from right Wmunu=" <<  Rneighbor[ix][iy][0].Wmunu[i][alpha][beta] << endl;
			}
		      positionU = (i+3*(beta+(4*(ix + ((nx+1)*iy)))));
		      Rneighbor[ix][iy][0].u[i][beta] = packageU[positionU];
		      Rneighbor[ix][iy][1].u[i][beta] = packageU2[positionU];
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
	      for(i=0; i<=DATA->rk_order; i++)
		{
		  positionPi_b = i+3*(ix + ((nx+1)*iy));
		  packagePi_b[positionPi_b] = arena[ix][iy][neta-1].pi_b[i];
		  packagePi_b2[positionPi_b] = arena[ix][iy][neta-2].pi_b[i];
		  for(beta=0; beta<4; beta++) // dissipative part
		    {
		      for(alpha=0; alpha<5; alpha++)
			{
			  positionDis = beta+4*(i+3*(alpha+(5*(ix + ((nx+1)*iy)))));
			  packageW[positionDis] = arena[ix][iy][neta-1].Wmunu[i][alpha][beta];
			  packageW2[positionDis] = arena[ix][iy][neta-2].Wmunu[i][alpha][beta];
			  packagePi[positionDis] = arena[ix][iy][neta-1].Pimunu[i][alpha][beta];
			  packagePi2[positionDis] = arena[ix][iy][neta-2].Pimunu[i][alpha][beta];
			}
		      positionU = (i+3*(beta+(4*(ix + ((nx+1)*iy)))));
		      packageU[positionU] = arena[ix][iy][neta-1].u[i][beta];
		      packageU2[positionU] = arena[ix][iy][neta-2].u[i][beta];
		    }
		}
	    }
	}
      
      //cout << "sending to the right Wmunu=" <<  arena[1][1][neta-1].Wmunu[0][1][1] << endl;
      
      MPI::COMM_WORLD.Send(packageW,sizeOfDataDis,MPI::DOUBLE,to,130);
      MPI::COMM_WORLD.Send(packageW2,sizeOfDataDis,MPI::DOUBLE,to,140);
      MPI::COMM_WORLD.Send(packagePi,sizeOfDataDis,MPI::DOUBLE,to,150);
      MPI::COMM_WORLD.Send(packagePi2,sizeOfDataDis,MPI::DOUBLE,to,160);
      MPI::COMM_WORLD.Send(packageU,sizeOfDataU,MPI::DOUBLE,to,170);
      MPI::COMM_WORLD.Send(packageU2,sizeOfDataU,MPI::DOUBLE,to,180);
      MPI::COMM_WORLD.Send(packagePi_b,sizeOfDataPi_b,MPI::DOUBLE,to,190);
      MPI::COMM_WORLD.Send(packagePi_b2,sizeOfDataPi_b,MPI::DOUBLE,to,200);
      //cout << " **"<<rank<<"** "<<  " done sending to the right " << endl;
    }
  // receiving and unwrapping the package
  if ( rank != 0 )
    {  
      //cout << " ++"<<rank<<"++ "<<  " receiving from the left " << endl;
      MPI::COMM_WORLD.Recv(packageW,sizeOfDataDis,MPI::DOUBLE,from,130);
      MPI::COMM_WORLD.Recv(packageW2,sizeOfDataDis,MPI::DOUBLE,from,140);
      MPI::COMM_WORLD.Recv(packagePi,sizeOfDataDis,MPI::DOUBLE,from,150);
      MPI::COMM_WORLD.Recv(packagePi2,sizeOfDataDis,MPI::DOUBLE,from,160);
      MPI::COMM_WORLD.Recv(packageU,sizeOfDataU,MPI::DOUBLE,from,170);
      MPI::COMM_WORLD.Recv(packageU2,sizeOfDataU,MPI::DOUBLE,from,180);
      MPI::COMM_WORLD.Recv(packagePi_b,sizeOfDataPi_b,MPI::DOUBLE,from,190);
      MPI::COMM_WORLD.Recv(packagePi_b2,sizeOfDataPi_b,MPI::DOUBLE,from,200);
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(i=0; i<=DATA->rk_order; i++)
		{
		  positionPi_b = i+3*(ix + ((nx+1)*iy));
		  Lneighbor[ix][iy][0].pi_b[i] = packagePi_b[positionPi_b];
		  Lneighbor[ix][iy][1].pi_b[i] = packagePi_b2[positionPi_b];
		  for(beta=0; beta<4; beta++) // dissipative part
		    {
		      for(alpha=0; alpha<5; alpha++)
			{
			  positionDis = beta+4*(i+3*(alpha+(5*(ix + ((nx+1)*iy)))));
			  Lneighbor[ix][iy][0].Wmunu[i][alpha][beta] = packageW[positionDis];
			  Lneighbor[ix][iy][1].Wmunu[i][alpha][beta] = packageW2[positionDis];
			  Lneighbor[ix][iy][0].Pimunu[i][alpha][beta] = packagePi[positionDis];
			  Lneighbor[ix][iy][1].Pimunu[i][alpha][beta] = packagePi2[positionDis];
			}
		      positionU = (i+3*(beta+(4*(ix + ((nx+1)*iy)))));
		      Lneighbor[ix][iy][0].u[i][beta] = packageU[positionU];
		      Lneighbor[ix][iy][1].u[i][beta] = packageU2[positionU];
		    }
		}
	    }
	}
      
      // cout << " done receiving from the left on rank " << rank << endl;
    }
  delete[] packageW;
  delete[] packageW2;
  delete[] packageU;
  delete[] packageU2;
  delete[] packagePi;
  delete[] packagePi2;
  delete[] packagePi_b;
  delete[] packagePi_b2;
  // if(rank==0)
    //    cout << "on rank " <<rank<< " " <<  " MPISendReceiveW done " << endl;
}//end MPISendReceive



/* %%%%%%%%%%%%%%%%%%%%%%%%%%%  Advance Local T %%%%%%%%%%%%%%%%%% */


int Advance::AdvanceLocalT(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			   Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank)
{
  // this advances the ideal part
 static Grid grid_rk;
 static double **qirk, *qi, *rhs, **w_rhs;
//  double tau_now, tau_next;
//  int trk_flag, rk_order_m1, flag;
//  int i, alpha, mu;
 static int ind=0;

 ind++;
 if(ind == 1)
 {
  qirk = util->mtx_malloc(5,4);
  qi = util->vector_malloc(5);
  rhs = util->vector_malloc(5);
  w_rhs = util->mtx_malloc(4,4); // this is new in dissipative code
  grid_rk.TJb = util->cube_malloc(DATA->rk_order,5,4);
  grid_rk.u = util->mtx_malloc(DATA->rk_order,4);
 } 
 
 FirstRKStepT(tau, DATA, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, rk_flag, qi, rhs, w_rhs, qirk, &grid_rk, size, rank);

 return 1; /* if successful */
}/* AdvanceLocalT */


/* %%%%%%%%%%%%%%%%% Advance Local W %%%%%%%%%% */


int Advance::AdvanceLocalW(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			   Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank)
{
 static Grid grid_rk;
 static double **qirk, *qi, *rhs, **w_rhs;
//  double tau_now, tau_next;
//  int trk_flag, rk_order_m1;
 int flag;
//  int i, alpha, mu;
 static int ind=0;

 ind++;
 if(ind == 1)
 {
  qirk = util->mtx_malloc(5,4);
  qi = util->vector_malloc(5);
  rhs = util->vector_malloc(5);
  w_rhs = util->mtx_malloc(4,4);
  grid_rk.TJb = util->cube_malloc(DATA->rk_order,5,4);
  grid_rk.u = util->mtx_malloc(DATA->rk_order,4);
 } 

 flag = FirstRKStepW(tau, DATA, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, rk_flag, qi, rhs, w_rhs, qirk, &grid_rk, size, rank);
 /* flag = 1 means it was successful 
    flag = -1 means that Wmunu had to be reverted */

 return flag; 
}/* AdvanceLocalW */


/* %%%%%%%%%%%%%%%%%%%%%% First steps begins here %%%%%%%%%%%%%%%%%% */



int Advance::FirstRKStepT(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, 
double *qi, double *rhs, double **w_rhs, double **qirk, Grid *grid_rk, int size, int rank)
{ 
  // this advances the ideal part
 double tau_now, tau_next, tau_rk, dwmn;
//  double tempf, p_rhs, temp_mu, temps, tempd;
 int alpha, flag;
//  int mu, nu;
//  static int ind=0;
 
 tau_now = tau;
 tau_next = tau + (DATA->delta_tau);
 
 if(rk_flag == 0) {tau_rk = tau_now;}
 else if(rk_flag > 0) {tau_rk = tau_next;}
   	else {fprintf(stderr,"rk_flag out of range.\n");exit(0);}

/* TEST */
 if(rk_flag==2) fprintf(stderr, "FirstRKStepT: rk_flag = %d\n", rk_flag);
/* TEST */

/*
 Solve partial_a T^{a mu} = -partial_a W^{a mu}
 Update T^{mu nu}
*/

 MakeDeltaQI(tau_rk, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, qi, rhs, DATA, rk_flag, size, rank);

 /* MakeDelatQI gets
    qi = q0 if rk_flag = 0 or
    qi = q0 + k1 if rk_flag = 1 */
 
 //rhs[alpha] is what MakeDeltaQI outputs. It is the spatial derivative part of partial_a T^{a mu} (including geometric terms)

 for(alpha=0; alpha<5; alpha++)
   {
     qirk[alpha][0] = qi[alpha] + rhs[alpha];
     if(!finite(qirk[alpha][0]))
       {
	 fprintf(stderr, "qirk[%d][0] = %e is a nan.\n", alpha, qirk[alpha][0]);
	 fprintf(stderr, "qi[%d] = %e\n", alpha, qi[alpha]);
	 fprintf(stderr, "rhs[%d] = %e\n", alpha, rhs[alpha]);
       }

     // now MakeWSource returns partial_a W^{a mu} (including geometric terms) 
     dwmn = diss->MakeWSource(tau_rk, alpha, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, DATA, rk_flag, size, rank);

     /* dwmn is the only one with the minus sign */
     qirk[alpha][0] -= dwmn*(DATA->delta_tau);
     
     // set baryon density back to zero if viscous correction made it non-zero
     // remove/modify if rho_b!=0 - this is only to remove the viscous correction that can make rho_b negative which we do not want.
     if(alpha==4 && qirk[alpha][0]!=0)
       qirk[alpha][0]=0.;

     /* if rk_flag > 0, we now have q0 + k1 + k2. So add q0 and multiply by 1/2 */
     if(rk_flag > 0)
      {
       qirk[alpha][0] += (grid_pt->TJb[0][alpha][0])*tau_now;
       qirk[alpha][0] *= 0.5;
      }
   }
 
 flag=
   reconst->ReconstIt(grid_rk, 0, tau_next, qirk, grid_pt,
		      grid_pt->epsilon, grid_pt->rhob, DATA, rk_flag); 
 
 if(flag==0)
   {
     reconst->ReconstError("grid_rk", 0, rk_flag+1, qi, qirk, grid_pt);
     return 0;
   }/* flag == 0 */
 else if(flag != 0) 
   {
     UpdateTJbRK(grid_rk, grid_pt, DATA, rk_flag); 
     /* TJb[rk_flag+1] is filled */
   }
 
 return 1;
}/* FirstRKStepT */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/*
   Done with T 
   Start W
*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int Advance::FirstRKStepW(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, 
double *qi, double *rhs, double **w_rhs, double **qirk, Grid *grid_rk, int size, int rank)
{ 
 double tau_now, tau_next, tempf, p_rhs, temps;
//  double temp_mu, eps_ratio, tempd;
//  double tau_rk;
//  int alpha, flag;
 int mu, nu, revert_flag;
//  static int ind=0;

 tau_now = tau;
 tau_next = tau + (DATA->delta_tau);
 
//  if(rk_flag == 0) {tau_rk = tau_now;}
//  else if(rk_flag > 0) {tau_rk = tau_next;}

/*
 Solve partial_a (u^a W^{mu nu}) = 0
 Update W^{mu nu}
*/

/* u[1][mu] now has u1[mu] */

/* calculate delta uWmunu  */
/* need to use u[0][mu], remember rk_flag = 0 here */
/* with the KT flux */
/* solve partial_tau (u^0 W^{kl}) = -partial_i (u^i W^{kl})  */
 

/* Advance uWmunu */

if(rk_flag == 0)
{
   diss->Make_uWRHS(tau_now, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, w_rhs, DATA, rk_flag, size, rank);
   for(mu=1; mu<=3; mu++)
    {
     for(nu=1; nu<=3; nu++)
      {
       tempf = (grid_pt->Wmunu[rk_flag][mu][nu])*(grid_pt->u[rk_flag][0]);

       temps = diss->Make_uWSource(tau_now, grid_pt, mu, nu, DATA, rk_flag); 
       tempf += temps*(DATA->delta_tau);
       tempf += w_rhs[mu][nu];

       grid_pt->Wmunu[rk_flag+1][mu][nu] = tempf/(grid_pt->u[rk_flag+1][0]);
      }
    }
}/* rk_flag == 0 */
else if(rk_flag > 0)
{
   diss->Make_uWRHS(tau_next, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, w_rhs, DATA, rk_flag, size, rank);
   for(mu=1; mu<=3; mu++)
    {
     for(nu=1; nu<=3; nu++)
      {
       tempf = (grid_pt->Wmunu[0][mu][nu])*(grid_pt->u[0][0]);

       temps = diss->Make_uWSource(tau_next, grid_pt, mu, nu, DATA, rk_flag); 
       tempf += temps*(DATA->delta_tau);
       tempf += w_rhs[mu][nu];

       tempf += (grid_pt->Wmunu[rk_flag][mu][nu])*(grid_pt->u[rk_flag][0]);
       tempf *= 0.5;
       
       grid_pt->Wmunu[rk_flag+1][mu][nu] = tempf/(grid_pt->u[rk_flag+1][0]);
      }
    }
}/* rk_flag > 0 */



/* calculate delta u pi */

if(rk_flag == 0)
{
/* calculate delta u^0 pi */
   diss->Make_uPRHS(tau_now, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, &p_rhs, DATA, rk_flag, size, rank);
   
   tempf = (grid_pt->pi_b[rk_flag])*(grid_pt->u[rk_flag][0]);
   
   temps = diss->Make_uPiSource(tau_now, grid_pt, DATA, rk_flag);
   tempf += temps*(DATA->delta_tau);
   tempf += p_rhs;
   
   grid_pt->pi_b[rk_flag+1] = tempf/(grid_pt->u[rk_flag+1][0]);
   
}/* rk_flag == 0 */
else if(rk_flag > 0)
{
/* calculate delta u^0 pi */
   diss->Make_uPRHS(tau_next, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, &p_rhs, DATA, rk_flag, size, rank);
   
   tempf = (grid_pt->pi_b[0])*(grid_pt->u[0][0]);
   
   temps = diss->Make_uPiSource(tau_next, grid_pt, DATA, rk_flag);
   tempf += temps*(DATA->delta_tau);
   tempf += p_rhs;
  
   tempf += (grid_pt->pi_b[1])*(grid_pt->u[0][0]);
   tempf *= 0.5;

   grid_pt->pi_b[rk_flag+1] = tempf/(grid_pt->u[rk_flag+1][0]);

}/* rk_flag > 0 */


/* update Pimunu */
   for(mu=0; mu<4; mu++)
    {
     for(nu=0; nu<4; nu++)
     {
      grid_pt->Pimunu[rk_flag+1][mu][nu]  = (grid_pt->u[rk_flag+1][mu]);
      grid_pt->Pimunu[rk_flag+1][mu][nu] *= (grid_pt->u[rk_flag+1][nu]);
      grid_pt->Pimunu[rk_flag+1][mu][nu] += DATA->gmunu[mu][nu];
      grid_pt->Pimunu[rk_flag+1][mu][nu] *= (grid_pt->pi_b[rk_flag+1]);
     }/* nu */
    }/* mu */


/* re-make Wmunu[3][3] so that Wmunu[mu][nu] is traceless */
grid_pt->Wmunu[rk_flag+1][3][3] = (2.*( grid_pt->u[rk_flag+1][1]*grid_pt->u[rk_flag+1][2]*grid_pt->Wmunu[rk_flag+1][1][2] 
                            + grid_pt->u[rk_flag+1][1]*grid_pt->u[rk_flag+1][3]*grid_pt->Wmunu[rk_flag+1][1][3]
                            + grid_pt->u[rk_flag+1][2]*grid_pt->u[rk_flag+1][3]*grid_pt->Wmunu[rk_flag+1][2][3] )
                         -( grid_pt->u[rk_flag+1][0]*grid_pt->u[rk_flag+1][0] - grid_pt->u[rk_flag+1][1]*grid_pt->u[rk_flag+1][1] )*grid_pt->Wmunu[rk_flag+1][1][1] 
                         -( grid_pt->u[rk_flag+1][0]*grid_pt->u[rk_flag+1][0] - grid_pt->u[rk_flag+1][2]*grid_pt->u[rk_flag+1][2] )*grid_pt->Wmunu[rk_flag+1][2][2]
			)/( grid_pt->u[rk_flag+1][0]*grid_pt->u[rk_flag+1][0] - grid_pt->u[rk_flag+1][3]*grid_pt->u[rk_flag+1][3] ) ;

/* make Wmunu[i][0] */
   for(mu=1; mu<=3; mu++)
    {
     tempf = 0.0;
     for(nu=1; nu<=3; nu++)
      {
       tempf +=
       (grid_pt->Wmunu[rk_flag+1][mu][nu])*(grid_pt->u[rk_flag+1][nu]); 
      }
     grid_pt->Wmunu[rk_flag+1][mu][0] = tempf/(grid_pt->u[rk_flag+1][0]);
    }
   

   for(mu=1; mu<=3; mu++)
    {
     tempf = 0.0;
     for(nu=1; nu<=3; nu++)
      {
       tempf +=
       (grid_pt->Wmunu[rk_flag+1][nu][mu])*(grid_pt->u[rk_flag+1][nu]); 
      }
     grid_pt->Wmunu[rk_flag+1][0][mu] = tempf/(grid_pt->u[rk_flag+1][0]);
    }


/* make Wmunu[0][0] */
     
     tempf = 0.0;
     for(nu=1; nu<=3; nu++)
      {
       tempf +=
       (grid_pt->Wmunu[rk_flag+1][0][nu])*(grid_pt->u[rk_flag+1][nu]); 
      }
     grid_pt->Wmunu[rk_flag+1][0][0] = tempf/(grid_pt->u[rk_flag+1][0]);
 
//   TestW(tau, DATA, grid_pt, rk_flag);

/// If the energy density of the fluid element is smaller than 0.01GeV
/// reduce Wmunu using the QuestRevert algorithm
if (grid_pt->epsilon < 0.01/hbarc){
revert_flag = 
  QuestRevert(tau, 0, grid_pt, rk_flag, DATA, size, rank);
}
  

/* if reverted, this is 1, otherwise, 0 */
   grid_pt->revert_flag = revert_flag;

   if(revert_flag == 1)
    {
     return -1;
    }
   else
    {
     return 1;
    }
}/* FirstRKStepW */


/* First steps done */


void Advance::UpdateTJbRK(Grid *grid_rk, Grid *grid_pt, InitData *DATA, int rk_flag)
{
 int trk_flag, mu, alpha;
//  double tempf;

 trk_flag = rk_flag+1;

 grid_pt->epsilon_t = grid_rk->epsilon;
 grid_pt->p_t = grid_rk->p;
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

int Advance::QuestRevert(double tau, int add, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank)
{
 int mu, nu;
 double rho_shear, rho_bulk, factor, pisize, bulksize;
//  int revert_flag=0;
 
 rho_shear = 0.;
 rho_bulk = 0.;
 factor = 10.;


 pisize = 
  (grid_pt->Wmunu[rk_flag+1][0][0]*grid_pt->Wmunu[rk_flag+1][0][0]
  +grid_pt->Wmunu[rk_flag+1][1][1]*grid_pt->Wmunu[rk_flag+1][1][1]
  +grid_pt->Wmunu[rk_flag+1][2][2]*grid_pt->Wmunu[rk_flag+1][2][2]
  +grid_pt->Wmunu[rk_flag+1][3][3]*grid_pt->Wmunu[rk_flag+1][3][3]
  -2.*(
        grid_pt->Wmunu[rk_flag+1][0][1]*grid_pt->Wmunu[rk_flag+1][0][1]
       +grid_pt->Wmunu[rk_flag+1][0][2]*grid_pt->Wmunu[rk_flag+1][0][2]
       +grid_pt->Wmunu[rk_flag+1][0][3]*grid_pt->Wmunu[rk_flag+1][0][3]
       )
  +2.*(
        grid_pt->Wmunu[rk_flag+1][1][2]*grid_pt->Wmunu[rk_flag+1][1][2]
       +grid_pt->Wmunu[rk_flag+1][1][3]*grid_pt->Wmunu[rk_flag+1][1][3]
       +grid_pt->Wmunu[rk_flag+1][2][3]*grid_pt->Wmunu[rk_flag+1][2][3]
       ));

 bulksize = 3.*grid_pt->pi_b[rk_flag+1]*grid_pt->pi_b[rk_flag+1] ;
       
 rho_shear = sqrt(  pisize/( grid_pt->epsilon*grid_pt->epsilon + 3.*grid_pt->p*grid_pt->p )  )/factor ; 

 rho_bulk  = sqrt(  bulksize/( grid_pt->epsilon*grid_pt->epsilon + 3.*grid_pt->p*grid_pt->p ) )/factor ;
 
 /// Reducing the shear stress tensor 
 if(rho_shear>0.00000000000000000001) 
   {
     
     for(mu=0; mu<4; mu++)
       {
	 for(nu=0; nu<4; nu++)
	   {   	       
	     grid_pt->Wmunu[rk_flag+1][mu][nu] = (tanh(rho_shear)/rho_shear)*grid_pt->Wmunu[rk_flag+1][mu][nu];
	   }
       }
   }
  
  /// Reducing bulk viscous pressure 
  if( rho_bulk>0.00000000000000000001) 
   {
     grid_pt->pi_b[rk_flag+1] = (tanh(rho_bulk)/rho_bulk)*grid_pt->pi_b[rk_flag+1];
     
     for(mu=0; mu<4; mu++)
       {
	 for(nu=0; nu<4; nu++)
	   {   	       
           grid_pt->Pimunu[rk_flag+1][mu][nu] = tanh(rho_bulk)/rho_bulk*grid_pt->Pimunu[rk_flag+1][mu][nu];
	   }
       }   
       
   }

return 0;


}/* QuestRevert */


void Advance::TestW(double tau, InitData *DATA, Grid *grid_pt, int rk_flag)
{
 int mu, nu;
 double trace, transv, nufac;
//  double mufac;

 trace = -grid_pt->Wmunu[rk_flag+1][0][0];
 for(mu=1; mu<=3; mu++)
  {
   trace += grid_pt->Wmunu[rk_flag+1][mu][mu];
  }
 
 transv = 0.0;
 for(mu=0; mu<=3; mu++)
  {
   for(nu=0; nu<=3; nu++)
    {
     nufac = (nu == 0 ? -1.0 : 1.0);
     transv +=
     (grid_pt->Wmunu[rk_flag+1][mu][nu])*(grid_pt->u[rk_flag+1][nu])*nufac;
    }
  }

if( (fabs(trace) > 1.0e-10 || fabs(transv) > 1.0e-10) && (grid_pt->epsilon > 0.3))
{
 fprintf(stderr, "TestW: trace = %e\n", trace);
 fprintf(stderr, "transv = %e\n", transv);
 fprintf(stderr, "at (%d, %d, %d) and tau = %e.\n", 
 grid_pt->position[1], grid_pt->position[2], grid_pt->position[3], tau);
 fprintf(stderr, "epsilon = %e\n", grid_pt->epsilon);
 fprintf(stderr, "W[0][0] = %e\n", grid_pt->Wmunu[rk_flag+1][0][0]);
 fprintf(stderr, "W[1][1] = %e\n", grid_pt->Wmunu[rk_flag+1][1][1]);
 fprintf(stderr, "W[2][2] = %e\n", grid_pt->Wmunu[rk_flag+1][2][2]);
 fprintf(stderr, "W[3][3] = %e\n", grid_pt->Wmunu[rk_flag+1][3][3]);
 fprintf(stderr, "u[0] = %e\n", grid_pt->u[rk_flag+1][0]);
 fprintf(stderr, "u[1] = %e\n", grid_pt->u[rk_flag+1][1]);
 fprintf(stderr, "u[2] = %e\n", grid_pt->u[rk_flag+1][2]);
 fprintf(stderr, "u[3] = %e\n", grid_pt->u[rk_flag+1][3]);
 fprintf(stderr, "\n");
}
 
}/* TestW */


void Advance::ProjectSpin2WS(double tau, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank)
{
 double Delta[4][4], tD[4][4], trace, norm, mfac, nfac, sum;
//  double  f[4][4];
 int m, n, a, b;

/* projecting the spin 2 part for W[rk_flag+1] */
 for(m=0; m<4; m++)
  {
   for(n=0; n<4; n++)
    {
     Delta[m][n] = 
     DATA->gmunu[m][n] + (grid_pt->u[rk_flag][m])*(grid_pt->u[rk_flag][n]);
    }}
 
 for(a=0; a<4; a++)
  {
   for(b=0; b<4; b++)
    {
     tD[a][b] = 0.0;
       for(n=0; n<4; n++)
        {
         nfac = (n==0 ? -1.0 : 1.0);
         tD[a][b] += Delta[a][n]*Delta[n][b]*nfac;
	}}}
 
 sum = 0.0; 
 for(a=0; a<4; a++)
  {
   for(b=0; b<4; b++)
    {
     sum += fabs(tD[a][b] - Delta[a][b]);
    }}
  if(sum > 1.0e-6 && grid_pt->epsilon > 0.3)
   {
    fprintf(stderr, "delta d check = %3.16e\n", sum);
   }
 
 sum = 0.0; 
 for(m=0; m<4; m++)
  {
   mfac = (m==0 ? -1.0: 1.0);
   for(n=0; n<4; n++)
    {
     nfac = (n==0 ? -1.0: 1.0);
     sum += Delta[m][n]*Delta[n][m]*mfac*nfac;
    }}
  if(fabs(sum-3.0) >  1.0e-6 && grid_pt->epsilon > 0.3)
    {
      double x=grid_pt->position[1]*DATA->delta_x-DATA->x_size/2;
      double y=grid_pt->position[2]*DATA->delta_y-DATA->y_size/2;
      double eta;
      if(size>1)
	eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2+static_cast<double>(rank)/static_cast<double>(size)*DATA->eta_size;
      else
	eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2;
      
      fprintf(stderr, "d trace check = %3.16e\n", sum);
      fprintf(stderr, "epsilon = %e\n", grid_pt->epsilon);
      fprintf(stderr, "u[0] = %e\n", grid_pt->u[rk_flag][0]);
      fprintf(stderr, "u[1] = %e\n", grid_pt->u[rk_flag][1]);
      fprintf(stderr, "u[2] = %e\n", grid_pt->u[rk_flag][2]);
      fprintf(stderr, "u[3] = %e\n", grid_pt->u[rk_flag][3]);
      fprintf(stderr, "ProjectSpin2WS at %f %f %f %f: exiting.\n",tau, x, y, eta);
      fprintf(stderr, "\n");
      
      //      exit(0);
    }

 trace = -(grid_pt->Wmunu[rk_flag+1][0][0]);
 norm = grid_pt->u[rk_flag][0]*grid_pt->u[rk_flag][0];
 for(a=1; a<=3; a++)
  {
   trace += grid_pt->Wmunu[rk_flag+1][a][a];
   norm -= grid_pt->u[rk_flag][a]*grid_pt->u[rk_flag][a];
  }


 if(fabs(norm -1.0) > 1.0e-6 && grid_pt->epsilon > 0.3)
  {
      double x=grid_pt->position[1]*DATA->delta_x-DATA->x_size/2;
      double y=grid_pt->position[2]*DATA->delta_y-DATA->y_size/2;
      double eta;
      if(size>1)
	eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2+static_cast<double>(rank)/static_cast<double>(size)*DATA->eta_size;
      else
	eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2;
   fprintf(stderr, "ProjectSpin2: norm = %3.16e\n", norm);
   fprintf(stderr, "ProjectSpin2: epsilon = %e\n", grid_pt->epsilon);
   fprintf(stderr, "ProjectSpin2WS at %f %f %f %f: exiting.\n",tau, x, y, eta);
   //   exit(0);
  }

 
 for(a=0; a<4; a++)
  {
   for(b=0; b<4; b++)
    {
     grid_pt->Wmunu[rk_flag+1][a][b] -= Delta[a][b]*trace/3.0; 
    }}
 
 trace = -(grid_pt->Wmunu[rk_flag+1][0][0]);
 for(a=1; a<=3; a++)
  {
   trace += grid_pt->Wmunu[rk_flag+1][a][a];
  }
 if(fabs(trace) > 1.0e-6 && grid_pt->epsilon > 0.3)
  {
      double x=grid_pt->position[1]*DATA->delta_x-DATA->x_size/2;
      double y=grid_pt->position[2]*DATA->delta_y-DATA->y_size/2;
      double eta;
      if(size>1)
	eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2+static_cast<double>(rank)/static_cast<double>(size)*DATA->eta_size;
      else
	eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2;
   fprintf(stderr, "ProjectSpin2WS: final trace = 3.16%e\n", trace);
   fprintf(stderr, "u0 = %e\n", grid_pt->u[rk_flag][0]);
   fprintf(stderr, "u1 = %e\n", grid_pt->u[rk_flag][1]);
   fprintf(stderr, "u2 = %e\n", grid_pt->u[rk_flag][2]);
   fprintf(stderr, "u3 = %e\n", grid_pt->u[rk_flag][3]);
   fprintf(stderr, "ProjectSpin2WS at %f %f %f %f: exiting.\n",tau, x, y, eta);
   //exit(0);
  }

}/* ProjectSpin2WS */



void Advance::ProjectSpin2W(double tau, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank)
{
 double Delta[4][4], trace, norm;
//  double  f[4][4], mfac, nfac;
 int m, n, a, b;

 for(m=0; m<4; m++)
  {
   for(n=0; n<4; n++)
    {
     Delta[m][n] = 
     DATA->gmunu[m][n] + (grid_pt->u[rk_flag+1][m])*(grid_pt->u[rk_flag+1][n]);
    }}

 trace = -(grid_pt->Wmunu[rk_flag+1][0][0]);
 norm = grid_pt->u[rk_flag+1][0]*grid_pt->u[rk_flag+1][0];
 for(a=1; a<=3; a++)
  {
   trace += grid_pt->Wmunu[rk_flag+1][a][a];
   norm -= grid_pt->u[rk_flag+1][a]*grid_pt->u[rk_flag+1][a];
  }
 if(fabs(norm -1.0) > 1.0e-6 && grid_pt->epsilon > 0.3)
  {
    double x=grid_pt->position[1]*DATA->delta_x-DATA->x_size/2;
    double y=grid_pt->position[2]*DATA->delta_y-DATA->y_size/2;
    double eta;
    if(size>1)
      eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2+static_cast<double>(rank)/static_cast<double>(size)*DATA->eta_size;
    else
      eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2;
   fprintf(stderr, "ProjectSpin2W: norm = %3.16e\n", norm);
   fprintf(stderr, "ProjectSpin2W: epsilon = %e\n", grid_pt->epsilon);
   fprintf(stderr, "ProjectSpin2W at %f %f %f %f: exiting.\n",tau, x, y, eta);
   //exit(0);
  }

/* TEST */
/*
 for(a=0; a<4; a++)
  {
   for(b=0; b<4; b++)
    {
     f[a][b] = 0.0;
     for(m=0; m<4; m++)
      {
       mfac = (m==0 ? -1.0 : 1.0);
       for(n=0; n<4; n++)
        {
         nfac = (n==0 ? -1.0 : 1.0);
	 f[a][b] +=
	 Delta[a][m]*Delta[b][n]*(grid_pt->Wmunu[rk_flag+1][m][n])*mfac*nfac;
	}
      }
     f[a][b] -= Delta[a][b]*trace/3.0; 
    }
  }
*/
 
 for(a=0; a<4; a++)
  {
   for(b=0; b<4; b++)
    {
/* TEST */
//     grid_pt->Wmunu[rk_flag+1][a][b] = (f[a][b]+f[b][a])/2.0;
     grid_pt->Wmunu[rk_flag+1][a][b] -= Delta[a][b]*trace/3.0; 
    }}
 
 trace = -(grid_pt->Wmunu[rk_flag+1][0][0]);
 for(a=1; a<=3; a++)
  {
   trace += grid_pt->Wmunu[rk_flag+1][a][a];
  }
 if(fabs(trace) > 1.0e-6 && grid_pt->epsilon > 0.3)
  {
    double x=grid_pt->position[1]*DATA->delta_x-DATA->x_size/2;
    double y=grid_pt->position[2]*DATA->delta_y-DATA->y_size/2;
    double eta;
    if(size>1)
      eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2+static_cast<double>(rank)/static_cast<double>(size)*DATA->eta_size;
    else
      eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2;
     fprintf(stderr, "ProjectSpin2W: final trace = %e\n", trace);
   fprintf(stderr, "u0 = %e\n", grid_pt->u[rk_flag+1][0]);
   fprintf(stderr, "u1 = %e\n", grid_pt->u[rk_flag+1][1]);
   fprintf(stderr, "u2 = %e\n", grid_pt->u[rk_flag+1][2]);
   fprintf(stderr, "u3 = %e\n", grid_pt->u[rk_flag+1][3]);
   fprintf(stderr, "ProjectSpin2W at %f %f %f %f: exiting.\n",tau, x, y, eta);
   //exit(0);
  }

}/* ProjectSpin2W */



void Advance::MakeDeltaQI(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			 Grid *Lneighbor2, Grid *Rneighbor2, double *qi, double *rhs, 
			 InitData *DATA, int rk_flag, int size, int rank) 
{
//   double xjet, yjet, xtrigger, ytrigger;
  static double delta[4], sumf;
//   static double tau_fac[4];
  static int alpha, i;
//   static int rk_order_m1;
//   static int nmax[4], flag;
  static double **DFmmp;
  static NbrQs NbrCells;
  static BdryCells HalfwayCells;
  static int ind=0;
//   double x=grid_pt->position[1]*DATA->delta_x-DATA->x_size/2;
//   double y=grid_pt->position[2]*DATA->delta_y-DATA->y_size/2;
//   double eta;
//   if(size>1)
//     eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2+static_cast<double>(rank)/static_cast<double>(size)*DATA->eta_size;
//   else
//     eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2;
  
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
  
//   nmax[1] = DATA->nx;
//   nmax[2] = DATA->ny;
//   nmax[3] = DATA->neta-1;

/* \partial_tau (tau Ttautau) + \partial_eta Tetatau 
           + \partial_x (tau Txtau) + \partial_y (tau Tytau) + Tetaeta = 0 */
/* \partial_tau (tau Ttaueta) + \partial_eta Teteta 
           + \partial_x (tau Txeta) + \partial_y (tau Txeta) + Tetatau = 0 */
/* \partial_tau (tau Txtau) + \partial_eta Tetax + \partial_x tau T_xx
  + \partial_y tau Tyx = 0 */

//      tau_fac[1] = tau;
//      tau_fac[2] = tau;
//      tau_fac[3] = 1.0;
 
 for(alpha=0; alpha<5; alpha++) 
  {
   qi[alpha] = grid_pt->TJb[rk_flag][alpha][0]*tau;
  }/* get qi first */

 /* implement Kurganov-Tadmor scheme */
 GetQIs(tau, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, qi, &NbrCells, rk_flag, DATA, size, rank);
 
//  flag = 
   MakeQIHalfs(qi, &NbrCells, &HalfwayCells, grid_pt, DATA);
 
//  flag = 
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

//    int jetPosition=DATA->includeJet;
//    int includeTrigger = DATA->includeTrigger;
   
   //plug in rates here:
   
   //   if ( DATA->includeJet>=1 )
   //  {
   // ...     
   
   rhs[alpha] = sumf*(DATA->delta_tau);
   
  }/* alpha */
 
 return;
}/* MakeDeltaQI */


void Advance::GetQIs
(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2,
 double *qi, NbrQs *NbrCells, int rk_flag, InitData *DATA, int size, int rank)
{
 int alpha, i;
 double tempg;
//  double T00, K00, tempf;
 int nmax[4];
 // auxiliary grids to store the received grids from the neighbor processor
//  int from;
//  int to;

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
	    
	    NbrCells->qip2[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip2[alpha][i] *= tempg;
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
	    
	    NbrCells->qim2[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim2[alpha][i] *= tempg;
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
	    
	    NbrCells->qip2[alpha][i] = grid_pt->nbr_p_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qip2[alpha][i] *= tempg;  
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
	    
	    NbrCells->qim2[alpha][i] = grid_pt->nbr_m_1[i]->TJb[rk_flag][alpha][0];
	    NbrCells->qim2[alpha][i] *= tempg;
	    
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


int Advance::MakeQIHalfs(double *qi, NbrQs *NbrCells, BdryCells *HalfwayCells, 
			 Grid *grid_pt, InitData *DATA)
{
 int alpha, direc;
//  int flag;
//  int nmax[4];
 double fphL, fphR, fmhL, fmhR;
 double gphL, gphR, gmhL, gmhR;
//  double x, y, eta;
 double tempf;

//   nmax[1] = DATA->nx;
//   nmax[2] = DATA->ny;
//   nmax[3] = DATA->neta-1;
  
  for(alpha=0; alpha<5; alpha++)
  {
   for(direc=1; direc<=3; direc++)
    {
      gphL = qi[alpha];
      fphL = 0.5*minmod->minmod_dx(NbrCells->qip1[alpha][direc], qi[alpha],
                           NbrCells->qim1[alpha][direc], DATA);
      
      gphR = NbrCells->qip1[alpha][direc];
      fphR = -0.5*minmod->minmod_dx(NbrCells->qip2[alpha][direc],
                            NbrCells->qip1[alpha][direc], qi[alpha], DATA);
      
      gmhL = NbrCells->qim1[alpha][direc];
      fmhL = 0.5*minmod->minmod_dx(qi[alpha], NbrCells->qim1[alpha][direc], 
                           NbrCells->qim2[alpha][direc], DATA);
      
      gmhR = qi[alpha];
      fmhR = -0.5*minmod->minmod_dx(NbrCells->qip1[alpha][direc], qi[alpha],
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


int Advance::ConstHalfwayCells(double tau, BdryCells *HalfwayCells, double *qi, Grid *grid_pt,
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


void Advance::MakeKTCurrents(double tau, double **DFmmp, Grid *grid_pt, 
			     BdryCells *HalfwayCells, int rk_flag)
{
 int i, alpha;
 double FiphL[5][4], FiphR[5][4], FimhL[5][4], FimhR[5][4];
 double Fiph[5][4], Fimh[5][4];
//  double delta[4];
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
 fprintf(stderr, "MakeKTCurrents: exiting.\n");
   
 //exit(0);
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


void Advance::MakeMaxSpeedAs(double tau, BdryCells *HalfwayCells, double aiph[], double aimh[], int rk_flag)
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




double Advance::MaxSpeed (double tau, int direc, Grid *grid_p, int rk_flag)
{
 double f, den, num;
 double rhob, utau, utau2, ux2, ux, p, eps, h; 
//  double  uy, ueta;
//  double deriv_p_rho, deriv_p_eps, deriv_p_h, rho_p_rho, h_p_h;
 double ut2mux2, ut, pe, rpr;
//  double rho;
 
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
	 fprintf(stderr, "MaxSpeed: exiting.\n");
   
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
   f =1.;
   //exit(0);
  }

 if(direc == 3) f /= tau;
 
 return f;
}/* MaxSpeed */


void Advance::InitNbrQs(NbrQs *NbrCells)
{
 (NbrCells->qip1) = util->mtx_malloc(5,4);
 (NbrCells->qip2) = util->mtx_malloc(5,4);
 (NbrCells->qim1) = util->mtx_malloc(5,4);
 (NbrCells->qim2) = util->mtx_malloc(5,4);
}/* InitNbrQs */


void Advance::InitTempGrids(BdryCells *HalfwayCells, int rk_order)
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


