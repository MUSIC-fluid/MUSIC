#include "util.h"
#include "data.h"
#include "grid.h"
#include "eos.h"
#include "evolve.h"
#include "advance.h"
//#include "ideal.h"

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

int Advance::AdvanceIt(double tau, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int rk_flag, int size, int rank){ 
  int ix, iy, ieta, flag;
  //cout << "0 AdvanceIt Wmunu=" << (Lneighbor[1][1][0]).Wmunu[rk_flag][1][1] << endl;
  //cout << "rank = " << rank << ", size = " << size << endl;

  if(DATA->Initial_profile != 7){
    for(ix=0; ix<=DATA->nx; ix++){
      for(iy=0; iy<=DATA->ny; iy++){
	for(ieta=0; ieta<DATA->neta; ieta++){
	  flag = AdvanceLocalT(tau, DATA, &(arena[ix][iy][ieta]), &(Lneighbor[ix][iy][0]), &(Rneighbor[ix][iy][0]), 
			       &(Lneighbor[ix][iy][1]), &(Rneighbor[ix][iy][1]), rk_flag, size, rank);
	  if(flag == 0){ 
	    return 0; 
	  }
	}/* ieta */
      }/*iy */
    }/* ix */
  }

  MPISendReceiveT(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);

  if(DATA->Initial_profile != 7){
    if(DATA->viscosity_flag == 1){
      for(ix=0; ix<=DATA->nx; ix++){
	for(iy=0; iy<=DATA->ny; iy++){
	  for(ieta=0; ieta<DATA->neta; ieta++){
	    flag = AdvanceLocalW(tau, DATA, &(arena[ix][iy][ieta]), &(Lneighbor[ix][iy][0]), &(Rneighbor[ix][iy][0]), 
				 &(Lneighbor[ix][iy][1]), &(Rneighbor[ix][iy][1]), rk_flag, size, rank);
	    if(flag==0){ 
	      return 0; 
	    }
	  }/* ieta */
	}/*iy */
      }/* ix */
    }/* if viscosity flag is set */
  }
  
  MPISendReceiveW(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);
  
  if(DATA->viscosity_flag == 1 && DATA->fluctuatingHydroFlag == 2){//linearized thermal fluctuations
    for(ix=0; ix<=DATA->nx; ix++){
      for(iy=0; iy<=DATA->ny; iy++){
	for(ieta=0; ieta<DATA->neta; ieta++){
	  flag = AdvanceLocalDeltaTAndW(tau, DATA, &(arena[ix][iy][ieta]), &(Lneighbor[ix][iy][0]), &(Rneighbor[ix][iy][0]),
					&(Lneighbor[ix][iy][1]), &(Rneighbor[ix][iy][1]), rk_flag, size, rank);
	}
	if(flag==0){ 
	  return 0; 
	}/* ieta */
      }/*iy */
    }/* ix */
  }

  if(DATA->fluctuatingHydroFlag == 1){
    MPISendReceiveXi(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);
  }
  if(DATA->fluctuatingHydroFlag == 2){
    MPISendReceiveXi2(DATA, arena, Lneighbor, Rneighbor, size, rank, rk_flag);
    //cout << "MPISendReceiveXi2 ended successfully..." << endl;
  }
  
  return 1;
}/* AdvanceIt */

void Advance::MPISendReceiveT(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag)
{
  // this sends and receives information from neighboring cells in the next processor in the eta direction 
  // and stores it in Lneighbor and Rneighbor (unless the processor is really at the edge of the total grid
  
  int ix, iy, ieta, nx, ny, neta, i, alpha, iflag;
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
  if(rank != 0){
      //      cout << " sending to the left on rank " << rank << endl;
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<5; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){ // rk_order 2 means rk_flag = 0, 1, 2
	    //if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
	    //cout << "arena[ix][iy][0].TJb[i][alpha][0]=" << arena[ix][iy][0].TJb[i][alpha][0] << endl;
	    position = (i+3*(alpha+(5*(ix + (nx*iy)))));
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
  if(rank != size-1){  
    MPI::COMM_WORLD.Recv(package,sizeOfData,MPI::DOUBLE,from,10+rk_flag);
    MPI::COMM_WORLD.Recv(package2,sizeOfData,MPI::DOUBLE,from,20+rk_flag);
    //cout << " receiving from the right on rank " << rank << endl;
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<5; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){
	    position = (i+3*(alpha+(5*(ix + (nx*iy)))));
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
  if(rank != size-1){
    //cout << " **"<<rank<<"** "<<  " sending to the right " << endl;
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<5; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){
	    position = (i+3*(alpha+(5*(ix + (nx*iy)))));
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
  if(rank != 0){  
    //cout << " ++"<<rank<<"++ "<<  " receiving from the left " << endl;
    MPI::COMM_WORLD.Recv(package,sizeOfData,MPI::DOUBLE,from,30+rk_flag);
    MPI::COMM_WORLD.Recv(package2,sizeOfData,MPI::DOUBLE,from,40+rk_flag);
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<5; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){
	    position = (i+3*(alpha+(5*(ix + (nx*iy)))));
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

void Advance::MPISendReceiveW(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag){
  // this sends and receives information from neighboring cells in the next processor in the eta direction 
  // and stores it in Lneighbor and Rneighbor (unless the processor is really at the edge of the total grid
  // Modified to send and receive Xi also. -CFY
  
  int ix, iy, ieta, nx, ny, neta, i, alpha, beta, iflag;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;
  
  int sizeOfDataDis = 5*3*(nx+1)*(ny+1)*4; // size of data package for W's and Pi's
  int sizeOfDataU = 4*3*(nx+1)*(ny+1); // size of data package for u's
  //int sizeOfDataPi_b = 3*(nx+1)*(ny+1); // size of data package for pi_b's
  
  int position, positionDis, positionU;
  //int positionPi_b;

  double *packageW; // dissipative parts
  //double *packagePi;
  double *packageW2; 
  //double *packagePi2;
  double *packageU; 
  double *packageU2;
  //double *packagePi_b;
  //double *packagePi_b2;

  packageW = new double[sizeOfDataDis];
  packageW2 = new double[sizeOfDataDis];
  //packagePi = new double[sizeOfDataDis];
  //packagePi2 = new double[sizeOfDataDis];
  packageU = new double[sizeOfDataU];
  packageU2 = new double[sizeOfDataU];
  //packagePi_b = new double[sizeOfDataPi_b];
  //packagePi_b2 = new double[sizeOfDataPi_b];

  // receive from the right / send to the left
  int from = rank+1;
  int to = rank-1;
  // packing the package to send
  if(rank != 0){
    //cout << " sending to the left on rank " << rank << endl;
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<5; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){ // rk_flag iterator
	    //if (ix ==nx/2 && iy==ny/2 && alpha == 0 && i == 0 )
	    //cout << "arena[ix][iy][0].TJb[i][alpha][0]=" << arena[ix][iy][0].TJb[i][alpha][0] << endl;
	    //positionPi_b = i+3*(ix + (nx*iy));
	    //packagePi_b[positionPi_b] = arena[ix][iy][0].pi_b[i];
	    //packagePi_b2[positionPi_b] = arena[ix][iy][1].pi_b[i];
	    if(alpha<4){
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionDis = beta+4*(i+3*(alpha+(5*(ix + (nx*iy)))));
		packageW[positionDis] = arena[ix][iy][0].Wmunu[i][alpha][beta];
		packageW2[positionDis] = arena[ix][iy][1].Wmunu[i][alpha][beta];
		//packagePi[positionDis] = arena[ix][iy][0].Pimunu[i][alpha][beta];
		//packagePi2[positionDis] = arena[ix][iy][1].Pimunu[i][alpha][beta];
		
		positionU = (i+3*(beta+(4*(ix + (nx*iy)))));
		packageU[positionU] = arena[ix][iy][0].u[i][beta];
		packageU2[positionU] = arena[ix][iy][1].u[i][beta];
	      }
	    }
	  }
	}
      }
    }
    MPI::COMM_WORLD.Send(packageW,sizeOfDataDis,MPI::DOUBLE,to,30);
    MPI::COMM_WORLD.Send(packageW2,sizeOfDataDis,MPI::DOUBLE,to,40);
    //MPI::COMM_WORLD.Send(packagePi,sizeOfDataDis,MPI::DOUBLE,to,50);
    //MPI::COMM_WORLD.Send(packagePi2,sizeOfDataDis,MPI::DOUBLE,to,60);
    MPI::COMM_WORLD.Send(packageU,sizeOfDataU,MPI::DOUBLE,to,70);
    MPI::COMM_WORLD.Send(packageU2,sizeOfDataU,MPI::DOUBLE,to,80);
    //MPI::COMM_WORLD.Send(packagePi_b,sizeOfDataPi_b,MPI::DOUBLE,to,90);
    //MPI::COMM_WORLD.Send(packagePi_b2,sizeOfDataPi_b,MPI::DOUBLE,to,100);
    //cout << " done sending to the left on rank " << rank << endl;
  }
  // receiving and unwrapping the package
  if(rank != size-1){  
    MPI::COMM_WORLD.Recv(packageW,sizeOfDataDis,MPI::DOUBLE,from,30);
    MPI::COMM_WORLD.Recv(packageW2,sizeOfDataDis,MPI::DOUBLE,from,40);
    //MPI::COMM_WORLD.Recv(packagePi,sizeOfDataDis,MPI::DOUBLE,from,50);
    //MPI::COMM_WORLD.Recv(packagePi2,sizeOfDataDis,MPI::DOUBLE,from,60);
    MPI::COMM_WORLD.Recv(packageU,sizeOfDataU,MPI::DOUBLE,from,70);
    MPI::COMM_WORLD.Recv(packageU2,sizeOfDataU,MPI::DOUBLE,from,80);
    //MPI::COMM_WORLD.Recv(packagePi_b,sizeOfDataPi_b,MPI::DOUBLE,from,90);
    //MPI::COMM_WORLD.Recv(packagePi_b2,sizeOfDataPi_b,MPI::DOUBLE,from,100);
    //cout << " receiving from the right on rank " << rank << endl;
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<5; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){
	    //positionPi_b = i+3*(ix + (nx*iy));
	    //Rneighbor[ix][iy][0].pi_b[i] = packagePi_b[positionPi_b];
	    //Rneighbor[ix][iy][1].pi_b[i] = packagePi_b2[positionPi_b];
	    if(alpha<4){
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionDis = beta+4*(i+3*(alpha+(5*(ix + (nx*iy)))));
		Rneighbor[ix][iy][0].Wmunu[i][alpha][beta] = packageW[positionDis];
		Rneighbor[ix][iy][1].Wmunu[i][alpha][beta] = packageW2[positionDis];
		//Rneighbor[ix][iy][0].Pimunu[i][alpha][beta] = packagePi[positionDis];
		//Rneighbor[ix][iy][1].Pimunu[i][alpha][beta] = packagePi2[positionDis];
		//cout << "received from right Wmunu=" <<  Rneighbor[ix][iy][0].Wmunu[i][alpha][beta] << endl;
		positionU = (i+3*(beta+(4*(ix + (nx*iy)))));
		Rneighbor[ix][iy][0].u[i][beta] = packageU[positionU];
		Rneighbor[ix][iy][1].u[i][beta] = packageU2[positionU];
	      }
	    }
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
  if(rank != size-1){
    //cout << " **"<<rank<<"** "<<  " sending to the right " << endl;
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<5; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){
	    //positionPi_b = i+3*(ix + (nx*iy));
	    //packagePi_b[positionPi_b] = arena[ix][iy][neta-1].pi_b[i];
	    //packagePi_b2[positionPi_b] = arena[ix][iy][neta-2].pi_b[i];
	    if(alpha<4){
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionDis = beta+4*(i+3*(alpha+(5*(ix + (nx*iy)))));
		packageW[positionDis] = arena[ix][iy][neta-1].Wmunu[i][alpha][beta];
		packageW2[positionDis] = arena[ix][iy][neta-2].Wmunu[i][alpha][beta];
		//packagePi[positionDis] = arena[ix][iy][neta-1].Pimunu[i][alpha][beta];
		//packagePi2[positionDis] = arena[ix][iy][neta-2].Pimunu[i][alpha][beta];
		
		positionU = (i+3*(beta+(4*(ix + (nx*iy)))));
		packageU[positionU] = arena[ix][iy][neta-1].u[i][beta];
		packageU2[positionU] = arena[ix][iy][neta-2].u[i][beta];
	      }
	    }
	  }
	}
      }
    }
    //cout << "sending to the right Wmunu=" <<  arena[1][1][neta-1].Wmunu[0][1][1] << endl;
    
    MPI::COMM_WORLD.Send(packageW,sizeOfDataDis,MPI::DOUBLE,to,130);
    MPI::COMM_WORLD.Send(packageW2,sizeOfDataDis,MPI::DOUBLE,to,140);
    //MPI::COMM_WORLD.Send(packagePi,sizeOfDataDis,MPI::DOUBLE,to,150);
    //MPI::COMM_WORLD.Send(packagePi2,sizeOfDataDis,MPI::DOUBLE,to,160);
    MPI::COMM_WORLD.Send(packageU,sizeOfDataU,MPI::DOUBLE,to,170);
    MPI::COMM_WORLD.Send(packageU2,sizeOfDataU,MPI::DOUBLE,to,180);
    //MPI::COMM_WORLD.Send(packagePi_b,sizeOfDataPi_b,MPI::DOUBLE,to,190);
    //MPI::COMM_WORLD.Send(packagePi_b2,sizeOfDataPi_b,MPI::DOUBLE,to,200);
    //cout << " **"<<rank<<"** "<<  " done sending to the right " << endl;
  }
  // receiving and unwrapping the package
  if(rank != 0){  
    //cout << " ++"<<rank<<"++ "<<  " receiving from the left " << endl;
    MPI::COMM_WORLD.Recv(packageW,sizeOfDataDis,MPI::DOUBLE,from,130);
    MPI::COMM_WORLD.Recv(packageW2,sizeOfDataDis,MPI::DOUBLE,from,140);
    //MPI::COMM_WORLD.Recv(packagePi,sizeOfDataDis,MPI::DOUBLE,from,150);
    //MPI::COMM_WORLD.Recv(packagePi2,sizeOfDataDis,MPI::DOUBLE,from,160);
    MPI::COMM_WORLD.Recv(packageU,sizeOfDataU,MPI::DOUBLE,from,170);
    MPI::COMM_WORLD.Recv(packageU2,sizeOfDataU,MPI::DOUBLE,from,180);
    //MPI::COMM_WORLD.Recv(packagePi_b,sizeOfDataPi_b,MPI::DOUBLE,from,190);
    //MPI::COMM_WORLD.Recv(packagePi_b2,sizeOfDataPi_b,MPI::DOUBLE,from,200);
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<5; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){
	    //positionPi_b = i+3*(ix + (nx*iy));
	    //Lneighbor[ix][iy][0].pi_b[i] = packagePi_b[positionPi_b];
	    //Lneighbor[ix][iy][1].pi_b[i] = packagePi_b2[positionPi_b];
	    if(alpha<4){
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionDis = beta+4*(i+3*(alpha+(5*(ix + (nx*iy)))));
		Lneighbor[ix][iy][0].Wmunu[i][alpha][beta] = packageW[positionDis];
		Lneighbor[ix][iy][1].Wmunu[i][alpha][beta] = packageW2[positionDis];
		//Lneighbor[ix][iy][0].Pimunu[i][alpha][beta] = packagePi[positionDis];
		//Lneighbor[ix][iy][1].Pimunu[i][alpha][beta] = packagePi2[positionDis];
		
		positionU = (i+3*(beta+(4*(ix + (nx*iy)))));
		Lneighbor[ix][iy][0].u[i][beta] = packageU[positionU];
		Lneighbor[ix][iy][1].u[i][beta] = packageU2[positionU];
	      }
	    }
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
  //delete[] packagePi;
  //delete[] packagePi2;
  //delete[] packagePi_b;
  //delete[] packagePi_b2;
  // if(rank==0)
  //    cout << "on rank " <<rank<< " " <<  " MPISendReceiveW done " << endl;
}//end MPISendReceive

void Advance::MPISendReceiveXi(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag){
  // this sends and receives information from neighboring cells in the next processor in the eta direction 
  // and stores it in Lneighbor and Rneighbor (unless the processor is really at the edge of the total grid
  // Modified to send and receive Xi also. -CFY
  
  int ix, iy, ieta, nx, ny, neta, i, alpha, beta, iflag;
  nx = DATA->nx; ny = DATA->ny; neta = DATA->neta;
  
  int sizeOfDataXi = 4*4*3*(nx+1)*(ny+1); // Size of data for Xi. -CFY
  
  int position, positionXi;
  
  double *packageXi;
  
  packageXi = new double[sizeOfDataXi];
  
  // receive from the right / send to the left
  int from = rank+1;
  int to = rank-1;
  // packing the package to send
  if(rank != 0){
    //cout << " sending to the left on rank " << rank << endl;
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<4; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){ // rk_flag iterator
	    for(beta=0; beta<4; beta++){ // dissipative part
	      positionXi = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy))) ;
	      packageXi[positionXi] = arena[ix][iy][0].Xi[i][alpha][beta];
	    }
	  }
	}
      }
    }
    MPI::COMM_WORLD.Send(packageXi, sizeOfDataXi, MPI::DOUBLE, to, 10);
  }
  // receiving and unwrapping the package
  if(rank != size-1 ){  
    MPI::COMM_WORLD.Recv(packageXi, sizeOfDataXi, MPI::DOUBLE, from, 10);
    //cout << " receiving from the right on rank " << rank << endl;
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<4; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){
	    for(beta=0; beta<4; beta++){ // dissipative part
	      positionXi = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy)));
	      Rneighbor[ix][iy][0].Xi[i][alpha][beta] = packageXi[positionXi];
	    }
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
  if(rank != size-1){
    //cout << " **"<<rank<<"** "<<  " sending to the right " << endl;
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<4; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){
	    for(beta=0; beta<4; beta++){ // dissipative part
	      positionXi = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy)));
	      packageXi[positionXi] = arena[ix][iy][neta-1].Xi[i][alpha][beta];
	    }
	  }
	}
      }
    }
    //cout << "sending to the right Wmunu=" <<  arena[1][1][neta-1].Wmunu[0][1][1] << endl;
    
    MPI::COMM_WORLD.Send(packageXi, sizeOfDataXi, MPI::DOUBLE, to, 20);
    //cout << " **"<<rank<<"** "<<  " done sending to the right " << endl;
  }
  // receiving and unwrapping the package
  if(rank != 0){  
    //cout << " ++"<<rank<<"++ "<<  " receiving from the left " << endl;
    MPI::COMM_WORLD.Recv(packageXi, sizeOfDataXi, MPI::DOUBLE, from, 20);
    for(ix=0; ix<=nx; ix++){
      for(iy=0; iy<=ny; iy++){
	for(alpha=0; alpha<4; alpha++){
	  for(i=0; i<=DATA->rk_order; i++){
	    for(beta=0; beta<4; beta++){ // dissipative part
	      positionXi = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy)));
	      Lneighbor[ix][iy][0].Xi[i][alpha][beta] = packageXi[positionXi];
	    }
	  }
	}
      }
    }
    // cout << " done receiving from the left on rank " << rank << endl;
  }
  delete[] packageXi;
  // if(rank==0)
  //    cout << "on rank " <<rank<< " " <<  " MPISendReceiveW done " << endl;
}//end MPISendReceive

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%  Advance Local T %%%%%%%%%%%%%%%%%% */

void Advance::MPISendReceiveXi2(InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank, int rk_flag){
  // this sends and receives information from neighboring cells in the next processor in the eta direction 
  // and stores it in Lneighbor and Rneighbor (unless the processor is really at the edge of the total grid
  // Modified to send and receive Xi also. -CFY

  if(size != 1){

    int ix, iy, ieta, nx, ny, neta, i, alpha, beta, iflag;
    nx = DATA->nx; ny = DATA->ny; neta = DATA->neta;

    int sizeOfDataP = 3*(nx+1)*(ny+1); // Size of data for deltaP. -CFY
    int sizeOfDataU = 4*3*(nx+1)*(ny+1); // Size of data for deltaU. -CFY
    int sizeOfDataW = 4*4*3*(nx+1)*(ny+1); // Size of data for deltaT and deltaW. -CFY

    int position, positionP, positionU, positionW;

    double *packageP, *packageU, *packageW, *packageT, *packageXi;

    packageP = new double[sizeOfDataP];
    packageU = new double[sizeOfDataU];
    packageW = new double[sizeOfDataW];
    packageT = new double[sizeOfDataW];
    //packageXi = new double[sizeOfDataW];

    //// receive from the right / send to the left
    //int from = rank+1;
    //int to = rank-1;
    int from, to;
    if(rank==0){
      to = size-1;
    }
    else{
      to = rank-1;
    }
    if(rank==size-1){
      from = 0;
    }
    else{
      from = rank+1;
    }

    //// packing the package to send
    if(rank != 0){
      //cout << " sending to the left on rank " << rank << " to rank " << to << endl;
      //cout << "rank = " << rank << ", to = " << to << ", from = " << from << endl;
      //cout << "size = " << size << ", rank = " << rank << endl;
      //cout << MPI::COMM_WORLD.Get_size() << endl;
      for(i=0; i<=DATA->rk_order; i++){ // rk_flag iterator
	for(ix=0; ix<=nx; ix++){
	  for(iy=0; iy<=ny; iy++){
	    positionP = i + 3*(ix + nx*iy);
	    packageP[positionP] = arena[ix][iy][0].deltaP[i];
	    for(alpha=0; alpha<4; alpha++){
	      positionU = i + 3*(alpha + 4*(ix + nx*iy));
	      packageU[positionU] = arena[ix][iy][0].deltaU[i][alpha];
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionW = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy))) ;
		packageW[positionW] = arena[ix][iy][0].deltaW[i][alpha][beta];
		packageT[positionW] = arena[ix][iy][0].deltaT[i][alpha][beta];
		//packageXi[positionW] = arena[ix][iy][0].Xi[i][alpha][beta];
	      }
	    }
	  }
	}
      }
      //cout << "about to send on rank " << rank << endl;
      MPI::COMM_WORLD.Send(packageP, sizeOfDataP, MPI::DOUBLE, to, 10);
      MPI::COMM_WORLD.Send(packageU, sizeOfDataU, MPI::DOUBLE, to, 20);
      MPI::COMM_WORLD.Send(packageW, sizeOfDataW, MPI::DOUBLE, to, 30);
      MPI::COMM_WORLD.Send(packageT, sizeOfDataW, MPI::DOUBLE, to, 40);
      //MPI::COMM_WORLD.Send(packageXi, sizeOfDataW, MPI::DOUBLE, to, 50);
      //cerr << "finished sending to the left on rank " << rank << endl;
    }
    // receiving and unwrapping the package
    if(rank != size-1 ){
      //cout << " receiving from the right on rank " << rank << endl;
      MPI::COMM_WORLD.Recv(packageP, sizeOfDataW, MPI::DOUBLE, from, 10);
      MPI::COMM_WORLD.Recv(packageU, sizeOfDataW, MPI::DOUBLE, from, 20);
      MPI::COMM_WORLD.Recv(packageW, sizeOfDataW, MPI::DOUBLE, from, 30);
      MPI::COMM_WORLD.Recv(packageT, sizeOfDataW, MPI::DOUBLE, from, 40);
      //MPI::COMM_WORLD.Recv(packageXi, sizeOfDataW, MPI::DOUBLE, from, 50);
      for(i=0; i<=DATA->rk_order; i++){
	for(ix=0; ix<=nx; ix++){
	  for(iy=0; iy<=ny; iy++){
	    positionP = i + 3*(ix + nx*iy);
	    Rneighbor[ix][iy][0].deltaP[i] = packageP[positionP];
	    for(alpha=0; alpha<4; alpha++){
	      positionU = i + 3*(alpha + 4*(ix + nx*iy));
	      Rneighbor[ix][iy][0].deltaU[i][alpha] = packageU[positionU];
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionW = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy)));
		Rneighbor[ix][iy][0].deltaW[i][alpha][beta] = packageW[positionW];
		Rneighbor[ix][iy][0].deltaT[i][alpha][beta] = packageT[positionW];
		//Rneighbor[ix][iy][0].Xi[i][alpha][beta] = packageXi[positionW];
	      }
	    }
	  }
	}
      }
      //cerr << " done receiving from the right on rank " << rank << endl;
    }

    //// packing the package to send
    if(rank == 0){
      //cout << " sending to the left on rank " << rank << " to rank " << to << endl;
      //cout << "rank = " << rank << ", to = " << to << ", from = " << from << endl;
      //cout << "size = " << size << ", rank = " << rank << endl;
      //cout << MPI::COMM_WORLD.Get_size() << endl;
      for(i=0; i<=DATA->rk_order; i++){ // rk_flag iterator
	for(ix=0; ix<=nx; ix++){
	  for(iy=0; iy<=ny; iy++){
	    positionP = i + 3*(ix + nx*iy);
	    packageP[positionP] = arena[ix][iy][0].deltaP[i];
	    for(alpha=0; alpha<4; alpha++){
	      positionU = i + 3*(alpha + 4*(ix + nx*iy));
	      packageU[positionU] = arena[ix][iy][0].deltaU[i][alpha];
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionW = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy))) ;
		packageW[positionW] = arena[ix][iy][0].deltaW[i][alpha][beta];
		packageT[positionW] = arena[ix][iy][0].deltaT[i][alpha][beta];
		//packageXi[positionW] = arena[ix][iy][0].Xi[i][alpha][beta];
	      }
	    }
	  }
	}
      }
      //cout << "about to send on rank " << rank << endl;
      MPI::COMM_WORLD.Send(packageP, sizeOfDataP, MPI::DOUBLE, to, 10);
      MPI::COMM_WORLD.Send(packageU, sizeOfDataU, MPI::DOUBLE, to, 20);
      MPI::COMM_WORLD.Send(packageW, sizeOfDataW, MPI::DOUBLE, to, 30);
      MPI::COMM_WORLD.Send(packageT, sizeOfDataW, MPI::DOUBLE, to, 40);
      //MPI::COMM_WORLD.Send(packageXi, sizeOfDataW, MPI::DOUBLE, to, 50);
      //cerr << "finished sending to the left on rank " << rank << endl;
    }
    // receiving and unwrapping the package
    if(rank == size-1 ){
      //cout << " receiving from the right on rank " << rank << endl;
      MPI::COMM_WORLD.Recv(packageP, sizeOfDataW, MPI::DOUBLE, from, 10);
      MPI::COMM_WORLD.Recv(packageU, sizeOfDataW, MPI::DOUBLE, from, 20);
      MPI::COMM_WORLD.Recv(packageW, sizeOfDataW, MPI::DOUBLE, from, 30);
      MPI::COMM_WORLD.Recv(packageT, sizeOfDataW, MPI::DOUBLE, from, 40);
      //MPI::COMM_WORLD.Recv(packageXi, sizeOfDataW, MPI::DOUBLE, from, 50);
      for(i=0; i<=DATA->rk_order; i++){
	for(ix=0; ix<=nx; ix++){
	  for(iy=0; iy<=ny; iy++){
	    positionP = i + 3*(ix + nx*iy);
	    Rneighbor[ix][iy][0].deltaP[i] = packageP[positionP];
	    for(alpha=0; alpha<4; alpha++){
	      positionU = i + 3*(alpha + 4*(ix + nx*iy));
	      Rneighbor[ix][iy][0].deltaU[i][alpha] = packageU[positionU];
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionW = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy)));
		Rneighbor[ix][iy][0].deltaW[i][alpha][beta] = packageW[positionW];
		Rneighbor[ix][iy][0].deltaT[i][alpha][beta] = packageT[positionW];
		//Rneighbor[ix][iy][0].Xi[i][alpha][beta] = packageXi[positionW];
	      }
	    }
	  }
	}
      }
      //cerr << " done receiving from the right on rank " << rank << endl;
    }

    //// receive from the left / send to the right
    //from = rank-1;
    //to = rank+1;
    if(rank==size-1){
      to = 0;
    }
    else{
      to = rank+1;
    }
    if(rank==0){
      from = size-1;
    }
    else{
      from = rank-1;
    }

    // packing the package to send
    if(rank != size-1){
      //cout << " **"<<rank<<"** "<<  " sending to the right " << endl;
      for(i=0; i<=DATA->rk_order; i++){
	for(ix=0; ix<=nx; ix++){
	  for(iy=0; iy<=ny; iy++){
	    positionP = i + 3*(ix + nx*iy);
	    packageP[positionP] = arena[ix][iy][neta-1].deltaP[i];
	    for(alpha=0; alpha<4; alpha++){
	      positionU = i + 3*(alpha + 4*(ix + nx*iy));
	      packageU[positionU] = arena[ix][iy][neta-1].deltaU[i][alpha];
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionW = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy)));
		packageW[positionW] = arena[ix][iy][neta-1].deltaW[i][alpha][beta];
		packageT[positionW] = arena[ix][iy][neta-1].deltaT[i][alpha][beta];
		//packageXi[positionW] = arena[ix][iy][neta-1].Xi[i][alpha][beta];
	      }
	    }
	  }
	}
      }
      //cerr << "sending to the right Wmunu=" <<  arena[1][1][neta-1].Wmunu[0][1][1] << endl;

      MPI::COMM_WORLD.Send(packageP, sizeOfDataP, MPI::DOUBLE, to, 110);
      MPI::COMM_WORLD.Send(packageU, sizeOfDataU, MPI::DOUBLE, to, 120);
      MPI::COMM_WORLD.Send(packageW, sizeOfDataW, MPI::DOUBLE, to, 130);
      MPI::COMM_WORLD.Send(packageT, sizeOfDataW, MPI::DOUBLE, to, 140);
      //MPI::COMM_WORLD.Send(packageXi, sizeOfDataW, MPI::DOUBLE, to, 150);
      //cerr << " **"<<rank<<"** "<<  " done sending to the right " << endl;
    }
    // receiving and unwrapping the package
    if(rank != 0){  
      //cout << " ++"<<rank<<"++ "<<  " receiving from the left " << endl;
      MPI::COMM_WORLD.Recv(packageP, sizeOfDataP, MPI::DOUBLE, from, 110);
      MPI::COMM_WORLD.Recv(packageU, sizeOfDataU, MPI::DOUBLE, from, 120);
      MPI::COMM_WORLD.Recv(packageW, sizeOfDataW, MPI::DOUBLE, from, 130);
      MPI::COMM_WORLD.Recv(packageT, sizeOfDataW, MPI::DOUBLE, from, 140);
      //MPI::COMM_WORLD.Recv(packageXi, sizeOfDataW, MPI::DOUBLE, from, 150);
      for(i=0; i<=DATA->rk_order; i++){
	for(ix=0; ix<=nx; ix++){
	  for(iy=0; iy<=ny; iy++){
	    positionP = i + 3*(ix + nx*iy);
	    Lneighbor[ix][iy][0].deltaP[i] = packageP[positionP];
	    for(alpha=0; alpha<4; alpha++){
	      positionU = i + 3*(alpha + 4*(ix + nx*iy));
	      Lneighbor[ix][iy][0].deltaU[i][alpha] = packageU[positionU];
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionW = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy)));
		Lneighbor[ix][iy][0].deltaW[i][alpha][beta] = packageW[positionW];
		Lneighbor[ix][iy][0].deltaT[i][alpha][beta] = packageT[positionW];
		//Lneighbor[ix][iy][0].Xi[i][alpha][beta] = packageXi[positionW];
	      }
	    }
	  }
	}
      }
      //cerr << " done receiving from the left on rank " << rank << endl;
    }

    // packing the package to send
    if(rank == size-1){
      //cout << " **"<<rank<<"** "<<  " sending to the right " << endl;
      for(i=0; i<=DATA->rk_order; i++){
	for(ix=0; ix<=nx; ix++){
	  for(iy=0; iy<=ny; iy++){
	    positionP = i + 3*(ix + nx*iy);
	    packageP[positionP] = arena[ix][iy][neta-1].deltaP[i];
	    for(alpha=0; alpha<4; alpha++){
	      positionU = i + 3*(alpha + 4*(ix + nx*iy));
	      packageU[positionU] = arena[ix][iy][neta-1].deltaU[i][alpha];
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionW = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy)));
		packageW[positionW] = arena[ix][iy][neta-1].deltaW[i][alpha][beta];
		packageT[positionW] = arena[ix][iy][neta-1].deltaT[i][alpha][beta];
		//packageXi[positionW] = arena[ix][iy][neta-1].Xi[i][alpha][beta];
	      }
	    }
	  }
	}
      }
      //cout << "sending to the right Wmunu=" <<  arena[1][1][neta-1].Wmunu[0][1][1] << endl;

      MPI::COMM_WORLD.Send(packageP, sizeOfDataP, MPI::DOUBLE, to, 110);
      MPI::COMM_WORLD.Send(packageU, sizeOfDataU, MPI::DOUBLE, to, 120);
      MPI::COMM_WORLD.Send(packageW, sizeOfDataW, MPI::DOUBLE, to, 130);
      MPI::COMM_WORLD.Send(packageT, sizeOfDataW, MPI::DOUBLE, to, 140);
      //MPI::COMM_WORLD.Send(packageXi, sizeOfDataW, MPI::DOUBLE, to, 150);
      //cerr << " **"<<rank<<"** "<<  " done sending to the right " << endl;
    }
    // receiving and unwrapping the package
    if(rank == 0){  
      //cout << " ++"<<rank<<"++ "<<  " receiving from the left " << endl;
      MPI::COMM_WORLD.Recv(packageP, sizeOfDataP, MPI::DOUBLE, from, 110);
      MPI::COMM_WORLD.Recv(packageU, sizeOfDataU, MPI::DOUBLE, from, 120);
      MPI::COMM_WORLD.Recv(packageW, sizeOfDataW, MPI::DOUBLE, from, 130);
      MPI::COMM_WORLD.Recv(packageT, sizeOfDataW, MPI::DOUBLE, from, 140);
      //MPI::COMM_WORLD.Recv(packageXi, sizeOfDataW, MPI::DOUBLE, from, 150);
      for(i=0; i<=DATA->rk_order; i++){
	for(ix=0; ix<=nx; ix++){
	  for(iy=0; iy<=ny; iy++){
	    positionP = i + 3*(ix + nx*iy);
	    Lneighbor[ix][iy][0].deltaP[i] = packageP[positionP];
	    for(alpha=0; alpha<4; alpha++){
	      positionU = i + 3*(alpha + 4*(ix + nx*iy));
	      Lneighbor[ix][iy][0].deltaU[i][alpha] = packageU[positionU];
	      for(beta=0; beta<4; beta++){ // dissipative part
		positionW = beta + 4*(i + 3*(alpha + 4*(ix + nx*iy)));
		Lneighbor[ix][iy][0].deltaW[i][alpha][beta] = packageW[positionW];
		Lneighbor[ix][iy][0].deltaT[i][alpha][beta] = packageT[positionW];
		//Lneighbor[ix][iy][0].Xi[i][alpha][beta] = packageXi[positionW];
	      }
	    }
	  }
	}
      }
      //cerr << " done receiving from the left on rank " << rank << endl;
    }

    delete[] packageP;
    delete[] packageU;
    delete[] packageW;
    delete[] packageT;
    //delete[] packageXi;
  }
  // if(rank==0)
  //    cout << "on rank " <<rank<< " " <<  " MPISendReceiveW done " << endl;
}//end MPISendReceive

int Advance::AdvanceLocalT(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			   Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank)
{
  // this advances the ideal part
 static Grid grid_rk;
 static double **qirk, *qi, *rhs, **w_rhs;
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
  w_rhs = util->mtx_malloc(4,4);
  grid_rk.TJb = util->cube_malloc(DATA->rk_order,5,4);
  grid_rk.u = util->mtx_malloc(DATA->rk_order,4);
 } 

 flag = FirstRKStepW(tau, DATA, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, rk_flag, qi, rhs, w_rhs, qirk, &grid_rk, size, rank);
 /* flag = 1 means it was successful 
    flag = -1 means that Wmunu had to be reverted */

 return flag; 
}/* AdvanceLocalW */

int Advance::AdvanceLocalDeltaTAndW(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
				       Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank){
  if(DATA->rk_order != 2){
    cerr << "From AdvanceLocalDeltaTAndW: rk_order != 2 not allowed! Exiting..." << endl;
    exit(1);
  }

  double delta[4];
  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = DATA->delta_eta;

  int nmax[4];
  nmax[0] = 0; //Should not be used
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;

  //if(grid_pt->position[1] == nmax[1]/2 && grid_pt->position[2] == nmax[2]/2 && grid_pt->position[3] == nmax[3]/2 && rank == size/2-1){
  //cerr << "T0 = {" << grid_pt->deltaT[rk_flag][0][0] << ", " << grid_pt->deltaT[rk_flag][0][1] << ", " << grid_pt->deltaT[rk_flag][0][2] << ", " << grid_pt->deltaT[rk_flag][0][3] 
  //	 << "}" << endl;
  //cerr << "deltaW = {" << grid_pt->deltaW[rk_flag][1][0] << ", " << grid_pt->deltaW[rk_flag][1][1] << ", " << grid_pt->deltaW[rk_flag][1][2] << ", " << grid_pt->deltaW[rk_flag][1][3] 
  //	 << "}" << endl;
  //cerr << "dW = {" << getDMuDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0)
  //	 << ", " << getDMuDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 1)
  //	 << ", " << getDMuDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 2)
  //	 << ", " << getDMuDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 3) << "}" << endl;
  //cerr << "deltaW[0][3] = " << grid_pt->deltaW[rk_flag][0][3] << endl;
  //cerr << "deltaW[3][3] = " << grid_pt->deltaW[rk_flag][3][3] << endl;
  //}

  //Solve for the change in deltaW.

  //For the given grid cell, determine the noise tensor:
  double entrDens = eos->s_func(grid_pt->epsilon_t, grid_pt->p_t, grid_pt->rhob_t, DATA);
  double shearV;
  shearV = (DATA->shear_to_s)*entrDens;

  double ZCalc[4][4];
  if(grid_pt->S[0][0] == 0. && grid_pt->S[0][1] == 0. && grid_pt->S[0][2] == 0.
     && grid_pt->S[1][0] == 0. && grid_pt->S[1][1] == 0. && grid_pt->S[1][2] == 0.
     && grid_pt->S[2][0] == 0. && grid_pt->S[2][1] == 0. && grid_pt->S[2][2] == 0.){
    for(int zi=0; zi<4; zi++){
      for(int zj=0; zj<4; zj++){
	ZCalc[zi][zj] = 0.;
      }
    }
  }
  else{
    //Determine the 3x3 Z tensor in the fluid rest frame:
    double ZRest[3][3];
    double multFactor = sqrt(2.*grid_pt->T_t*shearV);
    for(int mu=0; mu<3; mu++){
      for(int alpha=0; alpha<3; alpha++){
	ZRest[mu][alpha] = multFactor*grid_pt->S[mu][alpha];
      }
    }
    
    //Now, boost into the grid's frame:
    double uBoost[4];
    for(int mu=0; mu<4; mu++){
      uBoost[mu] = grid_pt->u[rk_flag][mu];
    }
    util->LorentzBoost3x3Tensor(ZRest, uBoost, ZCalc);
  }
  
  //if(grid_pt->position[1] == nmax[1]/2 && grid_pt->position[2] == nmax[2]/2 && grid_pt->position[3] == nmax[3]/2 && rank == size/2-1){
  //cerr << "ZCalc = {" << ZCalc[1][0] << ", " << ZCalc[1][1] << ", " << ZCalc[1][2] << ", " << ZCalc[1][3] << "}" << endl;
  //}


  //To make the code run efficiently, these are computed outside of the upcoming loops:
  double dMuUMu = getDMuUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank);
  double dMuDeltaUMu = getDMuDeltaUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank);
  //double dDeltaUW[4], deltaWDU[4], deltaUDU[4];
  double DUW[4], deltaUDUW[4], DDeltaUW[4], DUDeltaW[4];// DUXi[4];
  for(int ii=0; ii<4; ii++){
    DUW[ii] = getDUWMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii);
    deltaUDUW[ii] = getDeltaUDUWMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii);
    DDeltaUW[ii] = getDDeltaUWMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii);
    DUDeltaW[ii] = getDeltaWDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii);
    //DUXi[ii] = getXiDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii);
  }
  //deltaUWDU = getDeltaUWDU(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank);

  //The coefficient \delta \eta:
  double deltaEta;
  if(grid_pt->deltaP[rk_flag] == 0.){
    deltaEta = 0.;
  }
  else{
    deltaEta = (shearV/eos->get_dpOverde2(grid_pt->epsilon_t, grid_pt->rhob_t))
      *grid_pt->deltaP[rk_flag]/(grid_pt->p_t + grid_pt->epsilon_t);
  }

  //Finally, tauPi:
  double s_den = eos->s_func(grid_pt->epsilon_t, grid_pt->p_t, grid_pt->rhob_t, DATA);
  double shear = (DATA->shear_to_s)*s_den;
  double tauPi = 3.0*shear/(grid_pt->epsilon_t + grid_pt->p_t);
  tauPi = max(tauPi, DATA->tau_pi);
  if(!finite(tauPi)) tauPi = DATA->tau_pi;

  //cout << "OK starting evolution of viscous part..." << endl;
  for(int ii=0; ii<4; ii++){
    for(int ij=0; ij<4; ij++){
      if(ii > ij){
	grid_pt->deltaW[rk_flag+1][ii][ij] = grid_pt->deltaW[rk_flag+1][ij][ii];
	//grid_pt->Xi[rk_flag+1][ii][ij] = grid_pt->Xi[rk_flag+1][ij][ii];
      }
      else{
	if(rk_flag==0){
	  //grid_pt->deltaW[rk_flag+1][ii][ij] = grid_pt->deltaW[rk_flag][ii][ij]
	  //-DATA->delta_tau*( grid_pt->deltaW[rk_flag][ii][ij]
	  //		     -shearV*getDeltaS(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)
	  //		     -ZCalc[ii][ij] )/(DATA->tau_pi*grid_pt->u[rk_flag][0]);
	  //If T < T_min, the fluctuations are decoupled from the averaged quantities:
	  if(grid_pt->T_t < DATA->fluctuatingTMin || grid_pt->position[1]<20 || grid_pt->position[1]>=nmax[1]-20
	     || grid_pt->position[2]<20 || grid_pt->position[2]>=nmax[2]-20){
	    //grid_pt->deltaW[rk_flag+1][ii][ij] = grid_pt->deltaW[rk_flag][ii][ij]
	    //+DATA->delta_tau*
	    //( -getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    //				   dMuUMu, dMuDeltaUMu, ii, ij) );
	    /*
	    grid_pt->Xi[rk_flag+1][ii][ij] = grid_pt->Xi[rk_flag][ii][ij]
	      +(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiXiMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij) 

		+getProjectorSourceXi(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUXi)

	    	-getDeltaXiDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    				   dMuUMu, dMuDeltaUMu, ii, ij) );

	    grid_pt->deltaW[rk_flag+1][ii][ij] = grid_pt->deltaW[rk_flag][ii][ij]
	      +(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij) 

		+getProjectorSource(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUW, deltaUDUW, DDeltaUW, DUDeltaW)

	    	-getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    				   dMuUMu, dMuDeltaUMu, ii, ij) );
	    */
	    grid_pt->deltaW[rk_flag+1][ii][ij] = grid_pt->deltaW[rk_flag][ii][ij]
	      +(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij) 
		
	    	+getProjectorSource(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUW, deltaUDUW, DDeltaUW, DUDeltaW)
		
	    	-getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    				   dMuUMu, dMuDeltaUMu, ii, ij) );
	  }
	  else{
	    /*
	    grid_pt->Xi[rk_flag+1][ii][ij] = grid_pt->Xi[rk_flag][ii][ij]
	      +(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiXiMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)

	    	-getDeltaXiDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    				   dMuUMu, dMuDeltaUMu, ii, ij)

		+getProjectorSourceXi(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUXi)
		
		+(1./tauPi)*ZCalc[ii][ij]);
	    	
	    grid_pt->deltaW[rk_flag+1][ii][ij] = grid_pt->deltaW[rk_flag][ii][ij]
	      +(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)

	    	-getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    				   dMuUMu, dMuDeltaUMu, ii, ij)

		+getProjectorSource(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUW, deltaUDUW, DDeltaUW, DUDeltaW));
	    */	
	    grid_pt->deltaW[rk_flag+1][ii][ij] = grid_pt->deltaW[rk_flag][ii][ij]
	      +(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)
		
	    	////-getDeltaUDWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)
		
	    	-getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    				   dMuUMu, dMuDeltaUMu, ii, ij)
		
	    	+getProjectorSource(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUW, deltaUDUW, DDeltaUW, DUDeltaW)
	    	
	    	//+(1./DATA->tau_pi)*ZCalc[ii][ij]);
	    	+(1./tauPi)*ZCalc[ii][ij]);
	    
	    //-getDeltaDeltaW0S0Delta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, dMuUMu, ii, ij)
	    //	
	    //	//-getDeltaW0S0DeltaDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, dMuUMu, ii, ij) );
	  }
	}
	if(rk_flag==1){
	  //grid_pt->deltaW[rk_flag+1][ii][ij] = 0.5*(grid_pt->deltaW[0][ii][ij] + grid_pt->deltaW[rk_flag][ii][ij])
	  //-0.5*DATA->delta_tau*( grid_pt->deltaW[rk_flag][ii][ij]
	  //			 -shearV*getDeltaS(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)
	  //			 -ZCalc[ii][ij] )/(DATA->tau_pi*grid_pt->u[rk_flag][0]);
	  if(grid_pt->T_t < DATA->fluctuatingTMin || grid_pt->position[1]<20 || grid_pt->position[1]>=nmax[1]-20
             || grid_pt->position[2]<20 || grid_pt->position[2]>=nmax[2]-20){
	    //grid_pt->deltaW[rk_flag+1][ii][ij] = 0.5*(grid_pt->deltaW[0][ii][ij] + grid_pt->deltaW[rk_flag][ii][ij])
	    //+0.5*DATA->delta_tau*
	    //( -getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    //				   dMuUMu, dMuDeltaUMu, ii, ij) );
	    /*
	    grid_pt->Xi[rk_flag+1][ii][ij] = 0.5*(grid_pt->Xi[0][ii][ij] + grid_pt->Xi[rk_flag][ii][ij])
	      +0.5*(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiXiMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)

		+getProjectorSourceXi(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUXi)

	    	-getDeltaXiDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
				 dMuUMu, dMuDeltaUMu, ii, ij) );

	    grid_pt->deltaW[rk_flag+1][ii][ij] = 0.5*(grid_pt->deltaW[0][ii][ij] + grid_pt->deltaW[rk_flag][ii][ij])
	      +0.5*(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)

		+getProjectorSource(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUW, deltaUDUW, DDeltaUW, DUDeltaW)

	    	-getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    				   dMuUMu, dMuDeltaUMu, ii, ij) );
	    */
	    grid_pt->deltaW[rk_flag+1][ii][ij] = 0.5*(grid_pt->deltaW[0][ii][ij] + grid_pt->deltaW[rk_flag][ii][ij])
	      +0.5*(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)
		
	    	+getProjectorSource(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUW, deltaUDUW, DDeltaUW, DUDeltaW)
		
	    	-getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    				   dMuUMu, dMuDeltaUMu, ii, ij) );
	  }
	  else{
	    /*
	    grid_pt->Xi[rk_flag+1][ii][ij] = 0.5*(grid_pt->Xi[0][ii][ij] + grid_pt->Xi[rk_flag][ii][ij])
	      +0.5*(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiXiMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)

		+getProjectorSourceXi(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUXi)
	
	  	-getDeltaXiDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	  				   dMuUMu, dMuDeltaUMu, ii, ij)

		+(1./tauPi)*ZCalc[ii][ij]);

	    grid_pt->deltaW[rk_flag+1][ii][ij] = 0.5*(grid_pt->deltaW[0][ii][ij] + grid_pt->deltaW[rk_flag][ii][ij])
	      +0.5*(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)

		+getProjectorSource(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUW, deltaUDUW, DDeltaUW, DUDeltaW)
	
	  	-getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	  				   dMuUMu, dMuDeltaUMu, ii, ij));
	    */
	    grid_pt->deltaW[rk_flag+1][ii][ij] = 0.5*(grid_pt->deltaW[0][ii][ij] + grid_pt->deltaW[rk_flag][ii][ij])
	      +0.5*(DATA->delta_tau/grid_pt->u[rk_flag][0])*
	      ( -getUiDiDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)
	    	
	    	//-getDeltaUDWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)
		//	
	    	+getProjectorSource(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, DUW, deltaUDUW, DDeltaUW, DUDeltaW)
		
	    	-getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, deltaEta, 
	    				   dMuUMu, dMuDeltaUMu, ii, ij)
		
	    	//+(1./DATA->tau_pi)*ZCalc[ii][ij]);
	    	+(1./tauPi)*ZCalc[ii][ij]);
	    ////	
	    ////	//-getDeltaDeltaW0S0Delta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, dMuUMu, ii, ij)
	    ////	
	    ////	//-getDeltaW0S0DeltaDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, shearV, dMuUMu, ii, ij) );
	  }
	}
      }

      if(grid_pt->deltaW[rk_flag+1][ii][ij] > 1.e+3 ){
	cerr << "At position (" << grid_pt->position[1] << ", " << grid_pt->position[2] << ", " << grid_pt->position[3]+rank*nmax[3] << "), rk_step= " << rk_flag << ":" << endl;
	cerr << "deltaW[" << ii << "][" << ij << "] has become big!" << endl;
	cerr << "deltaW = " << grid_pt->deltaW[rk_flag+1][ii][ij] << endl;
	cerr << "deltaW_i-1 = " << grid_pt->deltaW[rk_flag][ii][ij] << endl;
	cerr << "T = " << grid_pt->T_t*hbarc << endl;
	cerr << "UiDiDeltaW=" << getUiDiDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)
	     << ", deltaU*D W= " << getDeltaUDWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij)
	     << ", Delta deltaW... Delta=" << getDeltaDeltaWDeltaSDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 
									shearV, deltaEta, dMuUMu, dMuDeltaUMu, ii, ij)
	     << ", Z=" << (1./DATA->tau_pi)*ZCalc[ii][ij] << endl;
	  //<< ", delta Delta W... Delta=" << getDeltaDeltaW0S0Delta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 
	  //							      shearV, dMuUMu, ii, ij)
	  //<< ", Delta W delta Delta=" << getDeltaW0S0DeltaDelta(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 
	  //							   shearV, dMuUMu, ii, ij) << endl;
	cerr << "delta S=" << getDeltaS(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, dMuUMu, dMuDeltaUMu)
	     << ", S=" << getS(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, ij, dMuUMu)
	     << ", partial*u=" << dMuUMu << ", partial*deltaU=" << dMuDeltaUMu << endl;
	cerr << "eta=" << shearV << ", deltaEta=" << deltaEta << ", w_0=" << grid_pt->epsilon_t+grid_pt->p_t << endl;
	cerr << "dpde=" << eos->get_dpOverde2(grid_pt->epsilon_t, grid_pt->rhob_t) << endl;
	cerr << "u = (" << grid_pt->u[rk_flag][0] << ", " << grid_pt->u[rk_flag][1] << ", " << 
	  grid_pt->u[rk_flag][2] << ", " << grid_pt->u[rk_flag][3] << ")" << endl;
      }
    }

  }
  
  //We solve for the perturbation in TJb first by finding deltaT[rk_flag+1][0][nu], using the MacCormack method:
  if(rk_flag == 0){
    for(int nu=0; nu<4; nu++){
      grid_pt->deltaT[rk_flag+1][0][nu] = grid_pt->deltaT[rk_flag][0][nu]
	-getDMuDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, nu)*delta[0];
	//-getDMuXiMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, nu)*delta[0];
      //We use sweeps in all three spatial dimensions:
      if(DATA->Initial_profile != 7){
	for(int i=1; i<3; i++){
	  //Outflowing boundary conditions:
	  if(grid_pt->position[i] == nmax[i]){
	    //grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	    //grid_pt->deltaT[rk_flag+1][0][nu] += grid_pt->deltaT[rk_flag][i][nu]*delta[0]/delta[i];
	    grid_pt->deltaT[rk_flag+1][0][nu] += -(grid_pt->nbr_p_1[i]->deltaT[rk_flag][i][nu] - grid_pt->deltaT[rk_flag][i][nu])*delta[0]/delta[i];
	    //putthisbackgrid_pt->deltaT[rk_flag+1][0][nu] += -(grid_pt->deltaT[rk_flag][i][nu] - grid_pt->nbr_m_1[i]->deltaT[rk_flag][i][nu])*delta[0]/delta[i];
	  }
	  else{
	    grid_pt->deltaT[rk_flag+1][0][nu] += -(grid_pt->nbr_p_1[i]->deltaT[rk_flag][i][nu] - grid_pt->deltaT[rk_flag][i][nu])*delta[0]/delta[i];
	  }
	}

	if(grid_pt->position[3] == nmax[3]){
	  if(size == 1){
	    //grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	    grid_pt->deltaT[rk_flag+1][0][nu] += grid_pt->deltaT[rk_flag][3][nu]*delta[0]/(tau*delta[3]);
	    //putthisbackgrid_pt->deltaT[rk_flag+1][0][nu] += -(grid_pt->deltaT[rk_flag][3][nu]-grid_pt->nbr_m_1[3]->deltaT[rk_flag][3][nu])*delta[0]/(tau*delta[3]);
	  }
	  else{
	    if(rank == size-1){
	      //grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	      grid_pt->deltaT[rk_flag+1][0][nu] += grid_pt->deltaT[rk_flag][3][nu]*delta[0]/(tau*delta[3]);
	      //putthisbackgrid_pt->deltaT[rk_flag+1][0][nu] += -(grid_pt->deltaT[rk_flag][3][nu]-grid_pt->nbr_m_1[3]->deltaT[rk_flag][3][nu])*delta[0]/(tau*delta[3]);
	    }
	    else{
	      grid_pt->deltaT[rk_flag+1][0][nu] += -(Rneighbor->deltaT[rk_flag][3][nu]-grid_pt->deltaT[rk_flag][3][nu])*delta[0]/(tau*delta[3]);
	    }
	  }
	}
	else{
	  grid_pt->deltaT[rk_flag+1][0][nu] += -(grid_pt->nbr_p_1[3]->deltaT[rk_flag][3][nu]-grid_pt->deltaT[rk_flag][3][nu])*delta[0]/(tau*delta[3]);
	}
	//From the Christoffel symbols:
	grid_pt->deltaT[rk_flag+1][0][nu] += -(delta[0]/tau)*grid_pt->deltaT[rk_flag][0][nu];
	if(nu==0){
	  grid_pt->deltaT[rk_flag+1][0][nu] += -(delta[0]/tau)*grid_pt->deltaT[rk_flag][3][3];
	}
	if(nu==3){
	  grid_pt->deltaT[rk_flag+1][0][nu] += -(delta[0]/tau)*grid_pt->deltaT[rk_flag][0][3];
	}
      }
      //When Initial_profile == 7, the coordinates are Cartesian:
      else{
	for(int i=1; i<3; i++){
	  //Outflowing boundary conditions:
	  if(grid_pt->position[i] == nmax[i]){
	    //grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	    //Trying periodic boundary conditions:
	    grid_pt->deltaT[rk_flag+1][0][nu] += -(grid_pt->nbr_p_1[i]->deltaT[rk_flag][i][nu] - grid_pt->deltaT[rk_flag][i][nu])*delta[0]/delta[i];
	  }
	  else{
	    grid_pt->deltaT[rk_flag+1][0][nu] += -(grid_pt->nbr_p_1[i]->deltaT[rk_flag][i][nu] - grid_pt->deltaT[rk_flag][i][nu])*delta[0]/delta[i];
	  }
	}
	if(grid_pt->position[3] == nmax[3]){
	  if(size == 1){
	    grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	  }
	  else{
	    if(rank == size-1){
	      grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	    }
	    else{
	      grid_pt->deltaT[rk_flag+1][0][nu] += -(Rneighbor->deltaT[rk_flag][3][nu]-grid_pt->deltaT[rk_flag][3][nu])*delta[0]/delta[1];
	    }
	  }
	}
	else{
	  grid_pt->deltaT[rk_flag+1][0][nu] += -(grid_pt->nbr_p_1[3]->deltaT[rk_flag][3][nu]-grid_pt->deltaT[rk_flag][3][nu])*delta[0]/delta[1];
	}
	////From the Christoffel symbols:
	//grid_pt->deltaT[rk_flag+1][0][nu] += -(delta[0]/tau)*grid_pt->deltaT[rk_flag][0][nu];
	//if(nu==0){
	//grid_pt->deltaT[rk_flag+1][0][nu] += -(delta[0]/tau)*grid_pt->deltaT[rk_flag][3][3];
	//}
	//if(nu==3){
	//grid_pt->deltaT[rk_flag+1][0][nu] += -(delta[0]/tau)*grid_pt->deltaT[rk_flag][0][3];
	//}
      }
    }
  }

  if(rk_flag == 1){
    for(int nu=0; nu<4; nu++){
      grid_pt->deltaT[rk_flag+1][0][nu] = 0.5*(grid_pt->deltaT[0][0][nu] + grid_pt->deltaT[rk_flag][0][nu])
	-0.5*getDMuDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, nu)*delta[0];
      //-0.5*getDMuXiMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, nu)*delta[0];
      //We use sweeps in all three spatial dimensions:

      if(DATA->Initial_profile != 7){
	for(int i=1; i<3; i++){
	  //Outflowing boundary conditions:
	  if(grid_pt->position[i] == 0){
	    //grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	    //grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*grid_pt->deltaT[rk_flag][i][nu]*delta[0]/delta[i];
	    grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->deltaT[rk_flag][i][nu]-grid_pt->nbr_m_1[i]->deltaT[rk_flag][i][nu])*delta[0]/delta[i];
	    //putthisbackgrid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->nbr_p_1[i]->deltaT[rk_flag][i][nu]-grid_pt->deltaT[rk_flag][i][nu])*delta[0]/delta[i];
	  }
	  else{
	    grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->deltaT[rk_flag][i][nu]-grid_pt->nbr_m_1[i]->deltaT[rk_flag][i][nu])*delta[0]/delta[i];
	  }
	}

	if(grid_pt->position[3] == 0){
	  if(size == 1){
	    //grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	    grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*grid_pt->deltaT[rk_flag][3][nu]*delta[0]/((tau+delta[0])*delta[3]);
	    //putthisbackgrid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->nbr_p_1[3]->deltaT[rk_flag][3][nu]-grid_pt->deltaT[rk_flag][3][nu])*delta[0]/((tau+delta[0])*delta[3]);
	  }
	  else{
	    if(rank == 0){
	      //grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	      grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*grid_pt->deltaT[rk_flag][3][nu]*delta[0]/((tau+delta[0])*delta[3]);
	      //putthisbackgrid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->nbr_p_1[3]->deltaT[rk_flag][3][nu]-grid_pt->deltaT[rk_flag][3][nu])*delta[0]/((tau+delta[0])*delta[3]);
	    }
	    else{
	      grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->deltaT[rk_flag][3][nu]-Lneighbor->deltaT[rk_flag][3][nu])*delta[0]/((tau+delta[0])*delta[3]);
	    }
	  }
	}
	else{
	  grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->deltaT[rk_flag][3][nu]-grid_pt->nbr_m_1[3]->deltaT[rk_flag][3][nu])*delta[0]/((tau+delta[0])*delta[3]);
	}
	//From the Christoffel symbols:
	grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(delta[0]/(tau+delta[0]))*grid_pt->deltaT[rk_flag][0][nu];
	if(nu==0){
	  grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(delta[0]/(tau+delta[0]))*grid_pt->deltaT[rk_flag][3][3];
	}
	if(nu==3){
	  grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(delta[0]/(tau+delta[0]))*grid_pt->deltaT[rk_flag][0][3];
	}
      }
      else{
	for(int i=1; i<3; i++){
	  //Outflowing boundary conditions:
	  if(grid_pt->position[i] == 0){
	    //grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	    //Trying periodic boundary conditions:
	    grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->deltaT[rk_flag][i][nu]-grid_pt->nbr_m_1[i]->deltaT[rk_flag][i][nu])*delta[0]/delta[i];
	  }
	  else{
	    grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->deltaT[rk_flag][i][nu]-grid_pt->nbr_m_1[i]->deltaT[rk_flag][i][nu])*delta[0]/delta[i];
	  }
	}

	if(grid_pt->position[3] == 0){
	  if(size == 1){
	    grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	  }
	  else{
	    if(rank == 0){
	      grid_pt->deltaT[rk_flag+1][0][nu] += 0.;
	    }
	    else{
	      grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->deltaT[rk_flag][3][nu]-Lneighbor->deltaT[rk_flag][3][nu])*delta[0]/delta[1];
	    }
	  }
	}
	else{
	  grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(grid_pt->deltaT[rk_flag][3][nu]-grid_pt->nbr_m_1[3]->deltaT[rk_flag][3][nu])*delta[0]/delta[1];
	}
	////From the Christoffel symbols:
	//grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(delta[0]/(tau+delta[0]))*grid_pt->deltaT[rk_flag][0][nu];
	//if(nu==0){
	//grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(delta[0]/(tau+delta[0]))*grid_pt->deltaT[rk_flag][3][3];
	//}
	//if(nu==3){
	//grid_pt->deltaT[rk_flag+1][0][nu] += -0.5*(delta[0]/(tau+delta[0]))*grid_pt->deltaT[rk_flag][0][3];
	//}
      }
    }
  }

  //Reconstruct the values for deltaP and deltaU:
  reconstructDeltaPAndDeltaU(tau, DATA, grid_pt, rk_flag);

  //return flag;
  return 1;

}

/* %%%%%%%%%%%%%%%%%%%%%% First steps begins here %%%%%%%%%%%%%%%%%% */

int Advance::FirstRKStepT(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, 
double *qi, double *rhs, double **w_rhs, double **qirk, Grid *grid_rk, int size, int rank)
{ 
  //cout << "FirstRKStepT begun..." << endl;

  // this advances the ideal part
 double tau_now, tau_next, tau_rk, tempf, p_rhs, temp_mu, temps, tempd, dwmn;
 int alpha, flag, mu, nu;
 static int ind=0;
 
 tau_now = tau;
 tau_next = tau + (DATA->delta_tau);
 
 if(rk_flag == 0) {tau_rk = tau_now;}
 else if(rk_flag > 0) {tau_rk = tau_next;}

/* TEST */
 if(rk_flag==2) fprintf(stderr, "FirstRKStepT: rk_flag = %d\n", rk_flag);
/* TEST */

 //First, advance Xi if required:
 if(DATA->fluctuatingHydroFlag == 1){
   //cout << "AdvanceXi should start..." << endl;
   AdvanceXi(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank);
 }
 
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

     // The effect of the stochastic term:
     double dXi = 0.;
     //Save time by only calculating this when the fluctuatingHydroFlag == 1. -CFY
     //if(DATA->fluctuatingHydroFlag == 1 && grid_pt->T > DATA->fluctuatingTMin){
     if(DATA->fluctuatingHydroFlag == 1){
       dXi = getDXi(tau, alpha, grid_pt, Lneighbor, Rneighbor, DATA, rk_flag, size, rank);
     }
     //cout << "dXi = " << dXi << ", dwmn = " << dwmn << endl;
 
     /* dwmn is the only one with the minus sign */
     qirk[alpha][0] -= (dwmn + dXi)*(DATA->delta_tau);

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
 double tau_now, tau_next, tempf, p_rhs, temp_mu, temps, eps_ratio, tempd, tau_rk;
 int alpha, flag, mu, nu, revert_flag;
 static int ind=0;

 tau_now = tau;
 tau_next = tau + (DATA->delta_tau);
 
 if(rk_flag == 0) {tau_rk = tau_now;}
 else if(rk_flag > 0) {tau_rk = tau_next;}

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

//if(rk_flag == 0)
//{
/* calculate delta u^0 pi */
  //diss->Make_uPRHS(tau_now, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, &p_rhs, DATA, rk_flag, size, rank);
   
  //tempf = (grid_pt->pi_b[rk_flag])*(grid_pt->u[rk_flag][0]);
   
  //temps = diss->Make_uPiSource(tau_now, grid_pt, DATA, rk_flag);
  //tempf += temps*(DATA->delta_tau);
  //tempf += p_rhs;
   
  //grid_pt->pi_b[rk_flag+1] = tempf/(grid_pt->u[rk_flag+1][0]);
  //temps = tempf = 0.;
   
//}/* rk_flag == 0 */
//else if(rk_flag > 0)
//{
/* calculate delta u^0 pi */
  //diss->Make_uPRHS(tau_next, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, &p_rhs, DATA, rk_flag, size, rank);
   
  //tempf = (grid_pt->pi_b[0])*(grid_pt->u[0][0]);
   
  //temps = diss->Make_uPiSource(tau_next, grid_pt, DATA, rk_flag);
  //tempf += temps*(DATA->delta_tau);
  //tempf += p_rhs;
  
  //tempf += (grid_pt->pi_b[1])*(grid_pt->u[0][0]);
  //tempf *= 0.5;

  //grid_pt->pi_b[rk_flag+1] = tempf/(grid_pt->u[rk_flag+1][0]);

  //temps = tempf = 0.;

//}/* rk_flag > 0 */


///* update Pimunu */
//for(mu=0; mu<4; mu++)
//{
//   for(nu=0; nu<4; nu++)
//   {
//    grid_pt->Pimunu[rk_flag+1][mu][nu]  = (grid_pt->u[rk_flag+1][mu]);
//    grid_pt->Pimunu[rk_flag+1][mu][nu] *= (grid_pt->u[rk_flag+1][nu]);
//    grid_pt->Pimunu[rk_flag+1][mu][nu] += DATA->gmunu[mu][nu];
//    grid_pt->Pimunu[rk_flag+1][mu][nu] *= (grid_pt->pi_b[rk_flag+1]);
//   }/* nu */
//  }/* mu */

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

   //if(isnan(grid_pt->Wmunu[2][0][1])) cout << "6 firstrksW Wmunu[2][0][1]" << grid_pt->Wmunu[2][0][1] << endl;

/* make Wmunu[mu][0] */
     
tempf = 0.0;
for(nu=1; nu<=3; nu++)
      {
       tempf +=
       (grid_pt->Wmunu[rk_flag+1][0][nu])*(grid_pt->u[rk_flag+1][nu]); 
      } 
grid_pt->Wmunu[rk_flag+1][0][0] = tempf/(grid_pt->u[rk_flag+1][0]);
  
//   TestW(tau, DATA, grid_pt, rk_flag);
   
revert_flag = 
  QuestRevert(tau, grid_pt, rk_flag, DATA, size, rank);

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
 double tempf;

 trk_flag = rk_flag+1;

 ////New (CFY):
 //grid_pt->prev_epsilon = grid_pt->epsilon;
 //grid_pt->prev_p = grid_pt->p;

 grid_pt->epsilon_t = grid_rk->epsilon;
 grid_pt->p_t = grid_rk->p;
 grid_pt->rhob_t = grid_rk->rhob;
 grid_pt->T_t = grid_rk->T;
  
 /* reconstructed grid_rk uses rk_flag 0 only */
 for(mu=0; mu<4; mu++)
   {
     tempf = grid_pt->u[trk_flag][mu] = grid_rk->u[0][mu];
     for(alpha=0; alpha<5; alpha++)
       {
	 grid_pt->TJb[trk_flag][alpha][mu] = grid_rk->TJb[0][alpha][mu];
       }/* alpha */
   }/* mu */
 
 ////At this point, update Xi, using the *_t values:
 //grid_pt->getXi(trk_flag, util, eos, DATA);
 
}/* UpdateTJbRK */


//int Advance::QuestRevert(double tau, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank)
//{
//double f11, f22, f33, g11, g22, g33, trace;
//int mu, nu, revert_flag;
//
//revert_flag = 0;
//
//trace = -grid_pt->Wmunu[rk_flag+1][0][0] - grid_pt->Xi[rk_flag+1][0][0];
//for(mu=1; mu<=3; mu++)
//{
// trace += grid_pt->Wmunu[rk_flag+1][mu][mu] + grid_pt->Xi[rk_flag+1][mu][mu];
//}
//
//f11 = grid_pt->TJb[rk_flag+1][1][1];
//g11 = fabs(grid_pt->Wmunu[rk_flag+1][1][1] + grid_pt->Xi[rk_flag+1][1][1]); // + grid_pt->Pimunu[rk_flag+1][1][1]);
// 
//f22 = grid_pt->TJb[rk_flag+1][2][2];
//g22 = fabs(grid_pt->Wmunu[rk_flag+1][2][2] + grid_pt->Xi[rk_flag+1][2][2]);// + grid_pt->Pimunu[rk_flag+1][2][2]);
// 
//f33 = grid_pt->TJb[rk_flag+1][3][3];
//g33 = fabs(grid_pt->Wmunu[rk_flag+1][3][3] + grid_pt->Xi[rk_flag+1][3][3]);// + grid_pt->Pimunu[rk_flag+1][3][3]);
//
//if( (f11 < g11) || (f22 < g22) || (f33 < g33) )
//{
//  //  cout << "reverting. epsilon=" <<  grid_pt->epsilon << endl;
//   double x=grid_pt->position[1]*DATA->delta_x-DATA->x_size/2;
//   double y=grid_pt->position[2]*DATA->delta_y-DATA->y_size/2;
//   double eta;
//   if(size>1)
//     eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2+static_cast<double>(rank)/static_cast<double>(size)*DATA->eta_size;
//   else
//     eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2;
//   double r;
//   r=sqrt(x*x+y*y);
//   if( grid_pt->epsilon > 10. ) // was 0.3 , changed by Schenke 2010/07/16  
//     {
//	 fprintf(stderr, 
//		 "QuestRevert: Viscous correction is bigger than T ");
//	 
//	 fprintf(stderr, "at (%d, %d, %d).\n", 
//		 grid_pt->position[1], grid_pt->position[2], grid_pt->position[3]);
//	 
//	 fprintf(stderr, "at (r=%f, eta=%f).\n",r, eta); 
//	 if (r<3.) fprintf(stderr, "WARNING: r=%f fm < 3 fm.\n",r); 
//	 
//	 grid_pt->trouble ++;
//	 //cout << grid_pt->trouble << endl;
//	 
//	 fprintf(stderr, 
//		 "Reporting because epsilon = %e is larger than 10/fm^4 \n", 
//		 grid_pt->epsilon);
//	 
//	 fprintf(stderr, "trace = %e\n", trace);
//	 if(DATA->zero_or_stop==2) // add this in name : 2 means that we set W^{munu}=0.9*T{munu}
//	   {
//	     fprintf(stderr, "Setting W^{munu}=0.95T^{munu}...\n");
//	   }
//	 if(DATA->zero_or_stop==1)
//	   {
//	     fprintf(stderr, "Reverting to the previous values...\n");
//	   }
//	 if(DATA->zero_or_stop==0)
//	   {
//	     fprintf(stderr, "Zeroing Wmunu...\n");
//	   }
//	 
//     }
//   
//   /* revert back */
//   /* reverting everything back to the previous value seems to not work
// that well. It leaves a core in the middle untouched. 
// must find a way to quench the viscous part but make the ideal part
// still proceed. However, this tends to make some part of the total T_ii
// negative. Need to reconcile those two. Let's try zero-ing.
// Zero-ing does stabilize. However, it looks like that it may make 
// the viscous effect too small. Also zeroing produces spikes a lot more than reverting.
// However, reverting will produce too large W^\mu\nu's at freeze-out that mess up
// the correction delta_f to the distribution function (it becomes way too large).
// */
//
// if(DATA->zero_or_stop==1)
//  {
//    //grid_pt->pi_b[rk_flag+1] = grid_pt->pi_b[rk_flag];
//   for(mu=0; mu<4; mu++)
//     {
//	 for(nu=0; nu<4; nu++)
//	   {
//	     //grid_pt->Pimunu[rk_flag+1][mu][nu] = grid_pt->Pimunu[rk_flag][mu][nu];
//	     grid_pt->Wmunu[rk_flag+1][mu][nu] = grid_pt->Wmunu[rk_flag][mu][nu]; 
//	     grid_pt->Xi[rk_flag+1][mu][nu] = grid_pt->Xi[rk_flag][mu][nu]; 
//	   }
//     } /* mu nu */
//  }
// else if(DATA->zero_or_stop==0)
//  {
//    //grid_pt->pi_b[rk_flag+1] = 0.0;
//    for(mu=0; mu<4; mu++)
//	{
//	  for(nu=0; nu<4; nu++)
//	    {
//	      //grid_pt->Pimunu[rk_flag+1][mu][nu] = 0.0;
//	      grid_pt->Wmunu[rk_flag+1][mu][nu] = 0.0;
//	      grid_pt->Xi[rk_flag+1][mu][nu] = 0.0;
//	    }
//	} /* mu nu */
//  }
// else if(DATA->zero_or_stop==2) 
//  {
//    //grid_pt->pi_b[rk_flag+1] = 0.95*grid_pt->pi_b[rk_flag];
//    for(mu=0; mu<4; mu++)
//	{
//	  for(nu=0; nu<4; nu++)
//	    {
//	      //grid_pt->Pimunu[rk_flag+1][mu][nu] = 0.95*grid_pt->Pimunu[rk_flag][mu][nu];
//	      grid_pt->Wmunu[rk_flag+1][mu][nu] = 0.95*grid_pt->Wmunu[rk_flag][mu][nu];
//	      grid_pt->Xi[rk_flag+1][mu][nu] = 0.95*grid_pt->Xi[rk_flag][mu][nu];
//	    }
//	} /* mu nu */
//  }
// 
// revert_flag = 1;
//}/* if something is wrong, revert back to the original values */
//return revert_flag;
//}/* QuestRevert */

int Advance::QuestRevert(double tau, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank)
{
  double f11, f22, f33, g11, g22, g33, trace;
  int mu, nu, revert_flag;
  
  revert_flag = 0;
  
  trace = -grid_pt->Wmunu[rk_flag+1][0][0];
  for(mu=1; mu<=3; mu++)
    {
      trace += grid_pt->Wmunu[rk_flag+1][mu][mu];
    }
  f11 = grid_pt->TJb[rk_flag+1][1][1];
  g11 = fabs(grid_pt->Wmunu[rk_flag+1][1][1]); // + grid_pt->Pimunu[rk_flag+1][1][1]);
  
  f22 = grid_pt->TJb[rk_flag+1][2][2];
  g22 = fabs(grid_pt->Wmunu[rk_flag+1][2][2]);// + grid_pt->Pimunu[rk_flag+1][2][2]);
  
  f33 = grid_pt->TJb[rk_flag+1][3][3];
  g33 = fabs(grid_pt->Wmunu[rk_flag+1][3][3]);// + grid_pt->Pimunu[rk_flag+1][3][3]);
  
  if( (f11 < g11) || (f22 < g22) || (f33 < g33) )
    {
      //  cout << "reverting. epsilon=" <<  grid_pt->epsilon << endl;
      double x=grid_pt->position[1]*DATA->delta_x-DATA->x_size/2;
      double y=grid_pt->position[2]*DATA->delta_y-DATA->y_size/2;
      double eta;
      if(size>1)
	eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2+static_cast<double>(rank)/static_cast<double>(size)*DATA->eta_size;
      else
	eta=grid_pt->position[3]*DATA->delta_eta-DATA->eta_size/2;
      double r;
      r=sqrt(x*x+y*y);
      if( grid_pt->epsilon > 10. ) // was 0.3 , changed by Schenke 2010/07/16  
	{
	  fprintf(stderr, 
		  "QuestRevert: Viscous correction is bigger than T ");
	  
	  fprintf(stderr, "at (%d, %d, %d).\n", 
		  grid_pt->position[1], grid_pt->position[2], grid_pt->position[3]);
	  
	  //	 fprintf(stderr, "at (r=%f, eta=%f).\n",r, eta); 
	  if (r<3.) fprintf(stderr, "WARNING: r=%f fm < 3 fm.\n",r); 
	  
	  grid_pt->trouble ++;
	  //cout << grid_pt->trouble << endl;
	  
	  fprintf(stderr, 
		  "Reporting because epsilon = %e is larger than 10/fm^4 \n", 
		  grid_pt->epsilon);
	  
	  fprintf(stderr, "trace = %e\n", trace);
	  if(DATA->zero_or_stop==2) // add this in name : 2 means that we set W^{munu}=0.9*T{munu}
	    {
	      fprintf(stderr, "Setting W^{munu}=0.95T^{munu}...\n");
	    }
	  if(DATA->zero_or_stop==1)
	    {
	      fprintf(stderr, "Reverting to the previous values...\n");
	    }
	  if(DATA->zero_or_stop==0)
	    {
	      fprintf(stderr, "Zeroing Wmunu...\n");
	    }
	  
	}
      
      /* revert back */
      /* reverting everything back to the previous value seems to not work
	 that well. It leaves a core in the middle untouched. 
	 must find a way to quench the viscous part but make the ideal part
	 still proceed. However, this tends to make some part of the total T_ii
	 negative. Need to reconcile those two. Let's try zero-ing.
	 Zero-ing does stabilize. However, it looks like that it may make 
	 the viscous effect too small. Also zeroing produces spikes a lot more than reverting.
	 However, reverting will produce too large W^\mu\nu's at freeze-out that mess up
	 the correction delta_f to the distribution function (it becomes way too large).
      */
      
      if(DATA->zero_or_stop==1)
	{
	  //grid_pt->pi_b[rk_flag+1] = grid_pt->pi_b[rk_flag];
	  for(mu=0; mu<4; mu++)
	    {
	      for(nu=0; nu<4; nu++)
		{
		  //grid_pt->Pimunu[rk_flag+1][mu][nu] = grid_pt->Pimunu[rk_flag][mu][nu];
		  grid_pt->Wmunu[rk_flag+1][mu][nu] = grid_pt->Wmunu[rk_flag][mu][nu]; 
		}
	    } /* mu nu */
	}
      else if(DATA->zero_or_stop==0)
	{
	  //grid_pt->pi_b[rk_flag+1] = 0.0;
	  for(mu=0; mu<4; mu++)
	    {
	      for(nu=0; nu<4; nu++)
		{
		  //grid_pt->Pimunu[rk_flag+1][mu][nu] = 0.0;
		  grid_pt->Wmunu[rk_flag+1][mu][nu] = 0.0;
		}
	    } /* mu nu */
	}
      else if(DATA->zero_or_stop==2) 
	{
	  //grid_pt->pi_b[rk_flag+1] = 0.95*grid_pt->pi_b[rk_flag];
	  for(mu=0; mu<4; mu++)
	    {
	      for(nu=0; nu<4; nu++)
		{
		  //grid_pt->Pimunu[rk_flag+1][mu][nu] = 0.95*grid_pt->Pimunu[rk_flag][mu][nu];
		  grid_pt->Wmunu[rk_flag+1][mu][nu] = 0.95*grid_pt->Wmunu[rk_flag][mu][nu];
		}
	    } /* mu nu */
	}
      
      revert_flag = 1;
    }/* if something is wrong, revert back to the original values */
  return revert_flag;
}/* QuestRevert */


void Advance::TestW(double tau, InitData *DATA, Grid *grid_pt, int rk_flag)
{
 int mu, nu;
 double trace, transv, mufac, nufac;

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
 double Delta[4][4], tD[4][4], f[4][4], trace, norm, mfac, nfac, sum;
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
 double Delta[4][4], f[4][4], trace, norm, mfac, nfac;
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
  
  for(alpha=0; alpha<5; alpha++){
    qi[alpha] = grid_pt->TJb[rk_flag][alpha][0]*tau;
  }/* get qi first */
  
  /* implement Kurganov-Tadmor scheme */
  GetQIs(tau, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, qi, &NbrCells, rk_flag, DATA, size, rank);
  
  flag = MakeQIHalfs(qi, &NbrCells, &HalfwayCells, grid_pt, DATA);
  
  flag = ConstHalfwayCells(tau, &HalfwayCells, qi, grid_pt, DATA, rk_flag, size, rank);
  
  MakeKTCurrents(tau, DFmmp, grid_pt, &HalfwayCells, rk_flag);
  
  for(alpha=0; alpha<5; alpha++){
    sumf = 0.0; 
    for(i=1; i<=3; i++){
      sumf += DFmmp[alpha][i]/delta[i];
    }/* i */
    
    if(alpha==0){
      sumf -= grid_pt->TJb[rk_flag][3][3];
    }
    else if(alpha==3){
      sumf -= grid_pt->TJb[rk_flag][3][0];
    }
    
    //   cout << tau << " " << x << " " << " " << y << " " << eta << endl;
    
    int jetPosition=DATA->includeJet;
    int includeTrigger = DATA->includeTrigger;
    
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


int Advance::MakeQIHalfs(double *qi, NbrQs *NbrCells, BdryCells *HalfwayCells, 
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
	  
	  tempf = DFmmp[alpha][i] = (Fimh[alpha][i] - Fiph[alpha][i]);
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

//Advances Xi to the current rk_order == 2. -CFY
void Advance::AdvanceXi(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank){
  //cout << "AdvanceXi started..." << endl;

  if(rk_flag != 0 && rk_flag != 1){
    cout << "Error: AdvanceXi called when rk_flag = " << rk_flag << ", this Runge-Kutta order is not yet instantiated! Exiting..." << endl;
    return;
  }

  double delta[4];
  int nmax[4];

  double Xi_i, Xi_ip1[4], Xi_im1[4];
  double uDXi, noiseTerm;
  
  nmax[0] = 0; //Should not be used. -CFY
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;
  
  double tauNow;
  if(rk_flag == 0){ tauNow = tau;}
  else{ tauNow = tau + DATA->delta_tau;}

  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = tauNow*DATA->delta_eta;
  
  //For the given grid cell, determine the noise tensor:
  double entrDens = eos->s_func(grid_pt->epsilon_t, grid_pt->p_t, grid_pt->rhob_t, DATA);
  //cout << "entrDens = " << entrDens << " when epsilon_t = " << grid_pt->epsilon_t << " and p_t = " << grid_pt->p_t << endl;
  double shearV = (DATA->shear_to_s)*entrDens;
  //cout << "shear viscosity = " << shearV << endl;
  
  //Determine the 3x3 Z tensor in the fluid rest frame:
  double ZRest[3][3];
  double multFactor = sqrt(2.*grid_pt->T_t*shearV);
  //cout << "T = " << grid_pt->T_t << endl;
  //cout << "multFactor = " << multFactor << endl;
  //cout << "ZRest = (" ;
  for(int mu=0; mu<3; mu++){
    for(int alpha=0; alpha<3; alpha++){
      ZRest[mu][alpha] = multFactor*grid_pt->S[mu][alpha];
      //cout << ZRest[mu][alpha] << " ";
      //ZRest[mu][alpha] = 2.*grid_pt->T_t*shearV*grid_pt->S[mu][alpha];
      //if(mu == alpha){ ZRest[mu][alpha] *= 4./3.; }
    }
    //cout << endl;
  }

  //Now, boost into the grid's frame:
  double uBoost[4];
  for(int mu=0; mu<4; mu++){
    uBoost[mu] = grid_pt->u[rk_flag][mu];
  }
  double ZCalc[4][4];

  //cout << "u = (" << uBoost[0] << ", " << uBoost[1] << ", " << uBoost[2] << ", " << uBoost[3] << ")" << endl;
  util->LorentzBoost3x3Tensor(ZRest, uBoost, ZCalc);
  //cout << "ZCalc = {" << ZCalc[0][0] << ", " << ZCalc[0][1] << ", " << ZCalc[0][2] << ", " << ZCalc[0][3] << endl;
  //cout << "         " << ZCalc[1][0] << ", " << ZCalc[1][1] << ", " << ZCalc[1][2] << ", " << ZCalc[1][3] << endl;
  //cout << "         " << ZCalc[2][0] << ", " << ZCalc[2][1] << ", " << ZCalc[2][2] << ", " << ZCalc[2][3] << endl;
  //cout << "         " << ZCalc[3][0] << ", " << ZCalc[3][1] << ", " << ZCalc[3][2] << ", " << ZCalc[3][3] << "}" << endl;


  //The loop over the indices of Xi:
  for(int mu=0; mu<4; mu++){
    for(int nu=0; nu<4; nu++){
      if(mu > nu){
	grid_pt->Xi[rk_flag+1][mu][nu] = grid_pt->Xi[rk_flag+1][nu][mu]; //Require Xi to be symmetric, and use the solution already obtained.
      }
      else{
	Xi_i = grid_pt->Xi[rk_flag][mu][nu];
	Xi_ip1[0] = Xi_im1[0] = 0.; //No time derivatives are involved here.
	for(int i=1; i<3; i++){
	  if(grid_pt->position[i] == nmax[i]){
	    Xi_ip1[i] = Xi_i ;
	    Xi_im1[i] = grid_pt->nbr_m_1[i]->Xi[rk_flag][mu][nu];
	  } 
	  else if(grid_pt->position[i] == 0){
	    Xi_ip1[i] = grid_pt->nbr_p_1[i]->Xi[rk_flag][mu][nu];
	    Xi_im1[i] = Xi_i;
	  }
	  else{
	    Xi_ip1[i] = grid_pt->nbr_p_1[i]->Xi[rk_flag][mu][nu];
	    Xi_im1[i] = grid_pt->nbr_m_1[i]->Xi[rk_flag][mu][nu];
	  }
	}/* i */

	if(grid_pt->position[3] == nmax[3]){
	  if(rank == size-1){
	    Xi_ip1[3] = Xi_i;
	    Xi_im1[3] = grid_pt->nbr_m_1[3]->Xi[rk_flag][mu][nu];
	  }
	  else{
	    Xi_ip1[3] = Rneighbor->Xi[rk_flag][mu][nu];
	    Xi_im1[3] = grid_pt->nbr_m_1[3]->Xi[rk_flag][mu][nu];
	  }
	}
	else if(grid_pt->position[3] == 0){
	  if(rank == 0){
	    Xi_ip1[3] = grid_pt->nbr_p_1[3]->Xi[rk_flag][mu][nu];
	    Xi_im1[3] = Xi_i;
	  }
	  else{
	    Xi_ip1[3] = grid_pt->nbr_p_1[3]->Xi[rk_flag][mu][nu];
	    Xi_im1[3] = Lneighbor->Xi[rk_flag][mu][nu];
	  }
	}
	else{
	  Xi_ip1[3] = grid_pt->nbr_p_1[3]->Xi[rk_flag][mu][nu];
	  Xi_im1[3] = grid_pt->nbr_m_1[3]->Xi[rk_flag][mu][nu];
	}
        
	uDXi = 0.; double DXi;
	for(int i=1; i<4; i++){
	  DXi = minmod->minmod_dx(Xi_ip1[i], Xi_i, Xi_im1[i], DATA)/delta[i];
	  if(i==3 && mu==0){
	    DXi += grid_pt->Xi[rk_flag][3][nu]/tauNow;
	  }
	  if(i==3 && nu==0){
	    DXi += grid_pt->Xi[rk_flag][mu][3]/tauNow;
	  }
	  if(i==3 && mu==3){
	    DXi += grid_pt->Xi[rk_flag][0][nu]/tauNow;
	  }
	  if(i==3 && nu==3){
	    DXi += grid_pt->Xi[rk_flag][mu][0]/tauNow;
	  }
	  //Finally, add to uDXi:
	  uDXi += grid_pt->u[rk_flag][i]*DXi;
	}
	
	//With uDXi determined, the stochastic term in the equation must also be determined:
	noiseTerm = (Xi_i-ZCalc[mu][nu])/DATA->tau_pi;
	
	//No reconstruction is necessary; if uDXi and noiseTerm and real numbers, the determination of Xi[rk_flag+1][mu][nu] is possible:

	if(rk_flag==0){
	  grid_pt->Xi[1][mu][nu] = Xi_i -(uDXi + noiseTerm)*delta[0]/grid_pt->u[rk_flag][0] ;
	}
	if(rk_flag==1){
	  grid_pt->Xi[2][mu][nu] = 0.5*(grid_pt->Xi[0][mu][nu] + Xi_i -(uDXi + noiseTerm)*delta[0]/grid_pt->u[rk_flag][0] );
	}
	if(rk_flag!=0 && rk_flag!=1){
	  cout << "AdvanceXi: rk_flag = " << rk_flag << "!" << endl;
	}
      }
    }
  }  

  //cout << "Xi[" << rk_flag << "] = {" << grid_pt->Xi[rk_flag][0][0] << ", " << grid_pt->Xi[rk_flag][0][1] << ", " << grid_pt->Xi[rk_flag][0][2] 
  //   << ", " << grid_pt->Xi[rk_flag][0][3] << "," << endl;
  //cout << "      " << grid_pt->Xi[rk_flag][1][0] << ", " << grid_pt->Xi[rk_flag][1][1] << ", " << grid_pt->Xi[rk_flag][1][2] << ", " 
  //   << grid_pt->Xi[rk_flag][1][3] << "," << endl;   
  //cout << "      " << grid_pt->Xi[rk_flag][2][0] << ", " << grid_pt->Xi[rk_flag][2][1] << ", " << grid_pt->Xi[rk_flag][2][2] << ", " 
  //   << grid_pt->Xi[rk_flag][2][3] << "," << endl;
  //cout << "      " << grid_pt->Xi[rk_flag][3][0] << ", " << grid_pt->Xi[rk_flag][3][1] << ", " << grid_pt->Xi[rk_flag][3][2] << ", " 
  //   << grid_pt->Xi[rk_flag][3][3] << "}" << endl;

  //cout << "Xi[" << rk_flag+1 << "] = {" << grid_pt->Xi[rk_flag+1][0][0] << ", " << grid_pt->Xi[rk_flag+1][0][1] << ", " << grid_pt->Xi[rk_flag+1][0][2] 
  //   << ", " << grid_pt->Xi[rk_flag+1][0][3] << "," << endl;
  //cout << "      " << grid_pt->Xi[rk_flag+1][1][0] << ", " << grid_pt->Xi[rk_flag+1][1][1] << ", " << grid_pt->Xi[rk_flag+1][1][2] << ", " 
  //   << grid_pt->Xi[rk_flag+1][1][3] << "," << endl;   
  //cout << "      " << grid_pt->Xi[rk_flag+1][2][0] << ", " << grid_pt->Xi[rk_flag+1][2][1] << ", " << grid_pt->Xi[rk_flag+1][2][2] << ", " 
  //   << grid_pt->Xi[rk_flag+1][2][3] << "," << endl;
  //cout << "      " << grid_pt->Xi[rk_flag+1][3][0] << ", " << grid_pt->Xi[rk_flag+1][3][1] << ", " << grid_pt->Xi[rk_flag+1][3][2] << ", " 
  //   << grid_pt->Xi[rk_flag+1][3][3] << "}" << endl;   
  
  return;
  
}

//Returns \tau \partial_\mu \Xi^{\mu \alpha}, despite the name. -CFY
double Advance::getDXi(double tau, int alpha, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		       InitData *DATA, int rk_flag, int size, int rank){
  double delta[4];
  int nmax[4];

  double dXidt, dXidxi, Xi_i, Xi_ip1, Xi_im1;
  double dXi = 0.;
  
  if(DATA->turn_on_shear){
    
    nmax[0] = 0; //Should not be used. -CFY
    nmax[1] = DATA->nx;
    nmax[2] = DATA->ny;
    nmax[3] = DATA->neta-1;
    
    delta[0] = DATA->delta_tau;
    delta[1] = DATA->delta_x;
    delta[2] = DATA->delta_y;
    delta[3] = DATA->delta_eta;
    
    if(rk_flag==0){
      dXidt = (grid_pt->Xi[rk_flag][0][alpha] - grid_pt->prev_Xi0[alpha])/delta[0];
    }
    else if(rk_flag > 0){
      dXidt = (grid_pt->Xi[rk_flag][0][alpha] - grid_pt->Xi[0][0][alpha])/delta[0];
    }
    
    dXi += dXidt;
    
    for(int i=1; i<=2; i++){
      Xi_i = grid_pt->Xi[rk_flag][i][alpha];
      
      if(grid_pt->position[i] == nmax[i]){
	Xi_ip1 = Xi_i ;
	Xi_im1 = grid_pt->nbr_m_1[i]->Xi[rk_flag][i][alpha];
      } 
      else if(grid_pt->position[i] == 0){
	Xi_ip1 = grid_pt->nbr_p_1[i]->Xi[rk_flag][i][alpha];
	Xi_im1 = Xi_i;
      }
      else{
	Xi_ip1 = grid_pt->nbr_p_1[i]->Xi[rk_flag][i][alpha];
	Xi_im1 = grid_pt->nbr_m_1[i]->Xi[rk_flag][i][alpha];
      }
      
      dXi += minmod->minmod_dx(Xi_ip1, Xi_i, Xi_im1, DATA)/delta[i]; 
    }/* i */

    Xi_i = grid_pt->Xi[rk_flag][3][alpha];
 
    if(grid_pt->position[3] == nmax[3])
      {
	if(rank == size-1){
	    Xi_ip1 = Xi_i;
	    Xi_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][alpha];
	}
	else{
	  Xi_ip1 = Rneighbor->Xi[rk_flag][3][alpha];
	  Xi_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][alpha];
	}
      }
    else if(grid_pt->position[3] == 0){
      if(rank == 0){
	Xi_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][alpha];
	Xi_im1 = Xi_i;
      }
      else{
	Xi_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][alpha];
	Xi_im1 = Lneighbor->Xi[rk_flag][3][alpha];
      }
    }
    else{
	Xi_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][alpha];
	Xi_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][alpha];
      }

    dXi += minmod->minmod_dx(Xi_ip1, Xi_i, Xi_im1, DATA)/(delta[3]*tau); 

    //\partial_m (tau Xi^mn) = Xi^0n + tau \partial_m Xi^mn

    //tau \partial_m Xi^mn
    dXi *= tau;

    //Xi^0n
    dXi += grid_pt->Xi[rk_flag][0][alpha];

    //Sources due to coordinate transform this is added to partial_m Xi^mn */

    if(alpha == 0){
      dXi += grid_pt->Xi[rk_flag][3][3];
    }

    if(alpha == 3){
      dXi += grid_pt->Xi[rk_flag][0][3];
    }

    if(isnan(dXi)){
      cout << "Problem in Advance::getDXi!" << endl;
    }
  }

  //cout << "dXi = " << dXi << endl;

  return dXi ;
}

//Uses MacCormack's method:
double Advance::getDMuDeltaWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int nu){
  if(DATA->rk_order != 2){
    cerr << "rk_order must be 2! Exiting..." << endl;
    exit(1);
  }
  
  double dMuDeltaWMuNu = 0.;
  double delta[4];
  int nmax[4];

  nmax[0] = 0; //Should not be used. -CFY
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;

  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = DATA->delta_eta;

  double w_i, w_ip1, w_im1;

  //if(rk_flag == 0){
  //dMuDeltaWMuNu += (grid_pt->deltaW[0][0][nu] - grid_pt->prev_deltaW[0][nu])/delta[0];
  //}
  //else{
  ////if(rk_flag > 1){
  ////Bug fix?
  ////if(rk_flag == 1){
  //  dMuDeltaWMuNu += (grid_pt->deltaW[rk_flag][0][nu] - grid_pt->deltaW[rk_flag-1][0][nu])/delta[0];
  //  //}
  //}

  //This is why the subroutine should only be called after the rk_flag+1-th step has been done:
  dMuDeltaWMuNu += (grid_pt->deltaW[rk_flag+1][0][nu] - grid_pt->deltaW[rk_flag][0][nu])/delta[0];

  if(DATA->Initial_profile != 7){
    for(int i=1; i<=2; i++){
      w_i = grid_pt->deltaW[rk_flag][i][nu];
      
      if(grid_pt->position[i] == nmax[i]){
	//putthisbackw_ip1 = w_i ;
	//w_ip1 = 0.;
	w_ip1 = grid_pt->nbr_p_1[i]->deltaW[rk_flag][i][nu];
	w_im1 = grid_pt->nbr_m_1[i]->deltaW[rk_flag][i][nu];
      }
      else if(grid_pt->position[i] == 0){
	w_ip1 = grid_pt->nbr_p_1[i]->deltaW[rk_flag][i][nu];
	//puthisbackw_im1 = w_i;
	w_im1 = grid_pt->nbr_m_1[i]->deltaW[rk_flag][i][nu];
	//w_im1 = 0.;
      }
      else{
	w_ip1 = grid_pt->nbr_p_1[i]->deltaW[rk_flag][i][nu];
	w_im1 = grid_pt->nbr_m_1[i]->deltaW[rk_flag][i][nu];
      }
      
      if(rk_flag == 0){
	dMuDeltaWMuNu += (w_ip1-w_i)/delta[i];
      }
      if(rk_flag == 1){
	dMuDeltaWMuNu += (w_i-w_im1)/delta[i];
      }
    }

    w_i = grid_pt->deltaW[rk_flag][3][nu];

    if(grid_pt->position[3] == nmax[3]){
      if(size == 1){
	//w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
	w_ip1 = w_i;
	//w_ip1 = 0.;
	w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
      }
      else{
	if(rank == size-1){
	  w_ip1 = w_i;
	  //w_ip1 = 0.;
	  w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
	}
	else{
	  w_ip1 = Rneighbor->deltaW[rk_flag][3][nu];
	  w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
	}
      }
    }
    else if(grid_pt->position[3] == 0){
      if(size == 1){
	w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
	//w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
	w_im1 = w_i;
	//w_im1 = 0.;
      }
      else{
	if(rank == 0){
	  //Bug fix? -CFY
	  w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
	  w_im1 = w_i;
	  //w_im1 = 0.;
	}
	else{
	  w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
	  w_im1 = Lneighbor->deltaW[rk_flag][3][nu];
	}
      }
    }
    else{
      w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
      w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
    }
    
    if(rk_flag == 0){
      dMuDeltaWMuNu += (w_ip1-w_i)/(delta[3]*tau);
    }
    if(rk_flag == 1){
      dMuDeltaWMuNu += (w_i-w_im1)/(delta[3]*(tau+delta[0]));
    }
    
    //The terms due to the coordinate transform:
    if(rk_flag == 0){
      dMuDeltaWMuNu += grid_pt->deltaW[rk_flag][0][nu]/tau;
      if(nu == 0){
	dMuDeltaWMuNu += grid_pt->deltaW[rk_flag][3][3]/tau;
      }
      
      if(nu == 3){
	dMuDeltaWMuNu += grid_pt->deltaW[rk_flag][3][0]/tau;
      }
    }
    if(rk_flag == 1){
      dMuDeltaWMuNu += grid_pt->deltaW[rk_flag][0][nu]/(tau+delta[0]);
      if(nu == 0){
	dMuDeltaWMuNu += grid_pt->deltaW[rk_flag][3][3]/(tau+delta[0]);
      }
      
      if(nu == 3){
	dMuDeltaWMuNu += grid_pt->deltaW[rk_flag][3][0]/(tau+delta[0]);
      }
    }
  }
  else{
    for(int i=1; i<=2; i++){
      w_i = grid_pt->deltaW[rk_flag][i][nu];
      
      //Trying periodic boundary conditions:
      //if(grid_pt->position[i] == nmax[i]){
      //w_ip1 = w_i ;
      //w_im1 = grid_pt->nbr_m_1[i]->deltaW[rk_flag][i][nu];
      //}
      //else if(grid_pt->position[i] == 0){
      //w_ip1 = grid_pt->nbr_p_1[i]->deltaW[rk_flag][i][nu];
      //w_im1 = w_i;
      //}
      //else{
      w_ip1 = grid_pt->nbr_p_1[i]->deltaW[rk_flag][i][nu];
      w_im1 = grid_pt->nbr_m_1[i]->deltaW[rk_flag][i][nu];
      //}
      
      if(rk_flag == 0){
	dMuDeltaWMuNu += (w_ip1-w_i)/delta[i];
      }
      if(rk_flag == 1){
	dMuDeltaWMuNu += (w_i-w_im1)/delta[i];
      }
    }
    
    w_i = grid_pt->deltaW[rk_flag][3][nu];
    
    if(grid_pt->position[3] == nmax[3]){
      if(size == 1){
	//w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
	w_ip1 = w_i;
	w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
      }
      else{
	if(rank == size-1){
	  w_ip1 = w_i;
	  w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
	}
	else{
	  w_ip1 = Rneighbor->deltaW[rk_flag][3][nu];
	  w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
	}
      }
    }
    else if(grid_pt->position[3] == 0){
      if(size == 1){
	w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
	//w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
	w_im1 = w_i;
      }
      else{
	if(rank == 0){
	  //Bug fix? -CFY
	  w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
	  w_im1 = w_i;
	}
	else{
	  w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
	  w_im1 = Lneighbor->deltaW[rk_flag][3][nu];
	}
      }
    }
    else{
      w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
      w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
    }
    
    if(rk_flag == 0){
      dMuDeltaWMuNu += (w_ip1-w_i)/delta[1];
    }
    if(rk_flag == 1){
      dMuDeltaWMuNu += (w_i-w_im1)/delta[1];
    }
  }
  
  return dMuDeltaWMuNu;

}

//Uses MacCormack's method:
double Advance::getDMuXiMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int nu){
  if(DATA->rk_order != 2){
    cerr << "rk_order must be 2! Exiting..." << endl;
    exit(1);
  }
  
  double dMuXiMuNu = 0.;
  double delta[4];
  int nmax[4];

  nmax[0] = 0; //Should not be used. -CFY
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;

  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = DATA->delta_eta;

  double w_i, w_ip1, w_im1;

  //if(rk_flag == 0){
  //dMuDeltaWMuNu += (grid_pt->deltaW[0][0][nu] - grid_pt->prev_deltaW[0][nu])/delta[0];
  //}
  //else{
  ////if(rk_flag > 1){
  ////Bug fix?
  ////if(rk_flag == 1){
  //  dMuDeltaWMuNu += (grid_pt->deltaW[rk_flag][0][nu] - grid_pt->deltaW[rk_flag-1][0][nu])/delta[0];
  //  //}
  //}

  //This is why the subroutine should only be called after the rk_flag+1-th step has been done:
  dMuXiMuNu += (grid_pt->Xi[rk_flag+1][0][nu] - grid_pt->Xi[rk_flag][0][nu])/delta[0];

  if(DATA->Initial_profile != 7){
    for(int i=1; i<=2; i++){
      w_i = grid_pt->Xi[rk_flag][i][nu];
      
      if(grid_pt->position[i] == nmax[i]){
	w_ip1 = w_i ;
	//w_ip1 = 0.;
	w_im1 = grid_pt->nbr_m_1[i]->Xi[rk_flag][i][nu];
      }
      else if(grid_pt->position[i] == 0){
	w_ip1 = grid_pt->nbr_p_1[i]->Xi[rk_flag][i][nu];
	w_im1 = w_i;
	//w_im1 = 0.;
      }
      else{
	w_ip1 = grid_pt->nbr_p_1[i]->Xi[rk_flag][i][nu];
	w_im1 = grid_pt->nbr_m_1[i]->Xi[rk_flag][i][nu];
      }
      
      if(rk_flag == 0){
	dMuXiMuNu += (w_ip1-w_i)/delta[i];
      }
      if(rk_flag == 1){
	dMuXiMuNu += (w_i-w_im1)/delta[i];
      }
    }

    w_i = grid_pt->Xi[rk_flag][3][nu];

    if(grid_pt->position[3] == nmax[3]){
      if(size == 1){
	//w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
	w_ip1 = w_i;
	//w_ip1 = 0.;
	w_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][nu];
      }
      else{
	if(rank == size-1){
	  w_ip1 = w_i;
	  //w_ip1 = 0.;
	  w_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][nu];
	}
	else{
	  w_ip1 = Rneighbor->Xi[rk_flag][3][nu];
	  w_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][nu];
	}
      }
    }
    else if(grid_pt->position[3] == 0){
      if(size == 1){
	w_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][nu];
	//w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
	w_im1 = w_i;
	//w_im1 = 0.;
      }
      else{
	if(rank == 0){
	  //Bug fix? -CFY
	  w_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][nu];
	  w_im1 = w_i;
	  //w_im1 = 0.;
	}
	else{
	  w_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][nu];
	  w_im1 = Lneighbor->Xi[rk_flag][3][nu];
	}
      }
    }
    else{
      w_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][nu];
      w_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][nu];
    }
    
    if(rk_flag == 0){
      dMuXiMuNu += (w_ip1-w_i)/(delta[3]*tau);
    }
    if(rk_flag == 1){
      dMuXiMuNu += (w_i-w_im1)/(delta[3]*(tau+delta[0]));
    }
    
    //The terms due to the coordinate transform:
    if(rk_flag == 0){
      dMuXiMuNu += grid_pt->Xi[rk_flag][0][nu]/tau;
      if(nu == 0){
	dMuXiMuNu += grid_pt->Xi[rk_flag][3][3]/tau;
      }
      
      if(nu == 3){
	dMuXiMuNu += grid_pt->Xi[rk_flag][3][0]/tau;
      }
    }
    if(rk_flag == 1){
      dMuXiMuNu += grid_pt->Xi[rk_flag][0][nu]/(tau+delta[0]);
      if(nu == 0){
	dMuXiMuNu += grid_pt->Xi[rk_flag][3][3]/(tau+delta[0]);
      }
      
      if(nu == 3){
	dMuXiMuNu += grid_pt->Xi[rk_flag][3][0]/(tau+delta[0]);
      }
    }
  }
  else{
    for(int i=1; i<=2; i++){
      w_i = grid_pt->Xi[rk_flag][i][nu];
      
      //Trying periodic boundary conditions:
      //if(grid_pt->position[i] == nmax[i]){
      //w_ip1 = w_i ;
      //w_im1 = grid_pt->nbr_m_1[i]->deltaW[rk_flag][i][nu];
      //}
      //else if(grid_pt->position[i] == 0){
      //w_ip1 = grid_pt->nbr_p_1[i]->deltaW[rk_flag][i][nu];
      //w_im1 = w_i;
      //}
      //else{
      w_ip1 = grid_pt->nbr_p_1[i]->Xi[rk_flag][i][nu];
      w_im1 = grid_pt->nbr_m_1[i]->Xi[rk_flag][i][nu];
      //}
      
      if(rk_flag == 0){
	dMuXiMuNu += (w_ip1-w_i)/delta[i];
      }
      if(rk_flag == 1){
	dMuXiMuNu += (w_i-w_im1)/delta[i];
      }
    }
    
    w_i = grid_pt->Xi[rk_flag][3][nu];
    
    if(grid_pt->position[3] == nmax[3]){
      if(size == 1){
	//w_ip1 = grid_pt->nbr_p_1[3]->deltaW[rk_flag][3][nu];
	w_ip1 = w_i;
	w_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][nu];
      }
      else{
	if(rank == size-1){
	  w_ip1 = w_i;
	  w_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][nu];
	}
	else{
	  w_ip1 = Rneighbor->Xi[rk_flag][3][nu];
	  w_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][nu];
	}
      }
    }
    else if(grid_pt->position[3] == 0){
      if(size == 1){
	w_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][nu];
	//w_im1 = grid_pt->nbr_m_1[3]->deltaW[rk_flag][3][nu];
	w_im1 = w_i;
      }
      else{
	if(rank == 0){
	  //Bug fix? -CFY
	  w_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][nu];
	  w_im1 = w_i;
	}
	else{
	  w_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][nu];
	  w_im1 = Lneighbor->Xi[rk_flag][3][nu];
	}
      }
    }
    else{
      w_ip1 = grid_pt->nbr_p_1[3]->Xi[rk_flag][3][nu];
      w_im1 = grid_pt->nbr_m_1[3]->Xi[rk_flag][3][nu];
    }
    
    if(rk_flag == 0){
      dMuXiMuNu += (w_ip1-w_i)/delta[1];
    }
    if(rk_flag == 1){
      dMuXiMuNu += (w_i-w_im1)/delta[1];
    }
  }
  
  return dMuXiMuNu;

}

//Computes D^\mu u^\nu, using rk_flag for the MacCormack method:
double Advance::getDMuUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu){
  if(DATA->rk_order != 2){
    cerr << "rk_order != 2! Exiting..." << endl;
    exit(1);
  }

  double delta[4];
  int nmax[4];

  nmax[0] = 0; //Should not be used. -CFY
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;

  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = DATA->delta_eta;

  double dMuUNu = 0.;
  if(mu==0){
    if(rk_flag==0){
      dMuUNu = (grid_pt->u[0][nu] - grid_pt->prev_u[0][nu])/delta[0];
    }
    else{
      dMuUNu = (grid_pt->u[rk_flag][nu] - grid_pt->u[rk_flag-1][nu])/delta[0];
    }
  }

  if(DATA->Initial_profile != 7){
    if(mu==1 || mu==2){
      if(rk_flag==0){
	if(grid_pt->position[mu]==nmax[mu]){
	  //dMuUNu = grid_pt->u[rk_flag][nu]/delta[mu];
	  dMuUNu = -(grid_pt->nbr_p_1[mu]->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/delta[mu];
	  //dMuUNu = 0.;
	  //putthisbackdMuUNu = -(grid_pt->u[rk_flag][nu] - grid_pt->nbr_m_1[mu]->u[rk_flag][nu])/delta[mu];
	}
	else{
	  dMuUNu = -(grid_pt->nbr_p_1[mu]->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/delta[mu];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[mu]==0){
	  //dMuUNu = -grid_pt->u[rk_flag][nu]/delta[mu];
	  dMuUNu = -(grid_pt->u[rk_flag][nu] - grid_pt->nbr_m_1[mu]->u[rk_flag][nu])/delta[mu];
	  //dMuUNu = 0.;
	  //putthisbackdMuUNu = -(grid_pt->nbr_p_1[mu]->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/delta[mu];
	}
	else{
	  dMuUNu = -(grid_pt->u[rk_flag][nu] - grid_pt->nbr_m_1[mu]->u[rk_flag][nu])/delta[mu];
	}
      }
    }

    if(mu==3){
      if(rk_flag==0){
	if(grid_pt->position[3]==nmax[3]){
	  if(size == 1 || rank == size-1){
	    //dMuUNu = grid_pt->u[rk_flag][nu]/(tau*delta[3]);
	    //dMuUNu = -(grid_pt->nbr_p_1[3]->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/(tau*delta[3]);
	    //dMuUNu = 0.;
	    dMuUNu = -(grid_pt->u[rk_flag][nu] - grid_pt->nbr_m_1[3]->u[rk_flag][nu])/(tau*delta[3]);
	  }
	  else{
	    dMuUNu = -(Rneighbor->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/(tau*delta[3]);
	  }
	}
	else{
	  dMuUNu = -(grid_pt->nbr_p_1[3]->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/(tau*delta[3]);
	}
	if(nu==0){
	  dMuUNu += -grid_pt->u[rk_flag][3]/tau;
	}
	if(nu==3){
	  dMuUNu += -grid_pt->u[rk_flag][0]/tau;
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[3]==0){
	  if(size == 1 || rank == 0){
	    //dMuUNu = -grid_pt->u[rk_flag][nu]/((tau+delta[0])*delta[3]);
	    ////dMuUNu = -(grid_pt->u[rk_flag][nu] - grid_pt->nbr_m_1[3]->u[rk_flag][nu])/((tau+delta[0])*delta[3]);
	    //dMuUNu = 0.;
	    dMuUNu = -(grid_pt->nbr_p_1[3]->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/((tau+delta[0])*delta[3]);
	  }
	  else{
	    dMuUNu = -(grid_pt->u[rk_flag][nu] - Lneighbor->u[rk_flag][nu])/((tau+delta[0])*delta[3]);
	  }
	}
	else{
	  dMuUNu = -(grid_pt->u[rk_flag][nu] - grid_pt->nbr_m_1[3]->u[rk_flag][nu])/((tau+delta[0])*delta[3]);
	}
	if(nu==0){
	  dMuUNu += -grid_pt->u[rk_flag][3]/(tau+delta[0]);
	}
	if(nu==3){
	  dMuUNu += -grid_pt->u[rk_flag][0]/(tau+delta[0]);
	}
      }
    }
  }
  else{
    if(mu==1 || mu==2){
      if(rk_flag==0){
	if(grid_pt->position[mu]==nmax[mu]){
	  //Trying periodic boundary conditions:
	  dMuUNu = -(grid_pt->nbr_p_1[mu]->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/delta[mu];
	  //dMuUNu = 0.;
	}
	else{
	  dMuUNu = -(grid_pt->nbr_p_1[mu]->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/delta[mu];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[mu]==0){
	  //Trying periodic boundary conditions:
	  dMuUNu = -(grid_pt->u[rk_flag][nu] - grid_pt->nbr_m_1[mu]->u[rk_flag][nu])/delta[mu];
	  //dMuUNu = 0.;
	}
	else{
	  dMuUNu = -(grid_pt->u[rk_flag][nu] - grid_pt->nbr_m_1[mu]->u[rk_flag][nu])/delta[mu];
	}
      }
    }

    if(mu==3){
      if(rk_flag==0){
	if(grid_pt->position[3]==nmax[3]){
	  if(size == 1 || rank == size-1){
	    //dMuUNu = -(grid_pt->nbr_p_1[3]->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/(tau*delta[3]);
	    dMuUNu = 0.;
	  }
	  else{
	    dMuUNu = -(Rneighbor->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/delta[1];
	  }
	}
	else{
	  dMuUNu = -(grid_pt->nbr_p_1[3]->u[rk_flag][nu] - grid_pt->u[rk_flag][nu])/delta[1];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[3]==0){
	  if(size == 1 || rank == 0){
	    //dMuUNu = -(grid_pt->u[rk_flag][nu] - grid_pt->nbr_m_1[3]->u[rk_flag][nu])/((tau+delta[0])*delta[3]);
	    dMuUNu = 0.;
	  }
	  else{
	    dMuUNu = -(grid_pt->u[rk_flag][nu] - Lneighbor->u[rk_flag][nu])/delta[1];
	  }
	}
	else{
	  dMuUNu = -(grid_pt->u[rk_flag][nu] - grid_pt->nbr_m_1[3]->u[rk_flag][nu])/delta[1];
	}
      }
    }
  }

  return dMuUNu;
}

//Computes \partial^\mu \delta u^\nu, using rk_flag for the MacCormack method:
double Advance::getDMuDeltaUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu){
  if(DATA->rk_order != 2){
    cerr << "rk_order != 2! Exiting..." << endl;
    exit(1);
  }

  double delta[4];
  int nmax[4];

  nmax[0] = 0; //Should not be used. -CFY
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;

  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = DATA->delta_eta;

  double dMuDeltaUNu = 0.;
  if(mu==0){
    if(rk_flag==0){
      dMuDeltaUNu = (grid_pt->deltaU[0][nu] - grid_pt->prev_deltaU[nu])/delta[0];
    }
    else{
      dMuDeltaUNu = (grid_pt->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag-1][nu])/delta[0];
    }
  }

  if(DATA->Initial_profile != 7){
    if(mu==1 || mu==2){
      if(rk_flag==0){
	if(grid_pt->position[mu]==nmax[mu]){
	  //dMuDeltaUNu = 0.;
	  dMuDeltaUNu = -(grid_pt->nbr_p_1[mu]->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/delta[mu];
	  //dMuDeltaUNu = grid_pt->deltaU[rk_flag][nu]/delta[mu];
	  //putthisbackdMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - grid_pt->nbr_m_1[mu]->deltaU[rk_flag][nu])/delta[mu];
	}
	else{
	  dMuDeltaUNu = -(grid_pt->nbr_p_1[mu]->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/delta[mu];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[mu]==0){
	  //dMuDeltaUNu = 0.;
	  dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - grid_pt->nbr_m_1[mu]->deltaU[rk_flag][nu])/delta[mu];
	  //dMuDeltaUNu = -grid_pt->deltaU[rk_flag][nu]/delta[mu];
	  //putthisbackdMuDeltaUNu = -(grid_pt->nbr_p_1[mu]->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/delta[mu];
	}
	else{
	  dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - grid_pt->nbr_m_1[mu]->deltaU[rk_flag][nu])/delta[mu];
	}
      }
    }

    if(mu==3){
      if(rk_flag==0){
	if(grid_pt->position[3]==nmax[3]){
	  if(size == 1 || rank == size-1){
	    //dMuDeltaUNu = grid_pt->deltaU[rk_flag][nu]/(tau*delta[3]);
	    ////dMuDeltaUNu = -(grid_pt->nbr_p_1[3]->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/(tau*delta[3]);
	    //dMuDeltaUNu = 0.;
	    dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - grid_pt->nbr_m_1[3]->deltaU[rk_flag][nu])/(tau*delta[3]);
	  }
	  else{
	    dMuDeltaUNu = -(Rneighbor->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/(tau*delta[3]);
	  }
	}
	else{
	  dMuDeltaUNu = -(grid_pt->nbr_p_1[3]->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/(tau*delta[3]);
	}
	if(nu==0){
	  dMuDeltaUNu += -grid_pt->deltaU[rk_flag][3]/tau;
	}
	if(nu==3){
	  dMuDeltaUNu += -grid_pt->deltaU[rk_flag][0]/tau;
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[3]==0){
	  if(size == 1 || rank == 0){
	    //dMuDeltaUNu = -grid_pt->deltaU[rk_flag][nu]/((tau+delta[0])*delta[3]);
	    ////dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - grid_pt->nbr_m_1[3]->deltaU[rk_flag][nu])/((tau+delta[0])*delta[3]);
	    //dMuDeltaUNu = 0.;
	    dMuDeltaUNu = -(grid_pt->nbr_p_1[3]->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/((tau+delta[0])*delta[3]);
	  }
	  else{
	    dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - Lneighbor->deltaU[rk_flag][nu])/((tau+delta[0])*delta[3]);
	  }
	}
	else{
	  dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - grid_pt->nbr_m_1[3]->deltaU[rk_flag][nu])/((tau+delta[0])*delta[3]);
	}
	if(nu==0){
	  dMuDeltaUNu += -grid_pt->deltaU[rk_flag][3]/(tau+delta[0]);
	}
	if(nu==3){
	  dMuDeltaUNu += -grid_pt->deltaU[rk_flag][0]/(tau+delta[0]);
	}
      }
    }
  }
  else{
    if(mu==1 || mu==2){
      if(rk_flag==0){
	if(grid_pt->position[mu]==nmax[mu]){
	  //dMuDeltaUNu = 0.;
	  //Trying periodic boundary conditions:
	  dMuDeltaUNu = -(grid_pt->nbr_p_1[mu]->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/delta[mu];
	}
	else{
	  dMuDeltaUNu = -(grid_pt->nbr_p_1[mu]->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/delta[mu];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[mu]==0){
	  //dMuDeltaUNu = 0.;
	  //Trying periodic boundary conditions:
	  dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - grid_pt->nbr_m_1[mu]->deltaU[rk_flag][nu])/delta[mu];
	}
	else{
	  dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - grid_pt->nbr_m_1[mu]->deltaU[rk_flag][nu])/delta[mu];
	}
      }
    }

    if(mu==3){
      if(rk_flag==0){
	if(grid_pt->position[3]==nmax[3]){
	  if(size == 1 || rank == size-1){
	    //dMuDeltaUNu = -(grid_pt->nbr_p_1[3]->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/(tau*delta[3]);
	    dMuDeltaUNu = 0.;
	  }
	  else{
	    dMuDeltaUNu = -(Rneighbor->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/delta[1];
	  }
	}
	else{
	  dMuDeltaUNu = -(grid_pt->nbr_p_1[3]->deltaU[rk_flag][nu] - grid_pt->deltaU[rk_flag][nu])/delta[1];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[3]==0){
	  if(size == 1 || rank == 0){
	    //dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - grid_pt->nbr_m_1[3]->deltaU[rk_flag][nu])/((tau+delta[0])*delta[3]);
	    dMuDeltaUNu = 0.;
	  }
	  else{
	    dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - Lneighbor->deltaU[rk_flag][nu])/delta[1];
	  }
	}
	else{
	  dMuDeltaUNu = -(grid_pt->deltaU[rk_flag][nu] - grid_pt->nbr_m_1[3]->deltaU[rk_flag][nu])/delta[1];
	}
      }
    }
  }

  return dMuDeltaUNu;
}

//Computes D_\mu \delta u^\mu, using rk_flag for the MacCormack method:
double Advance::getDMuUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank){
  if(DATA->rk_order != 2){
    cerr << "rk_order != 2! Exiting..." << endl;
    exit(1);
  }

  double delta[4];
  int nmax[4];

  nmax[0] = 0; //Should not be used. -CFY
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;

  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = DATA->delta_eta;

  double dMuUMu;
  if(rk_flag==0){
    dMuUMu = (grid_pt->u[0][0] - grid_pt->prev_u[0][0])/delta[0];
  }
  else{
    dMuUMu = (grid_pt->u[rk_flag][0] - grid_pt->u[rk_flag-1][0])/delta[0];
  }

  if(DATA->Initial_profile != 7){
    for(int mu=1; mu<3; mu++){
      if(rk_flag==0){
	if(grid_pt->position[mu]==nmax[mu]){
	  dMuUMu += (grid_pt->nbr_p_1[mu]->u[rk_flag][mu] - grid_pt->u[rk_flag][mu])/delta[mu];
	  //dMuUMu += -grid_pt->u[rk_flag][mu]/delta[mu];
	  //dMuUMu += 0.;
	  //putthisbackdMuUMu += (grid_pt->u[rk_flag][mu] - grid_pt->nbr_m_1[mu]->u[rk_flag][mu])/delta[mu];
	}
	else{
	  dMuUMu += (grid_pt->nbr_p_1[mu]->u[rk_flag][mu] - grid_pt->u[rk_flag][mu])/delta[mu];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[mu]==0){
	  dMuUMu += (grid_pt->u[rk_flag][mu] - grid_pt->nbr_m_1[mu]->u[rk_flag][mu])/delta[mu];
	  //dMuUMu += grid_pt->u[rk_flag][mu]/delta[mu];
	  //dMuUMu += 0.;
	  //putthisbackdMuUMu += (grid_pt->nbr_p_1[mu]->u[rk_flag][mu] - grid_pt->u[rk_flag][mu])/delta[mu];
	}
	else{
	  dMuUMu += (grid_pt->u[rk_flag][mu] - grid_pt->nbr_m_1[mu]->u[rk_flag][mu])/delta[mu];
	}
      }
    }

    if(rk_flag==0){
      if(grid_pt->position[3]==nmax[3]){
	if(size == 1 || rank == size-1){
	  ////dMuUMu += (grid_pt->nbr_p_1[3]->u[rk_flag][3] - grid_pt->u[rk_flag][3])/(tau*delta[3]);
	  //dMuUMu += -grid_pt->u[rk_flag][3]/(tau*delta[3]);
	  //dMuUMu += 0.;
	  dMuUMu += (grid_pt->u[rk_flag][3] - grid_pt->nbr_m_1[3]->u[rk_flag][3])/(tau*delta[3]);
	}
	else{
	  dMuUMu += (Rneighbor->u[rk_flag][3] - grid_pt->u[rk_flag][3])/(tau*delta[3]);
	}
      }
      else{
	dMuUMu += (grid_pt->nbr_p_1[3]->u[rk_flag][3] - grid_pt->u[rk_flag][3])/(tau*delta[3]);
      }
      dMuUMu += grid_pt->u[rk_flag][0]/tau;
    }
    if(rk_flag==1){
      if(grid_pt->position[3]==0){
	if(size == 1 || rank == 0){
	  //dMuUMu += (grid_pt->u[rk_flag][3] - grid_pt->nbr_m_1[3]->u[rk_flag][3])/((tau+delta[0])*delta[3]);
	  //dMuUMu += grid_pt->u[rk_flag][3]/((tau+delta[0])*delta[3]);
	  //dMuUMu += 0.;
	  dMuUMu += (grid_pt->nbr_p_1[3]->u[rk_flag][3] - grid_pt->u[rk_flag][3])/((tau+delta[0])*delta[3]);
	}
	else{
	  dMuUMu += (grid_pt->u[rk_flag][3] - Lneighbor->u[rk_flag][3])/((tau+delta[0])*delta[3]);
	}
      }
      else{
	dMuUMu += (grid_pt->u[rk_flag][3] - grid_pt->nbr_m_1[3]->u[rk_flag][3])/((tau+delta[0])*delta[3]);
      }
      dMuUMu += grid_pt->u[rk_flag][0]/(tau+delta[0]);
    }
  }
  else{
    for(int mu=1; mu<3; mu++){
      if(rk_flag==0){
	if(grid_pt->position[mu]==nmax[mu]){
	  //Trying periodic boundary conditions:
	  dMuUMu += (grid_pt->nbr_p_1[mu]->u[rk_flag][mu] - grid_pt->u[rk_flag][mu])/delta[mu];
	  //dMuUMu += 0.;
	}
	else{
	  dMuUMu += (grid_pt->nbr_p_1[mu]->u[rk_flag][mu] - grid_pt->u[rk_flag][mu])/delta[mu];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[mu]==0){
	  //Trying periodic boundary conditions:
	  dMuUMu += (grid_pt->u[rk_flag][mu] - grid_pt->nbr_m_1[mu]->u[rk_flag][mu])/delta[mu];
	  //dMuUMu += 0.;
	}
	else{
	  dMuUMu += (grid_pt->u[rk_flag][mu] - grid_pt->nbr_m_1[mu]->u[rk_flag][mu])/delta[mu];
	}
      }
    }

    if(rk_flag==0){
      if(grid_pt->position[3]==nmax[3]){
	if(size == 1 || rank == size-1){
	  //dMuUMu += (grid_pt->nbr_p_1[3]->u[rk_flag][3] - grid_pt->u[rk_flag][3])/(tau*delta[3]);
	  dMuUMu += 0.;
	}
	else{
	  dMuUMu += (Rneighbor->u[rk_flag][3] - grid_pt->u[rk_flag][3])/delta[1];
	}
      }
      else{
	dMuUMu += (grid_pt->nbr_p_1[3]->u[rk_flag][3] - grid_pt->u[rk_flag][3])/delta[1];
      }
    }
    if(rk_flag==1){
      if(grid_pt->position[3]==0){
	if(size == 1 || rank == 0){
	  //dMuUMu += (grid_pt->u[rk_flag][3] - grid_pt->nbr_m_1[3]->u[rk_flag][3])/((tau+delta[0])*delta[3]);
	  dMuUMu += 0.;
	}
	else{
	  dMuUMu += (grid_pt->u[rk_flag][3] - Lneighbor->u[rk_flag][3])/delta[1];
	}
      }
      else{
	dMuUMu += (grid_pt->u[rk_flag][3] - grid_pt->nbr_m_1[3]->u[rk_flag][3])/delta[1];
      }
    }
  }

  return dMuUMu;
}

//Computes D_\mu \delta u^\mu, using rk_flag for the MacCormack method:
double Advance::getDMuDeltaUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank){
  if(DATA->rk_order != 2){
    cerr << "rk_order != 2! Exiting..." << endl;
    exit(1);
  }

  double delta[4];
  int nmax[4];

  nmax[0] = 0; //Should not be used. -CFY
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;

  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = DATA->delta_eta;

  double dMuDeltaUMu;
  if(rk_flag==0){
    dMuDeltaUMu = (grid_pt->deltaU[0][0] - grid_pt->prev_deltaU[0])/delta[0];
  }
  else{
    dMuDeltaUMu = (grid_pt->deltaU[rk_flag][0] - grid_pt->deltaU[rk_flag-1][0])/delta[0];
  }

  if(DATA->Initial_profile != 7){
    for(int mu=1; mu<3; mu++){
      if(rk_flag==0){
	if(grid_pt->position[mu]==nmax[mu]){
	  //dMuDeltaUMu += 0.;
	  dMuDeltaUMu += (grid_pt->nbr_p_1[mu]->deltaU[rk_flag][mu] - grid_pt->deltaU[rk_flag][mu])/delta[mu];
	  //dMuDeltaUMu += -grid_pt->deltaU[rk_flag][mu]/delta[mu];
	  //putthisbackdMuDeltaUMu += (grid_pt->deltaU[rk_flag][mu] - grid_pt->nbr_m_1[mu]->deltaU[rk_flag][mu])/delta[mu];
	}
	else{
	  dMuDeltaUMu += (grid_pt->nbr_p_1[mu]->deltaU[rk_flag][mu] - grid_pt->deltaU[rk_flag][mu])/delta[mu];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[mu]==0){
	  //dMuDeltaUMu += 0.;
	  dMuDeltaUMu += (grid_pt->deltaU[rk_flag][mu] - grid_pt->nbr_m_1[mu]->deltaU[rk_flag][mu])/delta[mu];
	  //dMuDeltaUMu += grid_pt->deltaU[rk_flag][mu]/delta[mu];
	  //putthisbackdMuDeltaUMu += (grid_pt->nbr_p_1[mu]->deltaU[rk_flag][mu] - grid_pt->deltaU[rk_flag][mu])/delta[mu];
	}
	else{
	  dMuDeltaUMu += (grid_pt->deltaU[rk_flag][mu] - grid_pt->nbr_m_1[mu]->deltaU[rk_flag][mu])/delta[mu];
	}
      }
    }

    if(rk_flag==0){
      if(grid_pt->position[3]==nmax[3]){
	if(size == 1 || rank == size-1){
	  //dMuDeltaUMu += -grid_pt->deltaU[rk_flag][3]/(tau*delta[3]);
	  ////dMuDeltaUMu += (grid_pt->nbr_p_1[3]->deltaU[rk_flag][3] - grid_pt->deltaU[rk_flag][3])/(tau*delta[3]);
	  //dMuDeltaUMu += 0.;
	  dMuDeltaUMu += (grid_pt->deltaU[rk_flag][3] - grid_pt->nbr_m_1[3]->deltaU[rk_flag][3])/(tau*delta[3]);
	}
	else{
	  dMuDeltaUMu += (Rneighbor->deltaU[rk_flag][3] - grid_pt->deltaU[rk_flag][3])/(tau*delta[3]);
	}
      }
      else{
	dMuDeltaUMu += (grid_pt->nbr_p_1[3]->deltaU[rk_flag][3] - grid_pt->deltaU[rk_flag][3])/(tau*delta[3]);
      }
      dMuDeltaUMu += grid_pt->deltaU[rk_flag][0]/tau;
    }
    if(rk_flag==1){
      if(grid_pt->position[3]==0){
	if(size == 1 || rank == 0){
	  //dMuDeltaUMu += grid_pt->deltaU[rk_flag][3]/((tau+delta[0])*delta[3]);
	  ////dMuDeltaUMu += (grid_pt->deltaU[rk_flag][3] - grid_pt->nbr_m_1[3]->deltaU[rk_flag][3])/((tau+delta[0])*delta[3]);
	  //dMuDeltaUMu += 0.;
	  dMuDeltaUMu += (grid_pt->nbr_p_1[3]->deltaU[rk_flag][3] - grid_pt->deltaU[rk_flag][3])/((tau+delta[0])*delta[3]);
	}
	else{
	  dMuDeltaUMu += (grid_pt->deltaU[rk_flag][3] - Lneighbor->deltaU[rk_flag][3])/((tau+delta[0])*delta[3]);
	}
      }
      else{
	dMuDeltaUMu += (grid_pt->deltaU[rk_flag][3] - grid_pt->nbr_m_1[3]->deltaU[rk_flag][3])/((tau+delta[0])*delta[3]);
      }
      dMuDeltaUMu += grid_pt->deltaU[rk_flag][0]/(tau+delta[0]);
    }
  }
  else{
    for(int mu=1; mu<3; mu++){
      if(rk_flag==0){
	if(grid_pt->position[mu]==nmax[mu]){
	  //dMuDeltaUMu += 0.;
	  //Trying periodic boundary conditions:
	  dMuDeltaUMu += (grid_pt->nbr_p_1[mu]->deltaU[rk_flag][mu] - grid_pt->deltaU[rk_flag][mu])/delta[mu];
	}
	else{
	  dMuDeltaUMu += (grid_pt->nbr_p_1[mu]->deltaU[rk_flag][mu] - grid_pt->deltaU[rk_flag][mu])/delta[mu];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[mu]==0){
	  //dMuDeltaUMu += 0.;
	  //Trying periodic boundary conditions:
	  dMuDeltaUMu += (grid_pt->deltaU[rk_flag][mu] - grid_pt->nbr_m_1[mu]->deltaU[rk_flag][mu])/delta[mu];
	}
	else{
	  dMuDeltaUMu += (grid_pt->deltaU[rk_flag][mu] - grid_pt->nbr_m_1[mu]->deltaU[rk_flag][mu])/delta[mu];
	}
      }
    }

    if(rk_flag==0){
      if(grid_pt->position[3]==nmax[3]){
	if(size == 1 || rank == size-1){
	  //dMuDeltaUMu += (grid_pt->nbr_p_1[3]->deltaU[rk_flag][3] - grid_pt->deltaU[rk_flag][3])/(tau*delta[3]);
	  dMuDeltaUMu += 0.;
	}
	else{
	  dMuDeltaUMu += (Rneighbor->deltaU[rk_flag][3] - grid_pt->deltaU[rk_flag][3])/delta[1];
	}
      }
      else{
	dMuDeltaUMu += (grid_pt->nbr_p_1[3]->deltaU[rk_flag][3] - grid_pt->deltaU[rk_flag][3])/delta[1];
      }
    }
    if(rk_flag==1){
      if(grid_pt->position[3]==0){
	if(size == 1 || rank == 0){
	  //dMuDeltaUMu += (grid_pt->deltaU[rk_flag][3] - grid_pt->nbr_m_1[3]->deltaU[rk_flag][3])/((tau+delta[0])*delta[3]);
	  dMuDeltaUMu += 0.;
	}
	else{
	  dMuDeltaUMu += (grid_pt->deltaU[rk_flag][3] - Lneighbor->deltaU[rk_flag][3])/delta[1];
	}
      }
      else{
	dMuDeltaUMu += (grid_pt->deltaU[rk_flag][3] - grid_pt->nbr_m_1[3]->deltaU[rk_flag][3])/delta[1];
      }
    }
  }

  return dMuDeltaUMu;
}

double Advance::getDeltaMuUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu){
  double deltaMuUNu = getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu)
    -grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][0]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0, nu);
  for(int ii=1; ii<4; ii++){
    deltaMuUNu += grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][ii]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, nu);
  }

  return deltaMuUNu;
}

double Advance::getDeltaMuDeltaUNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu){
  double deltaMuDeltaUNu = getDMuDeltaUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu)
    -grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][0]*getDMuDeltaUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0, nu);
  for(int ii=1; ii<4; ii++){
    deltaMuDeltaUNu += grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][ii]*getDMuDeltaUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, nu);
  }

  return deltaMuDeltaUNu;
}

double Advance::getDeltaMuNu(Grid *grid_pt, int rk_flag, int mu, int nu){

  double DeltaMuNu = -grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][nu];
  if(mu==0 && nu==0){
    DeltaMuNu += 1.;
  }
  if(mu != 0 && mu==nu){
    DeltaMuNu -= 1.;
  }

  return DeltaMuNu;

}

//Here are the subroutines for determining the material derivatives needed for the projector source terms:
double Advance::getDUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu){
  double DUMu = grid_pt->u[rk_flag][0]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0, mu);
  for(int ii=1; ii<4; ii++){
    DUMu += -grid_pt->u[rk_flag][ii]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, mu);
  }

  return DUMu;
}

double Advance::getDeltaUDUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu){
  double deltaUDUMu = grid_pt->deltaU[rk_flag][0]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0, mu);
  for(int ii=1; ii<4; ii++){
    deltaUDUMu += -grid_pt->deltaU[rk_flag][ii]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, mu);
  }

  return deltaUDUMu;
}

double Advance::getDDeltaUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu){
  double DDeltaUMu = grid_pt->u[rk_flag][0]*getDMuDeltaUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0, mu);
  for(int ii=1; ii<4; ii++){
    DDeltaUMu += -grid_pt->u[rk_flag][ii]*getDMuDeltaUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, mu);
  }

  return DDeltaUMu;
}

double Advance::getDWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu){
  double DWMuNu = grid_pt->u[rk_flag][0]*getDAlphaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, 0);
  for(int ii=1; ii<4; ii++){
    DWMuNu += grid_pt->u[rk_flag][ii]*getDAlphaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, ii);
  }

  return DWMuNu;
}

double Advance::getDUWMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int  size, int rank, int mu){
  double DUWMu = getDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0)*grid_pt->Wmunu[rk_flag][0][mu];
  for(int ii=1; ii<4; ii++){
    DUWMu += -getDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii)*grid_pt->Wmunu[rk_flag][ii][mu];
  }

  return DUWMu;
}

double Advance::getDeltaUDUWMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu){
  double deltaUDUWMu = getDeltaUDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0)
    *grid_pt->Wmunu[rk_flag][0][mu];
  for(int ii=1; ii<4; ii++){
    deltaUDUWMu += -getDeltaUDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii)
      *grid_pt->Wmunu[rk_flag][ii][mu];
  }

  return deltaUDUWMu;
}

double Advance::getDDeltaUWMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu){

  double DDeltaUWMu = getDDeltaUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0)*grid_pt->Wmunu[rk_flag][0][mu];
  for(int ii=1; ii<4; ii++){
    DDeltaUWMu += -getDDeltaUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii)*grid_pt->Wmunu[rk_flag][ii][mu];
  }

  return DDeltaUWMu;
}

double Advance::getDeltaWDUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu){
  double DeltaWDUMu = getDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0)*grid_pt->deltaW[rk_flag][0][mu];
  for(int ii=1; ii<4; ii++){
    DeltaWDUMu += -getDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii)*grid_pt->deltaW[rk_flag][ii][mu];
  }

  return DeltaWDUMu;
}

double Advance::getXiDUMu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu){
  double XiDUMu = getDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0)*grid_pt->Xi[rk_flag][0][mu];
  for(int ii=1; ii<4; ii++){
    XiDUMu += -getDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii)*grid_pt->Xi[rk_flag][ii][mu];
  }

  return XiDUMu;
}

double Advance::getDeltaUWDU(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank){
  double DeltaUWDU = grid_pt->deltaU[rk_flag][0]*grid_pt->Wmunu[rk_flag][0][0]*getDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0);
  for(int ii=1; ii<4; ii++){
    DeltaUWDU += -grid_pt->deltaU[rk_flag][ii]*grid_pt->Wmunu[rk_flag][ii][0]*getDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0)
      -grid_pt->deltaU[rk_flag][0]*grid_pt->Wmunu[rk_flag][0][ii]*getDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii);
    for(int ij=1; ij<4; ij++){
      DeltaUWDU += grid_pt->deltaU[rk_flag][ii]*grid_pt->Wmunu[rk_flag][ii][ij]*getDUMu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ij);
    }
  }

  return DeltaUWDU;
}

double Advance::getProjectorSource(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, double *DUW, double *deltaUDUW, double *DDeltaUW, double *DUDeltaW){

  double pS = -grid_pt->deltaU[rk_flag][mu]*DUW[nu]
    -grid_pt->u[rk_flag][mu]*(deltaUDUW[nu] + DDeltaUW[nu] + DUDeltaW[nu])
    -grid_pt->deltaU[rk_flag][nu]*DUW[mu]
    -grid_pt->u[rk_flag][nu]*(deltaUDUW[mu] + DDeltaUW[mu] + DUDeltaW[mu])
    -grid_pt->deltaU[rk_flag][0]*getDAlphaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, 0);
  for(int ii=1; ii<4; ii++){
    pS += -grid_pt->deltaU[rk_flag][ii]*getDAlphaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, ii);
  }

  return pS;
}

double Advance::getProjectorSourceXi(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, double *DUXi){

  double pS = -grid_pt->u[rk_flag][mu]*DUXi[nu]
    -grid_pt->u[rk_flag][nu]*DUXi[mu];

  return pS;
}

double Advance::getDeltaDeltaMuNu(Grid *grid_pt, int rk_flag, int mu, int nu){

  double DeltaDeltaMuNu = -grid_pt->deltaU[rk_flag][mu]*grid_pt->u[rk_flag][nu]-grid_pt->u[rk_flag][mu]*grid_pt->deltaU[rk_flag][nu];

  return DeltaDeltaMuNu;

}

double Advance::getS(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, double dMuUMu){

  double S = -2.*dMuUMu*getDeltaMuNu(grid_pt, rk_flag, mu, nu)/3.
    +getDeltaMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu)
    +getDeltaMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, nu, mu);

  return S;
}

double Advance::getDeltaS(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, double dMuUMu, double dMuDeltaUMu){

  double S = -2.*dMuDeltaUMu*getDeltaMuNu(grid_pt, rk_flag, mu, nu)/3.

    -2.*dMuUMu*getDeltaDeltaMuNu(grid_pt, rk_flag, mu, nu)/3.

    +getDeltaMuDeltaUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu)

    +getDeltaMuDeltaUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, nu, mu)

    - grid_pt->deltaU[rk_flag][mu]*grid_pt->u[rk_flag][0]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0, nu)
    - grid_pt->deltaU[rk_flag][nu]*grid_pt->u[rk_flag][0]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0, mu)
    - grid_pt->u[rk_flag][mu]*grid_pt->deltaU[rk_flag][0]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0, nu)
    - grid_pt->u[rk_flag][nu]*grid_pt->deltaU[rk_flag][0]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, 0, mu);

  //CFY:
  for(int ii=1; ii<4; ii++){
    S += grid_pt->deltaU[rk_flag][mu]*grid_pt->u[rk_flag][ii]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, nu)
      + grid_pt->deltaU[rk_flag][nu]*grid_pt->u[rk_flag][ii]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, mu)
      + grid_pt->u[rk_flag][mu]*grid_pt->deltaU[rk_flag][ii]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, nu)
      + grid_pt->u[rk_flag][nu]*grid_pt->deltaU[rk_flag][ii]*getDMuUNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, ii, mu);
  }

  return S;
}

//Computes D_\alpha W^\mu\nu, using rk_flag for the MacCormack method:
double Advance::getDAlphaWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, int alpha){
  if(DATA->rk_order != 2){
    cerr << "rk_order != 2! Exiting..." << endl;
    exit(1);
  }

  double delta[4];
  int nmax[4];

  nmax[0] = 0; //Should not be used. -CFY
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;

  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = DATA->delta_eta;

  double dAlphaWMuNu = 0.;
  if(alpha==0){
    if(rk_flag==0){
      dAlphaWMuNu = (grid_pt->Wmunu[0][mu][nu] - grid_pt->prevWmunu[0][mu][nu])/delta[0];
    }
    else{
      dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag-1][mu][nu])/delta[0];
    }
  }

  if(DATA->Initial_profile != 7){
    if(alpha==1 || alpha==2){
      if(rk_flag==0){
	if(grid_pt->position[alpha]==nmax[alpha]){
	  dAlphaWMuNu = (grid_pt->nbr_p_1[alpha]->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaWMuNu = -grid_pt->Wmunu[rk_flag][mu][nu]/delta[alpha];
	  //putthisbackdAlphaWMuNu = 0.;
	}
	else{
	  dAlphaWMuNu = (grid_pt->nbr_p_1[alpha]->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag][mu][nu])/delta[alpha];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[alpha]==0){
	  dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->Wmunu[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaWMuNu = grid_pt->Wmunu[rk_flag][mu][nu]/delta[alpha];
	  //putthisbackdAlphaWMuNu = 0.;
	}
	else{
	  dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->Wmunu[rk_flag][mu][nu])/delta[alpha];
	}
      }
    }

    if(alpha==3){
      if(rk_flag==0){
	if(grid_pt->position[3]==nmax[3]){
	  if(size == 1 || rank == size-1){
	    ////dAlphaWMuNu = (grid_pt->nbr_p_1[3]->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag][mu][nu])/(tau*delta[3]);
	    //dAlphaWMuNu = -grid_pt->Wmunu[rk_flag][mu][nu]/(tau*delta[3]);
	    dAlphaWMuNu = 0.;
	  }
	  else{
	    dAlphaWMuNu = (Rneighbor->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag][mu][nu])/(tau*delta[3]);
	  }
	}
	else{
	  dAlphaWMuNu = (grid_pt->nbr_p_1[3]->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag][mu][nu])/(tau*delta[3]);
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[3]==0){
	  if(size == 1 || rank == 0){
	    ////dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->Wmunu[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	    //dAlphaWMuNu = grid_pt->Wmunu[rk_flag][mu][nu]/((tau+delta[0])*delta[3]);
	    dAlphaWMuNu = 0.;
	  }
	  else{
	    dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - Lneighbor->Wmunu[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	  }
	}
	else{
	  dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->Wmunu[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	}
      }
    }
    
    //Finally, the terms from the Christoffel symbols must be added:
    if(alpha==3){
      if(mu==0){
	if(rk_flag==0){
	  dAlphaWMuNu += grid_pt->Wmunu[rk_flag][3][nu]/tau;
	}
	if(rk_flag==1){
	  dAlphaWMuNu += grid_pt->Wmunu[rk_flag][3][nu]/(tau+delta[0]);
	}
      }
      if(mu==3){
	if(rk_flag==0){
	  dAlphaWMuNu += grid_pt->Wmunu[rk_flag][0][nu]/tau;
	}
	if(rk_flag==1){
	  dAlphaWMuNu += grid_pt->Wmunu[rk_flag][0][nu]/(tau+delta[0]);
	}
      }
      if(nu==0){
	if(rk_flag==0){
	  dAlphaWMuNu += grid_pt->Wmunu[rk_flag][mu][3]/tau;
	}
	if(rk_flag==1){
	  dAlphaWMuNu += grid_pt->Wmunu[rk_flag][mu][3]/(tau+delta[0]);
	}
      }
      if(nu==3){
	if(rk_flag==0){
	  dAlphaWMuNu += grid_pt->Wmunu[rk_flag][mu][0]/tau;
	}
	if(rk_flag==1){
	  dAlphaWMuNu += grid_pt->Wmunu[rk_flag][mu][0]/(tau+delta[0]);
	}
      }
    }
  }
  else{
    if(alpha==1 || alpha==2){
      if(rk_flag==0){
	if(grid_pt->position[alpha]==nmax[alpha]){
	  //Trying periodic boundary conditions:
	  dAlphaWMuNu = (grid_pt->nbr_p_1[alpha]->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaWMuNu = 0.;
	}
	else{
	  dAlphaWMuNu = (grid_pt->nbr_p_1[alpha]->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag][mu][nu])/delta[alpha];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[alpha]==0){
	  //Trying periodic boundary conditions:
	  dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->Wmunu[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaWMuNu = 0.;
	}
	else{
	  dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->Wmunu[rk_flag][mu][nu])/delta[alpha];
	}
      }
    }

    if(alpha==3){
      if(rk_flag==0){
	if(grid_pt->position[3]==nmax[3]){
	  if(size == 1 || rank == size-1){
	    //dAlphaWMuNu = (grid_pt->nbr_p_1[3]->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag][mu][nu])/(tau*delta[3]);
	    dAlphaWMuNu = 0.;
	  }
	  else{
	    dAlphaWMuNu = (Rneighbor->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag][mu][nu])/delta[1];
	  }
	}
	else{
	  dAlphaWMuNu = (grid_pt->nbr_p_1[3]->Wmunu[rk_flag][mu][nu] - grid_pt->Wmunu[rk_flag][mu][nu])/delta[1];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[3]==0){
	  if(size == 1 || rank == 0){
	    //dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->Wmunu[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	    dAlphaWMuNu = 0.;
	  }
	  else{
	    dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - Lneighbor->Wmunu[rk_flag][mu][nu])/delta[1];
	  }
	}
	else{
	  dAlphaWMuNu = (grid_pt->Wmunu[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->Wmunu[rk_flag][mu][nu])/delta[1];
	}
      }
    }
  }

  return dAlphaWMuNu;
}

//Computes (\delta u \cdot D)W^\mu \nu:
double Advance::getDeltaUDWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu){
  double deltaUDWMuNu = 0.;
  for(int ii=0; ii<4; ii++){
    deltaUDWMuNu += grid_pt->deltaU[rk_flag][ii]*getDAlphaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, ii);
  }
  return deltaUDWMuNu;
}

//Computes D_\alpha \delta W^\mu\nu, using rk_flag for the MacCormack method:
double Advance::getDAlphaDeltaWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, int alpha){
  if(DATA->rk_order != 2){
    cerr << "rk_order != 2! Exiting..." << endl;
    exit(1);
  }

  double delta[4];
  int nmax[4];

  nmax[0] = 0; //Should not be used. -CFY
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;

  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = DATA->delta_eta;

  double dAlphaDeltaWMuNu = 0.;
  if(alpha==0){
    if(rk_flag==0){
      dAlphaDeltaWMuNu = (grid_pt->deltaW[0][mu][nu] - grid_pt->prev_deltaW[mu][nu])/delta[0];
    }
    else{
      dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag-1][mu][nu])/delta[0];
    }
  }

  //warning, Warning, WARNING!!!!
  //if(DATA->Initial_profile != 7){
  if(DATA->Initial_profile != 7){
    if(alpha==1 || alpha==2){
      if(rk_flag==0){
	if(grid_pt->position[alpha]==nmax[alpha]){
	  dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[alpha]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaDeltaWMuNu = -grid_pt->deltaW[rk_flag][mu][nu]/delta[alpha];
	  //putthisbackdAlphaDeltaWMuNu = 0.;
	}
	else{
	  dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[alpha]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/delta[alpha];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[alpha]==0){
	  dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->deltaW[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaDeltaWMuNu = grid_pt->deltaW[rk_flag][mu][nu]/delta[alpha];
	  //putthisbackdAlphaDeltaWMuNu = 0.;
	}
	else{
	  dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->deltaW[rk_flag][mu][nu])/delta[alpha];
	}
      }
    }

    if(alpha==3){
      if(rk_flag==0){
	if(grid_pt->position[3]==nmax[3]){
	  if(size == 1 || rank == size-1){
	    ////dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[3]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/(tau*delta[3]);
	    //dAlphaDeltaWMuNu = -grid_pt->deltaW[rk_flag][mu][nu]/(tau*delta[3]);
	    dAlphaDeltaWMuNu = 0.;
	  }
	  else{
	    dAlphaDeltaWMuNu = (Rneighbor->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/(tau*delta[3]);
	  }
	}
	else{
	  dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[3]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/(tau*delta[3]);
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[3]==0){
	  if(size == 1 || rank == 0){
	    ////dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->deltaW[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	    //dAlphaDeltaWMuNu = grid_pt->deltaW[rk_flag][mu][nu]/((tau+delta[0])*delta[3]);
	    dAlphaDeltaWMuNu = 0.;
	  }
	  else{
	    dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - Lneighbor->deltaW[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	  }
	}
	else{
	  dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->deltaW[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	}
      }
    }
    
    //Finally, the terms from the Christoffel symbols must be added:
    if(alpha==3){
      if(mu==0){
	if(rk_flag==0){
	  dAlphaDeltaWMuNu += grid_pt->deltaW[rk_flag][3][nu]/tau;
	}
	if(rk_flag==1){
	  dAlphaDeltaWMuNu += grid_pt->deltaW[rk_flag][3][nu]/(tau+delta[0]);
	}
      }
      if(mu==3){
	if(rk_flag==0){
	  dAlphaDeltaWMuNu += grid_pt->deltaW[rk_flag][0][nu]/tau;
	}
	if(rk_flag==1){
	  dAlphaDeltaWMuNu += grid_pt->deltaW[rk_flag][0][nu]/(tau+delta[0]);
	}
      }
      if(nu==0){
	if(rk_flag==0){
	  dAlphaDeltaWMuNu += grid_pt->deltaW[rk_flag][mu][3]/tau;
	}
	if(rk_flag==1){
	  dAlphaDeltaWMuNu += grid_pt->deltaW[rk_flag][mu][3]/(tau+delta[0]);
	}
      }
      if(nu==3){
	if(rk_flag==0){
	  dAlphaDeltaWMuNu += grid_pt->deltaW[rk_flag][mu][0]/tau;
	}
	if(rk_flag==1){
	  dAlphaDeltaWMuNu += grid_pt->deltaW[rk_flag][mu][0]/(tau+delta[0]);
	}
      }
    }
  }
  else{
    if(alpha==1 || alpha==2){
      if(rk_flag==0){
	if(grid_pt->position[alpha]==nmax[alpha]){
	  //Trying periodic boundary conditions:
	  dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[alpha]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaDeltaWMuNu = 0.;
	}
	else{
	  dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[alpha]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/delta[alpha];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[alpha]==0){
	  //Trying periodic boundary conditions:
	  dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->deltaW[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaDeltaWMuNu = 0.;
	}
	else{
	  dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->deltaW[rk_flag][mu][nu])/delta[alpha];
	}
      }
    }

    if(alpha==3){
      if(rk_flag==0){
	if(grid_pt->position[3]==nmax[3]){
	  if(size == 1 || rank == size-1){
	    //dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[3]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/(tau*delta[3]);
	    dAlphaDeltaWMuNu = 0.;
	  }
	  else{
	    dAlphaDeltaWMuNu = (Rneighbor->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/delta[1];
	  }
	}
	else{
	  dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[3]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/delta[1];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[3]==0){
	  if(size == 1 || rank == 0){
	    //dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->deltaW[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	    dAlphaDeltaWMuNu = 0.;
	  }
	  else{
	    dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - Lneighbor->deltaW[rk_flag][mu][nu])/delta[1];
	  }
	}
	else{
	  dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->deltaW[rk_flag][mu][nu])/delta[1];
	}
      }
    }
  }

  return dAlphaDeltaWMuNu;
}

//Computes D_\alpha \Xi^\mu\nu, using rk_flag for the MacCormack method:
double Advance::getDAlphaXiMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu, int alpha){
  if(DATA->rk_order != 2){
    cerr << "rk_order != 2! Exiting..." << endl;
    exit(1);
  }

  double delta[4];
  int nmax[4];

  nmax[0] = 0; //Should not be used. -CFY
  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta-1;

  delta[0] = DATA->delta_tau;
  delta[1] = DATA->delta_x;
  delta[2] = DATA->delta_y;
  delta[3] = DATA->delta_eta;

  double dAlphaXiMuNu = 0.;
  if(alpha==0){
    if(rk_flag==0){
      dAlphaXiMuNu = (grid_pt->Xi[0][mu][nu] - grid_pt->prev_Xi[mu][nu])/delta[0];
    }
    else{
      dAlphaXiMuNu = (grid_pt->Xi[rk_flag][mu][nu] - grid_pt->Xi[rk_flag-1][mu][nu])/delta[0];
    }
  }

  if(DATA->Initial_profile != 7){
    if(alpha==1 || alpha==2){
      if(rk_flag==0){
	if(grid_pt->position[alpha]==nmax[alpha]){
	  ////dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[alpha]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaDeltaWMuNu = -grid_pt->deltaW[rk_flag][mu][nu]/delta[alpha];
	  dAlphaXiMuNu = 0.;
	}
	else{
	  dAlphaXiMuNu = (grid_pt->nbr_p_1[alpha]->Xi[rk_flag][mu][nu] - grid_pt->Xi[rk_flag][mu][nu])/delta[alpha];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[alpha]==0){
	  ////dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->deltaW[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaDeltaWMuNu = grid_pt->deltaW[rk_flag][mu][nu]/delta[alpha];
	  dAlphaXiMuNu = 0.;
	}
	else{
	  dAlphaXiMuNu = (grid_pt->Xi[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->Xi[rk_flag][mu][nu])/delta[alpha];
	}
      }
    }

    if(alpha==3){
      if(rk_flag==0){
	if(grid_pt->position[3]==nmax[3]){
	  if(size == 1 || rank == size-1){
	    ////dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[3]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/(tau*delta[3]);
	    //dAlphaDeltaWMuNu = -grid_pt->deltaW[rk_flag][mu][nu]/(tau*delta[3]);
	    dAlphaXiMuNu = 0.;
	  }
	  else{
	    dAlphaXiMuNu = (Rneighbor->Xi[rk_flag][mu][nu] - grid_pt->Xi[rk_flag][mu][nu])/(tau*delta[3]);
	  }
	}
	else{
	  dAlphaXiMuNu = (grid_pt->nbr_p_1[3]->Xi[rk_flag][mu][nu] - grid_pt->Xi[rk_flag][mu][nu])/(tau*delta[3]);
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[3]==0){
	  if(size == 1 || rank == 0){
	    ////dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->deltaW[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	    //dAlphaDeltaWMuNu = grid_pt->deltaW[rk_flag][mu][nu]/((tau+delta[0])*delta[3]);
	    dAlphaXiMuNu = 0.;
	  }
	  else{
	    dAlphaXiMuNu = (grid_pt->Xi[rk_flag][mu][nu] - Lneighbor->Xi[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	  }
	}
	else{
	  dAlphaXiMuNu = (grid_pt->Xi[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->Xi[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	}
      }
    }
    
    //Finally, the terms from the Christoffel symbols must be added:
    if(alpha==3){
      if(mu==0){
	if(rk_flag==0){
	  dAlphaXiMuNu += grid_pt->Xi[rk_flag][3][nu]/tau;
	}
	if(rk_flag==1){
	  dAlphaXiMuNu += grid_pt->Xi[rk_flag][3][nu]/(tau+delta[0]);
	}
      }
      if(mu==3){
	if(rk_flag==0){
	  dAlphaXiMuNu += grid_pt->Xi[rk_flag][0][nu]/tau;
	}
	if(rk_flag==1){
	  dAlphaXiMuNu += grid_pt->Xi[rk_flag][0][nu]/(tau+delta[0]);
	}
      }
      if(nu==0){
	if(rk_flag==0){
	  dAlphaXiMuNu += grid_pt->Xi[rk_flag][mu][3]/tau;
	}
	if(rk_flag==1){
	  dAlphaXiMuNu += grid_pt->Xi[rk_flag][mu][3]/(tau+delta[0]);
	}
      }
      if(nu==3){
	if(rk_flag==0){
	  dAlphaXiMuNu += grid_pt->Xi[rk_flag][mu][0]/tau;
	}
	if(rk_flag==1){
	  dAlphaXiMuNu += grid_pt->Xi[rk_flag][mu][0]/(tau+delta[0]);
	}
      }
    }
  }
  else{
    if(alpha==1 || alpha==2){
      if(rk_flag==0){
	if(grid_pt->position[alpha]==nmax[alpha]){
	  //Trying periodic boundary conditions:
	  dAlphaXiMuNu = (grid_pt->nbr_p_1[alpha]->Xi[rk_flag][mu][nu] - grid_pt->Xi[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaDeltaWMuNu = 0.;
	}
	else{
	  dAlphaXiMuNu = (grid_pt->nbr_p_1[alpha]->Xi[rk_flag][mu][nu] - grid_pt->Xi[rk_flag][mu][nu])/delta[alpha];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[alpha]==0){
	  //Trying periodic boundary conditions:
	  dAlphaXiMuNu = (grid_pt->Xi[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->Xi[rk_flag][mu][nu])/delta[alpha];
	  //dAlphaDeltaWMuNu = 0.;
	}
	else{
	  dAlphaXiMuNu = (grid_pt->Xi[rk_flag][mu][nu] - grid_pt->nbr_m_1[alpha]->Xi[rk_flag][mu][nu])/delta[alpha];
	}
      }
    }

    if(alpha==3){
      if(rk_flag==0){
	if(grid_pt->position[3]==nmax[3]){
	  if(size == 1 || rank == size-1){
	    //dAlphaDeltaWMuNu = (grid_pt->nbr_p_1[3]->deltaW[rk_flag][mu][nu] - grid_pt->deltaW[rk_flag][mu][nu])/(tau*delta[3]);
	    dAlphaXiMuNu = 0.;
	  }
	  else{
	    dAlphaXiMuNu = (Rneighbor->Xi[rk_flag][mu][nu] - grid_pt->Xi[rk_flag][mu][nu])/delta[1];
	  }
	}
	else{
	  dAlphaXiMuNu = (grid_pt->nbr_p_1[3]->Xi[rk_flag][mu][nu] - grid_pt->Xi[rk_flag][mu][nu])/delta[1];
	}
      }
      if(rk_flag==1){
	if(grid_pt->position[3]==0){
	  if(size == 1 || rank == 0){
	    //dAlphaDeltaWMuNu = (grid_pt->deltaW[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->deltaW[rk_flag][mu][nu])/((tau+delta[0])*delta[3]);
	    dAlphaXiMuNu = 0.;
	  }
	  else{
	    dAlphaXiMuNu = (grid_pt->Xi[rk_flag][mu][nu] - Lneighbor->Xi[rk_flag][mu][nu])/delta[1];
	  }
	}
	else{
	  dAlphaXiMuNu = (grid_pt->Xi[rk_flag][mu][nu] - grid_pt->nbr_m_1[3]->Xi[rk_flag][mu][nu])/delta[1];
	}
      }
    }
  }

  return dAlphaXiMuNu;
}

double Advance::getUiDiDeltaWMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu){
  double uiDiDeltaWMuNu = 0.;
  for(int ii=1; ii<4; ii++){
    uiDiDeltaWMuNu += grid_pt->u[rk_flag][ii]*getDAlphaDeltaWMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, ii);
  }

  return uiDiDeltaWMuNu;
}

double Advance::getUiDiXiMuNu(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, int size, int rank, int mu, int nu){
  double uiDiXiMuNu = 0.;
  for(int ii=1; ii<4; ii++){
    uiDiXiMuNu += grid_pt->u[rk_flag][ii]*getDAlphaXiMuNu(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, ii);
  }

  return uiDiXiMuNu;
}

double Advance::getDeltaDeltaWDeltaSDelta(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, 
					  int size, int rank, double eta, double deltaEta, double dMuUMu, double dMuDeltaUMu, 
					  int mu, int nu){
  //First, DeltaWTimesDeltaNu, a vector:
  double DeltaWDeltaSTimesDeltaNu[4];
  //Then, DeltaTimesDeltaWTimesDeltaMuNu, a single number:
  double DeltaTimesDeltaWDeltaSTimesDeltaMuNu = 0.;
  for(int ii=0; ii<4; ii++){
    DeltaWDeltaSTimesDeltaNu[ii] = 0.;
  }
  
  //tauPi:
  double s_den = eos->s_func(grid_pt->epsilon_t, grid_pt->p_t, grid_pt->rhob_t, DATA);
  double shear = (DATA->shear_to_s)*s_den;
  double tauPi = 3.0*shear/(grid_pt->epsilon_t + grid_pt->p_t);
  tauPi = max(tauPi, DATA->tau_pi);
  if(!finite(tauPi)) tauPi = DATA->tau_pi;

  if(grid_pt->T_t < DATA->fluctuatingTMin){
    //DeltaTimesDeltaWDeltaSTimesDeltaMuNu = grid_pt->deltaW[rk_flag][mu][nu]/DATA->tau_pi;
    DeltaTimesDeltaWDeltaSTimesDeltaMuNu = grid_pt->deltaW[rk_flag][mu][nu]/tauPi
      //warning, Warning, WARNING!!!!
      + (4./3.)*dMuDeltaUMu*grid_pt->Wmunu[rk_flag][mu][nu]
      + (4./3.)*dMuUMu*grid_pt->deltaW[rk_flag][mu][nu];
  }
  else{
    //DeltaTimesDeltaWDeltaSTimesDeltaMuNu = (grid_pt->deltaW[rk_flag][mu][nu]
    //					    - eta*getDeltaS(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, dMuUMu, dMuDeltaUMu)
    //					    - deltaEta*getS(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, dMuUMu) 
    //					    )/DATA->tau_pi
    //+ (4./3.)*dMuDeltaUMu*grid_pt->Wmunu[rk_flag][mu][nu]
    //+ (4./3.)*dMuUMu*grid_pt->deltaW[rk_flag][mu][nu];
    DeltaTimesDeltaWDeltaSTimesDeltaMuNu = (grid_pt->deltaW[rk_flag][mu][nu]
					    - eta*getDeltaS(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, dMuUMu, dMuDeltaUMu)
					    - deltaEta*getS(tau, DATA, grid_pt, Lneighbor, Rneighbor, rk_flag, size, rank, mu, nu, dMuUMu) 
					    )/tauPi
      + (4./3.)*dMuDeltaUMu*grid_pt->Wmunu[rk_flag][mu][nu]
      + (4./3.)*dMuUMu*grid_pt->deltaW[rk_flag][mu][nu];
  }

  return DeltaTimesDeltaWDeltaSTimesDeltaMuNu;
}

double Advance::getDeltaXiDelta(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, int rk_flag, 
				int size, int rank, double eta, double deltaEta, double dMuUMu, double dMuDeltaUMu, 
				int mu, int nu){
  double XiTimesDeltaNu[4];
  //Then, DeltaTimesDeltaWTimesDeltaMuNu, a single number:
  double DeltaTimesXiTimesDeltaMuNu = 0.;
  for(int ii=0; ii<4; ii++){
    XiTimesDeltaNu[ii] = 0.;
  }
  
  //tauPi:
  double s_den = eos->s_func(grid_pt->epsilon_t, grid_pt->p_t, grid_pt->rhob_t, DATA);
  double shear = (DATA->shear_to_s)*s_den;
  double tauPi = 3.0*shear/(grid_pt->epsilon_t + grid_pt->p_t);
  tauPi = max(tauPi, DATA->tau_pi);
  if(!finite(tauPi)) tauPi = DATA->tau_pi;

  if(grid_pt->T_t < DATA->fluctuatingTMin){
    DeltaTimesXiTimesDeltaMuNu = grid_pt->Xi[rk_flag][mu][nu]/tauPi
      //warning, Warning, WARNING!!!!
      + (4./3.)*dMuUMu*grid_pt->Xi[rk_flag][mu][nu];
  }
  else{
    DeltaTimesXiTimesDeltaMuNu = grid_pt->Xi[rk_flag][mu][nu]/tauPi
      + (4./3.)*dMuUMu*grid_pt->Xi[rk_flag][mu][nu];
  }

  return DeltaTimesXiTimesDeltaMuNu;
}

void Advance::reconstructDeltaPAndDeltaU(double tau, InitData *DATA, Grid *grid_pt, int rk_flag){
  if(grid_pt->deltaT[rk_flag+1][0][0] == 0. && grid_pt->deltaT[rk_flag+1][0][1] == 0.
     && grid_pt->deltaT[rk_flag+1][0][2] == 0. && grid_pt->deltaT[rk_flag+1][0][3] == 0.){
    grid_pt->deltaP[rk_flag+1] = 0.;
    for(int ii=0; ii<4; ii++){
      grid_pt->deltaU[rk_flag+1][ii] = 0.;
      for(int ij=0; ij<4; ij++){
	grid_pt->deltaT[rk_flag+1][ii][ij] = 0.;
      }
    }
  }
  else{
    double dedp;
    if(grid_pt->epsilon_t < 0.03){
      dedp = 1./eos->get_dpOverde2(0.03, grid_pt->rhob_t);
    }
    else{
      dedp = 1./eos->get_dpOverde2(grid_pt->epsilon_t, grid_pt->rhob_t);
    }
    double u0, u1, u2, u3;
    u1 = grid_pt->u[rk_flag+1][1];
    if(u1 > 1.5) u1 = 1.5;
    if(u1 < -1.5) u1 = -1.5;
    u2 = min(1.5, grid_pt->u[rk_flag+1][2]);
    if(u2 > 1.5) u2 = 1.5;
    if(u2 < -1.5) u2 = -1.5;
    u3 = min(1.5, grid_pt->u[rk_flag+1][3]);
    if(u3 > 1.5) u3 = 1.5;
    if(u3 < -1.5) u3 = -1.5;
    u0 = sqrt(1. + u1*u1 + u2*u2 + u3*u3);

    double A = (1. + dedp)*u0;
    if(A>100.){
      cerr << "At (" << grid_pt->position[1] << ", " << grid_pt->position[2] << ", " 
	   << grid_pt->position[3] << "), A = " << A << ", why?" << endl;
      cerr << "u = (" << u0 << ", " << u1 << ", " << u2 << ", " << u3 << ")." << endl;
      cerr << "dedp = " << dedp << "." << endl;
    }
    double B;
    if(grid_pt->epsilon_t < 0.03){
      B = 0.03 + eos->p_func(0.03, 0.);
    }
    else{
      B = grid_pt->epsilon_t + grid_pt->p_t;
    }

    double dT00 = grid_pt->deltaT[rk_flag+1][0][0];
    double uMudT0Mu = u0*grid_pt->deltaT[rk_flag+1][0][0]
      -u1*grid_pt->deltaT[rk_flag+1][0][1]
      -u2*grid_pt->deltaT[rk_flag+1][0][2]
      -u3*grid_pt->deltaT[rk_flag+1][0][3];
    
    double dX[4];

    dX[0] = (dT00 - 2.*u0*uMudT0Mu)/((1.-dedp)*u0*u0 - 1.);
    double deltaU0 = (dT00 - ((1.+dedp)*u0*u0 - 1.)*dX[0])/(2.*B*u0);
    dX[1] = (grid_pt->deltaT[rk_flag+1][0][1] - dX[0]*(1.+dedp)*u0*u1 - B*deltaU0*u1)/(B*u0);
    dX[2] = (grid_pt->deltaT[rk_flag+1][0][2] - dX[0]*(1.+dedp)*u0*u2 - B*deltaU0*u2)/(B*u0);
    dX[3] = (grid_pt->deltaT[rk_flag+1][0][3] - dX[0]*(1.+dedp)*u0*u3 - B*deltaU0*u3)/(B*u0);

    //The values of deltaT from these values of the variables:
    double dT0[4];
    dT0[0] = dX[0]*(-1. + A*u0)
      + 2.*B*u0*deltaU0;
    dT0[1] = dX[0]*A*u1
      + B*(u0*dX[1] + u1*deltaU0);
    dT0[2] = dX[0]*A*u2
      + B*(u0*dX[2] + u2*deltaU0);
    dT0[3] = dX[0]*A*u3
      + B*(u0*dX[3] + u3*deltaU0);
    
    double error = fabs(dT0[0] - grid_pt->deltaT[rk_flag+1][0][0])
      + fabs(dT0[1] - grid_pt->deltaT[rk_flag+1][0][1])
      + fabs(dT0[2] - grid_pt->deltaT[rk_flag+1][0][2])
      + fabs(dT0[3] - grid_pt->deltaT[rk_flag+1][0][3]);

    if(error > 1.e-8){
      cerr << "error = " << error << endl;
      cerr << "A = " << A << ", B = " << B << endl;
    }

    //Update the values of deltaP and deltaU:
    grid_pt->deltaP[rk_flag+1] = dX[0];
    grid_pt->deltaU[rk_flag+1][0] = deltaU0;
    grid_pt->deltaU[rk_flag+1][1] = dX[1];
    grid_pt->deltaU[rk_flag+1][2] = dX[2];
    grid_pt->deltaU[rk_flag+1][3] = dX[3];

    //Finally, reconstruct deltaT:
    double u0Vec[4]; u0Vec[0] = u0; u0Vec[1] =u1; u0Vec[2] =u2; u0Vec[3] =u3;
    for(int i=1; i<4; i++){
      grid_pt->deltaT[rk_flag+1][i][0] = grid_pt->deltaT[rk_flag+1][0][i];
      for(int j=1; j<4; j++){
	grid_pt->deltaT[rk_flag+1][i][j] = (1.+dedp)*dX[0]*u0Vec[i]*u0Vec[j]
	  + B*(u0Vec[i]*dX[j] + u0Vec[j]*dX[i]);
	if(i==j){
	  grid_pt->deltaT[rk_flag+1][i][j] += dX[0];
	}
      }
    }
  }

  return;
}
