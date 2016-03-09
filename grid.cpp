#include "util.h"
#include "grid.h"
#include <sstream>
#include <string>

using namespace std;

Grid *Grid::grid_v_malloc(int n1)
{
 Grid *d1_ptr;
// int i;

    /* pointer to the n1 array */
    d1_ptr = (Grid *) malloc (sizeof(Grid )*n1);

 return d1_ptr;
}/* grid_v_malloc */


Grid **Grid::grid_m_malloc(int n1, int n2)
{
    int i;// ,j;
    Grid **d1_ptr, *tmp_ptr;

    tmp_ptr = (Grid *)malloc(sizeof(Grid)*n1*n2);
    d1_ptr = (Grid **) malloc (sizeof(Grid *)*n1);

    for(i=0; i<n1; i++) 
     {
      d1_ptr[i] = &(tmp_ptr[i*n2]);
     }
    
return d1_ptr;
}/* grid_m_malloc */


Grid ***Grid::grid_c_malloc(int n1, int n2, int n3)
{
    int i,j,inc;//k,
    Grid ***d1_ptr, *tmp_ptr;

    tmp_ptr = (Grid *) malloc(sizeof(Grid)*n1*n2*n3);

    /* pointer to the n1*n2*n3 memory */

    d1_ptr = (Grid ***) malloc (sizeof(Grid **)*n1);

    for(i=0; i<n1; i++) 
     {
      d1_ptr[i] = (Grid **) malloc (sizeof(Grid *)*n2);
     } 
    
    for(i=0; i<n1; i++)
    {
     for(j=0; j<n2; j++) 
      {
       inc = n2*n3*i + n3*j;
       d1_ptr[i][j] = &(tmp_ptr[inc]);
      }
    }
    
    return d1_ptr;
}/* grid_c_malloc */

// output for e.g. MARTINI on a tau-x-y-z-grid: T, QGP fraction, flow velocities
// note that this is really z, not eta.
// also note that MARTINI's Hydro:xmax etc. are 0.5 times the values in the mpihydro input file, i.e.,
// MARTINI has the range -xmax to +xmax and mpihydro has -X_grid_size_in_fm/2 to +X_grid_size_in_fm/2 !! 
// only output one quarter of the symmetric solution
void Grid::OutputEvolutionDataXYZ(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
{
  Util *util;
  util = new Util();

  int position;
  int ix, iy, ieta, nx, ny, neta;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;

  /// MPI send and receive

  int sizeOfData = (nx+1)*(ny+1)*(neta);
  int to, from;
  
  double *eps;
  double *rhob;
  double *utau;
  double *ux;
  double *uy;
  double *ueta;
  double *Txx;
  double *Txy;
  double *Tyy;

  eps = (double *)malloc(sizeof(double)*sizeOfData);
  rhob = (double *)malloc(sizeof(double)*sizeOfData);
  utau = (double *)malloc(sizeof(double)*sizeOfData);
  ux = (double *)malloc(sizeof(double)*sizeOfData);
  uy = (double *)malloc(sizeof(double)*sizeOfData);
  ueta = (double *)malloc(sizeof(double)*sizeOfData);
  Txx = (double *)malloc(sizeof(double)*sizeOfData);
  Txy = (double *)malloc(sizeof(double)*sizeOfData);
  Tyy = (double *)malloc(sizeof(double)*sizeOfData);

  if (rank>0) //send all to rank 0
    {
      to = 0;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  position = ieta+(neta*(ix + (nx*iy)));
		  eps[position] = arena[ix][iy][ieta].epsilon;
		  rhob[position] = arena[ix][iy][ieta].rhob;
		  utau[position] = arena[ix][iy][ieta].u[0][0];
		  ux[position] = arena[ix][iy][ieta].u[0][1];
		  uy[position] = arena[ix][iy][ieta].u[0][2];
		  ueta[position] = arena[ix][iy][ieta].u[0][3];
		  if(DATA->fluctuatingHydroFlag == 2){
		    //Txx[position] = arena[ix][iy][ieta].TJb[0][1][1] + arena[ix][iy][ieta].deltaT[0][1][1];
		    //Txy[position] = arena[ix][iy][ieta].TJb[0][1][2] + arena[ix][iy][ieta].deltaT[0][1][2];
		    //Tyy[position] = arena[ix][iy][ieta].TJb[0][2][2] + arena[ix][iy][ieta].deltaT[0][2][2];
		    Txx[position] = arena[ix][iy][ieta].TJb[0][1][1] + arena[ix][iy][ieta].Wmunu[0][1][1]
		      + arena[ix][iy][ieta].deltaT[0][1][1] + arena[ix][iy][ieta].deltaW[0][1][1];
		    Txy[position] = arena[ix][iy][ieta].TJb[0][1][2] + arena[ix][iy][ieta].Wmunu[0][1][2]
		      + arena[ix][iy][ieta].deltaT[0][1][2] + arena[ix][iy][ieta].deltaW[0][1][2];
		    Tyy[position] = arena[ix][iy][ieta].TJb[0][2][2] + arena[ix][iy][ieta].Wmunu[0][2][2]
		      + arena[ix][iy][ieta].deltaT[0][2][2] + arena[ix][iy][ieta].deltaW[0][2][2];
		  }
		  else{
		    //Txx[position] = arena[ix][iy][ieta].TJb[0][1][1];
		    //Txy[position] = arena[ix][iy][ieta].TJb[0][1][2];
		    //Tyy[position] = arena[ix][iy][ieta].TJb[0][2][2];
		    Txx[position] = arena[ix][iy][ieta].TJb[0][1][1] + arena[ix][iy][ieta].Wmunu[0][1][1];
		    Txy[position] = arena[ix][iy][ieta].TJb[0][1][2] + arena[ix][iy][ieta].Wmunu[0][1][2];
		    Tyy[position] = arena[ix][iy][ieta].TJb[0][2][2] + arena[ix][iy][ieta].Wmunu[0][2][2];
		  }
		}
	    }
	}
      //      cout << "eps[5]=" <<eps[5]<< endl;
    
      MPI::COMM_WORLD.Send(eps,sizeOfData,MPI::DOUBLE,to,1);
      MPI::COMM_WORLD.Send(rhob,sizeOfData,MPI::DOUBLE,to,2);
      MPI::COMM_WORLD.Send(utau,sizeOfData,MPI::DOUBLE,to,3);
      MPI::COMM_WORLD.Send(ux,sizeOfData,MPI::DOUBLE,to,4);
      MPI::COMM_WORLD.Send(uy,sizeOfData,MPI::DOUBLE,to,5);
      MPI::COMM_WORLD.Send(ueta,sizeOfData,MPI::DOUBLE,to,6);
      MPI::COMM_WORLD.Send(Txx,sizeOfData,MPI::DOUBLE,to,7);
      MPI::COMM_WORLD.Send(Txy,sizeOfData,MPI::DOUBLE,to,8);
      MPI::COMM_WORLD.Send(Tyy,sizeOfData,MPI::DOUBLE,to,9);
    }
  
  if (rank==0) 
    {
      double ***epsFrom;
      double ***rhobFrom;
      double ***utauFrom;
      double ***uxFrom;
      double ***uyFrom;
      double ***uetaFrom;
      double ***TxxFrom;
      double ***TxyFrom;
      double ***TyyFrom;
      
      epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      utauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      TxxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      TxyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      TyyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  epsFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon;
		  rhobFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
		  utauFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][0];
		  uxFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][1];
		  uyFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][2];
		  uetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][3];
		  if(DATA->fluctuatingHydroFlag == 2){
		    //TxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][1] + arena[ix][iy][ieta].deltaT[0][1][1];
		    //TxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][2] + arena[ix][iy][ieta].deltaT[0][1][2];
		    //TyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][2][2] + arena[ix][iy][ieta].deltaT[0][2][2];
		    TxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][1] + arena[ix][iy][ieta].Wmunu[0][1][1]
		      + arena[ix][iy][ieta].deltaT[0][1][1] + arena[ix][iy][ieta].deltaW[0][1][1];
		    TxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][2] + arena[ix][iy][ieta].Wmunu[0][1][2]
		      + arena[ix][iy][ieta].deltaT[0][1][2] + arena[ix][iy][ieta].deltaW[0][1][2];
		    TyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][2][2] + arena[ix][iy][ieta].Wmunu[0][2][2]
		      + arena[ix][iy][ieta].deltaT[0][2][2] + arena[ix][iy][ieta].deltaW[0][2][2];
		  }
		  else{
		    //TxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][1];
		    //TxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][2];
		    //TyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][2][2];
		    TxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][1] + arena[ix][iy][ieta].Wmunu[0][1][1];
		    TxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][2] + arena[ix][iy][ieta].Wmunu[0][1][2];
		    TyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][2][2] + arena[ix][iy][ieta].Wmunu[0][2][2];
		  }
		}
	    }
	}
	

      for (int irank=1; irank<size; irank++)
	{
	  from = irank;
	  MPI::COMM_WORLD.Recv(eps,sizeOfData,MPI::DOUBLE,from,1);
	  MPI::COMM_WORLD.Recv(rhob,sizeOfData,MPI::DOUBLE,from,2);
	  MPI::COMM_WORLD.Recv(utau,sizeOfData,MPI::DOUBLE,from,3);
	  MPI::COMM_WORLD.Recv(ux,sizeOfData,MPI::DOUBLE,from,4);
	  MPI::COMM_WORLD.Recv(uy,sizeOfData,MPI::DOUBLE,from,5);
	  MPI::COMM_WORLD.Recv(ueta,sizeOfData,MPI::DOUBLE,from,6);
	  MPI::COMM_WORLD.Recv(Txx,sizeOfData,MPI::DOUBLE,from,7);
	  MPI::COMM_WORLD.Recv(Txy,sizeOfData,MPI::DOUBLE,from,8);
	  MPI::COMM_WORLD.Recv(Tyy,sizeOfData,MPI::DOUBLE,from,9);

	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		      position = ieta+(neta*(ix + (nx*iy)));
		      epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		      rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		      utauFrom[ix][iy][ieta+irank*neta] = utau[position];
		      uxFrom[ix][iy][ieta+irank*neta] = ux[position];
		      uyFrom[ix][iy][ieta+irank*neta] = uy[position];
		      uetaFrom[ix][iy][ieta+irank*neta] = ueta[position];
		      TxxFrom[ix][iy][ieta+irank*neta] = Txx[position];
		      TxyFrom[ix][iy][ieta+irank*neta] = Txy[position];
		      TyyFrom[ix][iy][ieta+irank*neta] = Tyy[position];
		    }
		}
	    }
	  //	  cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;

      stringstream numIn;
      numIn << DATA->seed;
      FILE *out_file;
      string out_name;
      out_name = "evolution_"+numIn.str()+".dat";
      FILE *ecc_file;
      string ecc_name;
      ecc_name = "momentumEccentricity_"+numIn.str()+".dat";
      //char* out_name = "evolution.dat";
      //out_file = fopen(out_name, "a");
      out_file = fopen(out_name.c_str(), "a");
      fprintf(out_file,"");
      ecc_file = fopen(ecc_name.c_str(), "a");
      fprintf(ecc_file,"");
      int iz, nz;
      double T, x, y, z, eta, delta_z, z_size, eps, etafrac;
      double ux, uy, ueta, utau, epsilon, rhob, QGPfrac;
      double eta_lower, eps_lower, eps_higher, rhob_lower, rhob_higher;
      double utau_lower, utau_higher, ux_lower, ux_higher, uy_lower, uy_higher, ueta_lower, ueta_higher, u0, uz;
      nz = 160;
      delta_z = 0.25;
      z_size = 40.; // z_size is fixed to 40 fm: corresponds to Hydro:zmax=20 in MARTINI!
      for(iz=0; iz<nz; iz++)
	{
	  // get temperature value and others by interpolation
	  z = iz*(delta_z) - (z_size)/2.0;
	  eta = asinh(z/tau);
	  //cout << "eta_size=" << (DATA->eta_size) << endl;
	  ieta = floor((eta + (DATA->eta_size)/2.0)/(DATA->delta_eta));
	  //cout << "ieta = " << ieta << endl;
	  eta_lower = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
	  etafrac = (eta-eta_lower)/DATA->delta_eta;
	  
	  for(iy=DATA->ny/2; iy<DATA->ny; iy++) // only do positive y
	    {
	      y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	      for(ix=DATA->nx/2; ix<DATA->nx; ix++) // only do positive x
		{
		  x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
		  eps_lower = epsFrom[ix][iy][ieta];
		  rhob_lower = rhobFrom[ix][iy][ieta];
		  utau_lower = utauFrom[ix][iy][ieta];
		  ux_lower = uxFrom[ix][iy][ieta];
		  uy_lower = uyFrom[ix][iy][ieta];
		  ueta_lower = uetaFrom[ix][iy][ieta];
		  
		  if (ieta+1 < DATA->neta*size)
		    {
		      eps_higher = epsFrom[ix][iy][ieta+1];
		      rhob_higher = rhobFrom[ix][iy][ieta+1];
		      utau_higher = utauFrom[ix][iy][ieta+1];
		      ux_higher = uxFrom[ix][iy][ieta+1];
		      uy_higher = uyFrom[ix][iy][ieta+1];
		      ueta_higher = uetaFrom[ix][iy][ieta+1];
		    }
		  else
		    {
		      eps_higher = eps_lower;
		      rhob_higher = rhob_lower;
		      utau_higher = utau_lower;
		      ux_higher = ux_lower;
		      uy_higher = uy_lower;
		      ueta_higher = ueta_lower;
		    }
		  
		  eps = eps_lower * (1.-etafrac) + (etafrac)*eps_higher;
		  rhob = rhob_lower * (1.-etafrac) + (etafrac)*rhob_higher;
		  utau = utau_lower * (1.-etafrac) + (etafrac)*utau_higher;
		  ux = ux_lower * (1.-etafrac) + (etafrac)*ux_higher;
		  uy = uy_lower * (1.-etafrac) + (etafrac)*uy_higher;
		  ueta = ueta_lower * (1.-etafrac) + (etafrac)*ueta_higher;
		  
		  if (DATA->whichEOS==1)
		    {
		      T = eos->interpolate(eps, rhob, 0);
		      QGPfrac=(eps*hbarc-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
		      //cout << T << " " << QGPfrac << endl;
		      if (QGPfrac>1.) QGPfrac = 1;
		      else if (QGPfrac<0.) QGPfrac=0.;
		    }
		  else if (DATA->whichEOS==2)
		    {
		      T = eos->interpolate2(eps, rhob, 1);
		      QGPfrac = eos->interpolate2(eps, rhob, 3);
		    }
		  
		  //Now these are the flow velocities as e,g, MARTINI needs them
		  u0 = ueta*sinh(eta)+utau*cosh(eta); // = gamma factor
		  ux = ux/u0;
		  uy = uy/u0;
		  uz = ueta*cosh(eta)+utau*sinh(eta);
		  uz /= u0;
		  
		  /* 	     epsilon = arena[ix][iy][ieta].epsilon; */
		  /* 	     rhob = arena[ix][iy][ieta].rhob; */
		  /* 	     utau = arena[ix][iy][ieta].u[0][0]; */
		  /* 	     ux = arena[ix][iy][ieta].u[0][1]; */
		  /* 	     uy = arena[ix][iy][ieta].u[0][2]; */
		  /* 	     ueta = arena[ix][iy][ieta].u[0][3]; */
		  
		  
		  // output with coordinates for testing:
		  //fprintf(out_file,"%e %e %e %e %e %e %e %e %e %e\n", 
		  //  tau, x, y, z, eta, T*hbarc, QGPfrac, ux, uy, uz);
		  //cout << "(tau,x,y,z)=(" << tau << ", " << x << ", " << y << ", " << z << ")" << endl;

		  //fprintf(out_file,"%e %e %e %e \n", x, y, z, T*hbarc);
		 
		  // exclude the actual coordinates from the output to save space:
		  fprintf(out_file,"%e %e %e %e %e\n", 
		  	  T*hbarc, QGPfrac, ux, uy, uz);
		  
		}/* ix */
	    }/* iy */
	}/* iz */

      //At this step, compute the momentum eccentricity appropriate for comparing with the v_2 at
      //the LHC for |y| < 1.0:
      int EtaRange = (int)(1./DATA->delta_eta);
      double weightedTxx = 0.;
      double weightedTxy = 0.;
      double weightedTyy = 0.;
      for(int ieta = neta*size/2-EtaRange; ieta <= neta*size/2+EtaRange; ieta++){
	//cerr << "ieta = " << ieta << endl;
	for(int ix=0; ix<nx; ix++){
	  for(int iy=0; iy<ny; iy++){
	    weightedTxx += TxxFrom[ix][iy][ieta];
	    weightedTxy += TxyFrom[ix][iy][ieta];
	    weightedTyy += TyyFrom[ix][iy][ieta];
	  }
	}
      }

      double e1 = (weightedTxx-weightedTyy)/(weightedTxx+weightedTyy);
      double e2 = 2.*weightedTxy/(weightedTxx+weightedTyy);

      fprintf(ecc_file,"%e %e %e\n", tau, e1, e2);

      //cerr << "<<T^xx-T^yy>>/<<T^xx+T^yy>> = " << (weightedTxx-weightedTyy)/(weightedTxx+weightedTyy) << endl;
      //cerr << "<<2*T^xy>>/<<T^xx+T^yy>> = " << 2.*weightedTxy/(weightedTxx+weightedTyy) << endl;

      fclose(out_file);
      fclose(ecc_file);

      util->cube_free(epsFrom,nx+1,ny+1,size*neta);
      util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
      util->cube_free(utauFrom,nx+1,ny+1,size*neta);
      util->cube_free(uxFrom,nx+1,ny+1,size*neta);
      util->cube_free(uyFrom,nx+1,ny+1,size*neta);
      util->cube_free(uetaFrom,nx+1,ny+1,size*neta);
      util->cube_free(TxxFrom,nx+1,ny+1,size*neta);
      util->cube_free(TxyFrom,nx+1,ny+1,size*neta);
      util->cube_free(TyyFrom,nx+1,ny+1,size*neta);
	  }
  free(eps);
  free(rhob);
  free(utau);
  free(ux);
  free(uy);
  free(ueta);
  free(Txx);
  free(Txy);
  free(Tyy);
  delete(util);
}/* OutputEvolutionDataXYZ */

// Output for MARTINI when the option Coordinates is set to "taueta".
// Because fluctuations in initial conditions are taken into account here,
// the range is now all x, y, and eta, not restricted to positive values.
// However, to save space, x, y, and eta are not written into the output, and
// the interpolation in MARTINI must match properly with the output generated
//here. -CFY 11/12/2010
void Grid::OutputEvolutionDataXYEta(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
{
  Util *util;
  util = new Util();
  
  int position;
  int ix, iy, ieta, nx, ny, neta;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;
  /// MPI send and receive
  
  int sizeOfData = (nx+1)*(ny+1)*(neta);
  int to, from;
  
  double *eps;
  double *rhob;
  double *utau;
  double *ux;
  double *uy;
  double *ueta;
  
  eps = (double *)malloc(sizeof(double)*sizeOfData);
  rhob = (double *)malloc(sizeof(double)*sizeOfData);
  utau = (double *)malloc(sizeof(double)*sizeOfData);
  ux = (double *)malloc(sizeof(double)*sizeOfData);
  uy = (double *)malloc(sizeof(double)*sizeOfData);
  ueta = (double *)malloc(sizeof(double)*sizeOfData);
  
  if (rank>0) //send all to rank 0
    {
      to = 0;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  position = ieta+(neta*(ix + (nx*iy)));
		  eps[position] = arena[ix][iy][ieta].epsilon;
		  rhob[position] = arena[ix][iy][ieta].rhob;
		  utau[position] = arena[ix][iy][ieta].u[0][0];
		  ux[position] = arena[ix][iy][ieta].u[0][1];
		  uy[position] = arena[ix][iy][ieta].u[0][2];
		  ueta[position] = arena[ix][iy][ieta].u[0][3];
		}
	    }
	}
      //      cout << "eps[5]=" <<eps[5]<< endl;
      
      MPI::COMM_WORLD.Send(eps,sizeOfData,MPI::DOUBLE,to,1);
      MPI::COMM_WORLD.Send(rhob,sizeOfData,MPI::DOUBLE,to,2);
      MPI::COMM_WORLD.Send(utau,sizeOfData,MPI::DOUBLE,to,3);
      MPI::COMM_WORLD.Send(ux,sizeOfData,MPI::DOUBLE,to,4);
      MPI::COMM_WORLD.Send(uy,sizeOfData,MPI::DOUBLE,to,5);
      MPI::COMM_WORLD.Send(ueta,sizeOfData,MPI::DOUBLE,to,6);
    }
  
  if (rank==0)
    {
      double ***epsFrom;
      double ***rhobFrom;
      double ***utauFrom;
      double ***uxFrom;
      double ***uyFrom;
      double ***uetaFrom;
      
      epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      utauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  epsFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon;
		  rhobFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
		  utauFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][0];
		  uxFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][1];
		  uyFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][2];
		  uetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][3];
		}
	    }
	}
      
      
      for (int irank=1; irank<size; irank++)
	{
	  from = irank;
	  MPI::COMM_WORLD.Recv(eps,sizeOfData,MPI::DOUBLE,from,1);
	  MPI::COMM_WORLD.Recv(rhob,sizeOfData,MPI::DOUBLE,from,2);
	  MPI::COMM_WORLD.Recv(utau,sizeOfData,MPI::DOUBLE,from,3);
	  MPI::COMM_WORLD.Recv(ux,sizeOfData,MPI::DOUBLE,from,4);
	  MPI::COMM_WORLD.Recv(uy,sizeOfData,MPI::DOUBLE,from,5);
	  MPI::COMM_WORLD.Recv(ueta,sizeOfData,MPI::DOUBLE,from,6);
	    
	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		      position = ieta+(neta*(ix + (nx*iy)));
		      epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		      rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		      utauFrom[ix][iy][ieta+irank*neta] = utau[position];
		      uxFrom[ix][iy][ieta+irank*neta] = ux[position];
		      uyFrom[ix][iy][ieta+irank*neta] = uy[position];
		      uetaFrom[ix][iy][ieta+irank*neta] = ueta[position];
		    }
		}
	    }
	  //      cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;
      
      stringstream numIn;
      numIn << DATA->seed;
      FILE *out_file;
      string out_name = "evolution_xyeta_"+numIn.str()+".dat";
      //char* out_name = "evolution_xyeta.dat";
      //out_file = fopen(out_name, "a");
      out_file = fopen(out_name.c_str(), "a");
      fprintf(out_file,"");
      //Although it is a little confusing, it is easiest to use (tau, x, y, eta) coordinates
      //and save at these points vx, vy, and vz. -CFY 11/16/2010
      double T1, u01, ux1, uy1, uz1, ueta1, utau1, epsilon1, rhob1, QGPfrac1;
      double eta;
      
      //No interpolation is necessary here!
      for(ieta=0; ieta<(DATA->neta)*size; ieta++){
	eta = ((double)ieta)*(DATA->delta_eta)-(DATA->eta_size)/2.0;
	for(iy=0; iy<=DATA->ny; iy++) //All y
	  {
	    for(ix=0; ix<=DATA->nx; ix++) // All x
	      {
		epsilon1 = epsFrom[ix][iy][ieta];
		rhob1 = rhobFrom[ix][iy][ieta];
		utau1 = utauFrom[ix][iy][ieta];
		ux1 = uxFrom[ix][iy][ieta];
		uy1 = uyFrom[ix][iy][ieta];
		ueta1 = uetaFrom[ix][iy][ieta];
		
		u01 = ueta1*sinh(eta)+utau1*cosh(eta); // = gamma factor
		ux1 = ux1/u01;
		uy1 = uy1/u01;
		uz1 = ueta1*cosh(eta)+utau1*sinh(eta);
		uz1 /= u01;
		
		if (DATA->whichEOS==1)
		  {
		    T1 = eos->interpolate(epsilon1, rhob1, 0);
		    QGPfrac1=(epsilon1*hbarc-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
		    //cout << T << " " << QGPfrac << endl;
		    if (QGPfrac1>1.) QGPfrac1 = 1;
		    else if (QGPfrac1<0.) QGPfrac1=0.;
		  }
		else if (DATA->whichEOS==2)
		  {
		    T1 = eos->interpolate2(epsilon1, rhob1, 1);
		    QGPfrac1 = eos->interpolate2(epsilon1, rhob1, 3);
		  }
		
		// exclude the actual coordinates from the output to save space:
		fprintf(out_file,"%e %e %e %e %e\n", T1*hbarc, QGPfrac1, ux1, uy1, uz1);
		
	      }/* ix */
	  }/* iy */
      }/* ieta */
      fclose(out_file);
      util->cube_free(epsFrom,nx+1,ny+1,size*neta);
      util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
      util->cube_free(utauFrom,nx+1,ny+1,size*neta);
      util->cube_free(uxFrom,nx+1,ny+1,size*neta);
      util->cube_free(uyFrom,nx+1,ny+1,size*neta);
      util->cube_free(uetaFrom,nx+1,ny+1,size*neta);
    }
  free(eps);
  free(rhob);
  free(utau);
  free(ux);
  free(uy);
  free(ueta);
  delete(util);
}/* OutputEvolutionDataXYEta */

void Grid::OutputPlotDataXYZ(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
{
  Util *util;
  util = new Util();
  double s;
  int position;
  int ix, iy, ieta, nx, ny, neta;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;

  /// MPI send and receive

  int sizeOfData = (nx+1)*(ny+1)*(neta);
  int to, from;
  
  double *eps;
  double *rhob;
  double *utau;
  double *ux;
  double *uy;
  double *ueta;
  double totalS = 0.;

  eps = (double *)malloc(sizeof(double)*sizeOfData);
  rhob = (double *)malloc(sizeof(double)*sizeOfData);
  utau = (double *)malloc(sizeof(double)*sizeOfData);
  ux = (double *)malloc(sizeof(double)*sizeOfData);
  uy = (double *)malloc(sizeof(double)*sizeOfData);
  ueta = (double *)malloc(sizeof(double)*sizeOfData);

  //cerr << "OK in OutputPlotDataXYZ after allocation" << endl;

  if (rank>0) //send all to rank 0
    {
      to = 0;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  position = ieta+(neta*(ix + (nx*iy)));
		  eps[position] = arena[ix][iy][ieta].epsilon;
		  rhob[position] = arena[ix][iy][ieta].rhob;
		  utau[position] = arena[ix][iy][ieta].u[0][0];
		  ux[position] = arena[ix][iy][ieta].u[0][1];
		  uy[position] = arena[ix][iy][ieta].u[0][2];
		  ueta[position] = arena[ix][iy][ieta].u[0][3];
		}
	    }
	}
      //      cout << "eps[5]=" <<eps[5]<< endl;
    
      MPI::COMM_WORLD.Send(eps,sizeOfData,MPI::DOUBLE,to,1);
      MPI::COMM_WORLD.Send(rhob,sizeOfData,MPI::DOUBLE,to,2);
      MPI::COMM_WORLD.Send(utau,sizeOfData,MPI::DOUBLE,to,3);
      MPI::COMM_WORLD.Send(ux,sizeOfData,MPI::DOUBLE,to,4);
      MPI::COMM_WORLD.Send(uy,sizeOfData,MPI::DOUBLE,to,5);
      MPI::COMM_WORLD.Send(ueta,sizeOfData,MPI::DOUBLE,to,6);
    }

  //cerr << "OK in OutputPlotDataXYZ after MPI Send" << endl;
  
  if (rank==0) 
    {
      double ***epsFrom;
      double ***rhobFrom;
      double ***utauFrom;
      double ***uxFrom;
      double ***uyFrom;
      double ***uetaFrom;
      
      epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      utauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  epsFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon;
		  rhobFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
		  utauFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][0];
		  uxFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][1];
		  uyFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][2];
		  uetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][3];
		}
	    }
	}
	

      for (int irank=1; irank<size; irank++)
	{
	  from = irank;
	  MPI::COMM_WORLD.Recv(eps,sizeOfData,MPI::DOUBLE,from,1);
	  MPI::COMM_WORLD.Recv(rhob,sizeOfData,MPI::DOUBLE,from,2);
	  MPI::COMM_WORLD.Recv(utau,sizeOfData,MPI::DOUBLE,from,3);
	  MPI::COMM_WORLD.Recv(ux,sizeOfData,MPI::DOUBLE,from,4);
	  MPI::COMM_WORLD.Recv(uy,sizeOfData,MPI::DOUBLE,from,5);
	  MPI::COMM_WORLD.Recv(ueta,sizeOfData,MPI::DOUBLE,from,6);
	  
	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		      position = ieta+(neta*(ix + (nx*iy)));
		      epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		      rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		      utauFrom[ix][iy][ieta+irank*neta] = utau[position];
		      uxFrom[ix][iy][ieta+irank*neta] = ux[position];
		      uyFrom[ix][iy][ieta+irank*neta] = uy[position];
		      uetaFrom[ix][iy][ieta+irank*neta] = ueta[position];
		    }
		}
	    }
	  //	  cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cerr << " RECEIVED ALL DATA FROM OTHER ranks" << endl;

      double eps;

      stringstream numIn;
      numIn << DATA->seed;
      FILE *s_file;
      string s_name = "entropy_"+numIn.str()+".dat";
      //char* s_name = "entropy.dat";
      //s_file = fopen(s_name, "a");
      s_file = fopen(s_name.c_str(), "a");

      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta*size; ieta++)
		{
		  if (DATA->whichEOS==2)
		    {
		      eps = epsFrom[ix][iy][ieta];
		      s = eos->interpolate2(eps, 0., 2); //entropy density
		    }
		  totalS += utauFrom[ix][iy][ieta]*s*DATA->delta_x*DATA->delta_y*DATA->delta_eta*tau;
		}
	    }
	}

      //cerr << "totalS calculated..." << endl;

      fprintf(s_file,"%f %f\n",tau,totalS);

      FILE *out_file;
      string out_name = "contourPlot_"+numIn.str()+".dat";
      //char* out_name = "contourPlot.dat";
      //out_file = fopen(out_name, "a");
      out_file = fopen(out_name.c_str(), "a");
      fprintf(out_file,"");
      int iz, nz;
      double T, x, y, z, eta, delta_z, z_size, etafrac;
      double ux, uy, ueta, utau, epsilon, rhob, QGPfrac;
      double eta_lower, eps_lower, eps_higher, rhob_lower, rhob_higher;
      double utau_lower, utau_higher, ux_lower, ux_higher, uy_lower, uy_higher, ueta_lower, ueta_higher, u0, uz;
      nz = 160;
      z_size = 40.; 
      delta_z = static_cast<double>(z_size)/static_cast<double>(nz);
      //cerr << "OK after outputting variable definitions in OutputPlotDataXYZ..." << endl;
      for(iz=0; iz<nz; iz++)
	//for(iz=nz/2; iz<=nz/2; iz++)
	{
	  // get temperature value and others by interpolation
	  z = iz*(delta_z) - (z_size)/2.0;
	  eta = asinh(z/tau);
	  //cout << "eta_size=" << (DATA->eta_size) << endl;
	  ieta = floor((eta + (DATA->eta_size)/2.0)/(DATA->delta_eta));
	  eta_lower = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
	  etafrac = (eta-eta_lower)/DATA->delta_eta;
	  
	  for(iy=0; iy<DATA->ny; iy++) // do all y
	    {
	      y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	      for(ix=0; ix<DATA->nx; ix++) // do all x
		{
		  x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
		  eps_lower = epsFrom[ix][iy][ieta];
		  rhob_lower = rhobFrom[ix][iy][ieta];
		  utau_lower = utauFrom[ix][iy][ieta];
		  ux_lower = uxFrom[ix][iy][ieta];
		  uy_lower = uyFrom[ix][iy][ieta];
		  ueta_lower = uetaFrom[ix][iy][ieta];
		  
		  if (ieta+1 < DATA->neta*size)
		    {
		      eps_higher = epsFrom[ix][iy][ieta+1];
		      rhob_higher = rhobFrom[ix][iy][ieta+1];
		      utau_higher = utauFrom[ix][iy][ieta+1];
		      ux_higher = uxFrom[ix][iy][ieta+1];
		      uy_higher = uyFrom[ix][iy][ieta+1];
		      ueta_higher = uetaFrom[ix][iy][ieta+1];
		    }
		  else
		    {
		      eps_higher = eps_lower;
		      rhob_higher = rhob_lower;
		      utau_higher = utau_lower;
		      ux_higher = ux_lower;
		      uy_higher = uy_lower;
		      ueta_higher = ueta_lower;
		    }
		  
		  eps = eps_lower * (1.-etafrac) + (etafrac)*eps_higher;
		  rhob = rhob_lower * (1.-etafrac) + (etafrac)*rhob_higher;
		  utau = utau_lower * (1.-etafrac) + (etafrac)*utau_higher;
		  ux = ux_lower * (1.-etafrac) + (etafrac)*ux_higher;
		  uy = uy_lower * (1.-etafrac) + (etafrac)*uy_higher;
		  ueta = ueta_lower * (1.-etafrac) + (etafrac)*ueta_higher;
		  
		  //cout << "utau_lower = " << utau_lower << ", ux_lower = " << ux_lower << ", uy_lower = " << uy_lower << ", ueta_lower = " << ueta_lower << endl;
		  //cout << "utau_higher = " << utau_higher << ", ux_higher = " << ux_higher << ", uy_higher = " << uy_higher << ", ueta_higher = " << ueta_higher << endl;

		  if (DATA->whichEOS==1)
		    {
		      T = eos->interpolate(eps, rhob, 0);
		      QGPfrac=(eps*hbarc-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
		      //cout << T << " " << QGPfrac << endl;
		      if (QGPfrac>1.) QGPfrac = 1;
		      else if (QGPfrac<0.) QGPfrac=0.;
		    }
		  else if (DATA->whichEOS==2)
		    {
		      T = eos->interpolate2(eps, rhob, 1);
		      QGPfrac = eos->interpolate2(eps, rhob, 3);
		      s = eos->interpolate2(eps, rhob, 2); //entropy density
		    }
		  
		  //Now these are the flow velocities as e,g, MARTINI needs them
		  u0 = ueta*sinh(eta)+utau*cosh(eta); // = gamma factor
		  ux = ux/u0;
		  uy = uy/u0;
		  uz = ueta*cosh(eta)+utau*sinh(eta);
		  uz /= u0;
		  
		  /* 	     epsilon = arena[ix][iy][ieta].epsilon; */
		  /* 	     rhob = arena[ix][iy][ieta].rhob; */
		  /* 	     utau = arena[ix][iy][ieta].u[0][0]; */
		  /* 	     ux = arena[ix][iy][ieta].u[0][1]; */
		  /* 	     uy = arena[ix][iy][ieta].u[0][2]; */
		  /* 	     ueta = arena[ix][iy][ieta].u[0][3]; */
		  
		  
		  // output with coordinates for testing:
		  //fprintf(out_file,"%e %e %e %e %e %e %e %e %e %e\n", 
		  //  tau, x, y, z, eta, T*hbarc, QGPfrac, ux, uy, uz);
		  //cout << "(tau,x,y,z)=(" << tau << ", " << x << ", " << y << ", " << z << ")" << endl;

		  //		  fprintf(out_file,"%e %e %e %e \n", x, y, z, T*hbarc);
		  fprintf(out_file,"%e %e %e %e %e\n", 
		  	  T*hbarc, QGPfrac, ux, uy, uz);
				 
		  // exclude the actual coordinates from the output to save space:
		  //fprintf(out_file,"%e %e %e %e %e\n", 
		  //	  T*hbarc, QGPfrac, ux, uy, uz);
		  
		  
		}/* ix */
	    }/* iy */
	}/* iz */
      fprintf(out_file,"\n");
	
      fclose(out_file);
      fclose(s_file);

      cout << "S(" << tau << " fm)=" << totalS << endl;
      util->cube_free(epsFrom,nx+1,ny+1,size*neta);
      util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
      util->cube_free(utauFrom,nx+1,ny+1,size*neta);
      util->cube_free(uxFrom,nx+1,ny+1,size*neta);
      util->cube_free(uyFrom,nx+1,ny+1,size*neta);
      util->cube_free(uetaFrom,nx+1,ny+1,size*neta);
     }
  free(eps);
  free(rhob);
  free(utau);
  free(ux);
  free(uy);
  free(ueta);
  delete(util);
  //cerr << "OK finishing OutputPlotDataXYZ..." << endl;
}/* OutputEvolutionDataXYZ */

void Grid::OutputFluctuatingPlotDataXYEta(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
{
  Util *util;
  util = new Util();
  double s;
  int position;
  int ix, iy, ieta, nx, ny, neta;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;

  /// MPI send and receive

  int sizeOfData = (nx+1)*(ny+1)*(neta);
  int to, from;
  
  double *PTo;
  double *deltaPTo;
  double *ETo;
  double *deltaETo;
  double *UxTo;
  double *UyTo;
  double *UzTo;
  double *deltaUtTo;
  double *deltaUxTo;
  double *deltaUyTo;
  double *deltaUzTo;
  double *WttTo;
  double *WtxTo;
  double *WtyTo;
  double *WtzTo;
  double *WxxTo;
  double *WxyTo;
  double *WxzTo;
  double *WyyTo;
  double *WyzTo;
  double *WzzTo;

  PTo = (double *)malloc(sizeof(double)*sizeOfData);
  deltaPTo = (double *)malloc(sizeof(double)*sizeOfData);
  ETo = (double *)malloc(sizeof(double)*sizeOfData);
  deltaETo = (double *)malloc(sizeof(double)*sizeOfData);
  UxTo = (double *)malloc(sizeof(double)*sizeOfData);
  UyTo = (double *)malloc(sizeof(double)*sizeOfData);
  UzTo = (double *)malloc(sizeof(double)*sizeOfData);
  deltaUtTo = (double *)malloc(sizeof(double)*sizeOfData);
  deltaUxTo = (double *)malloc(sizeof(double)*sizeOfData);
  deltaUyTo = (double *)malloc(sizeof(double)*sizeOfData);
  deltaUzTo = (double *)malloc(sizeof(double)*sizeOfData);
  WttTo = (double *)malloc(sizeof(double)*sizeOfData);
  WtxTo = (double *)malloc(sizeof(double)*sizeOfData);
  WtyTo = (double *)malloc(sizeof(double)*sizeOfData);
  WtzTo = (double *)malloc(sizeof(double)*sizeOfData);
  WxxTo = (double *)malloc(sizeof(double)*sizeOfData);
  WxyTo = (double *)malloc(sizeof(double)*sizeOfData);
  WxzTo = (double *)malloc(sizeof(double)*sizeOfData);
  WyyTo = (double *)malloc(sizeof(double)*sizeOfData);
  WyzTo = (double *)malloc(sizeof(double)*sizeOfData);
  WzzTo = (double *)malloc(sizeof(double)*sizeOfData);

  //cerr << "OK in OutputFluctuatingPlotDataXYZ after allocation" << endl;

  if (rank>0) //send all to rank 0
    {
      to = 0;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  position = ieta+(neta*(ix + (nx*iy)));
		  PTo[position] = arena[ix][iy][ieta].p_t;
		  deltaPTo[position] = arena[ix][iy][ieta].deltaP[0];
		  ETo[position] = arena[ix][iy][ieta].epsilon_t;
		  deltaETo[position] = arena[ix][iy][ieta].deltaP[0]/eos->get_dpOverde2(arena[ix][iy][ieta].epsilon_t, arena[ix][iy][ieta].rhob_t);
		  UxTo[position] = arena[ix][iy][ieta].u[0][1];
		  UyTo[position] = arena[ix][iy][ieta].u[0][2];
		  UzTo[position] = arena[ix][iy][ieta].u[0][3];
		  deltaUtTo[position] = arena[ix][iy][ieta].deltaU[0][0];
		  deltaUxTo[position] = arena[ix][iy][ieta].deltaU[0][1];
		  deltaUyTo[position] = arena[ix][iy][ieta].deltaU[0][2];
		  deltaUzTo[position] = arena[ix][iy][ieta].deltaU[0][3];
		  if(DATA->fluctuatingHydroFlag == 2){
		    WttTo[position] = arena[ix][iy][ieta].Wmunu[0][0][0] + arena[ix][iy][ieta].deltaW[0][0][0];
		    WtxTo[position] = arena[ix][iy][ieta].Wmunu[0][0][1] + arena[ix][iy][ieta].deltaW[0][0][1];
		    WtyTo[position] = arena[ix][iy][ieta].Wmunu[0][0][2] + arena[ix][iy][ieta].deltaW[0][0][2];
		    WtzTo[position] = arena[ix][iy][ieta].Wmunu[0][0][3] + arena[ix][iy][ieta].deltaW[0][0][3];
		    WxxTo[position] = arena[ix][iy][ieta].Wmunu[0][1][1] + arena[ix][iy][ieta].deltaW[0][1][1];
		    WxyTo[position] = arena[ix][iy][ieta].Wmunu[0][1][2] + arena[ix][iy][ieta].deltaW[0][1][2];
		    WxzTo[position] = arena[ix][iy][ieta].Wmunu[0][1][3] + arena[ix][iy][ieta].deltaW[0][1][3];
		    WyyTo[position] = arena[ix][iy][ieta].Wmunu[0][2][2] + arena[ix][iy][ieta].deltaW[0][2][2];
		    WyzTo[position] = arena[ix][iy][ieta].Wmunu[0][2][3] + arena[ix][iy][ieta].deltaW[0][2][3];
		    WzzTo[position] = arena[ix][iy][ieta].Wmunu[0][3][3] + arena[ix][iy][ieta].deltaW[0][3][3];
		  }
		  else{
		    WttTo[position] = arena[ix][iy][ieta].Wmunu[0][0][0];
		    WtxTo[position] = arena[ix][iy][ieta].Wmunu[0][0][1];
		    WtyTo[position] = arena[ix][iy][ieta].Wmunu[0][0][2];
		    WtzTo[position] = arena[ix][iy][ieta].Wmunu[0][0][3];
		    WxxTo[position] = arena[ix][iy][ieta].Wmunu[0][1][1];
		    WxyTo[position] = arena[ix][iy][ieta].Wmunu[0][1][2];
		    WxzTo[position] = arena[ix][iy][ieta].Wmunu[0][1][3];
		    WyyTo[position] = arena[ix][iy][ieta].Wmunu[0][2][2];
		    WyzTo[position] = arena[ix][iy][ieta].Wmunu[0][2][3];
		    WzzTo[position] = arena[ix][iy][ieta].Wmunu[0][3][3];
		  }
		}
	    }
	}
    
      MPI::COMM_WORLD.Send(PTo,sizeOfData,MPI::DOUBLE,to,1);
      MPI::COMM_WORLD.Send(deltaPTo,sizeOfData,MPI::DOUBLE,to,2);
      MPI::COMM_WORLD.Send(ETo,sizeOfData,MPI::DOUBLE,to,3);
      MPI::COMM_WORLD.Send(deltaETo,sizeOfData,MPI::DOUBLE,to,4);
      MPI::COMM_WORLD.Send(UxTo,sizeOfData,MPI::DOUBLE,to,5);
      MPI::COMM_WORLD.Send(UyTo,sizeOfData,MPI::DOUBLE,to,6);
      MPI::COMM_WORLD.Send(UzTo,sizeOfData,MPI::DOUBLE,to,7);
      MPI::COMM_WORLD.Send(deltaUtTo,sizeOfData,MPI::DOUBLE,to,8);
      MPI::COMM_WORLD.Send(deltaUxTo,sizeOfData,MPI::DOUBLE,to,9);
      MPI::COMM_WORLD.Send(deltaUyTo,sizeOfData,MPI::DOUBLE,to,10);
      MPI::COMM_WORLD.Send(deltaUzTo,sizeOfData,MPI::DOUBLE,to,11);
      MPI::COMM_WORLD.Send(WttTo,sizeOfData,MPI::DOUBLE,to,12);
      MPI::COMM_WORLD.Send(WtxTo,sizeOfData,MPI::DOUBLE,to,13);
      MPI::COMM_WORLD.Send(WtyTo,sizeOfData,MPI::DOUBLE,to,14);
      MPI::COMM_WORLD.Send(WtzTo,sizeOfData,MPI::DOUBLE,to,15);
      MPI::COMM_WORLD.Send(WxxTo,sizeOfData,MPI::DOUBLE,to,16);
      MPI::COMM_WORLD.Send(WxyTo,sizeOfData,MPI::DOUBLE,to,17);
      MPI::COMM_WORLD.Send(WxzTo,sizeOfData,MPI::DOUBLE,to,18);
      MPI::COMM_WORLD.Send(WyyTo,sizeOfData,MPI::DOUBLE,to,19);
      MPI::COMM_WORLD.Send(WyzTo,sizeOfData,MPI::DOUBLE,to,20);
      MPI::COMM_WORLD.Send(WzzTo,sizeOfData,MPI::DOUBLE,to,21);
    }

  //cerr << "OK in OutputFluctuatingPlotDataXYZ after MPI Send" << endl;
  
  if (rank==0) 
    {
      double ***PFrom;
      double ***deltaPFrom;
      double ***EFrom;
      double ***deltaEFrom;
      double ***UxFrom;
      double ***UyFrom;
      double ***UzFrom;
      double ***deltaUtFrom;
      double ***deltaUxFrom;
      double ***deltaUyFrom;
      double ***deltaUzFrom;
      double ***WttFrom;
      double ***WtxFrom;
      double ***WtyFrom;
      double ***WtzFrom;
      double ***WxxFrom;
      double ***WxyFrom;
      double ***WxzFrom;
      double ***WyyFrom;
      double ***WyzFrom;
      double ***WzzFrom;

      PFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      deltaPFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      EFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      deltaEFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      UxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      UyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      UzFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      deltaUtFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      deltaUxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      deltaUyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      deltaUzFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WttFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WtxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WtyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WtzFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WxxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WxyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WxzFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WyyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WyzFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WzzFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  PFrom[ix][iy][ieta] = arena[ix][iy][ieta].p_t;
		  deltaPFrom[ix][iy][ieta] = arena[ix][iy][ieta].deltaP[0];
		  EFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon_t;
		  deltaEFrom[ix][iy][ieta] = arena[ix][iy][ieta].deltaP[0]/eos->get_dpOverde2(arena[ix][iy][ieta].epsilon_t, arena[ix][iy][ieta].rhob_t);
		  UxFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][1];
		  UyFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][2];
		  UzFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][3];
		  deltaUtFrom[ix][iy][ieta] = arena[ix][iy][ieta].deltaU[0][0];
		  deltaUxFrom[ix][iy][ieta] = arena[ix][iy][ieta].deltaU[0][1];
		  deltaUyFrom[ix][iy][ieta] = arena[ix][iy][ieta].deltaU[0][2];
		  deltaUzFrom[ix][iy][ieta] = arena[ix][iy][ieta].deltaU[0][3];
		  if(DATA->fluctuatingHydroFlag == 2){
		    WttFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][0] + arena[ix][iy][ieta].deltaW[0][0][0];
		    WtxFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][1] + arena[ix][iy][ieta].deltaW[0][0][1];
		    WtyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][2] + arena[ix][iy][ieta].deltaW[0][0][2];
		    WtzFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][3] + arena[ix][iy][ieta].deltaW[0][0][3];
		    WxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][1] + arena[ix][iy][ieta].deltaW[0][1][1];
		    WxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][2] + arena[ix][iy][ieta].deltaW[0][1][2];
		    WxzFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][3] + arena[ix][iy][ieta].deltaW[0][1][3];
		    WyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][2][2] + arena[ix][iy][ieta].deltaW[0][2][2];
		    WyzFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][2][3] + arena[ix][iy][ieta].deltaW[0][2][3];
		    WzzFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][3][3] + arena[ix][iy][ieta].deltaW[0][3][3];
		  }
		  else{
		    WttFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][0];
		    WtxFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][1];
		    WtyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][2];
		    WtzFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][3];
		    WxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][1];
		    WxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][2];
		    WxzFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][3];
		    WyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][2][2];
		    WyzFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][2][3];
		    WzzFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][3][3];
		  }
		}
	    }
	}

      for (int irank=1; irank<size; irank++)
	{
	  from = irank;
	  MPI::COMM_WORLD.Recv(PTo,sizeOfData,MPI::DOUBLE,from,1);
	  MPI::COMM_WORLD.Recv(deltaPTo,sizeOfData,MPI::DOUBLE,from,2);
	  MPI::COMM_WORLD.Recv(ETo,sizeOfData,MPI::DOUBLE,from,3);
	  MPI::COMM_WORLD.Recv(deltaETo,sizeOfData,MPI::DOUBLE,from,4);
	  MPI::COMM_WORLD.Recv(UxTo,sizeOfData,MPI::DOUBLE,from,5);
	  MPI::COMM_WORLD.Recv(UyTo,sizeOfData,MPI::DOUBLE,from,6);
	  MPI::COMM_WORLD.Recv(UzTo,sizeOfData,MPI::DOUBLE,from,7);
	  MPI::COMM_WORLD.Recv(deltaUtTo,sizeOfData,MPI::DOUBLE,from,8);
	  MPI::COMM_WORLD.Recv(deltaUxTo,sizeOfData,MPI::DOUBLE,from,9);
	  MPI::COMM_WORLD.Recv(deltaUyTo,sizeOfData,MPI::DOUBLE,from,10);
	  MPI::COMM_WORLD.Recv(deltaUzTo,sizeOfData,MPI::DOUBLE,from,11);
	  MPI::COMM_WORLD.Recv(WttTo,sizeOfData,MPI::DOUBLE,from,12);
	  MPI::COMM_WORLD.Recv(WtxTo,sizeOfData,MPI::DOUBLE,from,13);
	  MPI::COMM_WORLD.Recv(WtyTo,sizeOfData,MPI::DOUBLE,from,14);
	  MPI::COMM_WORLD.Recv(WtzTo,sizeOfData,MPI::DOUBLE,from,15);
	  MPI::COMM_WORLD.Recv(WxxTo,sizeOfData,MPI::DOUBLE,from,16);
	  MPI::COMM_WORLD.Recv(WxyTo,sizeOfData,MPI::DOUBLE,from,17);
	  MPI::COMM_WORLD.Recv(WxzTo,sizeOfData,MPI::DOUBLE,from,18);
	  MPI::COMM_WORLD.Recv(WyyTo,sizeOfData,MPI::DOUBLE,from,19);
	  MPI::COMM_WORLD.Recv(WyzTo,sizeOfData,MPI::DOUBLE,from,20);
	  MPI::COMM_WORLD.Recv(WzzTo,sizeOfData,MPI::DOUBLE,from,21);
	  
	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		      position = ieta+(neta*(ix + (nx*iy)));
		      PFrom[ix][iy][ieta+irank*neta] = PTo[position];
		      deltaPFrom[ix][iy][ieta+irank*neta] = deltaPTo[position];
		      EFrom[ix][iy][ieta+irank*neta] = ETo[position];
		      deltaEFrom[ix][iy][ieta+irank*neta] = deltaETo[position];
		      UxFrom[ix][iy][ieta+irank*neta] = UxTo[position];
		      UyFrom[ix][iy][ieta+irank*neta] = UyTo[position];
		      UzFrom[ix][iy][ieta+irank*neta] = UzTo[position];
		      deltaUtFrom[ix][iy][ieta+irank*neta] = deltaUtTo[position];
		      deltaUxFrom[ix][iy][ieta+irank*neta] = deltaUxTo[position];
		      deltaUyFrom[ix][iy][ieta+irank*neta] = deltaUyTo[position];
		      deltaUzFrom[ix][iy][ieta+irank*neta] = deltaUzTo[position];
		      WttFrom[ix][iy][ieta+irank*neta] = WttTo[position];
		      WtxFrom[ix][iy][ieta+irank*neta] = WtxTo[position];
		      WtyFrom[ix][iy][ieta+irank*neta] = WtyTo[position];
		      WtzFrom[ix][iy][ieta+irank*neta] = WtzTo[position];
		      WxxFrom[ix][iy][ieta+irank*neta] = WxxTo[position];
		      WxyFrom[ix][iy][ieta+irank*neta] = WxyTo[position];
		      WxzFrom[ix][iy][ieta+irank*neta] = WxzTo[position];
		      WyyFrom[ix][iy][ieta+irank*neta] = WyyTo[position];
		      WyzFrom[ix][iy][ieta+irank*neta] = WyzTo[position];
		      WzzFrom[ix][iy][ieta+irank*neta] = WzzTo[position];
		    }
		}
	    }
	}
      
      //cerr << " RECEIVED ALL DATA FROM OTHER ranks" << endl;

      double eps;

      stringstream numIn;
      numIn << DATA->seed;
      FILE *out_file;
      string out_name = "contourFluctuatingPlot_"+numIn.str()+".dat";
      //char* out_name = "contourPlot.dat";
      //out_file = fopen(out_name, "a");
      out_file = fopen(out_name.c_str(), "a");
      fprintf(out_file,"");

      for(ieta=0; ieta<neta*size; ieta++){
	for(iy=0; iy<DATA->ny; iy++){
	  for(ix=0; ix<DATA->nx; ix++){
	    fprintf(out_file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
		    PFrom[ix][iy][ieta], deltaPFrom[ix][iy][ieta], EFrom[ix][iy][ieta], deltaEFrom[ix][iy][ieta], 
		    UxFrom[ix][iy][ieta], UyFrom[ix][iy][ieta], UzFrom[ix][iy][ieta],
		    deltaUtFrom[ix][iy][ieta], deltaUxFrom[ix][iy][ieta], deltaUyFrom[ix][iy][ieta], deltaUzFrom[ix][iy][ieta],
		    WttFrom[ix][iy][ieta], WtxFrom[ix][iy][ieta], WtyFrom[ix][iy][ieta], WtzFrom[ix][iy][ieta],
		    WxxFrom[ix][iy][ieta], WxyFrom[ix][iy][ieta], WxzFrom[ix][iy][ieta],
		    WyyFrom[ix][iy][ieta], WyzFrom[ix][iy][ieta],
		    WzzFrom[ix][iy][ieta]);
	  }/* ix */
	}/* iy */
      }/* ieta */
      fprintf(out_file,"\n");
	
      fclose(out_file);

      util->cube_free(PFrom,nx+1,ny+1,size*neta);
      util->cube_free(deltaPFrom,nx+1,ny+1,size*neta);
      util->cube_free(EFrom,nx+1,ny+1,size*neta);
      util->cube_free(deltaEFrom,nx+1,ny+1,size*neta);
      util->cube_free(UxFrom,nx+1,ny+1,size*neta);
      util->cube_free(UyFrom,nx+1,ny+1,size*neta);
      util->cube_free(UzFrom,nx+1,ny+1,size*neta);
      util->cube_free(deltaUtFrom,nx+1,ny+1,size*neta);
      util->cube_free(deltaUxFrom,nx+1,ny+1,size*neta);
      util->cube_free(deltaUyFrom,nx+1,ny+1,size*neta);
      util->cube_free(deltaUzFrom,nx+1,ny+1,size*neta);
      util->cube_free(WttFrom,nx+1,ny+1,size*neta);
      util->cube_free(WtxFrom,nx+1,ny+1,size*neta);
      util->cube_free(WtyFrom,nx+1,ny+1,size*neta);
      util->cube_free(WtzFrom,nx+1,ny+1,size*neta);
      util->cube_free(WxxFrom,nx+1,ny+1,size*neta);
      util->cube_free(WxyFrom,nx+1,ny+1,size*neta);
      util->cube_free(WxzFrom,nx+1,ny+1,size*neta);
      util->cube_free(WyyFrom,nx+1,ny+1,size*neta);
      util->cube_free(WyzFrom,nx+1,ny+1,size*neta);
      util->cube_free(WzzFrom,nx+1,ny+1,size*neta);
    }
  free(PTo);
  free(deltaPTo);
  free(ETo);
  free(deltaETo);
  free(UxTo);
  free(UyTo);
  free(UzTo);
  free(deltaUtTo);
  free(deltaUxTo);
  free(deltaUyTo);
  free(deltaUzTo);
  free(WttTo);
  free(WtxTo);
  free(WtyTo);
  free(WtzTo);
  free(WxxTo);
  free(WxyTo);
  free(WxzTo);
  free(WyyTo);
  free(WyzTo);
  free(WzzTo);
  delete(util);
  //cerr << "OK finishing OutputFluctuatingPlotDataXYZ..." << endl;
}/* OutputEvolutionDataXYZ */

void Grid::OutputXY(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
{
  Util *util;
  util = new Util();

  int position;
  int ix, iy, ieta, nx, ny, neta;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;

  /// MPI send and receive

  int sizeOfData = (nx+1)*(ny+1)*(neta);
  int to, from;
  
  double *trouble;
  double *eps;
  double *Txx;
  double *Tyy;
  double *Txy;
  double *rhob;
  double *utau;
  double *ux;
  double *uy;
  double *ueta;

  trouble = (double *)malloc(sizeof(double)*sizeOfData);
  eps = (double *)malloc(sizeof(double)*sizeOfData);
  Txx = (double *)malloc(sizeof(double)*sizeOfData);
  Tyy = (double *)malloc(sizeof(double)*sizeOfData);
  Txy = (double *)malloc(sizeof(double)*sizeOfData);
  rhob = (double *)malloc(sizeof(double)*sizeOfData);
  utau = (double *)malloc(sizeof(double)*sizeOfData);
  ux = (double *)malloc(sizeof(double)*sizeOfData);
  uy = (double *)malloc(sizeof(double)*sizeOfData);
  ueta = (double *)malloc(sizeof(double)*sizeOfData);

  if (rank>0) //send all to rank 0
    {
      to = 0;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  position = ieta+(neta*(ix + (nx*iy)));
		  eps[position] = arena[ix][iy][ieta].epsilon;
		  trouble[position] = static_cast<double>(arena[ix][iy][ieta].trouble);
		  Txx[position] = arena[ix][iy][ieta].TJb[0][1][1];
		  Tyy[position] = arena[ix][iy][ieta].TJb[0][2][2];
		  Txy[position] = arena[ix][iy][ieta].TJb[0][1][2];
		  rhob[position] = arena[ix][iy][ieta].rhob;
		  utau[position] = arena[ix][iy][ieta].u[0][0];
		  ux[position] = arena[ix][iy][ieta].u[0][1];
		  uy[position] = arena[ix][iy][ieta].u[0][2];
		  ueta[position] = arena[ix][iy][ieta].u[0][3];
		}
	    }
	}
      //      cout << "eps[5]=" <<eps[5]<< endl;
    
      MPI::COMM_WORLD.Send(trouble,sizeOfData,MPI::DOUBLE,to,0);
      MPI::COMM_WORLD.Send(eps,sizeOfData,MPI::DOUBLE,to,1);
      MPI::COMM_WORLD.Send(rhob,sizeOfData,MPI::DOUBLE,to,2);
      MPI::COMM_WORLD.Send(utau,sizeOfData,MPI::DOUBLE,to,3);
      MPI::COMM_WORLD.Send(ux,sizeOfData,MPI::DOUBLE,to,4);
      MPI::COMM_WORLD.Send(uy,sizeOfData,MPI::DOUBLE,to,5);
      MPI::COMM_WORLD.Send(ueta,sizeOfData,MPI::DOUBLE,to,6);
      MPI::COMM_WORLD.Send(Txx,sizeOfData,MPI::DOUBLE,to,7);
      MPI::COMM_WORLD.Send(Tyy,sizeOfData,MPI::DOUBLE,to,8);
      MPI::COMM_WORLD.Send(Txy,sizeOfData,MPI::DOUBLE,to,9);
    }
  
  if (rank==0) 
    {
      double ***troubleFrom;
      double ***epsFrom;
      double ***TxxFrom;
      double ***TyyFrom;
      double ***TxyFrom;
      double ***rhobFrom;
      double ***utauFrom;
      double ***uxFrom;
      double ***uyFrom;
      double ***uetaFrom;
      
      troubleFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      TxxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      TyyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      TxyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      utauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  troubleFrom[ix][iy][ieta] = static_cast<double>(arena[ix][iy][ieta].trouble);
		  epsFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon;
		  TxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][1];
		  TyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][2][2];
		  TxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][2];
		  rhobFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
		  utauFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][0];
		  uxFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][1];
		  uyFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][2];
		  uetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][3];
		}
	    }
	}
	

      for (int irank=1; irank<size; irank++)
	{
	  from = irank;
	  MPI::COMM_WORLD.Recv(trouble,sizeOfData,MPI::DOUBLE,from,0);
	  MPI::COMM_WORLD.Recv(eps,sizeOfData,MPI::DOUBLE,from,1);
	  MPI::COMM_WORLD.Recv(rhob,sizeOfData,MPI::DOUBLE,from,2);
	  MPI::COMM_WORLD.Recv(utau,sizeOfData,MPI::DOUBLE,from,3);
	  MPI::COMM_WORLD.Recv(ux,sizeOfData,MPI::DOUBLE,from,4);
	  MPI::COMM_WORLD.Recv(uy,sizeOfData,MPI::DOUBLE,from,5);
	  MPI::COMM_WORLD.Recv(ueta,sizeOfData,MPI::DOUBLE,from,6);
	  MPI::COMM_WORLD.Recv(Txx,sizeOfData,MPI::DOUBLE,from,7);
	  MPI::COMM_WORLD.Recv(Tyy,sizeOfData,MPI::DOUBLE,from,8);
	  MPI::COMM_WORLD.Recv(Txy,sizeOfData,MPI::DOUBLE,from,9);
	  
	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		      position = ieta+(neta*(ix + (nx*iy)));
		      troubleFrom[ix][iy][ieta+irank*neta] = trouble[position];
		      epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		      rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		      utauFrom[ix][iy][ieta+irank*neta] = utau[position];
		      uxFrom[ix][iy][ieta+irank*neta] = ux[position];
		      uyFrom[ix][iy][ieta+irank*neta] = uy[position];
		      uetaFrom[ix][iy][ieta+irank*neta] = ueta[position];
		      TxxFrom[ix][iy][ieta+irank*neta] = Txx[position];
		      TyyFrom[ix][iy][ieta+irank*neta] = Tyy[position];
		      TxyFrom[ix][iy][ieta+irank*neta] = Txy[position];
		    }
		}
	    }
	  //	  cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;

      stringstream numIn;
      numIn << DATA->seed;
      FILE *out_file;
      string out_name = "e_x_y_profile_"+numIn.str()+".dat";
      //char* out_name = "e_x_y_profile.dat";
      //out_file = fopen(out_name, "a");
      out_file = fopen(out_name.c_str(), "a");
      fprintf(out_file,"");
      int iz, nz;
      double trouble,trouble_lower,trouble_higher,T, x, y, z, eta, delta_z, z_size, eps, Txx, Tyy, Txy,etafrac;
      double ux, uy, ueta, utau, epsilon, rhob, QGPfrac;
      double eta_lower, eps_lower, eps_higher, rhob_lower, rhob_higher, Txx_lower, Txx_higher, Tyy_lower, Tyy_higher, Txy_lower, Txy_higher;
      double utau_lower, utau_higher, ux_lower, ux_higher, uy_lower, uy_higher, ueta_lower, ueta_higher, u0, uz;
      nz = 80;
      delta_z = 0.5;
      z_size = 40.; // z_size is fixed to 40 fm: corresponds to Hydro:zmax=20 in MARTINI!
      //iz=nz/2+1;
      iz=nz/2;//nz-1;
      // get temperature value and others by interpolation
      z = iz*(delta_z) - (z_size)/2.0;
      eta = asinh(z/tau);
      //      eta+=5;
      cerr << "eta =" << eta << endl;
      //cout << "eta_size=" << (DATA->eta_size) << endl;
      ieta = floor((eta + (DATA->eta_size)/2.0)/(DATA->delta_eta));
      eta_lower = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
      etafrac = (eta-eta_lower)/DATA->delta_eta;
      
      for(iy=0; iy<DATA->ny; iy++) 
	{
	  y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	  for(ix=0; ix<DATA->nx; ix++)
	    {
	      x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
	      trouble = troubleFrom[ix][iy][ieta];
	      eps_lower = epsFrom[ix][iy][ieta];
	      Txx_lower = TxxFrom[ix][iy][ieta];
	      Tyy_lower = TyyFrom[ix][iy][ieta];
	      Txy_lower = TxyFrom[ix][iy][ieta];
	      rhob_lower = rhobFrom[ix][iy][ieta];
	      utau_lower = utauFrom[ix][iy][ieta];
	      ux_lower = uxFrom[ix][iy][ieta];
	      uy_lower = uyFrom[ix][iy][ieta];
	      ueta_lower = uetaFrom[ix][iy][ieta];
	      
	      if (ieta+1 < DATA->neta*size)
		{
		  eps_higher = epsFrom[ix][iy][ieta+1];
		  trouble_higher = troubleFrom[ix][iy][ieta+1];
		  Txx_higher = TxxFrom[ix][iy][ieta+1];
		  Tyy_higher = TyyFrom[ix][iy][ieta+1];
		  Txy_higher = TxyFrom[ix][iy][ieta+1];
		  rhob_higher = rhobFrom[ix][iy][ieta+1];
		  utau_higher = utauFrom[ix][iy][ieta+1];
		  ux_higher = uxFrom[ix][iy][ieta+1];
		  uy_higher = uyFrom[ix][iy][ieta+1];
		  ueta_higher = uetaFrom[ix][iy][ieta+1];
		}
	      else
		{
		  eps_higher = eps_lower;
		  trouble_higher = trouble_lower;
		  Txx_higher = Txx_lower;
		  Tyy_higher = Tyy_lower;
		  Txy_higher = Txy_lower;
		  rhob_higher = rhob_lower;
		  utau_higher = utau_lower;
		  ux_higher = ux_lower;
		  uy_higher = uy_lower;
		  ueta_higher = ueta_lower;
		}
	      
	      eps = eps_lower * (1.-etafrac) + (etafrac)*eps_higher;
	      trouble = trouble_lower * (1.-etafrac) + (etafrac)*trouble_higher;
	      Txx = Txx_lower * (1.-etafrac) + (etafrac)*Txx_higher;
	      Tyy = Tyy_lower * (1.-etafrac) + (etafrac)*Tyy_higher;
	      Txy = Txy_lower * (1.-etafrac) + (etafrac)*Txy_higher;
	      rhob = rhob_lower * (1.-etafrac) + (etafrac)*rhob_higher;
	      utau = utau_lower * (1.-etafrac) + (etafrac)*utau_higher;
	      ux = ux_lower * (1.-etafrac) + (etafrac)*ux_higher;
	      uy = uy_lower * (1.-etafrac) + (etafrac)*uy_higher;
	      ueta = ueta_lower * (1.-etafrac) + (etafrac)*ueta_higher;
	      
	      if (DATA->whichEOS==1)
		{
		  T = eos->interpolate(eps, rhob, 0);
		  QGPfrac=(eps*hbarc-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
		  //cout << T << " " << QGPfrac << endl;
		  if (QGPfrac>1.) QGPfrac = 1;
		  else if (QGPfrac<0.) QGPfrac=0.;
		}
	      else if (DATA->whichEOS==2)
		{
		  T = eos->interpolate2(eps, rhob, 1);
		  QGPfrac = eos->interpolate2(eps, rhob, 3);
		}
	      
	      //Now these are the flow velocities as e,g, MARTINI needs them
	      u0 = ueta*sinh(eta)+utau*cosh(eta); // = gamma factor
	      ux = ux/u0;
	      uy = uy/u0;
	      uz = ueta*cosh(eta)+utau*sinh(eta);
	      uz /= u0;
	      
	      /* 	     epsilon = arena[ix][iy][ieta].epsilon; */
	      /* 	     rhob = arena[ix][iy][ieta].rhob; */
	      /* 	     utau = arena[ix][iy][ieta].u[0][0]; */
	      /* 	     ux = arena[ix][iy][ieta].u[0][1]; */
	      /* 	     uy = arena[ix][iy][ieta].u[0][2]; */
	      /* 	     ueta = arena[ix][iy][ieta].u[0][3]; */
	      
	      
	      // output with coordinates for testing:
	      //fprintf(out_file,"%e %e %e %e %e %e %e %e %e %e\n", 
	      //  tau, x, y, z, eta, T*hbarc, QGPfrac, ux, uy, uz);
	      //cout << "(tau,x,y,z)=(" << tau << ", " << x << ", " << y << ", " << z << ")" << endl;
	      
	      
	      fprintf(out_file,"%e %e %e %e %e %e %e %e %e %e %e %e\n",x, y, 
		      T*hbarc, QGPfrac, eps, Txx, Tyy, Txy, ux, uy, uz, trouble);
	      
	      if (ix==DATA->nx-1)
		fprintf(out_file,"\n");
	    }/* ix */
	}/* iy */
      fclose(out_file);
      util->cube_free(troubleFrom,nx+1,ny+1,size*neta);
      util->cube_free(epsFrom,nx+1,ny+1,size*neta);
      util->cube_free(TxxFrom,nx+1,ny+1,size*neta);
      util->cube_free(TyyFrom,nx+1,ny+1,size*neta);
      util->cube_free(TxyFrom,nx+1,ny+1,size*neta);
      util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
      util->cube_free(utauFrom,nx+1,ny+1,size*neta);
      util->cube_free(uxFrom,nx+1,ny+1,size*neta);
      util->cube_free(uyFrom,nx+1,ny+1,size*neta);
      util->cube_free(uetaFrom,nx+1,ny+1,size*neta);
    }
  free(trouble);
  free(eps);
  free(Txx);
  free(Tyy);
  free(Txy);
  free(rhob);
  free(utau);
  free(ux);
  free(uy);
  free(ueta);
  delete(util);
  //cerr << "OutputXY ended successfully..." << endl;
}/* OutputXY */


void Grid::PrintArena(Grid ***arena, InitData *DATA, double tau)
{
 int ix, iy, ieta;
 double x, y, eta;
 double ux, uy, ueta, utau, epsilon, rhob;
 for(ix=0; ix<=DATA->nx; ix++)
  {
   x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
   for(iy=0; iy<=DATA->ny; iy++)
    {
     y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
     for(ieta=0; ieta<DATA->neta; ieta++)
      {
       eta = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
       epsilon = arena[ix][iy][ieta].epsilon;
       rhob = arena[ix][iy][ieta].rhob;
       utau = arena[ix][iy][ieta].u[0][0];
       ux = arena[ix][iy][ieta].u[0][1];
       uy = arena[ix][iy][ieta].u[0][2];
       ueta = arena[ix][iy][ieta].u[0][3];

       printf("%e  %e  %e  %e  %e  %e  %e  %e  %e\n", 
               tau, x, y, eta, epsilon, rhob, ux/utau, uy/utau, ueta/utau);
      }/* ieta */
    }/* iy */
  }/* ix */
}/* PrintArena */

void Grid::PrintEtaEpsilon(Grid ***arena, InitData *DATA, double tau, int size, int rank)
{
 int ix, iy, ieta;
 double x, y, eta;
 double ux, uy, ueta, utau, epsilon, rhob;
 FILE *d_file;
 //char* d_name;
 string d_name;
 stringstream numIn;
 numIn << DATA->seed;
 if(rank==0) d_name = "e_profile_"+numIn.str()+".dat";
 else d_name = "e_profile2_"+numIn.str()+".dat";
 //if(rank==0) d_name = "e_profile.dat";
 //else d_name = "e_profile2.dat";
 //d_file = fopen(d_name, "a");
 d_file = fopen(d_name.c_str(), "a");
 fprintf(d_file,"\n");
 x=0.;
 y=0.;
 ix = floor(DATA->x_size/2/DATA->delta_x);
 iy = floor(DATA->y_size/2/DATA->delta_y);

 for(ieta=0; ieta<DATA->neta; ieta++)
   {
     eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
     //fprintf(stderr,"%d  %d  %d \n", ix, iy, ieta);
     epsilon = arena[ix][iy][ieta].epsilon;
     fprintf(d_file,"%e  %e  %e \n", tau, eta, epsilon*hbarc);
   }
 //fprintf(stderr,"grid says arena[20][20][0].epsilon = %e \n", arena[20][20][0].epsilon);
 fclose(d_file);
}/* PrintEtaEpsilon */

void Grid::PrintxEpsilon(Grid ***arena, InitData *DATA, double tau, int size, int rank)
{
 int ix, iy, ieta;
 double x, y, eta;
 double ux, uy, ueta, utau, epsilon, rhob;
 if(size==1 || (size==2 && rank==1))
   {
     stringstream numIn;
     numIn << DATA->seed;
     FILE *d_file;
     //char* d_name = "e_x_profile.dat";
     //d_file = fopen(d_name, "a");
     string d_name = "e_x_profile_"+numIn.str()+".dat";
     //d_file = fopen(d_name, "a");
     d_file = fopen(d_name.c_str(), "a");
     fprintf(d_file,"\n");
     FILE *d2_file;
     //char* d2_name = "e_y_profile.dat";
     //d2_file = fopen(d2_name, "a");
     string d2_name = "e_y_profile_"+numIn.str()+".dat";
     d2_file = fopen(d2_name.c_str(), "a");
     fprintf(d2_file,"\n");
     
     if (size==1) ieta = floor(DATA->eta_size/2/DATA->delta_eta);
     else if (size==2) ieta = 0;
     //if (size==1) ieta = floor(DATA->eta_size/2/DATA->delta_eta+DATA->eta_size/4/DATA->delta_eta)-2;
     //else if (size==2) ieta = floor(DATA->eta_size/2-6);
     iy = floor(DATA->y_size/2/DATA->delta_y);
     
     for(ix=0; ix<=DATA->nx; ix++)
       {
	 x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
	 //fprintf(stderr,"%d  %d  %d \n", ix, iy, ieta);
	 epsilon = arena[ix][iy][ieta].epsilon;
	 fprintf(d_file,"%e  %e  %e \n", tau, x, epsilon*hbarc);
       }
     fclose(d_file);
     
     ix = floor(DATA->x_size/2/DATA->delta_x);
     
     for(iy=0; iy<=DATA->ny; iy++)
       {
	 y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	 //fprintf(stderr,"%d  %d  %d \n", ix, iy, ieta);
	 epsilon = arena[ix][iy][ieta].epsilon;
	 fprintf(d2_file,"%e  %e  %e \n", tau, y, epsilon*hbarc);
       }
     fclose(d2_file);
   }
}/* PrintxEpsilon */


void Grid::PrintGrid(Grid *grid_p, int rk_order)
{
 int r, m,n,i,j,a;

 fprintf(stderr, "\nStarting PrintGrid...\n");
 fprintf(stderr, "TJb\n");
 for(r=0; r<rk_order; r++)
  {
   for(a=0; a<5; a++)
    {
     for(i=0; i<4; i++)
      {
       fprintf(stderr, "TJb[%d][%d][%d] = %e\n", r, a, i, grid_p->TJb[r][a][i]);
      }/* i */
    }/* a */
  }/* r */ 
 fprintf(stderr, "\n");
 
 fprintf(stderr, "u_mu\n");
 for(r=0; r<rk_order; r++)
  {
     for(i=0; i<4; i++)
      {
       fprintf(stderr, "u[%d][%d] = %e\n", r, i, grid_p->u[r][i]);
      }/* i */
  }/* r */ 
 fprintf(stderr, "\n");

 fprintf(stderr, "p = %e\n", grid_p->p);
 fprintf(stderr, "epsilon = %e\n", grid_p->epsilon);
 fprintf(stderr, "rhob = %e\n", grid_p->rhob);
 fprintf(stderr, "Done PrintGrid\n\n");
}/* PrintGrid */

void Grid::PrintAxy2(InitData *DATA, Grid ***arena, double tau)
{
 int ix, iy, ieta;
 double x, y, fxx, fyy, gxx, gyy, f00, g00, hxx, hyy, h00, kxx, kyy, k00;

 ieta = DATA->neta-1; //if i want midrapidity when using 2 processors
 for(ix=0; ix<=DATA->nx; ix++)
  {
   x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
   for(iy=0; iy<=DATA->ny; iy++)
    {
     y = iy*(DATA->delta_x) - (DATA->y_size/2.0);  
/* epsilon is in fm^-4 */
     gxx = arena[ix][iy][ieta].TJb[0][1][1];
     fxx = arena[ix][iy][ieta].Wmunu[0][1][1];
     //hxx = arena[ix][iy][ieta].Pimunu[0][1][1];
     kxx = gxx + fxx;// + hxx;
     
     gyy = arena[ix][iy][ieta].TJb[0][2][2];
     fyy = arena[ix][iy][ieta].Wmunu[0][2][2];
     //hyy = arena[ix][iy][ieta].Pimunu[0][2][2];
     kyy = gyy + fyy;// + hyy;
    
     g00 = arena[ix][iy][ieta].TJb[0][0][0];
     f00 = arena[ix][iy][ieta].Wmunu[0][0][0];
     //h00 = arena[ix][iy][ieta].Pimunu[0][0][0];
     k00 = g00 + f00;// + h00;
     
     printf(
     "%e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n", 
     x, y, 
     gxx, fxx, hxx, kxx,
     gyy, fyy, hyy, kyy,
     g00, f00, h00, k00,
     arena[ix][iy][ieta].epsilon);
    }
    printf("\n");
  }
  return;
}/* PrintAxy2 */


void Grid::PrintAxy(InitData *DATA, Grid ***arena, double tau)
{
 FILE *d_file;
 stringstream numIn;
 numIn << DATA->seed;
 //char* d_name = "e_x_y_profile.dat";
 //d_file = fopen(d_name, "a");
 string d_name = "e_x_y_profile_"+numIn.str()+".dat";
 d_file = fopen(d_name.c_str(), "a");
 int ix, iy, ieta;
 double x, y;

 ieta = DATA->neta/2;
 for(ix=0; ix<=DATA->nx; ix++)
  {
   x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
   for(iy=0; iy<=DATA->ny; iy++)
    {
     y = iy*(DATA->delta_y) - (DATA->y_size/2.0);  
/* epsilon is in fm^-4 */
     fprintf(d_file, "%e  %e  %e  %e  %e %e %e\n", 
               x, 
	       y, 
	       arena[ix][iy][ieta].TJb[0][1][1], 
	       arena[ix][iy][ieta].TJb[0][2][2],
	       arena[ix][iy][ieta].epsilon,
	       arena[ix][iy][ieta].T,
	       arena[ix][iy][ieta].mu
	    );
     if(iy==DATA->ny)
       fprintf(d_file,"\n");
    }
  }
 fclose(d_file);
  return;
}/* PrintAxy */

void Grid::ComputeAnisotropy(InitData *DATA, Grid ***arena, EOS *eos, double tau)
{
  stringstream numIn;
  numIn << DATA->seed;
  FILE *v2_file;
  string v2_name = "aniso_"+numIn.str()+".dat";
  //char* v2_name = "aniso.dat";
  //v2_file = fopen(v2_name, "a");
  v2_file = fopen(v2_name.c_str(), "a");
  
  int ix, iy, ieta;
  double v2, eps, x, y, epsp, T;
  double numerator, denominator;
  numerator = 0.;
  denominator = 0.;
  double totalT = 0.;

  ieta = DATA->neta/2;
  for(ix=0; ix<=DATA->nx; ix++)
    {
      x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
      for(iy=0; iy<=DATA->ny; iy++)
	{
	  y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
	  eps = arena[ix][iy][ieta].epsilon;
	  if (DATA->whichEOS==1)
	    {
	      T = eos->interpolate(eps, 0, 0);
	     }
	   else if (DATA->whichEOS==2)
	     {
	       T = eos->interpolate2(eps, 0, 1);
	     }
	   numerator += (arena[ix][iy][ieta].TJb[0][1][1]-arena[ix][iy][ieta].TJb[0][2][2]);
	   denominator += (arena[ix][iy][ieta].TJb[0][1][1]+arena[ix][iy][ieta].TJb[0][2][2]);
	   totalT+=T;
	}
    }
  totalT = eos->interpolate2(arena[DATA->nx/2][DATA->ny/2][DATA->neta/2].epsilon, 0, 1);
  v2 = numerator/denominator;
  fprintf(v2_file, "%f %f %f \n", tau, v2, totalT); // also print temperature
  fclose(v2_file);
  return;
}/* ComputeV2 */

void Grid::ComputeEccentricity(InitData *DATA, Grid ***arena, double tau)
{
  stringstream numIn;
  numIn << DATA->seed;
  FILE *ecc_file;
  string ecc_name = "eccentricity_"+numIn.str()+".dat";
  //ecc_file = fopen(ecc_name, "a");
  //char* ecc_name = "eccentricity.dat";
  ecc_file = fopen(ecc_name.c_str(), "a");
  
  int ix, iy, ieta;
  double value, x, y, epsp;
  double numerator, denominator;
  numerator = 0.;
  denominator = 0.;
  
  ieta = DATA->neta/2;
  for(ix=0; ix<=DATA->nx; ix++)
    {
      x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
      for(iy=0; iy<=DATA->ny; iy++)
	{
	  y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
	  numerator += (arena[ix][iy][ieta].epsilon)*(y*y-x*x);
	  denominator += (arena[ix][iy][ieta].epsilon)*(y*y+x*x);
	}
    }
  value = numerator/denominator;
  fprintf(ecc_file, "%f %f \n", tau, value);
  fclose(ecc_file);
  return;
}/* ComputeEccentricity */

void Grid::ComputeEnergyConservation(InitData *DATA, Grid ***arena, double tau)
{
  stringstream numIn;
  numIn << DATA->seed;
  FILE *ecc_file;
  string ecc_name = "energyConservation_"+numIn.str()+".dat";
  //ecc_file = fopen(ecc_name, "a");
  //char* ecc_name = "energyConservation.dat";
  ecc_file = fopen(ecc_name.c_str(), "a");
  
  int ix, iy, ieta;
  double value, x, y, epsp, eta;
  double TtautauPart, TetaetaPart, TetaetaPrev;
  TtautauPart = 0.;
  TetaetaPart = 0.;
  TetaetaPrev = 0.;
  
  for(ieta=0; ieta<DATA->neta; ieta++)
    {
      for(ix=0; ix<=DATA->nx; ix++)
	{
	  for(iy=0; iy<=DATA->ny; iy++)
	    {
	      if (ieta==0 || ieta==DATA->neta-1)
		{
		  if (ix==0 || ix==DATA->nx)
		    {
		      if (iy==0 || iy==DATA->ny)
			{
			  TtautauPart += 0.125*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.125*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		      else
			{
			  TtautauPart += 0.25*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.25*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		    }
		  else
		    {
		      if (iy==0 || iy==DATA->ny)
			{
			  TtautauPart += 0.25*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.25*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		      else
			{
			  TtautauPart += 0.5*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.5*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		    }
		}
	      else
		{  if (ix==0 || ix==DATA->nx)
		    {
		      if (iy==0 || iy==DATA->ny)
			{
			  TtautauPart += 0.25*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.25*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		      else
			{
			  TtautauPart += 0.5*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.5*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		    }
		  else
		    {
		      if (iy==0 || iy==DATA->ny)
			{
			  TtautauPart += 0.5*(arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau));
			  TetaetaPart += 0.5*(DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33));
			}
		      else
			{
			  TtautauPart += arena[ix][iy][ieta].TJb[0][0][0]*(tau+DATA->delta_tau)-arena[ix][iy][ieta].prev_T00*(tau);
			  TetaetaPart += DATA->delta_tau*0.5*(arena[ix][iy][ieta].TJb[0][3][3]+arena[ix][iy][ieta].prev_T33);
			}
		    }
		}
	    }	  
	}
    }
  value = (TtautauPart+TetaetaPart)/TtautauPart*100.;
  cout << "energy changes by " << value << "% per time step" << endl;
  cout << "Ttautau=" << TtautauPart << endl;
  cout << "Tetaeta=" << TetaetaPart << endl;
  cout << "sum=" << TtautauPart+TetaetaPart << endl;
  fprintf(ecc_file, "%f %f \n", tau, value);
  fclose(ecc_file);
  return;
}/* ComputeEnergyConservation */


void Grid::PrintdEdEta(InitData *DATA, Grid ***arena)
{
 int ix, iy, ieta, i;
 double eta, f, g, h, l, k, x, y;

 ix = DATA->nx/2; 
 iy = DATA->ny/2; 
 
 for(ieta=0; ieta<DATA->neta; ieta++)
  {
//   x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
   eta = ieta*(DATA->delta_eta) - (DATA->eta_size/2.0);  
   f = arena[ix][iy][ieta].u[0][1];
   g = arena[ix][iy][ieta].u[0][2];
   h = arena[ix][iy][ieta].u[0][3];
   l = arena[ix][iy][ieta].u[0][0];
   k = arena[ix][iy][ieta].epsilon;
   /* everything is in fm. epsilon is in fm^-4 */
   printf("%e  %e  %e  %e  %e\n", eta, f/l, g/l, h/l, k*hbarc);
  }
}

//Here, S represents the noise term \xi, divided by \sqrt{2.*\eta*T}:
void Grid::initS(gsl_rng *rng, double sigma, InitData *DATA){
  //This subroutine is almost the same as updateS, except here, 
  //prev_S is set to be the same as the newly sampled S. The values 
  //stored in S and prev_S are simply the W_i's (the sampling of the 
  //Wiener process), and must be multiplied by the proper function of 
  //eta, zeta, and T:
  if(T < DATA->fluctuatingTMin){
    for(int alpha=0; alpha < 3; alpha++){
      for(int mu=0; mu < 3; mu++){
	S[alpha][mu] = 0.;
      }
    }
  }
  else{
    for(int alpha=0; alpha < 3; alpha++){
      for(int mu=0; mu < 3; mu++){
	if(alpha > mu){
	  S[alpha][mu] = S[mu][alpha];
	}
	else if(alpha == mu){
	  S[alpha][mu] = sqrt(2.)*gsl_ran_gaussian(rng, 1./sigma);
	}
	else{
	  S[alpha][mu] = gsl_ran_gaussian(rng, 1./sigma);
	}
      }
    }
    
    //Remove the trace of S, making the proper correlation function in noise:
    double trS = S[0][0] + S[1][1] + S[2][2];
    for(int alpha=0; alpha < 3; alpha++){
      S[alpha][alpha] -= trS/3.;
    }
    
    ////Write out the array of noises:
    //cout << "S = {" << S[0][0] << ", " << S[0][1] << ", " << S[0][2] << ", " << endl;
    //cout << "     " << S[1][0] << ", " << S[1][1] << ", " << S[1][2] << ", " << endl;
    //cout << "     " << S[2][0] << ", " << S[2][1] << ", " << S[2][2] << "} " << endl;
  }

  return;
}

//Assumes u[0][mu] and Temp[0] have been determined:
void Grid::initXi(Util *util, EOS *eos, InitData *DATA){
  //Because epsilon must be that of rk_flag, this routine must be called within 
  //or immediately after ReconstIt:
  double entrDens = eos->s_func(epsilon_t, p_t, rhob_t, DATA);
  double shearV = (DATA->shear_to_s)*entrDens;
  
  //Determine the 3x3 Z tensor in the fluid rest frame:
  double multFactor = sqrt(2.*T_t*shearV);
  double ZRest[3][3];
  for(int mu=0; mu<3; mu++){
    for(int alpha=0; alpha<3; alpha++){
      ZRest[mu][alpha] = multFactor*S[mu][alpha];
    }
  }
  
  ////Write out the array of noises:
  //cout << "Z = {" << ZRest[0][0] << ", " << ZRest[0][1] << ", " << ZRest[0][2] << ", " << endl;
  //cout << "     " << ZRest[1][0] << ", " << ZRest[1][1] << ", " << ZRest[1][2] << ", " << endl;
  //cout << "     " << ZRest[2][0] << ", " << ZRest[2][1] << ", " << ZRest[2][2] << "} " << endl;
  
  //Now, boost into the grid's frame:
  double uBoost[4];
  for(int mu=0; mu<4; mu++){
    uBoost[mu] = u[0][mu];
  }
  double ZCalc[4][4];
  
  util->LorentzBoost3x3Tensor(ZRest, uBoost, ZCalc);
  
  for(int mu=0; mu<4; mu++){
    for(int alpha=0; alpha<4; alpha++){
      //Xi[0][mu][alpha] = ZCalc[mu][alpha];
      //Instantiate Xi to be zero:
      Xi[0][mu][alpha] = 0.;
    }
  }
  
  //cout << "Xi = {" << Xi[0][0][0] << ", " << Xi[0][0][1] << ", " << Xi[0][0][2] << ", " << Xi[0][0][3] << "," << endl;
  //cout << "      " << Xi[0][1][0] << ", " << Xi[0][1][1] << ", " << Xi[0][1][2] << ", " << Xi[0][1][3] << "," << endl;
  //cout << "      " << Xi[0][2][0] << ", " << Xi[0][2][1] << ", " << Xi[0][2][2] << ", " << Xi[0][2][3] << "," << endl;
  //cout << "      " << Xi[0][3][0] << ", " << Xi[0][3][1] << ", " << Xi[0][3][2] << ", " << Xi[0][3][3] << "}" << endl;

  return;
}

void Grid::updateS(gsl_rng *rng, double sigma, InitData *DATA){
  if(T < DATA->fluctuatingTMin){
    for(int alpha=0; alpha < 3; alpha++){
      for(int mu=0; mu < 3; mu++){
	S[alpha][mu] = 0.;
      }
    }
  }
  else{
    for(int alpha=0; alpha < 3; alpha++){
      for(int mu=0; mu < 3; mu++){
	if(alpha > mu){
	  S[alpha][mu] = S[mu][alpha];
	}
	else if(alpha == mu){
	  S[alpha][mu] = sqrt(2.)*gsl_ran_gaussian(rng, 1./sigma);
	}
	else{
	  S[alpha][mu] = gsl_ran_gaussian(rng, 1./sigma);
	}
      }
    }
    
    //Remove the trace of S, giving the proper correlation function:
    double trS = S[0][0] + S[1][1] + S[2][2];
    for(int alpha=0; alpha<3; alpha++){
      S[alpha][alpha] -= trS/3.;
    }
  }

  ////Write out the array of noises:
  //cout << "S = {" << S[0][0] << ", " << S[0][1] << ", " << S[0][2] << ", " << endl;
  //cout << "     " << S[1][0] << ", " << S[1][1] << ", " << S[1][2] << ", " << endl;
  //cout << "     " << S[2][0] << ", " << S[2][1] << ", " << S[2][2] << "} " << endl;

  return;
  
}

double Grid::getPrevDeltaT(int mu, int nu, EOS *eos){

  double dedp = 1./eos->get_dpOverde2(prev_epsilon, rhob);
  double prev_deltaT = (prev_epsilon + prev_p)*(prev_deltaU[mu]*prev_u[0][nu] + prev_u[0][mu]*prev_deltaU[nu])
    +(1. + dedp)*prev_deltaP*prev_u[0][mu]*prev_u[0][nu];
  if(mu==nu && mu==0){
    prev_deltaT += -prev_deltaP;
  }
  if(mu==nu && mu!=0){
    prev_deltaT += prev_deltaP;
  }

  return prev_deltaT;
}
