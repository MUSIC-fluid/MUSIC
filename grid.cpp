#include "util.h"
#include "grid.h"

using namespace std;

Grid *Grid::grid_v_malloc(int n1)
{
  Grid *d1_ptr;
  int i;
  
  /* pointer to the n1 array */
  d1_ptr = new Grid[n1];
  //(Grid *) malloc (sizeof(Grid )*n1);
  
  return d1_ptr;
}/* grid_v_malloc */


Grid **Grid::grid_m_malloc(int n1, int n2)
{
    int i, j;
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
    int i,j,k,inc;
    Grid ***d1_ptr, *tmp_ptr;

    d1_ptr = new Grid **[n1];

    for(i=0; i<n1; i++) 
     {
      d1_ptr[i] = new Grid *[n2];
     } 
    
    for(i=0; i<n1; i++)
    {
     for(j=0; j<n2; j++) 
      {
       d1_ptr[i][j] = new Grid[n3];
      }
    }
    return d1_ptr;
}/* grid_c_malloc */

// Grid ***Grid::grid_c_malloc(int n1, int n2, int n3)
// {
//     int i,j,k,inc;
//     Grid ***d1_ptr, *tmp_ptr;

//     tmp_ptr = (Grid *) malloc(sizeof(Grid)*n1*n2*n3);

//     /* pointer to the n1*n2*n3 memory */

//     d1_ptr = (Grid ***) malloc (sizeof(Grid **)*n1);

//     for(i=0; i<n1; i++) 
//      {
//       d1_ptr[i] = (Grid **) malloc (sizeof(Grid *)*n2);
//      } 
    
//     for(i=0; i<n1; i++)
//     {
//      for(j=0; j<n2; j++) 
//       {
//        inc = n2*n3*i + n3*j;
//        d1_ptr[i][j] = &(tmp_ptr[inc]);
//       }
//     }
    
//     return d1_ptr;
// }/* grid_c_malloc */

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
	  //	  cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;

      FILE *out_file;
      char* out_name = "evolution.dat";
      out_file = fopen(out_name, "a");
      fprintf(out_file,"");
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
}/* OutputEvolutionDataXYZ */


void Grid::OutputEvolutionOSCAR(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
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
	  //	  cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;

      FILE *out_file;
      char* out_name = "OSCAR.dat";
      out_file = fopen(out_name, "a");
      fprintf(out_file,"");
      int iz, nz;
      double T, x, y, z, eta, delta_z, z_size, eps, etafrac;
      double ux, uy, ueta, utau, epsilon, rhob, QGPfrac;
      double eta_lower, eps_lower, eps_higher, rhob_lower, rhob_higher;
      double utau_lower, utau_higher, ux_lower, ux_higher, uy_lower, uy_higher, ueta_lower, ueta_higher, u0, uz;
      for(ieta=0; ieta<neta; ieta++)
	{
	  // get temperature value and others by interpolation
	  //	  z = iz*(delta_z) - (z_size)/2.0;
	  eta = -(DATA->eta_size)/2.0 + ieta*DATA->delta_eta;
	  //cout << "eta_size=" << (DATA->eta_size) << endl;

	  //eta_lower = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
	  //etafrac = (eta-eta_lower)/DATA->delta_eta;
	  
	  for(iy=0; iy<DATA->ny; iy++)
	    {
	      y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	      for(ix=0; ix<DATA->nx; ix++)
		{
		  x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
		  eps = epsFrom[ix][iy][ieta];
		  rhob = rhobFrom[ix][iy][ieta];
		  utau = utauFrom[ix][iy][ieta];
		  ux = uxFrom[ix][iy][ieta];
		  uy = uyFrom[ix][iy][ieta];
		  ueta = uetaFrom[ix][iy][ieta];
		  
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
		  
		  // it ix iy ieta e p T R_qgp vx vy y_L [2*C entries, D entries, T entries]
		  fprintf(out_file,"%d %d %d %d %e %e %e %e %e %e %e\n",
			  static_cast<int>(((tau-DATA->tau0)/static_cast<double>(DATA->delta_tau)+0.0000000000001)/10.), 
			  ix, 
			  iy, 
			  ieta, 
			  eps*hbarc,
			  0.,
		  	  T*hbarc, 
			  QGPfrac, 
			  ux, 
			  uy, 
			  uz); //this has to change if I want to comply with OSCAR... :(
		  
		}/* ix */
	    }/* iy */
	}/* iz */
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
}/* OutputEvolutionDataXYZ */



void Grid::OutputPlotDataXYZ(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
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
	  //	  cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;

      FILE *out_file;
      char* out_name = "contourPlot.dat";
      out_file = fopen(out_name, "a");
      fprintf(out_file,"");
      int iz, nz;
      double T, x, y, z, eta, delta_z, z_size, eps, etafrac;
      double ux, uy, ueta, utau, epsilon, rhob, QGPfrac;
      double eta_lower, eps_lower, eps_higher, rhob_lower, rhob_higher;
      double utau_lower, utau_higher, ux_lower, ux_higher, uy_lower, uy_higher, ueta_lower, ueta_higher, u0, uz;
      nz = 160;
      z_size = 40.; 
      delta_z = static_cast<double>(z_size)/static_cast<double>(nz);
      for(iz=0; iz<nz; iz++)
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
}/* OutputEvolutionDataXYZ */


void Grid::getAverageTandPlasmaEvolution(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
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

  eps = new double[sizeOfData];
  rhob = new double[sizeOfData];

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
		}
	    }
	}
      //      cout << "eps[5]=" <<eps[5]<< endl;
    
      MPI::COMM_WORLD.Send(eps,sizeOfData,MPI::DOUBLE,to,1);
      MPI::COMM_WORLD.Send(rhob,sizeOfData,MPI::DOUBLE,to,2);
    }
  
  if (rank==0) 
    {
      double ***epsFrom;
      double ***rhobFrom;
      
      epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  epsFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon;
		  rhobFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
		}
	    }
	}
	

      for (int irank=1; irank<size; irank++)
	{
	  from = irank;
	  MPI::COMM_WORLD.Recv(eps,sizeOfData,MPI::DOUBLE,from,1);
	  MPI::COMM_WORLD.Recv(rhob,sizeOfData,MPI::DOUBLE,from,2);
	  
	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		      position = ieta+(neta*(ix + (nx*iy)));
		      epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		      rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		    }
		}
	    }
	  //	  cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;

      neta*=size;
      ofstream out_file;
      string out_name = "avgT.dat";
      out_file.open(out_name.c_str(), ios::out | ios::app );
      ofstream out_file2;
      string out_name2 = "plasmaEvolutionTime.dat";
      out_file2.open(out_name2.c_str(), ios::out | ios::app );

      int iz, nCells=0, nCells2=0;
      double T, totalT=0, totalT2=0, x, y, z, eta, delta_z, z_size, eps, etafrac;
      double ux, uy, ueta, utau, epsilon, rhob, QGPfrac;
      double eta_lower, eps_lower, eps_higher, rhob_lower, rhob_higher;
      for(ieta=0; ieta<neta; ieta++)
	{
	  for(iy=0; iy<DATA->ny; iy++) // do all y
	    {
	      for(ix=0; ix<DATA->nx; ix++) // do all x
		{
		  eps = epsFrom[ix][iy][ieta];
		  rhob = rhobFrom[ix][iy][ieta];
			  
		  if (DATA->whichEOS==1)
		    {
		      T = eos->interpolate(eps, rhob, 0);
		      QGPfrac=(eps*hbarc-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
		      //cout << T << " " << QGPfrac << endl;
		      if (QGPfrac>1.) QGPfrac = 1;
		      else if (QGPfrac<0.) QGPfrac=0.;
		    }
		  else if (DATA->whichEOS>1)
		    {
		      T = eos->interpolate2(eps, rhob, 1);
		      QGPfrac = eos->interpolate2(eps, rhob, 3);
		    }
		  
		  if(T>0.16/hbarc && QGPfrac==1)
		    {
		      totalT+=T;
		      nCells+=1;
		    }
	
		  if(T>0.16/hbarc)
		    {
		      totalT2+=T;
		      nCells2+=1;
		    }
		  
		}/* ix */
	    }/* iy */
	}/* iz */
	
      if(nCells==0 && DATA->plasmaEvolutionTime==0)
	{
	  DATA->plasmaEvolutionTime=tau-DATA->tau0;
	  out_file2 << setprecision(10) << tau-DATA->tau0;
     	}

      if(nCells2==0 && DATA->plasmaEvolutionTime2==0)
	{
	  DATA->plasmaEvolutionTime2=tau-DATA->tau0;
	  out_file2 << setprecision(10) << tau-DATA->tau0 << endl;
	}
      
      if (nCells>0)
	{
	  totalT/=nCells;
	  DATA->avgT+=totalT;
	  DATA->nSteps+=1;
	}

      if (nCells2>0)
	{
	  totalT2/=nCells2;
	  DATA->avgT2+=totalT2;
	  DATA->nSteps2+=1;
	}
      
      out_file << setprecision(10) << setw(10) << tau << " " << setw(10) <<  totalT*hbarc << " " << setw(10) << totalT2*hbarc << endl;
      
      out_file.close();
      out_file2.close();
      util->cube_free(epsFrom,nx+1,ny+1,size*neta);
      util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
    }

  delete [] eps;
  delete [] rhob;
  delete(util);
  
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
  double s;
  
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
  //  double *corr;
  
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
  //  corr = (double *)malloc(sizeof(double)*sizeOfData);

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
		  //corr[position] = arena[ix][iy][ieta].rhob;
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
      //    MPI::COMM_WORLD.Send(corr,sizeOfData,MPI::DOUBLE,to,10);
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
      //      double ***corrFrom;
      
      troubleFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      //corrFrom = util->cube_malloc(nx+1,ny+1,size*neta);
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
		  //corrFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
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
	  //	  MPI::COMM_WORLD.Recv(corr,sizeOfData,MPI::DOUBLE,from,10);
	  
	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		      position = ieta+(neta*(ix + (nx*iy)));
		      //	      corrFrom[ix][iy][ieta+irank*neta] = corr[position];
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
      FILE *ecc_file;
      char* ecc_name = "eccentricity.dat";
      ecc_file = fopen(ecc_name, "a");
      
      double value, rA, epsp;
      double numerator, denominator;
      double avcos, avsin, avcos3, avsin3;
      double Psi2, phiA, Psi3, avrSq, avr3;
      double eccentricity2, eccentricity3;

      FILE *out_file;
      char* out_name = "e_x_y_profile.dat";
      out_file = fopen(out_name, "a");
      fprintf(out_file,"");
      FILE *s_file;
      char* s_name = "entropy-eta.dat";
      s_file = fopen(s_name, "a");
      fprintf(s_file,"");
      FILE *out_file_2;
      char* out_name_2 = "e_x_y_profile_05.dat";
      out_file_2 = fopen(out_name_2, "a");
      fprintf(out_file_2,"");
      int iz, nz;
      double trouble,trouble_lower,trouble_higher,T, x, y, z, eta, delta_z, z_size, eps, Txx, Tyy, Txy,etafrac;
      //double corr_lower, corr_higher;
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
      //eta = 1.;
      //      eta+=5;
      //      cerr << "eta =" << eta << endl;
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
	      //corr_lower = corrFrom[ix][iy][ieta];
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
		  //corr_higher = corrFrom[ix][iy][ieta+1];
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
		  //corr_higher = corr_lower;
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
	      //corr = corr_lower * (1.-etafrac) + (etafrac)*corr_higher;
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
	      else if (DATA->whichEOS>1)
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
	      

	      //	      if(ix%5 == 0 && iy%5 == 0)
	      fprintf(out_file,"%e %e %e %e %e %e %e %e %e %e %e %e \n",tau, x, y, 
		      T*hbarc, QGPfrac, eps, Txx, Tyy, Txy, ux, uy, ueta);
	      
	      if (ix==DATA->nx-1)
	      //if (ix==DATA->nx-1 && iy%5==0)
		fprintf(out_file,"\n");
	    }/*  */
	}/* iy */
      fclose(out_file);

      // compute and print eccentricity and angles:
      numerator = 0.;
      denominator = 0.;
      avcos = 0.;
      avsin = 0.;
      avcos3 = 0.;
      avsin3 = 0.;
      avrSq=0.;
      avr3=0.;
      ieta = DATA->neta/2;
      for(ix=0; ix<=DATA->nx; ix++)
	{
	  x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
	  for(iy=0; iy<=DATA->ny; iy++)
	    {
	      y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
	      if (x>=0)
		{
		  phiA = atan(y/x);
		  if (x==0) 
		    {
		      if (y>=0) phiA=PI/2.;
		      else if (y<0) phiA=3.*PI/2.;
		    }
		}
	      else
		{
		  phiA = atan(y/x)+PI;
		}
	      rA = sqrt( x*x + y*y );
	      avrSq += rA*rA*(arena[ix][iy][ieta].epsilon); // compute average r^2
	      avr3 += rA*rA*rA*(arena[ix][iy][ieta].epsilon);
	      
	      avcos  += rA*rA*cos(2.*phiA)*(arena[ix][iy][ieta].epsilon);
	      avsin  += rA*rA*sin(2.*phiA)*(arena[ix][iy][ieta].epsilon);
	      avcos3 += rA*rA*rA*cos(3.*phiA)*(arena[ix][iy][ieta].epsilon);
	      avsin3 += rA*rA*rA*sin(3.*phiA)*(arena[ix][iy][ieta].epsilon);
	    }
	}
      Psi2 = (atan(avsin/avcos))/2.+PI/2.;
      Psi3 = (atan(avsin3/avcos3)+PI)/3.;
      eccentricity2 = sqrt(avcos*avcos+avsin*avsin)/avrSq;
      eccentricity3 = sqrt(avcos3*avcos3+avsin3*avsin3)/avr3;
      
      cout << "ecc2=" << eccentricity2 << endl;
      cout << "Psi2=" << Psi2 << endl;
      cout << "ecc3=" << eccentricity3 << endl;
      cout << "Psi3=" << Psi3 << endl;
      fprintf(ecc_file, "%f %f %f %f %f \n", tau, eccentricity2, Psi2, eccentricity3, Psi3);
      
      fclose(ecc_file);


      // compute entropy as a function of rapidity
      for(ieta=0; ieta<DATA->neta*size; ieta++) 
	{
	  s = 0.;
	  eta = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
	  for(iy=0; iy<DATA->ny; iy++) 
	    {
	      y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	      for(ix=0; ix<DATA->nx; ix++)
		{
		  eps = epsFrom[ix][iy][ieta];
		  rhob=0.;
		  if (DATA->whichEOS==1)
		    {
		      cout << " not implemented for EOS=1" << endl;
		    }
		  else if (DATA->whichEOS>1)
		    {
		      s += eos->interpolate2(eps, rhob, 2)*utauFrom[ix][iy][ieta]*tau;
		    }
		}/*  */
	    }/* iy */
	  fprintf(s_file,"%e %e %e \n", tau, eta, s);
	}//ieta
      fprintf(s_file,"\n");
      fclose(s_file);

//       eta = 0.5;
//       //      eta+=5;
//       cerr << "eta =" << eta << endl;
//       //cout << "eta_size=" << (DATA->eta_size) << endl;
//       ieta = floor((eta + (DATA->eta_size)/2.0)/(DATA->delta_eta));
//       eta_lower = ieta*(DATA->delta_eta) - (DATA->eta_size)/2.0;
//       etafrac = (eta-eta_lower)/DATA->delta_eta;
      
//       for(iy=0; iy<DATA->ny; iy++) 
// 	{
// 	  y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
// 	  for(ix=0; ix<DATA->nx; ix++)
// 	    {
// 	      x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
// 	      trouble = troubleFrom[ix][iy][ieta];
// 	      eps_lower = epsFrom[ix][iy][ieta];
// 	      Txx_lower = TxxFrom[ix][iy][ieta];
// 	      Tyy_lower = TyyFrom[ix][iy][ieta];
// 	      Txy_lower = TxyFrom[ix][iy][ieta];
// 	      rhob_lower = rhobFrom[ix][iy][ieta];
// 	      utau_lower = utauFrom[ix][iy][ieta];
// 	      ux_lower = uxFrom[ix][iy][ieta];
// 	      uy_lower = uyFrom[ix][iy][ieta];
// 	      ueta_lower = uetaFrom[ix][iy][ieta];
	      
// 	      if (ieta+1 < DATA->neta*size)
// 		{
// 		  eps_higher = epsFrom[ix][iy][ieta+1];
// 		  trouble_higher = troubleFrom[ix][iy][ieta+1];
// 		  Txx_higher = TxxFrom[ix][iy][ieta+1];
// 		  Tyy_higher = TyyFrom[ix][iy][ieta+1];
// 		  Txy_higher = TxyFrom[ix][iy][ieta+1];
// 		  rhob_higher = rhobFrom[ix][iy][ieta+1];
// 		  utau_higher = utauFrom[ix][iy][ieta+1];
// 		  ux_higher = uxFrom[ix][iy][ieta+1];
// 		  uy_higher = uyFrom[ix][iy][ieta+1];
// 		  ueta_higher = uetaFrom[ix][iy][ieta+1];
// 		}
// 	      else
// 		{
// 		  eps_higher = eps_lower;
// 		  trouble_higher = trouble_lower;
// 		  Txx_higher = Txx_lower;
// 		  Tyy_higher = Tyy_lower;
// 		  Txy_higher = Txy_lower;
// 		  rhob_higher = rhob_lower;
// 		  utau_higher = utau_lower;
// 		  ux_higher = ux_lower;
// 		  uy_higher = uy_lower;
// 		  ueta_higher = ueta_lower;
// 		}
	      
// 	      eps = eps_lower * (1.-etafrac) + (etafrac)*eps_higher;
// 	      trouble = trouble_lower * (1.-etafrac) + (etafrac)*trouble_higher;
// 	      Txx = Txx_lower * (1.-etafrac) + (etafrac)*Txx_higher;
// 	      Tyy = Tyy_lower * (1.-etafrac) + (etafrac)*Tyy_higher;
// 	      Txy = Txy_lower * (1.-etafrac) + (etafrac)*Txy_higher;
// 	      rhob = rhob_lower * (1.-etafrac) + (etafrac)*rhob_higher;
// 	      utau = utau_lower * (1.-etafrac) + (etafrac)*utau_higher;
// 	      ux = ux_lower * (1.-etafrac) + (etafrac)*ux_higher;
// 	      uy = uy_lower * (1.-etafrac) + (etafrac)*uy_higher;
// 	      ueta = ueta_lower * (1.-etafrac) + (etafrac)*ueta_higher;
	      
// 	      if (DATA->whichEOS==1)
// 		{
// 		  T = eos->interpolate(eps, rhob, 0);
// 		  QGPfrac=(eps*hbarc-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
// 		  //cout << T << " " << QGPfrac << endl;
// 		  if (QGPfrac>1.) QGPfrac = 1;
// 		  else if (QGPfrac<0.) QGPfrac=0.;
// 		}
// 	      else if (DATA->whichEOS>1)
// 		{
// 		  T = eos->interpolate2(eps, rhob, 1);
// 		  QGPfrac = eos->interpolate2(eps, rhob, 3);
// 		}
	      
// 	      //Now these are the flow velocities as e,g, MARTINI needs them
// 	      u0 = ueta*sinh(eta)+utau*cosh(eta); // = gamma factor
// 	      ux = ux/u0;
// 	      uy = uy/u0;
// 	      uz = ueta*cosh(eta)+utau*sinh(eta);
// 	      uz /= u0;
	      
// 	      /* 	     epsilon = arena[ix][iy][ieta].epsilon; */
// 	      /* 	     rhob = arena[ix][iy][ieta].rhob; */
// 	      /* 	     utau = arena[ix][iy][ieta].u[0][0]; */
// 	      /* 	     ux = arena[ix][iy][ieta].u[0][1]; */
// 	      /* 	     uy = arena[ix][iy][ieta].u[0][2]; */
// 	      /* 	     ueta = arena[ix][iy][ieta].u[0][3]; */
	      
	      
// 	      // output with coordinates for testing:
// 	      //fprintf(out_file,"%e %e %e %e %e %e %e %e %e %e\n", 
// 	      //  tau, x, y, z, eta, T*hbarc, QGPfrac, ux, uy, uz);
// 	      //cout << "(tau,x,y,z)=(" << tau << ", " << x << ", " << y << ", " << z << ")" << endl;
	      

// 	      if(ix%5 == 0 && iy%5 == 0)
// 		fprintf(out_file_2,"%e %e %e %e %e %e %e %e %e \n",tau, x, y, 
// 			T*hbarc, QGPfrac, eps, Txx, Tyy, Txy);
	      
// 	      if (ix==DATA->nx-1 && iy%5==0)
// 		fprintf(out_file_2,"\n");
// 	    }/* ix */
// 	}/* iy */
      fclose(out_file_2);

      util->cube_free(troubleFrom,nx+1,ny+1,size*neta);
      //     util->cube_free(corrFrom,nx+1,ny+1,size*neta);
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
  //  free(corr);
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
 string d_name;
 ofstream d_file;

 if(rank==0) d_name = "e_profile.dat";
 else d_name = "e_profile2.dat";
 d_file.open(d_name.c_str(), ios::out | ios::app );

 x=0.;
 y=0.;
 ix = floor(DATA->x_size/2/DATA->delta_x);
 iy = floor(DATA->y_size/2/DATA->delta_y);

 for(ieta=0; ieta<DATA->neta; ieta++)
   {
     eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
     //fprintf(stderr,"%d  %d  %d \n", ix, iy, ieta);
     epsilon = arena[ix][iy][ieta].epsilon;
     d_file << tau << " " << eta << " " << setprecision(8) <<epsilon*hbarc << endl;
   }

 d_file << endl;
 d_file.close();
}/* PrintEtaEpsilon */

void Grid::PrintxEpsilon(Grid ***arena, InitData *DATA, double tau, int size, int rank)
{
 int ix, iy, ieta;
 double x, y, eta;
 double ux, uy, ueta, utau, epsilon, rhob;
 if(size==1 || (size==2 && rank==1))
   {
     string d_name;
     ofstream d_file;
     string d2_name;
     ofstream d2_file;

     d_name = "e_x_profile.dat";
     d_file.open(d_name.c_str(), ios::out | ios::app );
     d2_name = "e_y_profile.dat";
     d2_file.open(d2_name.c_str(), ios::out | ios::app );

     if (size==1) ieta = floor(DATA->eta_size/2/DATA->delta_eta+0.00001);
     else if (size==2) ieta = 0;
     iy = floor(DATA->y_size/2/DATA->delta_y);
     
     for(ix=0; ix<=DATA->nx; ix++)
       {
	 x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
	 //fprintf(stderr,"%d  %d  %d \n", ix, iy, ieta);
	 epsilon = arena[ix][iy][ieta].epsilon;
	 d_file << setw(10) << tau << " " << setw(10) << x << " " << setw(10) << setprecision(8) <<epsilon*hbarc << endl;
       }
     d_file << endl; 
     d_file.close();
     
     ix = floor(DATA->x_size/2/DATA->delta_x);
     
     for(iy=0; iy<=DATA->ny; iy++)
       {
	 y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	 //fprintf(stderr,"%d  %d  %d \n", ix, iy, ieta);
	 epsilon = arena[ix][iy][ieta].epsilon;
	 d_file << setw(10) << tau << " " << setw(10) << y << " " << setw(10) << setprecision(8) <<epsilon*hbarc << endl;
       }
     d2_file << endl; 
     d2_file.close();
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
     hxx = arena[ix][iy][ieta].Pimunu[0][1][1];
     kxx = gxx + fxx + hxx;
     
     gyy = arena[ix][iy][ieta].TJb[0][2][2];
     fyy = arena[ix][iy][ieta].Wmunu[0][2][2];
     hyy = arena[ix][iy][ieta].Pimunu[0][2][2];
     kyy = gyy + fyy + hyy;
    
     g00 = arena[ix][iy][ieta].TJb[0][0][0];
     f00 = arena[ix][iy][ieta].Wmunu[0][0][0];
     h00 = arena[ix][iy][ieta].Pimunu[0][0][0];
     k00 = g00 + f00 + h00;
     
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
 char* d_name = "e_x_y_profile.dat";
 d_file = fopen(d_name, "a");
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

void Grid::ComputeAnisotropy(InitData *DATA, Grid ***arena, double tau)
{
  
  string d_name;
  ofstream d_file;
  d_name = "aniso.dat";
  d_file.open(d_name.c_str(), ios::out | ios::app );
  
  int ix, iy, ieta;
  double v2, x, y, epsp;
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
	  numerator += (arena[ix][iy][ieta].TJb[0][1][1]-arena[ix][iy][ieta].TJb[0][2][2]);
	  denominator += (arena[ix][iy][ieta].TJb[0][1][1]+arena[ix][iy][ieta].TJb[0][2][2]);
	}
    }
  v2 = numerator/denominator;
  d_file << setw(10) << tau << " " << setw(10) << setprecision(10) << v2 << endl;
  d_file.close();
  return;
}/* ComputeV2 */

void Grid::ComputeEccentricity(InitData *DATA, Grid ***arena, double tau)
{
  FILE *ecc_file;
  char* ecc_name = "eccentricity.dat";
  ecc_file = fopen(ecc_name, "a");
  
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
  FILE *ecc_file;
  char* ecc_name = "energyConservation.dat";
  ecc_file = fopen(ecc_name, "a");
  
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
