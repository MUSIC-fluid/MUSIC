#include "util.h"
#include "grid.h"

using namespace std;

Grid *Grid::grid_v_malloc(int n1)
{
  Grid *d1_ptr;
//   int i;
  
  /* pointer to the n1 array */
  d1_ptr = new Grid[n1];
  //(Grid *) malloc (sizeof(Grid )*n1);
  
  return d1_ptr;
}/* grid_v_malloc */


Grid **Grid::grid_m_malloc(int n1, int n2)
{
    int i;
//     int j;
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
    int i,j;
//     int k,inc;
    Grid ***d1_ptr;
//     Grid *tmp_ptr;

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
		  position = ieta+(neta*(ix + ((nx+1)*iy)));
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
		      position = ieta+(neta*(ix + ((nx+1)*iy)));
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
      const char* out_name = "evolution.dat";
      out_file = fopen(out_name, "a");
//       fprintf(out_file,"");
      int iz, nz;
      double T, z, eta, delta_z, z_size, eps, etafrac;
//       double x, y;
      double ux, uy, ueta, utau, rhob, QGPfrac;
//       double epsilon;
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
// 	      y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	      for(ix=DATA->nx/2; ix<DATA->nx; ix++) // only do positive x
		{
// 		  x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
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
		  
		  T = eos->get_temperature(eps, rhob);
		  //QGPfrac = eos->get_qgp_frac(eps, rhob);
              QGPfrac = 0.0;
		  
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
  double *p;
  double *Wtautau, *Wtaux, *Wtauy, *Wtaueta, *Wxx, *Wxy, *Wxeta, *Wyy, *Wyeta, *Wetaeta; // added by Maxime  
//  double *Ttautau, *Ttaux, *Ttauy, *Ttaueta, *Txx, *Txy, *Txeta, *Tyy, *Tyeta, *Tetaeta; // added by Maxime


  eps = (double *)malloc(sizeof(double)*sizeOfData);
  rhob = (double *)malloc(sizeof(double)*sizeOfData);
  utau = (double *)malloc(sizeof(double)*sizeOfData);
  ux = (double *)malloc(sizeof(double)*sizeOfData);
  uy = (double *)malloc(sizeof(double)*sizeOfData);
  ueta = (double *)malloc(sizeof(double)*sizeOfData);
  p = (double *)malloc(sizeof(double)*sizeOfData);

  // added by Maxime
  Wtautau = (double *)malloc(sizeof(double)*sizeOfData);
  Wtaux = (double *)malloc(sizeof(double)*sizeOfData);
  Wtauy = (double *)malloc(sizeof(double)*sizeOfData);
  Wtaueta = (double *)malloc(sizeof(double)*sizeOfData);
  Wxx = (double *)malloc(sizeof(double)*sizeOfData);
  Wxy = (double *)malloc(sizeof(double)*sizeOfData);
  Wxeta = (double *)malloc(sizeof(double)*sizeOfData);
  Wyy = (double *)malloc(sizeof(double)*sizeOfData); 
  Wyeta = (double *)malloc(sizeof(double)*sizeOfData);
  Wetaeta = (double *)malloc(sizeof(double)*sizeOfData);
  // end of "added by Maxime"

  // added by Maxime
//  Ttautau = (double *)malloc(sizeof(double)*sizeOfData);
//  Ttaux = (double *)malloc(sizeof(double)*sizeOfData);
//  Ttauy = (double *)malloc(sizeof(double)*sizeOfData);
//  Ttaueta = (double *)malloc(sizeof(double)*sizeOfData);
//  Txx = (double *)malloc(sizeof(double)*sizeOfData);
//  Txy = (double *)malloc(sizeof(double)*sizeOfData);
//  Txeta = (double *)malloc(sizeof(double)*sizeOfData);
//  Tyy = (double *)malloc(sizeof(double)*sizeOfData); 
//  Tyeta = (double *)malloc(sizeof(double)*sizeOfData);
//  Tetaeta = (double *)malloc(sizeof(double)*sizeOfData);
  // end of "added by Maxime"
  
  if (rank>0) //send all to rank 0
    {
      to = 0;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  position = ieta+(neta*(ix + ((nx+1)*iy)));
		  eps[position] = arena[ix][iy][ieta].epsilon;
		  rhob[position] = arena[ix][iy][ieta].rhob;
		  utau[position] = arena[ix][iy][ieta].u[0][0];
		  ux[position] = arena[ix][iy][ieta].u[0][1];
		  uy[position] = arena[ix][iy][ieta].u[0][2];
		  ueta[position] = arena[ix][iy][ieta].u[0][3];
		  p[position] = arena[ix][iy][ieta].p; 

		  // added by Maxime
		  Wtautau[position] = arena[ix][iy][ieta].Wmunu[0][0][0];
		  Wtaux[position] = arena[ix][iy][ieta].Wmunu[0][0][1];
		  Wtauy[position] = arena[ix][iy][ieta].Wmunu[0][0][2];
		  Wtaueta[position] = arena[ix][iy][ieta].Wmunu[0][0][3];
		  Wxx[position] = arena[ix][iy][ieta].Wmunu[0][1][1];
		  Wxy[position] = arena[ix][iy][ieta].Wmunu[0][1][2];
		  Wxeta[position] = arena[ix][iy][ieta].Wmunu[0][1][3];
		  Wyy[position] = arena[ix][iy][ieta].Wmunu[0][2][2];
		  Wyeta[position] = arena[ix][iy][ieta].Wmunu[0][2][3];
		  Wetaeta[position] = arena[ix][iy][ieta].Wmunu[0][3][3];
		  // end of "added by Maxime"

		  // added by Maxime
//		  Ttautau[position] = arena[ix][iy][ieta].TJb[0][0][0];
//		  Ttaux[position] = arena[ix][iy][ieta].TJb[0][0][1];
//		  Ttauy[position] = arena[ix][iy][ieta].TJb[0][0][2];
//		  Ttaueta[position] = arena[ix][iy][ieta].TJb[0][0][3];
//		  Txx[position] = arena[ix][iy][ieta].TJb[0][1][1];
//		  Txy[position] = arena[ix][iy][ieta].TJb[0][1][2];
//		  Txeta[position] = arena[ix][iy][ieta].TJb[0][1][3];
//		  Tyy[position] = arena[ix][iy][ieta].TJb[0][2][2];
//		  Tyeta[position] = arena[ix][iy][ieta].TJb[0][2][3];
//		  Tetaeta[position] = arena[ix][iy][ieta].TJb[0][3][3];
		  // end of "added by Maxime"
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

      MPI::COMM_WORLD.Send(p,sizeOfData,MPI::DOUBLE,to,17);

      // added by Maxime
      MPI::COMM_WORLD.Send(Wtautau,sizeOfData,MPI::DOUBLE,to,7);
      MPI::COMM_WORLD.Send(Wtaux,sizeOfData,MPI::DOUBLE,to,8);
      MPI::COMM_WORLD.Send(Wtauy,sizeOfData,MPI::DOUBLE,to,9);
      MPI::COMM_WORLD.Send(Wtaueta,sizeOfData,MPI::DOUBLE,to,10);
      MPI::COMM_WORLD.Send(Wxx,sizeOfData,MPI::DOUBLE,to,11);
      MPI::COMM_WORLD.Send(Wxy,sizeOfData,MPI::DOUBLE,to,12);
      MPI::COMM_WORLD.Send(Wxeta,sizeOfData,MPI::DOUBLE,to,13);
      MPI::COMM_WORLD.Send(Wyy,sizeOfData,MPI::DOUBLE,to,14);
      MPI::COMM_WORLD.Send(Wyeta,sizeOfData,MPI::DOUBLE,to,15);
      MPI::COMM_WORLD.Send(Wetaeta,sizeOfData,MPI::DOUBLE,to,16);
      // end of "added by Maxime"

      // added by Maxime
//      MPI::COMM_WORLD.Send(Ttautau,sizeOfData,MPI::DOUBLE,to,18);
//      MPI::COMM_WORLD.Send(Ttaux,sizeOfData,MPI::DOUBLE,to,19);
//      MPI::COMM_WORLD.Send(Ttauy,sizeOfData,MPI::DOUBLE,to,20);
//      MPI::COMM_WORLD.Send(Ttaueta,sizeOfData,MPI::DOUBLE,to,21);
//      MPI::COMM_WORLD.Send(Txx,sizeOfData,MPI::DOUBLE,to,22);
//      MPI::COMM_WORLD.Send(Txy,sizeOfData,MPI::DOUBLE,to,23);
//      MPI::COMM_WORLD.Send(Txeta,sizeOfData,MPI::DOUBLE,to,24);
//      MPI::COMM_WORLD.Send(Tyy,sizeOfData,MPI::DOUBLE,to,25);
//      MPI::COMM_WORLD.Send(Tyeta,sizeOfData,MPI::DOUBLE,to,26);
//      MPI::COMM_WORLD.Send(Tetaeta,sizeOfData,MPI::DOUBLE,to,27);
      // end of "added by Maxime"
    }
  
  if (rank==0)
    {
      double ***epsFrom;
      double ***rhobFrom;
      double ***utauFrom;
      double ***uxFrom;
      double ***uyFrom;
      double ***uetaFrom;
      double ***WtautauFrom, ***WtauxFrom, ***WtauyFrom, ***WtauetaFrom, ***WxxFrom; // added by Maxime
      double ***WxyFrom, ***WxetaFrom, ***WyyFrom, ***WyetaFrom, ***WetaetaFrom; // added by Maxime
//      double ***TtautauFrom, ***TtauxFrom, ***TtauyFrom, ***TtauetaFrom, ***TxxFrom; // added by Maxime
//      double ***TxyFrom, ***TxetaFrom, ***TyyFrom, ***TyetaFrom, ***TetaetaFrom; // added by Maxime
      double ***pFrom;
      
      epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      utauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      uetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      pFrom = util->cube_malloc(nx+1,ny+1,size*neta);

      // added by Maxime
      WtautauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WtauxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WtauyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WtauetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WxxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WxyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WxetaFrom = util->cube_malloc(nx+1,ny+1,size*neta); 
      WyyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WyetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
      WetaetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);  
      // end of "added by Maxime"   

       // added by Maxime
//      TtautauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
//      TtauxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
//      TtauyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
//      TtauetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
//      TxxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
//      TxyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
//      TxetaFrom = util->cube_malloc(nx+1,ny+1,size*neta); 
//      TyyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
//      TyetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
//      TetaetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);  
      // end of "added by Maxime"   
      
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
		  pFrom[ix][iy][ieta] = arena[ix][iy][ieta].p;
		  
                  // added by Maxime
                  WtautauFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][0]; 
                  WtauxFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][1];
                  WtauyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][2];
                  WtauetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][3];
                  WxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][1];
                  WxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][2];
                  WxetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][3];
                  WyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][2][2];
                  WyetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][2][3];
                  WetaetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][3][3];
	  		// end of "added by Maxime"

			// added by Maxime
//	  		TtautauFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][0][0]; 
//	  		TtauxFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][0][1];
//	  		TtauyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][0][2];
//	  		TtauetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][0][3];
//	  		TxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][1];
//	  		TxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][2];
//	  		TxetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][1][3];
//	  		TyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][2][2];
//	  		TyetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][2][3];
//	  		TetaetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].TJb[0][3][3];
	  		// end of "added by Maxime"
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
	  
	  MPI::COMM_WORLD.Recv(p,sizeOfData,MPI::DOUBLE,from,17);

	  // added by Maxime
	  MPI::COMM_WORLD.Recv(Wtautau,sizeOfData,MPI::DOUBLE,from,7);
	  MPI::COMM_WORLD.Recv(Wtaux,sizeOfData,MPI::DOUBLE,from,8);
 	  MPI::COMM_WORLD.Recv(Wtauy,sizeOfData,MPI::DOUBLE,from,9);
	  MPI::COMM_WORLD.Recv(Wtaueta,sizeOfData,MPI::DOUBLE,from,10);
	  MPI::COMM_WORLD.Recv(Wxx,sizeOfData,MPI::DOUBLE,from,11);
	  MPI::COMM_WORLD.Recv(Wxy,sizeOfData,MPI::DOUBLE,from,12);
	  MPI::COMM_WORLD.Recv(Wxeta,sizeOfData,MPI::DOUBLE,from,13);
	  MPI::COMM_WORLD.Recv(Wyy,sizeOfData,MPI::DOUBLE,from,14);
	  MPI::COMM_WORLD.Recv(Wyeta,sizeOfData,MPI::DOUBLE,from,15);
   	  MPI::COMM_WORLD.Recv(Wetaeta,sizeOfData,MPI::DOUBLE,from,16);
	  // end of "added by Maxime"

	  // added by Maxime
//	  MPI::COMM_WORLD.Recv(Ttautau,sizeOfData,MPI::DOUBLE,from,18);
//	  MPI::COMM_WORLD.Recv(Ttaux,sizeOfData,MPI::DOUBLE,from,19);
// 	  MPI::COMM_WORLD.Recv(Ttauy,sizeOfData,MPI::DOUBLE,from,20);
//	  MPI::COMM_WORLD.Recv(Ttaueta,sizeOfData,MPI::DOUBLE,from,21);
//	  MPI::COMM_WORLD.Recv(Txx,sizeOfData,MPI::DOUBLE,from,22);
//	  MPI::COMM_WORLD.Recv(Txy,sizeOfData,MPI::DOUBLE,from,23);
//	  MPI::COMM_WORLD.Recv(Txeta,sizeOfData,MPI::DOUBLE,from,24);
//	  MPI::COMM_WORLD.Recv(Tyy,sizeOfData,MPI::DOUBLE,from,25);
//	  MPI::COMM_WORLD.Recv(Tyeta,sizeOfData,MPI::DOUBLE,from,26);
//   	  MPI::COMM_WORLD.Recv(Tetaeta,sizeOfData,MPI::DOUBLE,from,27);
	  // end of "added by Maxime"
	  
	  for(ix=0; ix<=nx; ix++)
	    {
	      for(iy=0; iy<=ny; iy++)
		{
		  for(ieta=0; ieta<neta; ieta++)
		    {
		    position = ieta+(neta*(ix + ((nx+1)*iy)));
		    epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		    rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		    utauFrom[ix][iy][ieta+irank*neta] = utau[position];
		    uxFrom[ix][iy][ieta+irank*neta] = ux[position];
		    uyFrom[ix][iy][ieta+irank*neta] = uy[position];
		    uetaFrom[ix][iy][ieta+irank*neta] = ueta[position];
		    pFrom[ix][iy][ieta+irank*neta] = p[position];

		    // added by Maxime
		    WtautauFrom[ix][iy][ieta+irank*neta] = Wtautau[position];
		    WtauxFrom[ix][iy][ieta+irank*neta] = Wtaux[position];
  		    WtauyFrom[ix][iy][ieta+irank*neta] = Wtauy[position];
		    WtauetaFrom[ix][iy][ieta+irank*neta] = Wtaueta[position];
		    WxxFrom[ix][iy][ieta+irank*neta] = Wxx[position];
 		    WxyFrom[ix][iy][ieta+irank*neta] = Wxy[position];
		    WxetaFrom[ix][iy][ieta+irank*neta] = Wxeta[position];
		    WyyFrom[ix][iy][ieta+irank*neta] = Wyy[position];
		    WyetaFrom[ix][iy][ieta+irank*neta] = Wyeta[position];
		    WetaetaFrom[ix][iy][ieta+irank*neta] = Wetaeta[position];
		    // end of added by Maxime

		    // added by Maxime
//		    TtautauFrom[ix][iy][ieta+irank*neta] = Ttautau[position];
//		    TtauxFrom[ix][iy][ieta+irank*neta] = Ttaux[position];
//  		    TtauyFrom[ix][iy][ieta+irank*neta] = Ttauy[position];
//		    TtauetaFrom[ix][iy][ieta+irank*neta] = Ttaueta[position];
//		    TxxFrom[ix][iy][ieta+irank*neta] = Txx[position];
// 		    TxyFrom[ix][iy][ieta+irank*neta] = Txy[position];
//		    TxetaFrom[ix][iy][ieta+irank*neta] = Txeta[position];
//		    TyyFrom[ix][iy][ieta+irank*neta] = Tyy[position];
//		    TyetaFrom[ix][iy][ieta+irank*neta] = Tyeta[position];
//		    TetaetaFrom[ix][iy][ieta+irank*neta] = Tetaeta[position];
		    // end of added by Maxime
		    }
		}
	    }
	  //      cout << "epsFrom[0][0][irank*neta+5]=" <<epsFrom[0][0][irank*neta+5]<< endl;
	}
      
      //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;
  //set up file name
  const string out_name_xyeta = "evolution_xyeta.dat";
  const string out_name_W_xyeta = "evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
  string out_open_mode;
  FILE *out_file_xyeta;
  FILE *out_file_W_xyeta;

  //If it's the first timestep, overwrite the previous file
  if (tau == DATA->tau0) {
     out_open_mode = "w";
  }
  else {
     out_open_mode = "a";	
  }
  //If we output in binary, set the mode accordingly
  if (0 == DATA->outputBinaryEvolution) {
    out_open_mode += "b";
  }

  out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());
  out_file_W_xyeta = fopen(out_name_W_xyeta.c_str(), out_open_mode.c_str());
      
//      out_file_T = fopen(out_name_T,"a");
      //fprintf(out_file,"");
      //Although it is a little confusing, it is easiest to use (tau, x, y, eta) coordinates
      //and save at these points vx, vy, and vz. -CFY 11/16/2010
      double T1, u01, ux1, uy1, uz1, ueta1, utau1, epsilon1, rhob1, QGPfrac1;
      double eta;
//       double entropy; // added by Maxime
//       double pressure; // added by Maxime
      double Wtt, Wtx, Wty, Wtz, Wzz, Wxz, Wyz; // added by Maxime
//      double Ttt, Ttx, Tty, Ttz, Tzz, Txz, Tyz; // added by Maxime
      double Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxeta, Wyy, Wyeta, Wetaeta; // added by Maxime
//      double Ttautau, Ttaux, Ttauy, Ttaueta, Txx, Txy, Txeta, Tyy, Tyeta, Tetaeta; // added by Maxime
      double div_factor, pressure1;

      //No interpolation is necessary here!
      for(ieta=0; ieta<(DATA->neta)*size; ieta++){
	eta = ((double)ieta)*(DATA->delta_eta)-(DATA->eta_size)/2.0;
	for(iy=0; iy<=DATA->ny; iy++) //All y
	  {
	    for(ix=0; ix<=DATA->nx; ix++) // All x
	      {
		epsilon1 = epsFrom[ix][iy][ieta];
		pressure1 = pFrom[ix][iy][ieta];
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
	
		T1=eos->get_temperature(epsilon1,rhob1);
//		QGPfrac1=eos->get_qgp_frac(epsilon1,rhob1);
            QGPfrac1 = 0.0;
// 		entropy=eos->get_entropy(epsilon1, rhob1);
//		if(T > 0.12)
//		{
                  div_factor=(epsilon1+pressure1);
		  Wtautau = WtautauFrom[ix][iy][ieta]/(div_factor); // I need <del_mu u_nu> not Wmunu so I divide by the viscosity -Maxime
		  Wtaux = WtauxFrom[ix][iy][ieta]/(div_factor);
		  Wtauy = WtauyFrom[ix][iy][ieta]/(div_factor);
		  Wtaueta = WtauetaFrom[ix][iy][ieta]/(div_factor);
		  Wxx = WxxFrom[ix][iy][ieta]/(div_factor);
		  Wxy = WxyFrom[ix][iy][ieta]/(div_factor);
		  Wxeta = WxetaFrom[ix][iy][ieta]/(div_factor);
		  Wyy = WyyFrom[ix][iy][ieta]/(div_factor);
		  Wyeta = WyetaFrom[ix][iy][ieta]/(div_factor);
		  Wetaeta = WetaetaFrom[ix][iy][ieta]/(div_factor);

		  Wtt = pow(cosh(eta),2)*Wtautau + pow(tau*sinh(eta),2)*Wetaeta*pow(tau,-2) + 2*tau*cosh(eta)*sinh(eta)*Wtaueta*pow(tau,-1);
		  Wtx = cosh(eta)*Wtaux + tau*sinh(eta)*Wxeta*pow(tau,-1);
		  Wty = cosh(eta)*Wtauy + tau*sinh(eta)*Wyeta*pow(tau,-1);
		  Wtz = cosh(eta)*sinh(eta)*Wtautau + tau*( pow(cosh(eta),2) + pow(sinh(eta),2) )*Wtaueta*pow(tau,-1) + tau*tau*cosh(eta)*sinh(eta)*Wetaeta*pow(tau,-2);
		  Wxz = sinh(eta)*Wtaux + tau*cosh(eta)*Wxeta*pow(tau,-1);
		  Wyz = sinh(eta)*Wtauy + tau*cosh(eta)*Wyeta*pow(tau,-1);
		  Wzz = pow(sinh(eta),2)*Wtautau + pow(tau*cosh(eta),2)*Wetaeta*pow(tau,-2) + 2*tau*cosh(eta)*sinh(eta)*Wtaueta*pow(tau,-1);
//		}
//		else
//		{
//			Wtt = 0.0;
//			Wtx = 0.0;
//			Wty = 0.0;
//			Wtz = 0.0;
//			Wxx = 0.0;
//			Wyy = 0.0;
//			Wxy = 0.0;
//			Wxz = 0.0;
//			Wyz = 0.0;
//			Wzz = 0.0;
//		}
		

//		Ttt = pow(cosh(eta),2)*Ttautau + pow(tau*sinh(eta),2)*Tetaeta*pow(tau,-2) + 2*tau*cosh(eta)*sinh(eta)*Ttaueta*pow(tau,-1);
//		Ttx = cosh(eta)*Ttaux + tau*sinh(eta)*Txeta*pow(tau,-1);
//		Tty = cosh(eta)*Ttauy + tau*sinh(eta)*Tyeta*pow(tau,-1);
//		Ttz = cosh(eta)*sinh(eta)*Ttautau + tau*( pow(cosh(eta),2) + pow(sinh(eta),2) )*Ttaueta*pow(tau,-1) + tau*tau*cosh(eta)*sinh(eta)*Tetaeta*pow(tau,-2);
//		Txz = sinh(eta)*Ttaux + tau*cosh(eta)*Txeta*pow(tau,-1);
//		Tyz = sinh(eta)*Ttauy + tau*cosh(eta)*Tyeta*pow(tau,-1);
//		Tzz = pow(sinh(eta),2)*Ttautau + pow(tau*cosh(eta),2)*Tetaeta*pow(tau,-2) + 2*tau*cosh(eta)*sinh(eta)*Ttaueta*pow(tau,-1);

                // exclude the actual coordinates from the output to save space:
                if (0 == DATA->outputBinaryEvolution) {
		  fprintf(out_file_xyeta,"%e %e %e %e %e\n", T1*hbarc, QGPfrac1, ux1, uy1, uz1);
		  if (1 == DATA->viscosity_flag) {
		    fprintf(out_file_W_xyeta,"%e %e %e %e %e %e %e %e %e %e\n",Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz); 
		  }
		}
		else {
		  double array[]={T1*hbarc, QGPfrac1, ux1, uy1, uz1};
		  fwrite(array,sizeof(double),5,out_file_xyeta);
		  //Write Wmunu/shear only if viscosity is on --- no need to fill a file with zeros in the ideal case
		  if (1 == DATA->viscosity_flag) {
		  	double array2[]={Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz};
		  	fwrite(array2,sizeof(double),10,out_file_W_xyeta);
		  }
		}
		//fprintf(out_file_T,"%e %e %e %e %e %e %e %e %e %e\n",Ttt,Ttx,Tty,Ttz,Txx,Txy,Txz,Tyy,Tyz,Tzz);

	      }/* ix */
	  }/* iy */
      }/* ieta */
      fclose(out_file_xyeta);
      fclose(out_file_W_xyeta);
//      fclose(out_file_T);

      /*End of hydro output in tau,x,y,eta*/

      util->cube_free(epsFrom,nx+1,ny+1,size*neta);
      util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
      util->cube_free(utauFrom,nx+1,ny+1,size*neta);
      util->cube_free(uxFrom,nx+1,ny+1,size*neta);
      util->cube_free(uyFrom,nx+1,ny+1,size*neta);
      util->cube_free(uetaFrom,nx+1,ny+1,size*neta);
      util->cube_free(pFrom,nx+1,ny+1,size*neta);

      // added by Maxime
      util->cube_free(WtautauFrom,nx+1,ny+1,size*neta);
      util->cube_free(WtauxFrom,nx+1,ny+1,size*neta);
      util->cube_free(WtauyFrom,nx+1,ny+1,size*neta);
      util->cube_free(WtauetaFrom,nx+1,ny+1,size*neta);
      util->cube_free(WxxFrom,nx+1,ny+1,size*neta);
      util->cube_free(WxyFrom,nx+1,ny+1,size*neta);
      util->cube_free(WxetaFrom,nx+1,ny+1,size*neta);
      util->cube_free(WyyFrom,nx+1,ny+1,size*neta);
      util->cube_free(WyetaFrom,nx+1,ny+1,size*neta);
      util->cube_free(WetaetaFrom,nx+1,ny+1,size*neta);
      // end of added by Maxime

      // added by Maxime
//      util->cube_free(TtautauFrom,nx+1,ny+1,size*neta);
//      util->cube_free(TtauxFrom,nx+1,ny+1,size*neta);
//      util->cube_free(TtauyFrom,nx+1,ny+1,size*neta);
//      util->cube_free(TtauetaFrom,nx+1,ny+1,size*neta);
//      util->cube_free(TxxFrom,nx+1,ny+1,size*neta);
//      util->cube_free(TxyFrom,nx+1,ny+1,size*neta);
//      util->cube_free(TxetaFrom,nx+1,ny+1,size*neta);
//      util->cube_free(TyyFrom,nx+1,ny+1,size*neta);
//      util->cube_free(TyetaFrom,nx+1,ny+1,size*neta);
//      util->cube_free(TetaetaFrom,nx+1,ny+1,size*neta);
      // end of added by Maxime
    }
  free(eps);
  free(rhob);
  free(utau);
  free(ux);
  free(uy);
  free(ueta);
  free(p);

  // added by Maxime
  free(Wtautau);
  free(Wtaux);
  free(Wtauy);
  free(Wtaueta);
  free(Wxx);
  free(Wxy);
  free(Wxeta);
  free(Wyy);
  free(Wyeta);
  free(Wetaeta);
  // end of "added by Maxime"

  // added by Maxime
//  free(Ttautau);
//  free(Ttaux);
//  free(Ttauy);
//  free(Ttaueta);
//  free(Txx);
//  free(Txy);
//  free(Txeta);
//  free(Tyy);
//  free(Tyeta);
//  free(Tetaeta);
  // end of "added by Maxime"

  delete(util);
}/* OutputEvolutionDataXYEta */

// output hydro evolution history with muB information for RHIC BES
void Grid::OutputEvolutionDataXYEta_finite_muB(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
{
  Util *util;
  util = new Util();
  
  int position;
  int ix, iy, ieta, nx, ny, neta;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;

  // for MPI send and receive
  int sizeOfData = (nx+1)*(ny+1)*(neta);
  int to, from;
  
  double *eps, *rhob;
  double *p;
  double *utau, *ux, *uy, *ueta;
  double *Wtautau, *Wtaux, *Wtauy, *Wtaueta, *Wxx, *Wxy, *Wxeta, *Wyy, *Wyeta, *Wetaeta;

  eps = (double *)malloc(sizeof(double)*sizeOfData);
  rhob = (double *)malloc(sizeof(double)*sizeOfData);
  utau = (double *)malloc(sizeof(double)*sizeOfData);
  ux = (double *)malloc(sizeof(double)*sizeOfData);
  uy = (double *)malloc(sizeof(double)*sizeOfData);
  ueta = (double *)malloc(sizeof(double)*sizeOfData);
  p = (double *)malloc(sizeof(double)*sizeOfData);

  Wtautau = (double *)malloc(sizeof(double)*sizeOfData);
  Wtaux = (double *)malloc(sizeof(double)*sizeOfData);
  Wtauy = (double *)malloc(sizeof(double)*sizeOfData);
  Wtaueta = (double *)malloc(sizeof(double)*sizeOfData);
  Wxx = (double *)malloc(sizeof(double)*sizeOfData);
  Wxy = (double *)malloc(sizeof(double)*sizeOfData);
  Wxeta = (double *)malloc(sizeof(double)*sizeOfData);
  Wyy = (double *)malloc(sizeof(double)*sizeOfData); 
  Wyeta = (double *)malloc(sizeof(double)*sizeOfData);
  Wetaeta = (double *)malloc(sizeof(double)*sizeOfData);

  if (rank>0) //send all to rank 0
  {
    to = 0;
    for(ix=0; ix<=nx; ix++)
    {
      for(iy=0; iy<=ny; iy++)
      {
        for(ieta=0; ieta<neta; ieta++)
        {
          position = ieta+(neta*(ix + ((nx+1)*iy)));

          eps[position] = arena[ix][iy][ieta].epsilon;
          p[position] = arena[ix][iy][ieta].p; 
          rhob[position] = arena[ix][iy][ieta].rhob;

          utau[position] = arena[ix][iy][ieta].u[0][0];
          ux[position] = arena[ix][iy][ieta].u[0][1];
          uy[position] = arena[ix][iy][ieta].u[0][2];
          ueta[position] = arena[ix][iy][ieta].u[0][3];
          
          Wtautau[position] = arena[ix][iy][ieta].Wmunu[0][0][0];
          Wtaux[position] = arena[ix][iy][ieta].Wmunu[0][0][1];
          Wtauy[position] = arena[ix][iy][ieta].Wmunu[0][0][2];
          Wtaueta[position] = arena[ix][iy][ieta].Wmunu[0][0][3];
          Wxx[position] = arena[ix][iy][ieta].Wmunu[0][1][1];
          Wxy[position] = arena[ix][iy][ieta].Wmunu[0][1][2];
          Wxeta[position] = arena[ix][iy][ieta].Wmunu[0][1][3];
          Wyy[position] = arena[ix][iy][ieta].Wmunu[0][2][2];
          Wyeta[position] = arena[ix][iy][ieta].Wmunu[0][2][3];
          Wetaeta[position] = arena[ix][iy][ieta].Wmunu[0][3][3];
        }
      }
    }

    MPI::COMM_WORLD.Send(eps,sizeOfData,MPI::DOUBLE,to,1);
    MPI::COMM_WORLD.Send(rhob,sizeOfData,MPI::DOUBLE,to,2);
    MPI::COMM_WORLD.Send(p,sizeOfData,MPI::DOUBLE,to,17);

    MPI::COMM_WORLD.Send(utau,sizeOfData,MPI::DOUBLE,to,3);
    MPI::COMM_WORLD.Send(ux,sizeOfData,MPI::DOUBLE,to,4);
    MPI::COMM_WORLD.Send(uy,sizeOfData,MPI::DOUBLE,to,5);
    MPI::COMM_WORLD.Send(ueta,sizeOfData,MPI::DOUBLE,to,6);
    
    MPI::COMM_WORLD.Send(Wtautau,sizeOfData,MPI::DOUBLE,to,7);
    MPI::COMM_WORLD.Send(Wtaux,sizeOfData,MPI::DOUBLE,to,8);
    MPI::COMM_WORLD.Send(Wtauy,sizeOfData,MPI::DOUBLE,to,9);
    MPI::COMM_WORLD.Send(Wtaueta,sizeOfData,MPI::DOUBLE,to,10);
    MPI::COMM_WORLD.Send(Wxx,sizeOfData,MPI::DOUBLE,to,11);
    MPI::COMM_WORLD.Send(Wxy,sizeOfData,MPI::DOUBLE,to,12);
    MPI::COMM_WORLD.Send(Wxeta,sizeOfData,MPI::DOUBLE,to,13);
    MPI::COMM_WORLD.Send(Wyy,sizeOfData,MPI::DOUBLE,to,14);
    MPI::COMM_WORLD.Send(Wyeta,sizeOfData,MPI::DOUBLE,to,15);
    MPI::COMM_WORLD.Send(Wetaeta,sizeOfData,MPI::DOUBLE,to,16);
  }
  
  if (rank==0)
  {
    double ***epsFrom, ***rhobFrom;
    double ***pFrom;

    double ***utauFrom, ***uxFrom, ***uyFrom, ***uetaFrom;
    double ***WtautauFrom, ***WtauxFrom, ***WtauyFrom, ***WtauetaFrom, ***WxxFrom;
    double ***WxyFrom, ***WxetaFrom, ***WyyFrom, ***WyetaFrom, ***WetaetaFrom;
      
    epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    pFrom = util->cube_malloc(nx+1,ny+1,size*neta);

    utauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    uxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    uyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    uetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);

    WtautauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    WtauxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    WtauyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    WtauetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    WxxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    WxyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    WxetaFrom = util->cube_malloc(nx+1,ny+1,size*neta); 
    WyyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    WyetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
    WetaetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);  
      
    // record rank 0 first
    for(ix=0; ix<=nx; ix++)
    {
      for(iy=0; iy<=ny; iy++)
      {
        for(ieta=0; ieta<neta; ieta++)
	  {
          epsFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon;
          rhobFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
          pFrom[ix][iy][ieta] = arena[ix][iy][ieta].p;

          utauFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][0];
          uxFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][1];
          uyFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][2];
          uetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][3];
	  
          WtautauFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][0]; 
          WtauxFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][1];
          WtauyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][2];
          WtauetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][3];
          WxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][1];
          WxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][2];
          WxetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][3];
          WyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][2][2];
          WyetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][2][3];
          WetaetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][3][3];
	  }
	}
    }
      
    // then record other rank
    for (int irank=1; irank<size; irank++)
    {
      from = irank;
      MPI::COMM_WORLD.Recv(eps,sizeOfData,MPI::DOUBLE,from,1);
      MPI::COMM_WORLD.Recv(rhob,sizeOfData,MPI::DOUBLE,from,2);
      MPI::COMM_WORLD.Recv(p,sizeOfData,MPI::DOUBLE,from,17);

      MPI::COMM_WORLD.Recv(utau,sizeOfData,MPI::DOUBLE,from,3);
      MPI::COMM_WORLD.Recv(ux,sizeOfData,MPI::DOUBLE,from,4);
      MPI::COMM_WORLD.Recv(uy,sizeOfData,MPI::DOUBLE,from,5);
      MPI::COMM_WORLD.Recv(ueta,sizeOfData,MPI::DOUBLE,from,6);
      
      MPI::COMM_WORLD.Recv(Wtautau,sizeOfData,MPI::DOUBLE,from,7);
      MPI::COMM_WORLD.Recv(Wtaux,sizeOfData,MPI::DOUBLE,from,8);
      MPI::COMM_WORLD.Recv(Wtauy,sizeOfData,MPI::DOUBLE,from,9);
      MPI::COMM_WORLD.Recv(Wtaueta,sizeOfData,MPI::DOUBLE,from,10);
      MPI::COMM_WORLD.Recv(Wxx,sizeOfData,MPI::DOUBLE,from,11);
      MPI::COMM_WORLD.Recv(Wxy,sizeOfData,MPI::DOUBLE,from,12);
      MPI::COMM_WORLD.Recv(Wxeta,sizeOfData,MPI::DOUBLE,from,13);
      MPI::COMM_WORLD.Recv(Wyy,sizeOfData,MPI::DOUBLE,from,14);
      MPI::COMM_WORLD.Recv(Wyeta,sizeOfData,MPI::DOUBLE,from,15);
      MPI::COMM_WORLD.Recv(Wetaeta,sizeOfData,MPI::DOUBLE,from,16);
	  
      for(ix=0; ix<=nx; ix++)
      {
        for(iy=0; iy<=ny; iy++)
        {
          for(ieta=0; ieta<neta; ieta++)
          {
            position = ieta+(neta*(ix + ((nx+1)*iy)));

            epsFrom[ix][iy][ieta+irank*neta] = eps[position];
            rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
            pFrom[ix][iy][ieta+irank*neta] = p[position];

            utauFrom[ix][iy][ieta+irank*neta] = utau[position];
            uxFrom[ix][iy][ieta+irank*neta] = ux[position];
            uyFrom[ix][iy][ieta+irank*neta] = uy[position];
            uetaFrom[ix][iy][ieta+irank*neta] = ueta[position];
        
            WtautauFrom[ix][iy][ieta+irank*neta] = Wtautau[position];
            WtauxFrom[ix][iy][ieta+irank*neta] = Wtaux[position];
            WtauyFrom[ix][iy][ieta+irank*neta] = Wtauy[position];
            WtauetaFrom[ix][iy][ieta+irank*neta] = Wtaueta[position];
            WxxFrom[ix][iy][ieta+irank*neta] = Wxx[position];
            WxyFrom[ix][iy][ieta+irank*neta] = Wxy[position];
            WxetaFrom[ix][iy][ieta+irank*neta] = Wxeta[position];
            WyyFrom[ix][iy][ieta+irank*neta] = Wyy[position];
            WyetaFrom[ix][iy][ieta+irank*neta] = Wyeta[position];
            WetaetaFrom[ix][iy][ieta+irank*neta] = Wetaeta[position];
	    }
	  }
	}
    }
    //cout << " RECEIVED ALL DATA FROM OTHER ranks" << endl;

    //set up file name
    const string out_name_xyeta = "evolution_xyeta.dat";
    const string out_name_W_xyeta = "evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
    string out_open_mode;
    FILE *out_file_xyeta;
    FILE *out_file_W_xyeta;
    
    //If it's the first timestep, overwrite the previous file
    if (tau == DATA->tau0)
    {
      out_open_mode = "w";
    }
    else
    {
      out_open_mode = "a";	
    }
    //If we output in binary, set the mode accordingly
    if (0 == DATA->outputBinaryEvolution)
    {
      out_open_mode += "b";
    }

    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());
    out_file_W_xyeta = fopen(out_name_W_xyeta.c_str(), out_open_mode.c_str());
    
    // output hydro history
    for(ieta = 0; ieta < (DATA->neta)*size; ieta++)
    {
      double eta = ((double)ieta)*(DATA->delta_eta) - (DATA->eta_size)/2.0;
      for(iy = 0; iy < DATA->ny+1; iy++) //All y
	{
        for(ix = 0; ix<DATA->nx+1; ix++) // All x
        {
          double epsilon1 = epsFrom[ix][iy][ieta];
          double rhob1 = rhobFrom[ix][iy][ieta];

          double utau1 = utauFrom[ix][iy][ieta];
          double ux1 = uxFrom[ix][iy][ieta];
          double uy1 = uyFrom[ix][iy][ieta];
          double ueta1 = uetaFrom[ix][iy][ieta];

          // calculate lab frame velocity
          double u01 = ueta1*sinh(eta)+utau1*cosh(eta); // = gamma factor
          double vx1 = ux1/u01;
          double vy1 = uy1/u01;
          double vz1 = (ueta1*cosh(eta)+utau1*sinh(eta))/u01;
          
          double T1=eos->get_temperature(epsilon1, rhob1);
          double muB1 = eos->get_mu(epsilon1, rhob1);

          double Wtautau = WtautauFrom[ix][iy][ieta];
          double Wtaux = WtauxFrom[ix][iy][ieta];
          double Wtauy = WtauyFrom[ix][iy][ieta];
          double Wtaueta = WtauetaFrom[ix][iy][ieta];
          double Wxx = WxxFrom[ix][iy][ieta];
          double Wxy = WxyFrom[ix][iy][ieta];
          double Wxeta = WxetaFrom[ix][iy][ieta];
          double Wyy = WyyFrom[ix][iy][ieta];
          double Wyeta = WyetaFrom[ix][iy][ieta];
          double Wetaeta = WetaetaFrom[ix][iy][ieta];

          double Wtt = pow(cosh(eta),2)*Wtautau + pow(tau*sinh(eta),2)*Wetaeta*pow(tau,-2) + 2*tau*cosh(eta)*sinh(eta)*Wtaueta*pow(tau,-1);
          double Wtx = cosh(eta)*Wtaux + tau*sinh(eta)*Wxeta*pow(tau,-1);
          double Wty = cosh(eta)*Wtauy + tau*sinh(eta)*Wyeta*pow(tau,-1);
          double Wtz = cosh(eta)*sinh(eta)*Wtautau + tau*( pow(cosh(eta),2) + pow(sinh(eta),2) )*Wtaueta*pow(tau,-1) + tau*tau*cosh(eta)*sinh(eta)*Wetaeta*pow(tau,-2);
          double Wxz = sinh(eta)*Wtaux + tau*cosh(eta)*Wxeta*pow(tau,-1);
          double Wyz = sinh(eta)*Wtauy + tau*cosh(eta)*Wyeta*pow(tau,-1);
          double Wzz = pow(sinh(eta),2)*Wtautau + pow(tau*cosh(eta),2)*Wetaeta*pow(tau,-2) + 2*tau*cosh(eta)*sinh(eta)*Wtaueta*pow(tau,-1);
          
          if (DATA->outputBinaryEvolution == 0)
          {
            fprintf(out_file_xyeta, "%e %e %e %e %e\n", T1*hbarc, muB1*hbarc, vx1, vy1, vz1);
		if (DATA->viscosity_flag == 1)
            {
              fprintf(out_file_W_xyeta, "%e %e %e %e %e %e %e %e %e %e\n", 
                      Wtt, Wtx, Wty, Wtz, Wxx, Wxy, Wxz, Wyy, Wyz, Wzz); 
		}
	    }
	    else
          {
            double array[]={T1*hbarc, muB1*hbarc, vx1, vy1, vz1};
            fwrite(array, sizeof(double), 5, out_file_xyeta);
		if (DATA->viscosity_flag == 1)
            {
              double array2[]={Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz};
              fwrite(array2, sizeof(double), 10, out_file_W_xyeta);
		}
	    }
	  }/* ix */
	}/* iy */
    }/* ieta */
    fclose(out_file_xyeta);
    fclose(out_file_W_xyeta);

    /*End of hydro output in tau,x,y,eta*/
    // clean up
    util->cube_free(epsFrom,nx+1,ny+1,size*neta);
    util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
    util->cube_free(utauFrom,nx+1,ny+1,size*neta);
    util->cube_free(uxFrom,nx+1,ny+1,size*neta);
    util->cube_free(uyFrom,nx+1,ny+1,size*neta);
    util->cube_free(uetaFrom,nx+1,ny+1,size*neta);
    util->cube_free(pFrom,nx+1,ny+1,size*neta);
    
    util->cube_free(WtautauFrom,nx+1,ny+1,size*neta);
    util->cube_free(WtauxFrom,nx+1,ny+1,size*neta);
    util->cube_free(WtauyFrom,nx+1,ny+1,size*neta);
    util->cube_free(WtauetaFrom,nx+1,ny+1,size*neta);
    util->cube_free(WxxFrom,nx+1,ny+1,size*neta);
    util->cube_free(WxyFrom,nx+1,ny+1,size*neta);
    util->cube_free(WxetaFrom,nx+1,ny+1,size*neta);
    util->cube_free(WyyFrom,nx+1,ny+1,size*neta);
    util->cube_free(WyetaFrom,nx+1,ny+1,size*neta);
    util->cube_free(WetaetaFrom,nx+1,ny+1,size*neta);
  } // if(rank == 0)

  // clean up
  free(eps);
  free(rhob);
  free(utau);
  free(ux);
  free(uy);
  free(ueta);
  free(p);

  free(Wtautau);
  free(Wtaux);
  free(Wtauy);
  free(Wtaueta);
  free(Wxx);
  free(Wxy);
  free(Wxeta);
  free(Wyy);
  free(Wyeta);
  free(Wetaeta);
  delete(util);
}/* OutputEvolutionDataXYEta */

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
		  position = ieta+(neta*(ix + ((nx+1)*iy)));
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
		      position = ieta+(neta*(ix + ((nx+1)*iy)));
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
      const char* out_name = "OSCAR.dat";
      out_file = fopen(out_name, "a");
//       fprintf(out_file,"");
//       int iz, nz;
      double T, eta, eps;
//       double z, delta_z, z_size, etafrac;
//       double x, y;
      double ux, uy, ueta, utau, rhob, QGPfrac;
//       double epsilon;
//       double eta_lower, eps_lower, eps_higher, rhob_lower, rhob_higher;
//       double utau_lower, utau_higher, ux_lower, ux_higher, uy_lower, uy_higher, ueta_lower, ueta_higher; 
      double u0, uz;
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
// 	      y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	      for(ix=0; ix<DATA->nx; ix++)
		{
// 		  x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
		  eps = epsFrom[ix][iy][ieta];
		  rhob = rhobFrom[ix][iy][ieta];
		  utau = utauFrom[ix][iy][ieta];
		  ux = uxFrom[ix][iy][ieta];
		  uy = uyFrom[ix][iy][ieta];
		  ueta = uetaFrom[ix][iy][ieta];
		  
		  T = eos->get_temperature(eps, rhob);
		  //QGPfrac = eos->get_qgp_frac(eps, rhob);
              QGPfrac = 0.0;
		  
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
		  position = ieta+(neta*(ix + ((nx+1)*iy)));
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
		      position = ieta+(neta*(ix + ((nx+1)*iy)));
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
      const char* out_name = "contourPlot.dat";
      out_file = fopen(out_name, "a");
//       fprintf(out_file,"");
      int iz, nz;
      double T, z, eta, delta_z, z_size, eps, etafrac;
//       double x, y;
      double ux, uy, ueta, utau, rhob, QGPfrac;
//       double epsilon;
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
// 	      y = iy*(DATA->delta_y) - (DATA->y_size)/2.0;
	      for(ix=0; ix<DATA->nx; ix++) // do all x
		{
// 		  x = ix*(DATA->delta_x) - (DATA->x_size)/2.0;
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
		  
		  T = eos->get_temperature(eps, rhob);
		  //QGPfrac = eos->get_qgp_frac(eps, rhob);
              QGPfrac = 0.0;
		  
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


void Grid::Tmax_profile(Grid ***arena, InitData *DATA, EOS *eos, double tau, int size, int rank)
{
  Util *util;
  util = new Util();

  int position;
  int ix, iy, ieta, nx, ny, neta;
  nx = DATA->nx;
  ny = DATA->ny;
  neta = DATA->neta;

  int sizeOfData = (nx+1)*(ny+1)*(neta);
  int to, from;
  
  double *eps;
  double *rhob;

  eps = new double[sizeOfData];
  rhob = new double[sizeOfData];

   if (rank>0)
    {
      to = 0;
      for(ix=0; ix<=nx; ix++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ieta=0; ieta<neta; ieta++)
		{
		  position = ieta+(neta*(ix + ((nx+1)*iy)));
		  eps[position] = arena[ix][iy][ieta].epsilon;
		  rhob[position] = arena[ix][iy][ieta].rhob;
		}
	    }
	}
    
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
		      position = ieta+(neta*(ix + ((nx+1)*iy)));
		      epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		      rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		    }
		}
	    }
	}

      double e;

      double * e_max_eta = util->vector_malloc(size*neta);
      for(ieta=0; ieta<size*neta; ieta++)
	{
	  for(iy=0; iy<=ny; iy++)
	    {
	      for(ix=0; ix<=nx; ix++)
		{
		  e = epsFrom[ix][iy][ieta];
		  if (e > e_max_eta[ieta]) e_max_eta[ieta] = e;
		}
	    }
	}
      double * e_max_x = util->vector_malloc(nx+1);
      for(ix=0; ix<=nx; ix++)
      {
            for(ieta=0; ieta<size*neta; ieta++)
	    {
	        for(iy=0; iy<=ny; iy++)
	        {
		  e = epsFrom[ix][iy][ieta];
		  if (e > e_max_x[ix]) e_max_x[ix] = e;
		}
	    }
      }
      double * e_max_y = util->vector_malloc(ny+1);
      for(iy=0; iy<=ny; iy++)
      {
            for(ieta=0; ieta<size*neta; ieta++)
	    {
	        for(ix=0; ix<=nx; ix++)
	        {
		  e = epsFrom[ix][iy][ieta];
		  if (e > e_max_y[iy]) e_max_y[iy] = e;
		}
	    }
      }

	int Nfiles = 3;
	string strTmaxAll[] = {"T_max_eta.txt", "T_max_x.txt", "T_max_y.txt"};
	int NRec[] = {neta*size, nx+1, ny+1};
	double * arrRec;
	for (int ifile = 0; ifile < Nfiles; ifile++)
	{
		if (! ( util->fileExists( strTmaxAll[ifile] ) ) ) 
		{
			std::ofstream Tmax (strTmaxAll[ifile].c_str(), std::ofstream::out | std::fstream::app);
			Tmax
				<< "#tau0= " <<  DATA->tau0
				<< " nx= " << DATA->nx 
				<< " ny= " << DATA->ny 
				<< " neta= " << DATA->neta
				<< " x_size= " << DATA->x_size
				<< " y_size= " << DATA->y_size
				<< " eta_size= " << DATA->eta_size
			<< endl;
			Tmax.close();
		}
			
		std::ofstream Tmax (strTmaxAll[ifile].c_str(), std::ofstream::out | std::fstream::app);
		Tmax << std::setprecision (4) << setw(5);

		switch (ifile)
		{
			case 0: arrRec = e_max_eta; break;
			case 1: arrRec = e_max_x; break;
			case 2: arrRec = e_max_y; break;
		}
		for (int iRec = 0; iRec < NRec[ifile]; iRec++) 
		Tmax << eos->get_temperature(arrRec[iRec], 0.)*hbarc << "\t";
		Tmax << endl;
		Tmax.close();
	}
      
      delete [] e_max_x;
      delete [] e_max_y;
      delete [] e_max_eta;
      util->cube_free(epsFrom,nx+1,ny+1,size*neta);
      util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
    }

  delete [] eps;
  delete [] rhob;
  delete(util);
  
}


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
		  position = ieta+(neta*(ix + ((nx+1)*iy)));
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
		      position = ieta+(neta*(ix + ((nx+1)*iy)));
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

//       int iz;
      int nCells=0, nCells2=0;
//       double x, y, z, eta, delta_z, z_size, etafrac;
      double T, eps;
      double totalT=0, totalT2=0;
//       double ux, uy, ueta, utau, epsilon;
      double rhob, QGPfrac;
//       double eta_lower, eps_lower, eps_higher, rhob_lower, rhob_higher;
      for(ieta=0; ieta<neta; ieta++)
	{
	  for(iy=0; iy<DATA->ny; iy++) // do all y
	    {
	      for(ix=0; ix<DATA->nx; ix++) // do all x
		{
		  eps = epsFrom[ix][iy][ieta];
		  rhob = rhobFrom[ix][iy][ieta];
			  
		  T=eos->get_temperature(eps,rhob);
		  //QGPfrac=eos->get_qgp_frac(eps,rhob);
              QGPfrac = 0.0;
		  
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

//   delete [] eps;
//   delete [] rhob;
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
		  position = ieta+(neta*(ix + ((nx+1)*iy)));
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
		      position = ieta+(neta*(ix + ((nx+1)*iy)));
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
      const char* ecc_name = "eccentricity.dat";
      ecc_file = fopen(ecc_name, "a");
      
//       double value, epsp;
      double rA;
//       double numerator, denominator;
      double avcos, avsin, avcos3, avsin3;
      double Psi2, phiA, Psi3, avrSq, avr3;
      double eccentricity2, eccentricity3;

      FILE *out_file;
      const char* out_name = "e_x_y_profile.dat";
      out_file = fopen(out_name, "a");
//       fprintf(out_file,"");
      FILE *s_file;
      const char* s_name = "entropy-eta.dat";
      s_file = fopen(s_name, "a");
//       fprintf(s_file,"");
      FILE *out_file_2;
      const char* out_name_2 = "e_x_y_profile_05.dat";
      out_file_2 = fopen(out_name_2, "a");
//       fprintf(out_file_2,"");
      int iz, nz;
//       double trouble_lower;
      double T, x, y, z, eta, delta_z, z_size, eps, Txx, Tyy, Txy,etafrac;
//       double trouble, trouble_higher;
      //double corr_lower, corr_higher;
      double ux, uy, ueta, utau, rhob, QGPfrac;
//       double epsilon;
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
// 	      trouble = troubleFrom[ix][iy][ieta];
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
// 		  trouble_higher = troubleFrom[ix][iy][ieta+1];
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
// 		  trouble_higher = trouble_lower;
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
// 	      trouble = trouble_lower * (1.-etafrac) + (etafrac)*trouble_higher;
	      //corr = corr_lower * (1.-etafrac) + (etafrac)*corr_higher;
	      Txx = Txx_lower * (1.-etafrac) + (etafrac)*Txx_higher;
	      Tyy = Tyy_lower * (1.-etafrac) + (etafrac)*Tyy_higher;
	      Txy = Txy_lower * (1.-etafrac) + (etafrac)*Txy_higher;
	      rhob = rhob_lower * (1.-etafrac) + (etafrac)*rhob_higher;
	      utau = utau_lower * (1.-etafrac) + (etafrac)*utau_higher;
	      ux = ux_lower * (1.-etafrac) + (etafrac)*ux_higher;
	      uy = uy_lower * (1.-etafrac) + (etafrac)*uy_higher;
	      ueta = ueta_lower * (1.-etafrac) + (etafrac)*ueta_higher;
	      
	      T = eos->get_temperature(eps, rhob);
	      //QGPfrac = eos->get_qgp_frac(eps, rhob);
            QGPfrac = 0.0;
	      
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
//       numerator = 0.;
//       denominator = 0.;
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
//  double x, y;
 double eta;
//  double ux, uy, ueta, utau, rhob;
 double epsilon;
 string d_name;
 ofstream d_file;

 if(rank==0) d_name = "e_profile.dat";
 else d_name = "e_profile2.dat";
 d_file.open(d_name.c_str(), ios::out | ios::app );

//  x=0.;
//  y=0.;
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
 double x, y;
//  double ux, uy, ueta, utau, rhob, eta;
 double epsilon;
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
 int r, i, a;
//  int m, n, j;

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
 const char* d_name = "e_x_y_profile.dat";
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

void Grid::print_qmu_evolution(InitData *DATA, Grid ***arena, double tau, EOS *eos, int rank)
{
  string d_name;
  ofstream d_file;
  d_name = "qmu_evo.dat";
  d_file.open(d_name.c_str(), ios::out | ios::app );
  
  int ieta = DATA->neta/2;
  double eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
  int iy = (int)(DATA->ny/2);
  for(int ix=0; ix<=DATA->nx; ix++)
  {
    double x = ix*(DATA->delta_x) - (DATA->x_size/2.0);
    double y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
    double rhob = arena[ix][iy][ieta].rhob;
    double epsilon = arena[ix][iy][ieta].epsilon;
    double muB = eos->get_mu(epsilon, rhob);
    double kappa = 0.2*rhob/muB;
    d_file << scientific << setw(18) << setprecision(8)
           << tau << "  " << x << "  " << y << "  " << eta << "  " 
           << rhob << "  " ;
    double D_mu = arena[ix][iy][ieta].a[0][4];
    for(int nu = 0; nu < 4; nu++)
    {
        double q_nu = arena[ix][iy][ieta].Wmunu[0][4][nu];
        double partial_mu = arena[ix][iy][ieta].dUsup[0][4][nu];
        double flow_u = arena[ix][iy][ieta].u[0][nu];
        double NS_term = -kappa*(partial_mu + flow_u*D_mu);
        d_file << scientific << setw(18) << setprecision(8)
               << q_nu << "  " << NS_term << "  " ;
               //<< partial_mu << "  " << flow_u*D_mu << "  " << NS_term << "  " ;
    }
    d_file << endl;
  }
  int ix = (int)(DATA->ny/2);
  for(int iy=0; iy<=DATA->ny; iy++)
  {
    double x = ix*(DATA->delta_x) - (DATA->x_size/2.0);
    double y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
    double rhob = arena[ix][iy][ieta].rhob;
    double epsilon = arena[ix][iy][ieta].epsilon;
    double muB = eos->get_mu(epsilon, rhob);
    double kappa = 0.2*rhob/muB;
    d_file << scientific << setw(18) << setprecision(8)
           << tau << "  " << x << "  " << y << "  " << eta << "  "
           << rhob << "  " ;
    double D_mu = arena[ix][iy][ieta].a[0][4];
    for(int nu = 0; nu < 4; nu++)
    {
        double q_nu = arena[ix][iy][ieta].Wmunu[0][4][nu];
        double partial_mu = arena[ix][iy][ieta].dUsup[0][4][nu];
        double flow_u = arena[ix][iy][ieta].u[0][nu];
        double NS_term = -kappa*(partial_mu + flow_u*D_mu);
        d_file << scientific << setw(18) << setprecision(8)
               << q_nu << "  " << NS_term << "  " ;
               //<< partial_mu << "  " << flow_u*D_mu << "   " ;
    }
    d_file << endl;
  }
  d_file.close();
  
  return;
}/* print qmu evolution */

void Grid::print_rhob_evolution(InitData *DATA, Grid ***arena, double tau, EOS* eos, int rank)
{
  string d_name;
  ofstream d_file;
  d_name = "rhoB_evo.dat";
  d_file.open(d_name.c_str(), ios::out | ios::app );
  
  int ieta = DATA->neta/2;
  double eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
  for(int ix=0; ix<=DATA->nx; ix++)
  {
    double x = ix*(DATA->delta_x) - (DATA->x_size/2.0);
    int iy = (int)(DATA->ny/2);
    double y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
    double rhob = arena[ix][iy][ieta].rhob;
    double ed = arena[ix][iy][ieta].epsilon;
    double mub = eos->get_mu(ed, rhob);
    double temperature = eos->get_temperature(ed, rhob);
    d_file << scientific << setw(18) << setprecision(8)
           << tau << "  " << x << "  " << y << "  " << eta << "  "
           << ed << "  " << rhob << "  " 
           << mub << "  " << temperature << "  "
           << mub/temperature << endl;
  }
  for(int iy=0; iy<=DATA->ny; iy++)
  {
    double y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
    int ix = (int)(DATA->nx/2);
    double x = ix*(DATA->delta_x) - (DATA->x_size/2.0);
    double rhob = arena[ix][iy][ieta].rhob;
    double ed = arena[ix][iy][ieta].epsilon;
    double mub = eos->get_mu(ed, rhob);
    double temperature = eos->get_temperature(ed, rhob);
    d_file << scientific << setw(18) << setprecision(8)
           << tau << "  " << x << "  " << y << "  "  << eta << "  "
           << ed << "  " << rhob << "  "
           << mub << "  " << temperature << "  "
           << mub/temperature << endl;
  }
  d_file.close();
  
  return;
}/* print rhob evolution */

void Grid::print_fireball_evolution_on_phasediagram(InitData *DATA, Grid ***arena, double tau, EOS* eos, int rank)
{
  int n_skip = 25;

  int nT = 401;
  int nmu = 801;
  double T0 = 0.0;
  double T_max = 0.4;
  double dT = (T_max - T0)/(nT - 1);
  double mu0 = 0.0;
  double mu_max = 0.8;
  double dmu = (mu_max - mu0)/(nmu - 1);

  double eta_goal = 0.0;
  double e_dec = 0.1;   // GeV/fm^3
  double volume_element = DATA->delta_x * DATA->delta_y;

  if((int)((tau - DATA->tau0)/DATA->delta_tau) % n_skip == 0)
  {
     for(int ieta = 0; ieta < DATA->neta; ieta++)
     {
        double eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;
        if(fabs(eta - eta_goal) < 1e-5)  // mid rapidity
        {
            double** dVdeta = new double* [nT];
            for(int iT = 0; iT < nT; iT++)
            {
               dVdeta[iT] = new double [nmu];
               for(int imu = 0; imu < nmu; imu++)
                  dVdeta[iT][imu] = 0.0;
            }

            ostringstream filename;
            filename << "fireball_evo_phasediagram_tau_" << tau << ".dat";
            ofstream d_file;
            d_file.open(filename.str().c_str());

            for(int ix=0; ix<=DATA->nx; ix++)
            {
                for(int iy = 0; iy <= DATA->ny; iy++)
                {
                    double rhob = arena[ix][iy][ieta].rhob;
                    double ed = arena[ix][iy][ieta].epsilon;
                    if(ed*hbarc > e_dec)
                    {
                       double mub = eos->get_mu(ed, rhob)*hbarc;   // GeV
                       double temperature = eos->get_temperature(ed, rhob)*hbarc;  // GeV

                       int T_idx = (int)((temperature - T0)/dT);
                       int mu_idx = (int)((mub - mu0)/dmu);
                       
                       if(T_idx >= 0 && T_idx < nT)
                          if(mu_idx >=0 && mu_idx < nmu)
                             dVdeta[T_idx][mu_idx] += volume_element;
                    }
               }
            }

            for(int iT = 0; iT < nT; iT++)
            {
               for(int imu = 0; imu < nmu; imu++)
               {
                   d_file << scientific << setw(18) << setprecision(8)
                          << dVdeta[iT][imu] << "   ";
               }
               d_file << endl;
            }
            d_file.close();

            for(int iT = 0; iT < nT; iT++)
                delete [] dVdeta[iT];
            delete [] dVdeta;
        }
     }
  }
  return;
}

void Grid::ComputeAnisotropy(InitData *DATA, Grid ***arena, double tau)
{
  
  string d_name;
  ofstream d_file;
  d_name = "aniso.dat";
  d_file.open(d_name.c_str(), ios::out | ios::app );
  
  int ix, iy, ieta;
  double v2;
//   double x, y, epsp;
  double numerator, denominator;
  numerator = 0.;
  denominator = 0.;
  
  ieta = DATA->neta/2;
  for(ix=0; ix<=DATA->nx; ix++)
    {
//       x = ix*(DATA->delta_x) - (DATA->x_size/2.0);  
      for(iy=0; iy<=DATA->ny; iy++)
	{
// 	  y = iy*(DATA->delta_x) - (DATA->y_size/2.0);
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
  const char* ecc_name = "eccentricity.dat";
  ecc_file = fopen(ecc_name, "a");
  
  int ix, iy, ieta;
  double value, x, y;
//   double epsp;
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
  const char* ecc_name = "energyConservation.dat";
  ecc_file = fopen(ecc_name, "a");
  
  int ix, iy, ieta;
  double value;
//   double x, y, epsp, eta;
  double TtautauPart, TetaetaPart;//, TetaetaPrev;
  TtautauPart = 0.;
  TetaetaPart = 0.;
//   TetaetaPrev = 0.;
  
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
 int ix, iy, ieta;
//  int i;
 double eta, f, g, h, l, k;
//  double x, y;

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

void Grid::Output_hydro_information_header(InitData *DATA, EOS *eos) {

	string fname = "hydro_info_header_h";

	//Open output file
	ofstream outfile;
	outfile.open(fname.c_str());

	outfile << "const int MUSIC_real_nx=" << DATA->nx+1 << ";\n";
	outfile << "const int MUSIC_real_ny=" << DATA->ny+1 << ";\n";
	//DATA->neta is _not_ the actual number of cells in eta, it is the number of cells in eta _per_ processor
	outfile << "const int MUSIC_real_neta=" << DATA->neta*MPI::COMM_WORLD.Get_size() << ";\n";

	//double x_size; /* in fermi -x_size/2 < x < x_size/2 */
	//double y_size; /* in fermi, ditto */
	//double eta_size; /* ditto */

	outfile << "const double MUSIC_tau0=" << DATA->tau0 << ";\n";

	outfile << "const double MUSIC_dx=" << DATA->delta_x << ";\n";
	outfile << "const double MUSIC_dy=" << DATA->delta_y << ";\n";
	outfile << "const double MUSIC_deta=" << DATA->delta_eta << ";\n";
	outfile << "const double MUSIC_effective_dtau=" << DATA->output_evolution_every_N_timesteps*DATA->delta_tau << ";\n";

	outfile << "const bool MUSIC_with_shear_viscosity=" << ((DATA->viscosity_flag)&&(DATA->turn_on_shear)) << ";\n";
	outfile << "const bool MUSIC_with_bulk_viscosity=" << ((DATA->viscosity_flag)&&(DATA->turn_on_bulk)) << ";\n";

	outfile << "const double MUSIC_kinetic_FO_temperature_in_GeV=";
	if (DATA->useEpsFO) {
		outfile << eos->get_temperature(DATA->epsilonFreeze/hbarc,0.0)*hbarc;
	}
	else {
		outfile << DATA->TFO;
	}
	outfile << ";\n";

	outfile << "const bool MUSIC_outputBinaryEvolution=" << DATA->outputBinaryEvolution << ";\n";

	outfile.close();

}
