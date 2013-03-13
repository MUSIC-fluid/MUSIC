#include "util.h"
#include "data.h"
#include "grid.h"
#include "minmod.h"
#include "u_derivative.h"

using namespace std;

U_derivative::U_derivative()
{
  minmod = new Minmod;
}

// destructor
U_derivative::~U_derivative()
{
  delete minmod;
}

int U_derivative::MakedU(double tau, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int rk_flag, int size, int rank)
{
 int ix, iy, ieta, nx, ny, neta, flag, ic, mu, nu;
 double g, f, h, tfac, gfac, temp_nu;
 Grid *grid_pt;

 if(DATA->viscosity_flag == 0) return 1;

 nx = DATA->nx;
 ny = DATA->ny;
 neta = DATA->neta-1;

 cout << "";

 for(ix=0; ix<=nx; ix++)
  {
   for(iy=0; iy<=ny; iy++)
    {
     for(ieta=0; ieta<=neta; ieta++)
      {
	/* this calculates du/dx, du/dy, (du/deta)/tau */
   
	flag = MakeDSpatial(tau, DATA, &(arena[ix][iy][ieta]), &(Lneighbor[ix][iy][0]), &(Rneighbor[ix][iy][0]), 
			    &(Lneighbor[ix][iy][1]), &(Rneighbor[ix][iy][1]), rk_flag, size, rank); 

	/* this calculates du/dtau */
        flag = MakeDTau(tau, DATA, &(arena[ix][iy][ieta]), rk_flag); 
      
	grid_pt= &(arena[ix][iy][ieta]);
	
	/* TEST TEST */
// 	temp_nu = -(grid_pt->u[rk_flag][0])*(grid_pt->dUsup[rk_flag][0][0]);
// 	temp_nu += (grid_pt->u[rk_flag][1])*(grid_pt->dUsup[rk_flag][1][0]);
// 	temp_nu += (grid_pt->u[rk_flag][2])*(grid_pt->dUsup[rk_flag][2][0]);
// 	temp_nu += (grid_pt->u[rk_flag][3])*(grid_pt->dUsup[rk_flag][3][0]);
// 	if(fabs(temp_nu) > 0.01)
// 	 { 
// 	  fprintf(stderr, "partial u^2 = %e\n", temp_nu);
// 	  fprintf(stderr, "This can't happen.\n");
// 	  exit(0);
// 	 }
	/* Test TEST */
      }/* ieta */
    }/*iy */
  }/* ix */

 // cout << "first part done" << endl;

 // cout << "second part" << endl;

/* theta_u = partial_mu u^mu */
 for(ix=0; ix<=nx; ix++)
  {
   for(iy=0; iy<=ny; iy++)
    {
     for(ieta=0; ieta<=neta; ieta++)
      {
       g = 0.0;
       for(mu=0; mu<4; mu++)
        {
	 if(mu==0) gfac = -1.0;
	 else gfac = 1.0;

	 g += arena[ix][iy][ieta].dUsup[rk_flag][mu][mu]*gfac;
	}/* mu */
/* need to add the tau-eta coordinate source term */
       g += arena[ix][iy][ieta].u[rk_flag][0]/tau;
       arena[ix][iy][ieta].theta_u[rk_flag] = g;
      
      }/* ieta */
    }/*iy */
  }/* ix */

 // cout << "second part done " << endl;

 // cout << "third part" << endl;

/* a[mu] = u^nu partial_nu u_mu */
 for(ix=0; ix<=nx; ix++)
  {
   for(iy=0; iy<=ny; iy++)
    {
     for(ieta=0; ieta<=neta; ieta++)
      {
       for(mu=0; mu<4; mu++)
        {
         f = 0.0;
	 for(ic=0; ic<4; ic++)
          {
	   /* eta derivative already has 1/tau 
	   remember: du[m][n] = u^{m,n} = partial^n u^m 
	   so g_{ln} u^l\partial^n u^m*/

	   tfac = (ic==0 ? -1.0 : 1.0);
	  
           h =  (arena[ix][iy][ieta].u[rk_flag][ic]);
           /*
	   h += (arena[ix][iy][ieta].prev_u[rk_flag][ic]);
	   h *= 0.5;
	   */
	   h *= (arena[ix][iy][ieta].dUsup[rk_flag][mu][ic]);
	     
	   f += h*tfac;
	  }/* ic  */
	 arena[ix][iy][ieta].a[rk_flag][mu] = f;
	
	}/* mu */
      }/* ieta */
    }/*iy */
  }/* ix */
 // cout << "third part done" << endl;

 return 1; /* successful */
}/* MakedU */


int U_derivative::MakeDSpatial(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			       Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank)
{
 int nmax[4], m, n;
 double g, f, fp1, fm1, taufactor;
 double delta[4];
 
 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta-1;

 delta[1] = DATA->delta_x;
 delta[2] = DATA->delta_y;
 delta[3] = DATA->delta_eta;

/* partial_n u_m */
/* for u[i] */
 for(m=1; m<=3; m++)
  {
   for(n=1; n<=2; n++)
    {
     taufactor = 1.0;
     f = grid_pt->u[rk_flag][m];
     if(grid_pt->position[n] == nmax[n]) 
      {
       fp1 = grid_pt->u[rk_flag][m];
       fm1 = grid_pt->nbr_m_1[n]->u[rk_flag][m];
      }
     else if(grid_pt->position[n] == 0) 
      {
       fp1 = grid_pt->nbr_p_1[n]->u[rk_flag][m];
       fm1 = grid_pt->u[rk_flag][m];
      }
     else
      {
       fp1 = grid_pt->nbr_p_1[n]->u[rk_flag][m];
       fm1 = grid_pt->nbr_m_1[n]->u[rk_flag][m];
      }

     g = minmod->minmod_dx(fp1, f, fm1, DATA);
     g /= delta[n]*taufactor;
     grid_pt->dUsup[rk_flag][m][n] = g;
    
    }

   n=3;
   taufactor = tau;
   f = grid_pt->u[rk_flag][m];
   if(grid_pt->position[n] == nmax[n]) 
     {
       if (rank==size-1)
	 {
	   fp1 = grid_pt->u[rk_flag][m];
	   fm1 = grid_pt->nbr_m_1[n]->u[rk_flag][m];
	 }
       else
	 {
	   fp1 = Rneighbor->u[rk_flag][m];
	   fm1 = grid_pt->nbr_m_1[n]->u[rk_flag][m];
	 }
     }
   else if(grid_pt->position[n] == 0) 
     {
       if(rank==0)
	 {
	   fp1 = grid_pt->nbr_p_1[n]->u[rk_flag][m];
	   fm1 = grid_pt->u[rk_flag][m];
	 }
       else
	 {
	   fp1 = grid_pt->nbr_p_1[n]->u[rk_flag][m];
	   fm1 = Lneighbor->u[rk_flag][m];
	 } 
     }
   else
     {
       fp1 = grid_pt->nbr_p_1[n]->u[rk_flag][m];
       fm1 = grid_pt->nbr_m_1[n]->u[rk_flag][m];
     }

   g = minmod->minmod_dx(fp1, f, fm1, DATA);
   g /= delta[n]*taufactor;
   grid_pt->dUsup[rk_flag][m][n] = g;
   
  }
 

 /* for u[0], use u[0]u[0] = 1 + u[i]u[i] */
 /* u[0]_m = u[i]_m (u[i]/u[0]) */
 
 /* for u[0] */
 
 for(n=1; n<=3; n++)
   {
     f = 0.0;
     for(m=1; m<=3; m++)
       {
	 /* (partial_n u^m) u[m] */
	 f += (grid_pt->dUsup[rk_flag][m][n])*(grid_pt->u[rk_flag][m]);
       }
     f /= grid_pt->u[rk_flag][0];
     grid_pt->dUsup[rk_flag][0][n] = f;
   }
 
 return 1;
}/* MakeDSpatial */


int U_derivative::MakeDTau(double tau, InitData *DATA, Grid *grid_pt, int rk_flag)
{
 int m, n;
 double g, f;

/* this makes dU[m][0] = partial^tau u^m */
/* note the minus sign at the end because of g[0][0] = -1 */
/* rk_flag is 0, 2, 4, ... */

if(rk_flag == 0)
{
 for(m=1; m<=3; m++)
  {
   /* first order is more stable */
   f = (grid_pt->u[rk_flag][m]);
   f -= (grid_pt->prev_u[0][m]);
   f /= (DATA->delta_tau);
   
   grid_pt->dUsup[rk_flag][m][0] = -f; /* g00 = -1 */
  }/* m */
}/* rk_flag = 0 */
else if(rk_flag > 0)
{
 for(m=1; m<=3; m++)
  {
   /* first order */
   f = (grid_pt->u[rk_flag][m]); // this is from the prev full RK step 
   f -= (grid_pt->u[0][m]);
   f /= (DATA->delta_tau);
   grid_pt->dUsup[rk_flag][m][0] = -f; /* g00 = -1 */
  }/* m */
}

 /* I have now partial^tau u^i */
 /* I need to calculate (u^i partial^tau u^i) = u^0 partial^tau u^0 */
 /* u_0 d^0 u^0 + u_m d^0 u^m = 0 */
 /* -u^0 d^0 u^0 + u_m d^0 u^m = 0 */
 /* d^0 u^0 = u_m d^0 u^m/u^0 */

   f = 0.0;
   for(m=1; m<=3; m++)
    {
     /* (partial_0 u^m) u[m] */
     f += (grid_pt->dUsup[rk_flag][m][0])*(grid_pt->u[rk_flag][m]);
    }
   f /= grid_pt->u[rk_flag][0];
   grid_pt->dUsup[rk_flag][0][0] = f;
 
 return 1;
}/* MakeDTau */


