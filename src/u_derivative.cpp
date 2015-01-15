#include "util.h"
#include "data.h"
#include "grid.h"
#include "minmod.h"
#include "eos.h"
#include "u_derivative.h"

using namespace std;

// Sangyong Nov 18 2014
// added EOS in the argument 
U_derivative::U_derivative(EOS *eosIn, InitData* DATA_in)
{
   //Sangyong Nov 18 2014: added eos
   //eos = new EOS;
   eos = eosIn;
   minmod = new Minmod(DATA_in);
}

// destructor
U_derivative::~U_derivative()
{
   delete minmod;
   //Sangyong Nov 18 2014: added delete eos
   //delete eos;
}

int U_derivative::MakedU(double tau, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int rk_flag, int size, int rank)
{
   //int flag;
   //double temp_nu, nu;
   //Grid *grid_pt;

   // ideal hydro: no need to evaluate any flow derivatives
   if(DATA->viscosity_flag == 0)
       return 1;  

   int nx = DATA->nx;
   int ny = DATA->ny;
   int neta = DATA->neta-1;

   //cout << "";

   for(int ix=0; ix<=nx; ix++)
   {
      for(int iy=0; iy<=ny; iy++)
      {
         for(int ieta=0; ieta<=neta; ieta++)
         {
	      /* this calculates du/dx, du/dy, (du/deta)/tau */
            MakeDSpatial(tau, DATA, &(arena[ix][iy][ieta]), 
                         &(Lneighbor[ix][iy][0]), &(Rneighbor[ix][iy][0]), 
			       &(Lneighbor[ix][iy][1]), &(Rneighbor[ix][iy][1]), rk_flag, size, rank);
            /* this calculates du/dtau */
            MakeDTau(tau, DATA, &(arena[ix][iy][ieta]), rk_flag); 
      
            //grid_pt= &(arena[ix][iy][ieta]);
            /* TEST TEST */
            //temp_nu = -(grid_pt->u[rk_flag][0])*(grid_pt->dUsup[rk_flag][0][0]);
            //temp_nu += (grid_pt->u[rk_flag][1])*(grid_pt->dUsup[rk_flag][1][0]);
            //temp_nu += (grid_pt->u[rk_flag][2])*(grid_pt->dUsup[rk_flag][2][0]);
            //temp_nu += (grid_pt->u[rk_flag][3])*(grid_pt->dUsup[rk_flag][3][0]);
            //if(fabs(temp_nu) > 0.01)
            // { 
            //  fprintf(stderr, "partial u^2 = %e\n", temp_nu);
            //  fprintf(stderr, "This can't happen.\n");
            //  exit(0);
            // }
	      /* Test TEST */
         }/* ieta */
      }/*iy */
   }/* ix */

   //cout << "first part done" << endl;
   //cout << "second part" << endl;

   for(int ix=0; ix<=nx; ix++)
   {
      for(int iy=0; iy<=ny; iy++)
      {
         for(int ieta=0; ieta<=neta; ieta++)
         {
            double u_local[4];
            for(int i = 0; i < 4; i++)
               u_local[i] = arena[ix][iy][ieta].u[rk_flag][i];
            double dUsup_local[5][4];   // dUsup[m][n] = partial^n u^m
            for(int i = 0; i < 5; i++)
               for(int j = 0; j < 4; j++)
               {
                  dUsup_local[i][j] = arena[ix][iy][ieta].dUsup[rk_flag][i][j];
               }
            double partial_mu_u_supmu = 0.0;
            for(int mu=0; mu<4; mu++)
            {
               double gfac = (mu==0 ? -1.0 : 1.0);
               partial_mu_u_supmu += dUsup_local[mu][mu]*gfac; // for expansion rate: theta

               double u_supnu_partial_nu_u_supmu = 0.0;
	         for(int nu=0; nu<4; nu++)
               {
                  double tfac = (nu==0 ? -1.0 : 1.0);
                  u_supnu_partial_nu_u_supmu += tfac*u_local[nu]*dUsup_local[mu][nu];
               }
               arena[ix][iy][ieta].a[rk_flag][mu] = u_supnu_partial_nu_u_supmu;
	      }/* mu */
            // derivatives of chemical potentials for diffusion
            int mu = 4;
            double u_supnu_partial_nu_muoverT = 0.0;
	      for(int nu=0; nu<4; nu++)
            {
               double tfac = (nu==0 ? -1.0 : 1.0);
               u_supnu_partial_nu_muoverT += tfac*u_local[nu]*dUsup_local[mu][nu];
            }
            arena[ix][iy][ieta].a[rk_flag][mu] = u_supnu_partial_nu_muoverT;

            //double check = 0.0;
            //for(int nu = 0; nu < 4; nu++)
            //{
            //    double tfac = (nu==0 ? -1.0 : 1.0);
            //    check += tfac*u_local[nu]*(dUsup_local[nu][4] + u_local[nu]*u_supnu_partial_nu_muoverT);
            //}
            //cout << check << endl;

            // expansion rate need to add the tau-eta coordinate source term
            arena[ix][iy][ieta].theta_u[rk_flag] = partial_mu_u_supmu + u_local[0]/tau;
         }/* ieta */
      }/*iy */
   }/* ix */

   // calculate the velocity shear tensor sigma^{\mu\nu}
   // immigrate the code from dissipative.cpp
   for(int ix=0; ix<=nx; ix++)
   {
      for(int iy=0; iy<=ny; iy++)
      {
         for(int ieta=0; ieta<=neta; ieta++)
         {
            double dUsup_local[4][4];
            double u_local[4];
            double a_local[4];
            double theta_u_local = arena[ix][iy][ieta].theta_u[rk_flag];
            double sigma_local[4][4];
            for(int i = 0; i < 4; i++)
            {
               u_local[i] = arena[ix][iy][ieta].u[rk_flag][i];
               a_local[i] = arena[ix][iy][ieta].a[rk_flag][i];
               for(int j = 0; j < 4; j++)
               {
                  dUsup_local[i][j] = arena[ix][iy][ieta].dUsup[rk_flag][i][j];
               }
            }
            for(int a=0; a<4; a++)
            {
               for(int b=0; b<4; b++)
               {
                  sigma_local[a][b] = (
                      (dUsup_local[a][b] + dUsup_local[b][a])/2.
                      - (DATA->gmunu[a][b] + u_local[a]*u_local[b])*theta_u_local/3.
                      + u_local[0]/tau*DATA->gmunu[a][3]*DATA->gmunu[b][3]
                      - u_local[3]/tau/2.*(DATA->gmunu[a][3]*DATA->gmunu[b][0] 
                                           + DATA->gmunu[b][3]*DATA->gmunu[a][0])
                      + u_local[3]*u_local[0]/tau/2.*(DATA->gmunu[a][3]*u_local[b] 
                                                      + DATA->gmunu[b][3]*u_local[a])
                      - u_local[3]*u_local[3]/tau/2.*(DATA->gmunu[a][0]*u_local[b] 
                                                      + DATA->gmunu[b][0]*u_local[a])
                      + (u_local[a]*a_local[b] + u_local[b]*a_local[a])/2.);
                  arena[ix][iy][ieta].sigma[rk_flag][a][b] = sigma_local[a][b];
               }
            }
         }/* ieta */
      }/*iy */
   }/* ix */

   return 1; /* successful */
}/* MakedU */


int U_derivative::MakeDSpatial(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			       Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank)
{
 int nmax[4], m, n;
 double g, f, fp1, fm1, taufactor;
 double delta[4];
 //Sangyong Nov 18 2014: added these doubles
 double rhob, eps, muB, T;
 
 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta-1;

 delta[1] = DATA->delta_x;
 delta[2] = DATA->delta_y;
 delta[3] = DATA->delta_eta;

/* dUsup[m][n] = partial_n u_m */
/* for u[i] */
 for(m=1; m<=3; m++)
  {
   // partial_n u[m]
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

     g = minmod->minmod_dx(fp1, f, fm1);
     g /= delta[n]*taufactor;
     grid_pt->dUsup[rk_flag][m][n] = g;
    }// n=x,y
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

   g = minmod->minmod_dx(fp1, f, fm1);
   g /= delta[n]*taufactor;
   grid_pt->dUsup[rk_flag][m][n] = g;
   
  }// m=x,y,eta
 

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
 
// Sangyong Nov 18 2014
// Here we make derivatives of muB/T
// dUsup[rk_flag][4][n] = partial_n (muB/T)
// partial_x (muB/T) and partial_y (muB/T) first

   m = 4; // means (muB/T)
   for(n=1; n<=2; n++)
   {
     taufactor = 1.0;
     // f = grid_pt->rhob_t;
     if(rk_flag == 0)
     {
       rhob = grid_pt->rhob;
       eps = grid_pt->epsilon;
     }
     else
     {
       rhob = grid_pt->rhob_t;
       eps = grid_pt->epsilon_t;
     }
     muB = eos->get_mu(eps, rhob);
     T = eos->get_temperature(eps, rhob);
     f = muB/T; 
    
     if(grid_pt->position[n] == nmax[n]) 
     {
       // fp1 = grid_pt->rhob;
       if(rk_flag == 0)
       {
         rhob = grid_pt->rhob;
         eps = grid_pt->epsilon;
       }
       else
       {
         rhob = grid_pt->rhob_t;
         eps = grid_pt->epsilon_t;
       }
       muB = eos->get_mu(eps, rhob);
       T = eos->get_temperature(eps, rhob);
       fp1 = muB/T; 
       
       // fm1 = grid_pt->nbr_m_1[n]->rhob;
       if(rk_flag == 0)
       {
         rhob = grid_pt->nbr_m_1[n]->rhob;
         eps = grid_pt->nbr_m_1[n]->epsilon;
       }
       else
       {
         rhob = grid_pt->nbr_m_1[n]->rhob_t;
         eps = grid_pt->nbr_m_1[n]->epsilon_t;
       }
       muB = eos->get_mu(eps, rhob);
       T = eos->get_temperature(eps, rhob);
       fm1 = muB/T; 
     }
     else if(grid_pt->position[n] == 0) 
     {
       // fp1 = grid_pt->nbr_p_1[n]->rhob;
       if(rk_flag == 0)
       {
         rhob = grid_pt->nbr_p_1[n]->rhob;
         eps = grid_pt->nbr_p_1[n]->epsilon;
       }
       else
       {
         rhob = grid_pt->nbr_p_1[n]->rhob_t;
         eps = grid_pt->nbr_p_1[n]->epsilon_t;
       }
       muB = eos->get_mu(eps, rhob);
       T = eos->get_temperature(eps, rhob);
       fp1 = muB/T; 
       
       // fm1 = grid_pt->rhob;
       if(rk_flag == 0)
       {
         rhob = grid_pt->rhob;
         eps = grid_pt->epsilon;
       }
       else
       {
         rhob = grid_pt->rhob_t;
         eps = grid_pt->epsilon_t;
       }
       muB = eos->get_mu(eps, rhob);
       T = eos->get_temperature(eps, rhob);
       fm1 = muB/T; 
     }
     else
     {
       //fp1 = grid_pt->nbr_p_1[n]->rhob;
       if(rk_flag == 0)
       {
         rhob = grid_pt->nbr_p_1[n]->rhob;
         eps = grid_pt->nbr_p_1[n]->epsilon;
       }
       else
       {
         rhob = grid_pt->nbr_p_1[n]->rhob_t;
         eps = grid_pt->nbr_p_1[n]->epsilon_t;
       }
       muB = eos->get_mu(eps, rhob);
       T = eos->get_temperature(eps, rhob);
       fp1 = muB/T; 
       
       // fm1 = grid_pt->nbr_m_1[n]->rhob;
       if(rk_flag == 0)
       {
         rhob = grid_pt->nbr_m_1[n]->rhob;
         eps = grid_pt->nbr_m_1[n]->epsilon;
       }
       else
       {
         rhob = grid_pt->nbr_m_1[n]->rhob_t;
         eps = grid_pt->nbr_m_1[n]->epsilon_t;
       }
       muB = eos->get_mu(eps, rhob);
       T = eos->get_temperature(eps, rhob);
       fm1 = muB/T; 
     }

     g = minmod->minmod_dx(fp1, f, fm1);
     g /= delta[n]*taufactor;
     grid_pt->dUsup[rk_flag][m][n] = g;
   }// n=x,y
  
  // eta derivative
   n=3; // means eta
   taufactor = tau;
//   f = grid_pt->rhob_t;
   if(rk_flag == 0)
   {
     rhob = grid_pt->rhob;
     eps = grid_pt->epsilon;
   }
   else
   {
     rhob = grid_pt->rhob_t;
     eps = grid_pt->epsilon_t;
   }
   muB = eos->get_mu(eps, rhob);
   T = eos->get_temperature(eps, rhob);
   f = muB/T; 
   
   if(grid_pt->position[n] == nmax[n]) 
   {
     if (rank==size-1)
     {
      // fp1 = grid_pt->rhob;
        if(rk_flag == 0)
        {
          rhob = grid_pt->rhob;
          eps = grid_pt->epsilon;
        }
        else
        {
          rhob = grid_pt->rhob_t;
          eps = grid_pt->epsilon_t;
        }
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        fp1 = muB/T; 
     
      // fm1 = grid_pt->nbr_m_1[n]->rhob;
        if(rk_flag == 0)
        {
          rhob = grid_pt->nbr_m_1[n]->rhob;
          eps = grid_pt->nbr_m_1[n]->epsilon;
        }
        else
        {
          rhob = grid_pt->nbr_m_1[n]->rhob_t;
          eps = grid_pt->nbr_m_1[n]->epsilon_t;
        }
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        fm1 = muB/T; 
     }
     else
     {
      // fp1 = Rneighbor->rhob;
        if(rk_flag == 0)
        {
          rhob = Rneighbor->rhob;
          eps = Rneighbor->epsilon;
        }
        else
        {
          rhob = Rneighbor->rhob_t;
          eps = Rneighbor->epsilon_t;
        }
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        fp1 = muB/T; 
      
      // fm1 = grid_pt->nbr_m_1[n]->rhob;
        if(rk_flag == 0)
        {
          rhob = grid_pt->nbr_m_1[n]->rhob;
          eps = grid_pt->nbr_m_1[n]->epsilon;
        }
        else
        {
          rhob = grid_pt->nbr_m_1[n]->rhob_t;
          eps = grid_pt->nbr_m_1[n]->epsilon_t;
        }
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        fm1 = muB/T; 
     }
   }
   else if(grid_pt->position[n] == 0) 
   {
     if(rank==0)
     {
       // fp1 = grid_pt->nbr_p_1[n]->rhob;
         if(rk_flag == 0)
         {
           rhob = grid_pt->nbr_p_1[n]->rhob;
           eps = grid_pt->nbr_p_1[n]->epsilon;
         }
         else
         {
           rhob = grid_pt->nbr_p_1[n]->rhob_t;
           eps = grid_pt->nbr_p_1[n]->epsilon_t;
         }
         muB = eos->get_mu(eps, rhob);
         T = eos->get_temperature(eps, rhob);
         fp1 = muB/T; 
       
       // fm1 = grid_pt->rhob;
         if(rk_flag == 0)
         {
           rhob = grid_pt->rhob;
           eps = grid_pt->epsilon;
         }
         else
         {
           rhob = grid_pt->rhob_t;
           eps = grid_pt->epsilon_t;
         }
         muB = eos->get_mu(eps, rhob);
         T = eos->get_temperature(eps, rhob);
         fm1 = muB/T; 
     }
     else
     {
      // fp1 = grid_pt->nbr_p_1[n]->rhob;
        if(rk_flag == 0)
        {
          rhob = grid_pt->nbr_p_1[n]->rhob;
          eps = grid_pt->nbr_p_1[n]->epsilon;
        }
        else
        {
          rhob = grid_pt->nbr_p_1[n]->rhob_t;
          eps = grid_pt->nbr_p_1[n]->epsilon_t;
        }
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        fp1 = muB/T; 
      
      // fm1 = Lneighbor->rhob;
        if(rk_flag == 0)
        {
          rhob = Lneighbor->rhob;
          eps = Lneighbor->epsilon;
        }
        else
        {
          rhob = Lneighbor->rhob_t;
          eps = Lneighbor->epsilon_t;
        }
        muB = eos->get_mu(eps, rhob);
        T = eos->get_temperature(eps, rhob);
        fm1 = muB/T; 
     } 
   }
   else
   {
    // fp1 = grid_pt->nbr_p_1[n]->rhob;
      if(rk_flag == 0)
      {
        rhob = grid_pt->nbr_p_1[n]->rhob;
        eps = grid_pt->nbr_p_1[n]->epsilon;
      }
      else
      {
        rhob = grid_pt->nbr_p_1[n]->rhob_t;
        eps = grid_pt->nbr_p_1[n]->epsilon_t;
      }
      muB = eos->get_mu(eps, rhob);
      T = eos->get_temperature(eps, rhob);
      fp1 = muB/T; 
    
    // fm1 = grid_pt->nbr_m_1[n]->rhob;
      if(rk_flag == 0)
      {
        rhob = grid_pt->nbr_m_1[n]->rhob;
        eps = grid_pt->nbr_m_1[n]->epsilon;
      }
      else
      {
        rhob = grid_pt->nbr_m_1[n]->rhob_t;
        eps = grid_pt->nbr_m_1[n]->epsilon_t;
      }
      muB = eos->get_mu(eps, rhob);
      T = eos->get_temperature(eps, rhob);
      fm1 = muB/T; 
   }
   g = minmod->minmod_dx(fp1, f, fm1);
   g /= delta[n]*taufactor;
   grid_pt->dUsup[rk_flag][m][n] = g;

 return 1;
}/* MakeDSpatial */


int U_derivative::MakeDTau(double tau, InitData *DATA, Grid *grid_pt, int rk_flag)
{
 int m;
 double f;
 // Sangyong Nov 18 2014: added these doubles 
 double tildemu, tildemu_prev, rhob, eps, muB, T;

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

// Sangyong Nov 18 2014
// Here we make the time derivative of (muB/T)
if(rk_flag == 0)
{
 m = 4;  
   // first order is more stable 
   // backward derivative
   // current values
   // f = (grid_pt->rhob);
   rhob = grid_pt->rhob;
   eps = grid_pt->epsilon;
   muB = eos->get_mu(eps, rhob);
   T = eos->get_temperature(eps, rhob);
   tildemu = muB/T;

   // f -= (grid_pt->rhob_prev);
   rhob = grid_pt->prev_rhob;
   eps = grid_pt->prev_epsilon;
   muB = eos->get_mu(eps, rhob);
   T = eos->get_temperature(eps, rhob);
   tildemu_prev = muB/T;
   
   f = (tildemu - tildemu_prev)/(DATA->delta_tau);
   
   grid_pt->dUsup[rk_flag][m][0] = -f; /* g00 = -1 */
}/* rk_flag = 0 */
else if(rk_flag > 0)
{
 m = 4;  
   // first order 
   // forward derivative
//   f = (grid_pt->rhob_t); // this is from the prev full RK step 
//   f -= (grid_pt->rhob_prev);
//   f /= (DATA->delta_tau);
   
   rhob = grid_pt->rhob_t;
   eps = grid_pt->epsilon_t;
   muB = eos->get_mu(eps, rhob);
   T = eos->get_temperature(eps, rhob);
   tildemu = muB/T;

   // f -= (grid_pt->rhob_prev);
   rhob = grid_pt->rhob;
   eps = grid_pt->epsilon;
   muB = eos->get_mu(eps, rhob);
   T = eos->get_temperature(eps, rhob);
   tildemu_prev = muB/T;
   
   f = (tildemu - tildemu_prev)/(DATA->delta_tau);
   
   grid_pt->dUsup[rk_flag][m][0] = -f; /* g00 = -1 */
}
// Ends Sangyong's addition Nov 18 2014
 
 return 1;
}/* MakeDTau */


