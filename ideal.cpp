#include "util.h"
#include "data.h"
#include "grid.h"
#include "reconst.h"
#include "ideal.h"
#include "minmod.h"


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


void MakeDeltaQI(double tau, Grid *grid_pt, double *qi, double *rhs, 
InitData *DATA, int rk_flag) 
{
 static double delta[4], sumf, dwmn, tempf;
 static double tau_fac[4];
 static int alpha, i, rk_order_m1, nmax[4], flag;
 static double **DFmmp;
 static NbrQs NbrCells;
 static BdryCells HalfwayCells;
 static int ind=0;

 ind++;
 if(ind==1)
  {
   InitTempGrids(&HalfwayCells, DATA->rk_order); 
   InitNbrQs(&NbrCells);
   DFmmp = mtx_malloc(5,4);
  }

 delta[1] = DATA->delta_x;
 delta[2] = DATA->delta_y;
 delta[3] = DATA->delta_eta;

 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta;

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

/* implement Kurganov-Tadmor scheme for the IDEAL fluid part */
 GetQIs(tau, grid_pt, qi, &NbrCells, rk_flag, DATA);

 flag = 
 MakeQIHalfs(qi, &NbrCells, &HalfwayCells, grid_pt, DATA);

 flag = 
 ConstHalfwayCells(tau, &HalfwayCells, qi, grid_pt, DATA, rk_flag);

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
     sumf -= tempf;
    }
   else if(alpha==3)
    {
     sumf -= grid_pt->TJb[rk_flag][3][0];
     sumf -= tempf;
    }
   
   rhs[alpha] = sumf*(DATA->delta_tau);
  }/* alpha */
   
 return;
}/* MakeDeltaQI */


/* %%%%%%%%%%%%%%%%% Kurganov-Tadmor %%%%%%%%%%%%%%%%% */


void GetQIs
(double tau, Grid *grid_pt,
 double *qi, NbrQs *NbrCells, int rk_flag, InitData *DATA)
{
 int alpha, i;
 double tempg, T00, K00, tempf;
 int nmax[4];

 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta;

 tempg = tau;

 for(alpha=0; alpha<5; alpha++)
  {
    /* qs from the neighbors */
    /* implement outflow boundary condition - simply set the two outside
     * values the same as at the boundary */

     for(i=1;i<=3;i++)
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
  }/* alpha */

 return; 
}/* GetQIs */     


int MakeQIHalfs
(double *qi, NbrQs *NbrCells, BdryCells *HalfwayCells, 
 Grid *grid_pt, InitData *DATA)
{
 int alpha, direc, nmax[4], flag;
 double fphL, fphR, fmhL, fmhR;
 double gphL, gphR, gmhL, gmhR;
 double x, y, eta, tempf;

  nmax[1] = DATA->nx;
  nmax[2] = DATA->ny;
  nmax[3] = DATA->neta;
  
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


int ConstHalfwayCells
 (double tau, BdryCells *HalfwayCells, double *qi, Grid *grid_pt,
  InitData *DATA, int rk_flag)
{
 int direc, flag;
 double epsilon_init, rhob_init;

 epsilon_init = grid_pt->epsilon;
 rhob_init = grid_pt->rhob;

 for(direc=1; direc<=3; direc++)
 {
 /* for each direction, reconstruct half-way cells */

  flag=
  Reconst(&(HalfwayCells->grid_p_h_L[direc]), 
            direc, tau, HalfwayCells->qiphL, grid_pt,
	    epsilon_init, rhob_init, DATA, rk_flag); 
  if(flag==0) {
  ReconstError("grid_p_h_L", direc, rk_flag, qi, HalfwayCells->qiphL, grid_pt);
    return 0;
   }/* if Reconst returns error */

  flag=
  Reconst(&(HalfwayCells->grid_p_h_R[direc]), 
            direc, tau, HalfwayCells->qiphR, grid_pt,
	    epsilon_init, rhob_init, DATA, rk_flag); 
  if(flag==0) {
  ReconstError("grid_p_h_R", direc, rk_flag, qi, HalfwayCells->qiphR, grid_pt);
    return 0;
   }
  
  flag=
  Reconst(&(HalfwayCells->grid_m_h_L[direc]), 
            direc, tau, HalfwayCells->qimhL, grid_pt,
	    epsilon_init, rhob_init, DATA, rk_flag); 
  if(flag==0) {
  ReconstError("grid_m_h_L", direc, rk_flag, qi, HalfwayCells->qimhL, grid_pt);
    return 0;
   }
  
  flag=
  Reconst(&(HalfwayCells->grid_m_h_R[direc]), 
            direc, tau, HalfwayCells->qimhR, grid_pt,
	    epsilon_init, rhob_init, DATA, rk_flag); 
  if(flag==0) {
  ReconstError("grid_m_h_R", direc, rk_flag, qi, HalfwayCells->qimhR, grid_pt);
    return 0;
   }
 }/* direc */
 return 1;/* upon successful execution */
}/* ConstHalfwayCells */


void MakeKTCurrents
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



void MakeMaxSpeedAs
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




double MaxSpeed (double tau, int direc, Grid *grid_p, int rk_flag)
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
 
 pe = p_e_func(eps, rhob);
 rpr = rhob*p_rho_func(eps, rhob);

 den = -((1.0 + pe)*rpr*(-1.0 + Power(ut,2))) 
       + h*(pe + Power(ut,2) - pe*Power(ut,2));
 
 num = 
 sqrt(-((h*pe + rpr + pe*rpr)*(h*(pe*(-1.0 + ut2mux2) - ut2mux2) 
          + (1.0 + pe)*rpr*(-1.0 + ut2mux2)))) 
	  - h*(-1.0 + pe)*ut*ux - rpr*ut*ux - pe*rpr*ut*ux;
 
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
     fprintf(stderr, "SpeedMax = %e\n is smaller than v = %e.\n", f, ux/ut);
     fprintf(stderr, "Can't happen.\n");
     exit(0);
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



void InitNbrQs(NbrQs *NbrCells)
{
 (NbrCells->qip1) = mtx_malloc(5,4);
 (NbrCells->qip2) = mtx_malloc(5,4);
 (NbrCells->qim1) = mtx_malloc(5,4);
 (NbrCells->qim2) = mtx_malloc(5,4);
}/* InitNbrQs */


void InitTempGrids(BdryCells *HalfwayCells, int rk_order)
{
 int direc;

 (HalfwayCells->grid_p_h_L) = grid_v_malloc(4);
 (HalfwayCells->grid_p_h_R) = grid_v_malloc(4);
 (HalfwayCells->grid_m_h_L) = grid_v_malloc(4);
 (HalfwayCells->grid_m_h_R) = grid_v_malloc(4);

 HalfwayCells->qiphL = mtx_malloc(5,4);
 HalfwayCells->qiphR = mtx_malloc(5,4);
 HalfwayCells->qimhL = mtx_malloc(5,4);
 HalfwayCells->qimhR = mtx_malloc(5,4);
 
 for(direc=0; direc<4; direc++)
  {
    (HalfwayCells->grid_p_h_L)[direc].TJb = cube_malloc(rk_order,5,4);
    (HalfwayCells->grid_p_h_R)[direc].TJb = cube_malloc(rk_order,5,4);
    (HalfwayCells->grid_m_h_L)[direc].TJb = cube_malloc(rk_order,5,4);
    (HalfwayCells->grid_m_h_R)[direc].TJb = cube_malloc(rk_order,5,4);
    
    (HalfwayCells->grid_p_h_L)[direc].u = mtx_malloc(rk_order,4);
    (HalfwayCells->grid_p_h_R)[direc].u = mtx_malloc(rk_order,4);
    (HalfwayCells->grid_m_h_L)[direc].u = mtx_malloc(rk_order,4);
    (HalfwayCells->grid_m_h_R)[direc].u = mtx_malloc(rk_order,4);
   }
 return;
}/* InitTempGrids */


