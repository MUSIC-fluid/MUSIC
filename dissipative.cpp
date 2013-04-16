#include "util.h"
#include "grid.h"
#include "data.h"
#include "eos.h"
#include "dissipative.h"

using namespace std;

Diss::Diss(EOS *eosIn)
{
  eos = new EOS;
  eos = eosIn;
  minmod = new Minmod;
}

// destructor
Diss::~Diss()
{
  delete eos;
  delete minmod;
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Dissipative parts */

/* this is the only one that is being subtracted in the rhs */
double Diss::MakeWSource(double tau, int alpha, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			 Grid *Lneighbor2, Grid *Rneighbor2, InitData *DATA,int rk_flag, int size, int rank)
{
 double sgp1, sgm1, bgp1, bgm1;
 double delta[4], taufactor, tf;
 double shear_on, bulk_on, sf, bf, sg, bg;
 int i, nmax[4], mu, nu;

 if(DATA->turn_on_shear) shear_on = 1.0;
 else shear_on = 0.0;
 
 if(DATA->turn_on_bulk) bulk_on = 1.0;
 else bulk_on = 0.0;
 
 /* calculate d_m (tau W^{m,alpha}) + (geom source terms) */ 

 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta-1;
 
 delta[1] = DATA->delta_x;
 delta[2] = DATA->delta_y;
 delta[3] = DATA->delta_eta;

/* partial_tau W^tau alpha */

/* this is partial_tau evaluated at tau */
/* this is the first step. so rk_flag = 0 */
 
 sf = 0.0;
 bf = 0.0;

if(rk_flag==0)
{
/* first order is more stable */
 tf = (grid_pt->Wmunu[rk_flag][0][alpha]);
 tf -= (grid_pt->prevWmunu[0][0][alpha]);
 tf /= (DATA->delta_tau);
}
else if(rk_flag > 0)
{
/* first order since we don't know next values yet */
 tf = (grid_pt->Wmunu[rk_flag][0][alpha]);
 tf -= (grid_pt->Wmunu[0][0][alpha]);
 tf /= (DATA->delta_tau);
}
 
sf += tf;

/* bulk pressure term */

if(rk_flag==0)
{
/* second order since we know prev and pprev */
 tf  = (grid_pt->Pimunu[rk_flag][0][alpha]);
 tf -= (grid_pt->prevPimunu[rk_flag][0][alpha]);
 tf /= (DATA->delta_tau);
}
else if(rk_flag > 0)
{
/* first order since we don't know next values yet */
 tf  = (grid_pt->Pimunu[rk_flag][0][alpha]);
 tf -= (grid_pt->Pimunu[0][0][alpha]);
 tf /= (DATA->delta_tau);
}
 
 bf += tf;

 for(i=1; i<=2; i++) // x and y
  {
    taufactor = 1.0;
    
    sg = grid_pt->Wmunu[rk_flag][i][alpha];
    bg = grid_pt->Pimunu[rk_flag][i][alpha];
   
   if(grid_pt->position[i] == nmax[i])
    {
     sgp1 = sg;
     bgp1 = bg;
     sgm1 = grid_pt->nbr_m_1[i]->Wmunu[rk_flag][i][alpha];
     bgm1 = grid_pt->nbr_m_1[i]->Pimunu[rk_flag][i][alpha];
    } 
   else if(grid_pt->position[i] == 0)
    {
     sgp1 = grid_pt->nbr_p_1[i]->Wmunu[rk_flag][i][alpha];
     bgp1 = grid_pt->nbr_p_1[i]->Pimunu[rk_flag][i][alpha];
     sgm1 = sg;
     bgm1 = bg;
    }
   else
    {
     sgp1 = grid_pt->nbr_p_1[i]->Wmunu[rk_flag][i][alpha];
     bgp1 = grid_pt->nbr_p_1[i]->Pimunu[rk_flag][i][alpha];
     sgm1 = grid_pt->nbr_m_1[i]->Wmunu[rk_flag][i][alpha];
     bgm1 = grid_pt->nbr_m_1[i]->Pimunu[rk_flag][i][alpha];
    }
   
   sf += minmod->minmod_dx(sgp1, sg, sgm1, DATA)/delta[i]/taufactor; 
   bf += minmod->minmod_dx(bgp1, bg, bgm1, DATA)/delta[i]/taufactor; 
  }/* i */

 i=3;
 taufactor = tau;
 
 sg = grid_pt->Wmunu[rk_flag][i][alpha];
 bg = grid_pt->Pimunu[rk_flag][i][alpha];
 
 if(grid_pt->position[i] == nmax[i])
   {
     if(rank == size-1) // for the right most rank do boundary condition on the right
       {
	 sgp1 = sg;
	 bgp1 = bg;
	 sgm1 = grid_pt->nbr_m_1[i]->Wmunu[rk_flag][i][alpha];
	 bgm1 = grid_pt->nbr_m_1[i]->Pimunu[rk_flag][i][alpha];
       }
     else
       {
	 sgp1 = Rneighbor->Wmunu[rk_flag][i][alpha];
	 bgp1 = Rneighbor->Pimunu[rk_flag][i][alpha];
	 sgm1 = grid_pt->nbr_m_1[i]->Wmunu[rk_flag][i][alpha];
	 bgm1 = grid_pt->nbr_m_1[i]->Pimunu[rk_flag][i][alpha];
       }
   } 
 else if(grid_pt->position[i] == 0)
   {
     if(rank == 0) // for the left most rank do boundary condition on the left
       {
	 sgp1 = grid_pt->nbr_p_1[i]->Wmunu[rk_flag][i][alpha];
	 bgp1 = grid_pt->nbr_p_1[i]->Pimunu[rk_flag][i][alpha];
	 sgm1 = sg;
	 bgm1 = bg;
       }
     else
       {
	 sgp1 = grid_pt->nbr_p_1[i]->Wmunu[rk_flag][i][alpha];
	 bgp1 = grid_pt->nbr_p_1[i]->Pimunu[rk_flag][i][alpha];
	 sgm1 = Lneighbor->Wmunu[rk_flag][i][alpha];
	 bgm1 = Lneighbor->Pimunu[rk_flag][i][alpha];
       }
   }
 else
   {
     sgp1 = grid_pt->nbr_p_1[i]->Wmunu[rk_flag][i][alpha];
     bgp1 = grid_pt->nbr_p_1[i]->Pimunu[rk_flag][i][alpha];
     sgm1 = grid_pt->nbr_m_1[i]->Wmunu[rk_flag][i][alpha];
     bgm1 = grid_pt->nbr_m_1[i]->Pimunu[rk_flag][i][alpha];
   }
 
 sf += minmod->minmod_dx(sgp1, sg, sgm1, DATA)/delta[i]/taufactor; 
 bf += minmod->minmod_dx(bgp1, bg, bgm1, DATA)/delta[i]/taufactor; 

 /* partial_m (tau W^mn) = W^0n + tau partial_m W^mn */

 /* tau partial_m W^mn */
 sf *= tau;
 bf *= tau;

/* add W^0n */
 sf += grid_pt->Wmunu[rk_flag][0][alpha];
 bf += grid_pt->Pimunu[rk_flag][0][alpha];

/* sources due to coordinate transform this is added to partial_m W^mn */

 if(alpha == 0)
  {
   sf += grid_pt->Wmunu[rk_flag][3][3];
   bf += grid_pt->Pimunu[rk_flag][3][3];
  }
 
 if(alpha == 3)
  {
   sf += grid_pt->Wmunu[rk_flag][0][3];
   bf += grid_pt->Pimunu[rk_flag][0][3];
  }

 tf = (sf*shear_on + bf*bulk_on);

/*
 if(tf == 0.0 && tau > DATA->tau0+DATA->delta_tau) 
  {
   fprintf(stderr, "MakeWSource: tf = %e\n", tf);
   fprintf(stderr, "MakeWSource: sf = %e\n", sf);
   fprintf(stderr, "MakeWSource: bf = %e\n", bf);
   fprintf(stderr, "MakeWSource: shear_on = %e\n", shear_on);
   fprintf(stderr, "MakeWSource: bulk_on =  %e\n", bulk_on);
   fprintf(stderr, "MakeWSource: tau =  %e\n", tau);
   fprintf(stderr, "MakeWSource: rk_flag =  %d\n", rk_flag);
  }
*/

 if (isnan(tf)) cout << "sf=" << sf << " bf=" << bf << " sg=" << sg << " bg=" << bg << " Wmunu[" << rk_flag << "]=" << grid_pt->Wmunu[rk_flag][0][alpha]
		     << " Pimunu[" << rk_flag << "]=" << grid_pt->Pimunu[rk_flag][0][alpha] 
		     << " prevWmunu=" << grid_pt->prevWmunu[0][0][alpha] << endl;

 // cout << "tf=" << tf << endl;


 return tf;
}/* MakeWSource */

/* MakeWSource is for Tmunu */


double Diss::Make_uWSource(double tau, Grid *grid_pt, int mu, int nu, InitData *DATA, int rk_flag)
{
 double tempf, tempg, temps, tau_pi;
 double SW, s_den, shear, shear_to_s, T, epsilon, rhob, Ttr;
 int ic;
 Ttr = 0.18;
 
 epsilon = grid_pt->epsilon;
 rhob = grid_pt->rhob; 
	
 if(DATA->turn_on_shear == 0) return 0.0;

 if(DATA->T_dependent_shear_to_s == 1)
   {
     T=eos->get_temperature(epsilon,rhob)*hbarc;
     
     if(T < 0.18)
       {
	 shear_to_s=0.681-0.0594*T/Ttr-0.544*(T/Ttr)*(T/Ttr);
       }
     else
       {
	 shear_to_s=-0.289+0.288*T/Ttr+0.0818*(T/Ttr)*(T/Ttr);
       }
   }
 else
   {
     shear_to_s = DATA->shear_to_s;
   }
 
 //cout << "eta/s=" << shear_to_s << endl;
 
 s_den = eos->get_entropy(epsilon, rhob);
 shear = (shear_to_s)*s_den;
 tau_pi = 3.0*shear/(grid_pt->epsilon + grid_pt->p);

 //set viscosity to zero if energy density is very low (~1/10 of freeze out energy density)
 if (grid_pt->epsilon < 0.01/hbarc) shear = 0.;


 tau_pi = maxi(tau_pi, DATA->tau_pi);
 if(!finite(tau_pi)) tau_pi = DATA->tau_pi;


/*
 fprintf(stderr, "tau_pi = %e\n", tau_pi);
 fprintf(stderr, "h = %e\n", grid_pt->epsilon + grid_pt->p);
*/

/* This source has many terms */
/* everyting in the 1/(tau_pi) piece is here */
/* third step in the split-operator time evol 
   use Wmunu[rk_flag] and u[rk_flag] with rk_flag = 0 */
    
    tempf = 0.0;
    tempg = (1.0 + (4.0/3.0)*(tau_pi)*(grid_pt->theta_u[rk_flag]) );
    tempg *= (grid_pt->Wmunu[rk_flag][mu][nu]);
    tempf -= tempg;
 
   /* remember dUsup[m][n] = u^{m,n} = partial^n u^m */

/* this is a crude approximation. fix this. */
/* u+d+s+g e+p \approx 21 T^4 */
/* T = ((e+p)/21)^(1/4) */
/* s = (e+p)/T = (e+p)^(3/4)/(21)^(1/4) \approx (e+p)^(3/4)/2.14 */

/* inside the parenthesis */

    temps = 0.0;

    temps -=
    (grid_pt->dUsup[rk_flag][nu][mu] + grid_pt->dUsup[rk_flag][mu][nu]);

    tempg = (2.0*grid_pt->u[rk_flag][0]);
    tempg *= (DATA->gmunu[mu][3])*(DATA->gmunu[nu][3])/tau;
    temps -= tempg;
    
    tempg = (2.0/3.0)*(DATA->gmunu[mu][nu] 
                        + (grid_pt->u[rk_flag][mu])*(grid_pt->u[rk_flag][nu]));
    tempg *= (grid_pt->theta_u[rk_flag]);
    temps += tempg;
    
    temps +=
    (grid_pt->u[rk_flag][3]*(DATA->gmunu[mu][3])*(DATA->gmunu[nu][0]))/tau;
    temps +=
    (grid_pt->u[rk_flag][3]*(DATA->gmunu[nu][3])*(DATA->gmunu[mu][0]))/tau;


/* start: This is what's causing the most problems */
/* because it has the time derivative */

    temps -= ( grid_pt->u[rk_flag][mu]*grid_pt->a[rk_flag][nu] 
              +grid_pt->u[rk_flag][nu]*grid_pt->a[rk_flag][mu] );


/* end: This is what's causing the most problems */
    
    tempg = (grid_pt->u[rk_flag][3])*(grid_pt->u[rk_flag][0]);
    tempg *= ( (grid_pt->u[rk_flag][mu])*(DATA->gmunu[nu][3]) 
              +(grid_pt->u[rk_flag][nu])*(DATA->gmunu[mu][3]) );
    tempg /= tau;
    temps -= tempg;
    
    tempg = (grid_pt->u[rk_flag][3])*(grid_pt->u[rk_flag][3]);
    tempg *= ( (grid_pt->u[rk_flag][mu])*(DATA->gmunu[nu][0]) 
              +(grid_pt->u[rk_flag][nu])*(DATA->gmunu[mu][0]) );
    tempg /= tau;
    temps += tempg;


    tempf += temps*shear;
    SW = tempf/(tau_pi);
    
    return SW;
}/* Make_uWSource */


int Diss::Make_uWRHS(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			 Grid *Lneighbor2, Grid *Rneighbor2, double **w_rhs, InitData *DATA, int rk_flag, int size, int rank)
{
 int mu, nu, direc, nmax[4], ic;
 double f, fp1, fm1, fp2, fm2, ux, delta[4], sumf;
 double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
 double uWphR, uWphL, uWmhR, uWmhL, WphR, WphL, WmhR, WmhL;
 double HWph, HWmh, taufactor, HW, SW, ic_fac;
/*  HW[4][4][4], SW[4][4][4] */
 double tempf, tempg, shear, sum, shear_on;
 double s_den;

/* Kurganov-Tadmor for Wmunu */
/* implement 
  partial_tau (utau Wmn) + (1/tau)partial_eta (ueta Wmn) 
  + partial_x (ux Wmn) + partial_y (uy Wmn) + utau Wmn/tau = SW 
  or the right hand side of,
  partial_tau (utau Wmn) = 
  - (1/tau)partial_eta (ueta Wmn) - partial_x (ux Wmn) - partial_y (uy Wmn) 
  - utau Wmn/tau + SW 
  */

/* the local velocity is just u_x/u_tau, u_y/u_tau, u_eta/tau/u_tau */
/* KT flux is given by 
   H_{j+1/2} = (fRph + fLph)/2 - ax(uRph - uLph) 
   Here fRph = ux WmnRph and ax uRph = |ux/utau|_max utau Wmn */

/* This is the second step in the operator splitting. it uses
   rk_flag+1 as initial condition */

 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta-1;

 delta[1] = DATA->delta_x;
 delta[2] = DATA->delta_y;
 delta[3] = DATA->delta_eta;

 if(DATA->turn_on_shear) shear_on = 1.0;
 else shear_on = 0.0;

 for(mu=0; mu<4; mu++)
  {
   for(nu=0; nu<4; nu++)
    {
     sum = 0.0;
     for(direc=1; direc<=3; direc++)
      {
       if(direc==3) taufactor = tau;
       else taufactor = 1.0;

/* Get_uWmns */
       Get_uWmns(tau, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, mu, nu, direc, &g, &f, &gp1, &fp1, &gp2, &fp2, 
		 &gm1, &fm1, &gm2, &fm2, DATA, rk_flag, size, rank);
 
/*  MakeuWmnHalfs */
/* uWmn */
       uWphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f, DATA); 
       uWphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1, DATA);
       uWmhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1, DATA);
       uWmhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2, DATA);

/* just Wmn */
       WphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g, DATA); 
       WphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1, DATA);
       WmhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1, DATA);
       WmhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2, DATA);

       a = fabs(grid_pt->u[rk_flag][direc]);
       a /= grid_pt->u[rk_flag][0];


       if(direc<3) // x,y direction
	 {
	   if(grid_pt->position[direc] == 0)
	     {
	       am1 = a;
	     }
	   else
	     {
	       am1 = fabs(grid_pt->nbr_m_1[direc]->u[rk_flag][direc]);
	       am1 /= grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	     }
	   
	   if(grid_pt->position[direc] == nmax[direc])
	     {
	       ap1 = a;
	     }
	   else
	     {
	       ap1 = fabs(grid_pt->nbr_p_1[direc]->u[rk_flag][direc]);
	       ap1 /= grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	     }
       	 }
       else if (direc==3)
	 {
	   if(grid_pt->position[direc] == 0)
	     {
	       if (rank==0)
		 {
		   am1 = a;
		 }
	       else
		 {
		   am1 = fabs(Lneighbor->u[rk_flag][direc]);
		   am1 /= Lneighbor->u[rk_flag][0];
		 }
	     }
	   else
	     {
	       am1 = fabs(grid_pt->nbr_m_1[direc]->u[rk_flag][direc]);
	       am1 /= grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	     }
	   
	   if(grid_pt->position[direc] == nmax[direc])
	     {
	       if (rank==size-1)
		 {
		   ap1 = a;
		 }
	       else
		 {
		   ap1 = fabs(Rneighbor->u[rk_flag][direc]);
		   ap1 /= Rneighbor->u[rk_flag][0];
		 }		 
	     }
	   else
	     {
	       ap1 = fabs(grid_pt->nbr_p_1[direc]->u[rk_flag][direc]);
	       ap1 /= grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	     }
	 }

       ax = maxi(a, ap1);
       HWph = (uWphR + uWphL) - ax*(WphR - WphL);
       HWph *= 0.5;

       ax = maxi(a, am1); 
       HWmh = (uWmhR + uWmhL) - ax*(WmhR - WmhL);
       HWmh *= 0.5;

       HW = (HWph - HWmh)/delta[direc]/taufactor;

/* make partial_i (u^i Wmn) */
       sum += -HW;

      }/* direction */
    
/* add a source term -u^tau Wmn/tau due to the coordinate change to tau-eta */
    
    sum -= (grid_pt->u[rk_flag][0])*(grid_pt->Wmunu[rk_flag][mu][nu])/tau;
    sum += (grid_pt->theta_u[rk_flag])*(grid_pt->Wmunu[rk_flag][mu][nu]);
   
   /* this is from udW = d(uW) - Wdu = RHS */
   /* or d(uW) = udW + Wdu */
   /* this term is being added to the rhs so that -4/3 + 1 = -1/3 */

    
/* other source terms due to the coordinate change to tau-eta */
    tempf = 0.0;

    tempf += -(DATA->gmunu[3][mu])*(grid_pt->Wmunu[rk_flag][0][nu]);
    tempf += -(DATA->gmunu[3][nu])*(grid_pt->Wmunu[rk_flag][0][mu]);
    tempf +=  (DATA->gmunu[0][mu])*(grid_pt->Wmunu[rk_flag][3][nu]);
    tempf +=  (DATA->gmunu[0][nu])*(grid_pt->Wmunu[rk_flag][3][mu]);
  
    tempf +=
    (grid_pt->Wmunu[rk_flag][3][nu])*(grid_pt->u[rk_flag][mu])*(grid_pt->u[rk_flag][0]);
    tempf += 
    (grid_pt->Wmunu[rk_flag][3][mu])*(grid_pt->u[rk_flag][nu])*(grid_pt->u[rk_flag][0]);
    tempf -=
    (grid_pt->Wmunu[rk_flag][0][nu])*(grid_pt->u[rk_flag][mu])*(grid_pt->u[rk_flag][3]);
    tempf -=
    (grid_pt->Wmunu[rk_flag][0][mu])*(grid_pt->u[rk_flag][nu])*(grid_pt->u[rk_flag][3]);
    
    tempf *= (grid_pt->u[rk_flag][3]/tau);
    
    for(ic=0; ic<4; ic++)
    {
     ic_fac = (ic==0 ? -1.0 : 1.0);
     tempf +=
     (grid_pt->Wmunu[rk_flag][ic][nu])*(grid_pt->u[rk_flag][mu])*(grid_pt->a[rk_flag][ic])*ic_fac;
     tempf +=
     (grid_pt->Wmunu[rk_flag][ic][mu])*(grid_pt->u[rk_flag][nu])*(grid_pt->a[rk_flag][ic])*ic_fac;
    }
/* SYM TEST */
    tempg = 0.0;

    tempg += -(DATA->gmunu[3][nu])*(grid_pt->Wmunu[rk_flag][0][mu]);
    tempg += -(DATA->gmunu[3][mu])*(grid_pt->Wmunu[rk_flag][0][nu]);
    tempg +=  (DATA->gmunu[0][nu])*(grid_pt->Wmunu[rk_flag][3][mu]);
    tempg +=  (DATA->gmunu[0][mu])*(grid_pt->Wmunu[rk_flag][3][nu]);
  
    tempg +=
    (grid_pt->Wmunu[rk_flag][3][mu])*(grid_pt->u[rk_flag][nu])*(grid_pt->u[rk_flag][0]);
    tempg += 
    (grid_pt->Wmunu[rk_flag][3][nu])*(grid_pt->u[rk_flag][mu])*(grid_pt->u[rk_flag][0]);
    tempg -=
    (grid_pt->Wmunu[rk_flag][0][mu])*(grid_pt->u[rk_flag][nu])*(grid_pt->u[rk_flag][3]);
    tempg -=
    (grid_pt->Wmunu[rk_flag][0][nu])*(grid_pt->u[rk_flag][mu])*(grid_pt->u[rk_flag][3]);
    
    tempg *= (grid_pt->u[rk_flag][3]/tau);
    
    for(ic=0; ic<4; ic++)
    {
     ic_fac = (ic==0 ? -1.0 : 1.0);
     tempg +=
     (grid_pt->Wmunu[rk_flag][ic][mu])*(grid_pt->u[rk_flag][nu])*(grid_pt->a[rk_flag][ic])*ic_fac;
     tempg +=
     (grid_pt->Wmunu[rk_flag][ic][nu])*(grid_pt->u[rk_flag][mu])*(grid_pt->a[rk_flag][ic])*ic_fac;
    }
/* SYM TEST */

    sum += (tempf+tempg)/2.0;
     
     w_rhs[mu][nu] = sum*(DATA->delta_tau)*shear_on;
   //   cout << "w_rhs=" << w_rhs << ", sum=" << sum << ", tempf=" << tempf << ", HW=" << HW << endl; 
//      cout << "a=" << a << ", am1=" << am1 << ", ap1=" << ap1 << ", ax=" << ax << endl;
//      cout << "uWphR=" << uWphR << ", uWphL=" << uWphL << ", uWmhR=" << uWmhR << ", uWmhL=" << uWmhL << endl;
//      cout << "fp1=" << fp1 << ", fm1=" << fm1 << ", f=" << f << endl;
    }/* nu */
  }/* mu */

 return 1; /* if successful */
}/* Make_uWRHS */



//needs to be parallelized:
void Diss::Get_uWmns(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		     Grid *Lneighbor2, Grid *Rneighbor2, int mu, int nu, int direc,
		     double *g, double *f, double *gp1, double *fp1,
		     double *gp2, double *fp2, double *gm1, double *fm1, double *gm2, double *fm2, 
		     InitData *DATA, int rk_flag, int size, int rank) 
{
 double tf, tg, tgp1, tfp1, tgp2, tfp2, tgm1, tfm1, tgm2, tfm2;
 int nmax[4];

/* Get_uWmns */
/* this is the last step of split-operator evolution */
/* should use Wmunu[rk_flag+1] and u[rk_flag] */
/* recall that in this file, rk_flag = 0 always. */

 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta-1;

       tg = grid_pt->Wmunu[rk_flag][mu][nu];
       tf = tg*grid_pt->u[rk_flag][direc];
       tg *=   grid_pt->u[rk_flag][0];
       
       if (direc<3) // x and y direction
	 {
	   if(grid_pt->position[direc] == nmax[direc])
	     {
	       tgp2 = tg; 
	       tfp2 = tf;
	       
	       tgp1 = tg;
	       tfp1 = tf;
	 
	       tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
	     }
	   else if(grid_pt->position[direc] == (nmax[direc]-1))
	     {
	       tgp2 = tg;
	       tfp2 = tf;
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
	     }
	   else if(grid_pt->position[direc] == 0)
	     {
	       tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = tg;
	       tfm1 = tf;
	       
	       tgm2 = tg;
	       tfm2 = tf;
	     }
	   else if(grid_pt->position[direc] == 1)
	     {
	       tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = tg;
	       tfm2 = tf;
	     }
	   else 
	     {
	       tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
	     }
	 }


       if (direc==3) // eta direction needs extra care in mpi:
	 {
	   if(grid_pt->position[direc] == nmax[direc])
	     {
	       if(rank == size-1) // for the right most rank do boundary condition on the right
		 {
		   tgp2 = tg; 
		   tfp2 = tf;
		   
		   tgp1 = tg;
		   tfp1 = tf;
		   
		   tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
		 }
	       else // for all other ranks use values from neighboring CPUs
		 {
		   tgp2 = Rneighbor2->Wmunu[rk_flag][mu][nu];
		   tfp2 = tgp2*Rneighbor2->u[rk_flag][direc];
		   tgp2 *=     Rneighbor2->u[rk_flag][0];
		   
		   tgp1 = Rneighbor->Wmunu[rk_flag][mu][nu];
		   tfp1 = tgp1*Rneighbor->u[rk_flag][direc];
		   tgp1 *=     Rneighbor->u[rk_flag][0];
		     
		   tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
	// 	   tgm1 = Lneighbor->Wmunu[rk_flag][mu][nu];
// 		   tfm1 = tgm1*Lneighbor->u[rk_flag][direc];
// 		   tgm1 *=     Lneighbor->u[rk_flag][0];
		   
// 		   tgm2 = Lneighbor2->Wmunu[rk_flag][mu][nu];
// 		   tfm2 = tgm2*Lneighbor2->u[rk_flag][direc];
// 		   tgm2 *=     Lneighbor2->u[rk_flag][0];
		   
		 }
	     }
	   else if(grid_pt->position[direc] == (nmax[direc]-1))
	     {
	       if(rank == size-1) // for the right most rank do boundary condition on the right
		 {
		   tgp2 = tg;
		   tfp2 = tf;
		   
		   tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
		 }
	       else // for all other ranks use values from neighboring CPUs
		 {
		   tgp2 = Rneighbor->Wmunu[rk_flag][mu][nu];
		   tfp2 = tgp2*Rneighbor->u[rk_flag][direc];
		   tgp2 *=     Rneighbor->u[rk_flag][0];
		   
		   tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
		 
		 }
	     }
	   else if(grid_pt->position[direc] == 0)
	     {
	      if(rank == 0) // for the left most rank do boundary condition on the left
		{
		  tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
		  tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
		  tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
		  
		  tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
		  tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		  tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		  
		  tgm1 = tg;
		  tfm1 = tf;
		  
		  tgm2 = tg;
		  tfm2 = tf;
		}
	      
	      else // for all other ranks use values from neighboring CPUs
		{
		  tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
		  tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
		  tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
		  
		  tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
		  tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		  tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		  
		  tgm1 = Lneighbor->Wmunu[rk_flag][mu][nu];
		  tfm1 = tgm1*Lneighbor->u[rk_flag][direc];
		  tgm1 *=     Lneighbor->u[rk_flag][0];
		  
		  tgm2 = Lneighbor2->Wmunu[rk_flag][mu][nu];
		  tfm2 = tgm2*Lneighbor2->u[rk_flag][direc];
		  tgm2 *=     Lneighbor2->u[rk_flag][0];
		}
	     }
	   else if(grid_pt->position[direc] == 1)
	     {
	       if(rank == 0) // for the left most rank do boundary condition on the left
		 {
		   tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = tg;
		   tfm2 = tf;
		 }
	       else // for all other ranks use values from neighboring CPUs
		 {
		   tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = Lneighbor->Wmunu[rk_flag][mu][nu];
		   tfm2 = tgm2*Lneighbor->u[rk_flag][direc];
		   tgm2 *=     Lneighbor->u[rk_flag][0];
		 }
	     }
	   else // normal case (not at a boundary)
	     {
	       tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
	     }
	 }

       //cout << "tg=" << tg << ", tf=" << tf << ", tfp1=" << tfp1 << ", tfm1=" << tfm1 << endl;

       *g = tg;
       *f = tf;
       *gp1 = tgp1;
       *fp1 = tfp1;
       *gm1 = tgm1;
       *fm1 = tfm1;
       *gp2 = tgp2;
       *fp2 = tfp2;
       *gm2 = tgm2;
       *fm2 = tfm2;
       
       return;
}/* Get_uWmns */



int Diss::Make_uPRHS(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		     Grid *Lneighbor2, Grid *Rneighbor2, double *p_rhs, InitData *DATA, int rk_flag, int size, int rank)
{
 int mu, nu, direc, nmax[4], ic;
 double f, fp1, fm1, fp2, fm2, ux, delta[4];
 double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
 double uPiphR, uPiphL, uPimhR, uPimhL, PiphR, PiphL, PimhR, PimhL;
 double HPiph, HPimh, taufactor, HPi, SPi;
 double tempf, tempg, sum;
 double s_den, bulk_on;

/* Kurganov-Tadmor for Pi */
/* implement 
  partial_tau (utau Pi) + (1/tau)partial_eta (ueta Pi) 
  + partial_x (ux Pi) + partial_y (uy Pi) + utau Pi/tau = SP 
  or the right hand side of
  partial_tau (utau Pi) = -
  (1/tau)partial_eta (ueta Pi) - partial_x (ux Pi) - partial_y (uy Pi)
  - utau Pi/tau + SP 
  */

/* the local velocity is just u_x/u_tau, u_y/u_tau, u_eta/tau/u_tau */
/* KT flux is given by 
   H_{j+1/2} = (fRph + fLph)/2 - ax(uRph - uLph) 
   Here fRph = ux PiRph and ax uRph = |ux/utau|_max utau Pin */

/* This is the second step in the operator splitting. it uses
   rk_flag+1 as initial condition */

 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta-1;

 delta[1] = DATA->delta_x;
 delta[2] = DATA->delta_y;
 delta[3] = DATA->delta_eta;

 if(DATA->turn_on_bulk) bulk_on = 1.0;
 else bulk_on = 0.0;

     sum = 0.0;
     for(direc=1; direc<=3; direc++)
      {
       if(direc==3) taufactor = tau;
       else taufactor = 1.0;

/* Get_uPis */
       Get_uPis(tau, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, direc, &g, &f, &gp1, &fp1, &gp2, &fp2, 
		&gm1, &fm1, &gm2, &fm2, DATA, rk_flag, size, rank);

/*  Make upi Halfs */
/* uPi */
       uPiphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f, DATA); 
       uPiphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1, DATA);
       uPimhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1, DATA);
       uPimhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2, DATA);

/* just Pi */
       PiphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g, DATA); 
       PiphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1, DATA);
       PimhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1, DATA);
       PimhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2, DATA);

/* MakePimnCurrents following Kurganov-Tadmor */
    
       a = fabs(grid_pt->u[rk_flag][direc]);
       a /= grid_pt->u[rk_flag][0];
  
       if (direc<3) // x,y direction
	 {
	   if(grid_pt->position[direc] == 0)
	     {
	       am1 = a;
	     }
	   else
	     {
	       am1 = fabs(grid_pt->nbr_m_1[direc]->u[rk_flag][direc]);
	       am1 /= grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	     }
	   
	   if(grid_pt->position[direc] == nmax[direc])
	     {
	       ap1 = a;
	     }
	   else
	     {
	       ap1 = fabs(grid_pt->nbr_p_1[direc]->u[rk_flag][direc]);
	       ap1 /= grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	     } 
	 }
       else if (direc==3) // eta direction
	 {
	   if(grid_pt->position[direc] == 0)
	     {
	       if (rank==0) // for left most rank use boundary conditions
		 {
		   am1 = a;
		 }
	       else // use value from neighboring CPU for all other ranks 
		 {
		   am1 = fabs(Lneighbor->u[rk_flag][direc]);
		   am1 /= Lneighbor->u[rk_flag][0];
		 }
	     }
	   else
	     {
	       am1 = fabs(grid_pt->nbr_m_1[direc]->u[rk_flag][direc]);
	       am1 /= grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	     }
	     
	   
	   if(grid_pt->position[direc] == nmax[direc])
	     {
	       if (rank == size-1) // for right most CPU use boundary conditions
		 {
		   ap1 = a;
		 }
	       else // use value from neighboring CPU for all other ranks 
		 {
		   ap1 = fabs(Rneighbor->u[rk_flag][direc]);
		   ap1 /= Rneighbor->u[rk_flag][0];
		 }
	     }
	   else // usual case (not at a boundary)
	     {
	       ap1 = fabs(grid_pt->nbr_p_1[direc]->u[rk_flag][direc]);
	       ap1 /= grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	     } 
	 }
       ax = maxi(a, ap1);
       HPiph = (uPiphR + uPiphL) - ax*(PiphR - PiphL);
       HPiph *= 0.5;
       
       ax = maxi(a, am1); 
       HPimh = (uPimhR + uPimhL) - ax*(PimhR - PimhL);
       HPimh *= 0.5;
      
       HPi = (HPiph - HPimh)/delta[direc]/taufactor;

/* make partial_i (u^i Pi) */
       sum += -HPi;

      }/* direction */
       
/* add a source term due to the coordinate change to tau-eta */
     sum -= (grid_pt->pi_b[rk_flag])*(grid_pt->u[rk_flag][0])/tau;

     *p_rhs = sum*(DATA->delta_tau)*bulk_on;

 return 1; /* if successful */
}/* Make_uPRHS */


void Diss::Get_uPis(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		    Grid *Lneighbor2, Grid *Rneighbor2, int direc,
		    double *g, double *f, 
		    double *gp1, double *fp1, double *gp2, double *fp2, 
		    double *gm1, double *fm1, double *gm2, double *fm2, 
		    InitData *DATA, int rk_flag, int size, int rank) 
{
 double tf, tg, tgp1, tfp1, tgp2, tfp2, tgm1, tfm1, tgm2, tfm2;
 int nmax[4];

/* Get_uPis */
 nmax[1] = DATA->nx;
 nmax[2] = DATA->ny;
 nmax[3] = DATA->neta-1;

       tg = grid_pt->pi_b[rk_flag];
       tf = tg*grid_pt->u[rk_flag][direc];
       tg *=   grid_pt->u[rk_flag][0];

       if (direc<3) // x,y direction
	 {
	   if(grid_pt->position[direc] == nmax[direc])
	     {
	       tgp2 = tg; 
	       tfp2 = tf;
	       
	       tgp1 = tg;
	       tfp1 = tf;
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
	     }
	   else if(grid_pt->position[direc] == (nmax[direc]-1))
	     {
	       tgp2 = tg;
	       tfp2 = tf;
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
	     }
	   else if(grid_pt->position[direc] == 0)
	     {
	       tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = tg;
	       tfm1 = tf;
	       
	       tgm2 = tg;
	       tfm2 = tf;
	     }
	   else if(grid_pt->position[direc] == 1)
	     {
	       tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = tg;
	       tfm2 = tf;
	     }
	   else 
	     {
	       tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
	     }
	 }
       else if (direc==3) // eta direction needs special care in mpi:
	 {
	   if(grid_pt->position[direc] == nmax[direc])
	     {
	       if(rank==size-1) // for right most rank use boundary condition
		 {
		   tgp2 = tg; 
		   tfp2 = tf;
		   
		   tgp1 = tg;
		   tfp1 = tf;
		   
		   tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->pi_b[rk_flag];
		   tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
		 }
	       else // for all other ranks use values from neighboring CPUs
		 {
		   tgp2 = Rneighbor2->pi_b[rk_flag];
		   tfp2 = tgp2*Rneighbor2->u[rk_flag][direc];
		   tgp2 *=     Rneighbor2->u[rk_flag][0];
		   
		   tgp1 = Rneighbor->pi_b[rk_flag];
		   tfp1 = tgp1*Rneighbor->u[rk_flag][direc];
		   tgp1 *=     Rneighbor->u[rk_flag][0];
	

		   tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->pi_b[rk_flag];
		   tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
		 }
	     }
	   else if(grid_pt->position[direc] == (nmax[direc]-1))
	     {
	       if(rank==size-1) // for right most rank use boundary condition
		 {
		   tgp2 = tg;
		   tfp2 = tf;
		   
		   tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
		   tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->pi_b[rk_flag];
		   tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
		 }
	       else // for all other ranks use values from neighboring CPUs
		 {
		   tgp2 = Rneighbor->pi_b[rk_flag];
		   tfp2 = tgp2*Rneighbor->u[rk_flag][direc];
		   tgp2 *=     Rneighbor->u[rk_flag][0];
		
		   tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
		   tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->pi_b[rk_flag];
		   tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
		 }
	     }
	   else if(grid_pt->position[direc] == 0)
	     {
	       if(rank==0) // for left most rank use boundary condition
		 {
		   tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->pi_b[rk_flag];
		   tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
		   tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgm1 = tg;
		   tfm1 = tf;
		   
		   tgm2 = tg;
		   tfm2 = tf;
		 }
	       else // for all other ranks use values from neighboring CPUs
		 {
		   tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->pi_b[rk_flag];
		   tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
		   tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		  
		   tgm1 = Lneighbor->pi_b[rk_flag];
		   tfm1 = tgm1*Lneighbor->u[rk_flag][direc];
		   tgm1 *=     Lneighbor->u[rk_flag][0];
		   
		   tgm2 = Lneighbor2->pi_b[rk_flag];
		   tfm2 = tgm2*Lneighbor2->u[rk_flag][direc];
		   tgm2 *=     Lneighbor2->u[rk_flag][0];
		 }       
	     }
	   else if(grid_pt->position[direc] == 1)
	     {
	       if(rank==0) // for left most rank use boundary condition
		 {
		   tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->pi_b[rk_flag];
		   tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
		   tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
		   
		   tgm2 = tg;
		   tfm2 = tf;
		 }
	       else // for all other ranks use values from neighboring CPUs
		 {
		   tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->pi_b[rk_flag];
		   tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
		   tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
		   tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
		   
		   tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
		   tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
		   tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];

		   tgm2 = Lneighbor->pi_b[rk_flag];
		   tfm2 = tgm2*Lneighbor->u[rk_flag][direc];
		   tgm2 *=     Lneighbor->u[rk_flag][0];
		 }
	     }
	   else // usual case (not at a boundary)
	     {
	       tgp2 = grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp2 = tgp2*grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_1[direc]->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm2 = tgm2*grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_1[direc]->nbr_m_1[direc]->u[rk_flag][0];
	     }
	 }

       *g = tg;
       *f = tf;
       *gp1 = tgp1;
       *fp1 = tfp1;
       *gm1 = tgm1;
       *fm1 = tfm1;
       *gp2 = tgp2;
       *fp2 = tfp2;
       *gm2 = tgm2;
       *fm2 = tfm2;
	   
}/* Get_uPis */


double Diss::Make_uPiSource
(double tau, Grid *grid_pt, InitData *DATA, int rk_flag)
{
 double tempf;
 double s_den, bulk;

/* this is the first step in the split-operator evolution */
/* use rk_flag = 0 */
/* this is a crude approximation. fix this. */


    if(DATA->turn_on_bulk == 0) return 0.0;

    s_den = eos->get_entropy(grid_pt->epsilon, grid_pt->rhob);
    
    bulk = (DATA->bulk_to_s)*s_den;

    tempf = 0.0;
    tempf -= (grid_pt->pi_b[rk_flag]);
    tempf -= bulk*(grid_pt->theta_u[rk_flag]);
    tempf -= (1.0/3.0)*(DATA->tau_b_pi)*(grid_pt->theta_u[rk_flag])*(grid_pt->pi_b[rk_flag]);
   
    return tempf/(DATA->tau_b_pi);
}/* Make_uPiSource */


