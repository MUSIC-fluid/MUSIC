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
 int i, nmax[4];

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
   	else {fprintf(stderr,"rk_flag out of range.\n");exit(0);}
 
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


 if (isnan(tf)) cout << "sf=" << sf << " bf=" << bf << " sg=" << sg << " bg=" << bg << " Wmunu[" << rk_flag << "]=" << grid_pt->Wmunu[rk_flag][0][alpha]
		     << " Pimunu[" << rk_flag << "]=" << grid_pt->Pimunu[rk_flag][0][alpha] 
		     << " prevWmunu=" << grid_pt->prevWmunu[0][0][alpha] << endl;


 return tf;
}/* MakeWSource */

/* MakeWSource is for Tmunu */


double Diss::Make_uWSource(double tau, Grid *grid_pt, int mu, int nu, InitData *DATA, int rk_flag)
{
 double tempf, tau_pi;
//  double tempg, temps;
 double SW, s_den, shear, shear_to_s, T, epsilon, rhob, Ttr;
 int a, b;
 double sigma[4][4], gamma, ueta;
 double NS_term;
 Ttr = 0.18;  /// phase transition temperature


/// Useful variables to define
gamma = grid_pt->u[rk_flag][0];
ueta  = grid_pt->u[rk_flag][3];
epsilon = grid_pt->epsilon;
rhob = grid_pt->rhob;


 if(DATA->turn_on_shear == 0) return 0.0;

 if(DATA->T_dependent_shear_to_s == 1)
   {
     T=eos->get_temperature(epsilon,rhob)*hbarc;

     if(T < Ttr)
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


/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///                 Defining transport coefficients                        ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
 s_den = eos->get_entropy(epsilon, rhob);
 shear = (shear_to_s)*s_den;
 tau_pi = 5.0*shear/(grid_pt->epsilon + grid_pt->p);
  

 //tau_pi = maxi(tau_pi, DATA->tau_pi);
 if(!finite(tau_pi)) {tau_pi = DATA->delta_tau; cout << "tau_pi was infinite ..." << endl;}

 /// transport coefficient for nonlinear terms -- shear only terms -- 4Mar2013
 double transport_coefficient, transport_coefficient2, transport_coefficient3, transport_coefficient4;
 /// transport coefficients of a massless gas of single component particles
 transport_coefficient  = 9./70.*tau_pi/shear*(4./5.) ;
 transport_coefficient2 = 4./3.*tau_pi;
 transport_coefficient3 = 10./7.*tau_pi;
 transport_coefficient4 = 2.*tau_pi;

 /// transport coefficient for nonlinear terms -- coupling to bulk viscous pressure -- 4Mar2013
 double transport_coefficient_b, transport_coefficient2_b;
 /// transport coefficients not yet known -- fixed to zero
 transport_coefficient_b  = 6./5.*tau_pi ;
 transport_coefficient2_b = 0.;


/* This source has many terms */
/* everyting in the 1/(tau_pi) piece is here */
/* third step in the split-operator time evol 
   use Wmunu[rk_flag] and u[rk_flag] with rk_flag = 0 */

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///            Wmunu + transport_coefficient2*Wmunu*theta                  ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///

/// full term is
    tempf = -(1.0 + transport_coefficient2*(grid_pt->theta_u[rk_flag]) )*(grid_pt->Wmunu[rk_flag][mu][nu]);
/// FOR GUBSER ANALYTIC
//     tempf = -( transport_coefficient2*(grid_pt->theta_u[rk_flag]) )*(grid_pt->Wmunu[rk_flag][mu][nu]);


/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///             Navier-Stokes Term -- -2.*shear*sigma^munu                 ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///

// / remember: dUsup[m][n] = partial^n u^m  ///
// / remember:  a[n]  =  u^m*partial_m u^n  ///

    for( a=0;a<=3;a++ ){
    for( b=0;b<=3;b++ ){

    sigma[a][b] = ( grid_pt->dUsup[rk_flag][a][b] + grid_pt->dUsup[rk_flag][b][a] )/2.

                - ( DATA->gmunu[a][b] + (grid_pt->u[rk_flag][a])*(grid_pt->u[rk_flag][b]) )*(grid_pt->theta_u[rk_flag])/3.

                + gamma/tau*(DATA->gmunu[a][3])*(DATA->gmunu[b][3])

                - ueta/tau/2.*( (DATA->gmunu[a][3])*(DATA->gmunu[b][0]) + (DATA->gmunu[b][3])*(DATA->gmunu[a][0]) )

                + ueta*gamma/tau/2.*( (DATA->gmunu[a][3])*(grid_pt->u[rk_flag][b]) +( DATA->gmunu[b][3])*(grid_pt->u[rk_flag][a]) )

                - ueta*ueta/tau/2.*( (DATA->gmunu[a][0])*(grid_pt->u[rk_flag][b]) + (DATA->gmunu[b][0])*(grid_pt->u[rk_flag][a]) )

                + ( grid_pt->u[rk_flag][a]*grid_pt->a[rk_flag][b] + grid_pt->u[rk_flag][b]*grid_pt->a[rk_flag][a] )/2. ;

    }
    }

/// full Navier-Stokes term is
    NS_term = -2.*shear*sigma[mu][nu] ; // sign changes according to metric sign convention


/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///                             Vorticity Term                             ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
double omega[4][4];
double term1_Vorticity;
double Vorticity_term;

// / remember: dUsup[m][n] = partial^n u^m  ///
// / remember:  a[n]  =  u^m*partial_m u^n  ///

    for( a=0;a<=3;a++ ){
    for( b=0;b<=3;b++ ){

    omega[a][b] = ( grid_pt->dUsup[rk_flag][a][b] - grid_pt->dUsup[rk_flag][b][a] )/2.
                + ueta/tau/2.*( DATA->gmunu[a][0]*DATA->gmunu[b][3] - DATA->gmunu[b][0]*DATA->gmunu[a][3] )
                - ueta*gamma/tau/2.*( DATA->gmunu[a][3]*grid_pt->u[rk_flag][b] - DATA->gmunu[b][3]*grid_pt->u[rk_flag][a] )
                + ueta*ueta/tau/2.*( DATA->gmunu[a][0]*grid_pt->u[rk_flag][b] - DATA->gmunu[b][0]*grid_pt->u[rk_flag][a] )
                + ( grid_pt->u[rk_flag][a]*grid_pt->a[rk_flag][b] - grid_pt->u[rk_flag][b]*grid_pt->a[rk_flag][a] )/2. ;

    }}

   term1_Vorticity = (  -grid_pt->Wmunu[rk_flag][mu][0]*omega[nu][0] - grid_pt->Wmunu[rk_flag][nu][0]*omega[mu][0]
                       + grid_pt->Wmunu[rk_flag][mu][1]*omega[nu][1] + grid_pt->Wmunu[rk_flag][nu][1]*omega[mu][1]
                       + grid_pt->Wmunu[rk_flag][mu][2]*omega[nu][2] + grid_pt->Wmunu[rk_flag][nu][2]*omega[mu][2]
                       + grid_pt->Wmunu[rk_flag][mu][3]*omega[nu][3] + grid_pt->Wmunu[rk_flag][nu][3]*omega[mu][3] )/2.;

/// multiply term by its respective transport coefficient
   term1_Vorticity = transport_coefficient4*term1_Vorticity;

/// full term is
    Vorticity_term = term1_Vorticity ;


/// //////////////////////////////////////////////////////////////////////////// ///
/// //////////////////////////////////////////////////////////////////////////// ///
///                   Add nonlinear term in shear-stress tensor                  ///
///   transport_coefficient3*Delta(mu nu)(alpha beta)*Wmu gamma sigma nu gamma   ///
/// //////////////////////////////////////////////////////////////////////////// ///
/// //////////////////////////////////////////////////////////////////////////// ///
double Wsigma, Wsigma_term;
double term1_Wsigma, term2_Wsigma;

    Wsigma =
  (
   grid_pt->Wmunu[rk_flag][0][0]*sigma[0][0]
  +grid_pt->Wmunu[rk_flag][1][1]*sigma[1][1]
  +grid_pt->Wmunu[rk_flag][2][2]*sigma[2][2]
  +grid_pt->Wmunu[rk_flag][3][3]*sigma[3][3]

  -2.*(
        grid_pt->Wmunu[rk_flag][0][1]*sigma[0][1]
       +grid_pt->Wmunu[rk_flag][0][2]*sigma[0][2]
       +grid_pt->Wmunu[rk_flag][0][3]*sigma[0][3]
       )

  +2.*(
        grid_pt->Wmunu[rk_flag][1][2]*sigma[1][2]
       +grid_pt->Wmunu[rk_flag][1][3]*sigma[1][3]
       +grid_pt->Wmunu[rk_flag][2][3]*sigma[2][3]
       )
   );

   term1_Wsigma = (  - grid_pt->Wmunu[rk_flag][mu][0]*sigma[nu][0] - grid_pt->Wmunu[rk_flag][nu][0]*sigma[mu][0]
                     + grid_pt->Wmunu[rk_flag][mu][1]*sigma[nu][1] + grid_pt->Wmunu[rk_flag][nu][1]*sigma[mu][1]
                     + grid_pt->Wmunu[rk_flag][mu][2]*sigma[nu][2] + grid_pt->Wmunu[rk_flag][nu][2]*sigma[mu][2]
                     + grid_pt->Wmunu[rk_flag][mu][3]*sigma[nu][3] + grid_pt->Wmunu[rk_flag][nu][3]*sigma[mu][3] )/2.;

   term2_Wsigma = -(1./3.)*( DATA->gmunu[mu][nu] + grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][nu] )*Wsigma;

/// multiply term by its respective transport coefficient
   term1_Wsigma = transport_coefficient3*term1_Wsigma;
   term2_Wsigma = transport_coefficient3*term2_Wsigma;

/// full term is
   Wsigma_term = -term1_Wsigma - term2_Wsigma ;

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///               Add nonlinear term in shear-stress tensor                ///
///   transport_coefficient*Delta(mu nu)(alpha beta)*Wmu gamma Wnu gamma   ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
    double Wsquare, WW_term;
    double term1_WW, term2_WW;

    Wsquare =
  (
   grid_pt->Wmunu[rk_flag][0][0]*grid_pt->Wmunu[rk_flag][0][0]
  +grid_pt->Wmunu[rk_flag][1][1]*grid_pt->Wmunu[rk_flag][1][1]
  +grid_pt->Wmunu[rk_flag][2][2]*grid_pt->Wmunu[rk_flag][2][2]
  +grid_pt->Wmunu[rk_flag][3][3]*grid_pt->Wmunu[rk_flag][3][3]

  -2.*(
       grid_pt->Wmunu[rk_flag][0][1]*grid_pt->Wmunu[rk_flag][0][1]
       +grid_pt->Wmunu[rk_flag][0][2]*grid_pt->Wmunu[rk_flag][0][2]
       +grid_pt->Wmunu[rk_flag][0][3]*grid_pt->Wmunu[rk_flag][0][3]
       )

  +2.*(
       grid_pt->Wmunu[rk_flag][1][2]*grid_pt->Wmunu[rk_flag][1][2]
       +grid_pt->Wmunu[rk_flag][1][3]*grid_pt->Wmunu[rk_flag][1][3]
       +grid_pt->Wmunu[rk_flag][2][3]*grid_pt->Wmunu[rk_flag][2][3]
       )
   );

   term1_WW = - grid_pt->Wmunu[rk_flag][mu][0]*grid_pt->Wmunu[rk_flag][nu][0]
           + grid_pt->Wmunu[rk_flag][mu][1]*grid_pt->Wmunu[rk_flag][nu][1]
           + grid_pt->Wmunu[rk_flag][mu][2]*grid_pt->Wmunu[rk_flag][nu][2]
           + grid_pt->Wmunu[rk_flag][mu][3]*grid_pt->Wmunu[rk_flag][nu][3] ;

   term2_WW = -(1./3.)*( DATA->gmunu[mu][nu] + grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][nu] )*Wsquare;

/// multiply term by its respective transport coefficient
   term1_WW = term1_WW*transport_coefficient;
   term2_WW = term2_WW*transport_coefficient;

/// full term is
   WW_term = -term1_WW - term2_WW ; // sign changes according to metric sign convention


/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///               Add coupling to bulk viscous pressure                    ///
///              transport_coefficient_b*Bulk*sigma^mu nu                  ///
///               transport_coefficient2_b*Bulk*W^mu nu                    ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
    double Bulk_Sigma, Bulk_Sigma_term;
    double Bulk_W, Bulk_W_term;
    double Coupling_to_Bulk;

    Bulk_Sigma = grid_pt->pi_b[rk_flag]*sigma[mu][nu];

    Bulk_W = grid_pt->pi_b[rk_flag]*grid_pt->Wmunu[rk_flag][mu][nu];

/// multiply term by its respective transport coefficient
   Bulk_Sigma_term = Bulk_Sigma*transport_coefficient_b;
   Bulk_W_term     = Bulk_W*transport_coefficient2_b;

/// full term is
   Coupling_to_Bulk = -Bulk_Sigma_term + Bulk_W_term ;  // first term: sign changes according to metric sign convention


/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///

/// Comment the appropriate line in order to include the corresponding new term
WW_term=0.;
Vorticity_term = 0.;
Wsigma_term = 0.;


/// final answer is
  SW = ( NS_term + tempf + Vorticity_term + Wsigma_term + WW_term + Coupling_to_Bulk)/(tau_pi);



    return SW;
}/* Make_uWSource */


int Diss::Make_uWRHS(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			 Grid *Lneighbor2, Grid *Rneighbor2, double **w_rhs, InitData *DATA, int rk_flag, int size, int rank)
{
 int mu, nu, direc, nmax[4], ic;
 double f, fp1, fm1, fp2, fm2, delta[4];
//  double ux;
 double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
 double uWphR, uWphL, uWmhR, uWmhL, WphR, WphL, WmhR, WmhL;
 double HWph, HWmh, taufactor, HW, ic_fac;
//  double SW;
/*  HW[4][4][4], SW[4][4][4] */
 double tempf, sum, shear_on;
//  double tempg;

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

    sum += tempf;

     w_rhs[mu][nu] = sum*(DATA->delta_tau)*shear_on;

    }/* nu */
  }/* mu */

 return 1; /* if successful */
}/* Make_uWRHS */



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
 int direc, nmax[4];
 double f, fp1, fm1, fp2, fm2, delta[4];
//  double ux;
 double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
 double uPiphR, uPiphL, uPimhR, uPimhL, PiphR, PiphL, PimhR, PimhL;
 double HPiph, HPimh, taufactor, HPi;
 double sum;
 double bulk_on;

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
   	else {fprintf(stderr,"direc out of range.\n");exit(0);}

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
 double s_den, shear, bulk;
 double Bulk_Relax_time, transport_coeff1, transport_coeff2, transport_coeff1_s, transport_coeff2_s;
 double NS_term, BB_term;
 double Final_Answer;
 

/// Useful variables to define
double gamma, ueta, cs2;
gamma = grid_pt->u[rk_flag][0];
ueta  = grid_pt->u[rk_flag][3];


    if(DATA->turn_on_bulk == 0) return 0.0;

    /// defining bulk viscosity coefficient
    s_den = eos->get_entropy(grid_pt->epsilon, grid_pt->rhob);
    cs2 = eos->p_e_func(grid_pt->epsilon, grid_pt->rhob);    // cs2 is the velocity of sound squared
    shear = (DATA->shear_to_s)*s_den;                   // shear viscosity = constant * entropy density
    bulk = (DATA->bulk_to_s)*shear*(1./3.-cs2)*(1./3.-cs2);  // bulk viscosity = constant * shear viscosity * (1/3-cs2)**2
   // parameter DATA->bulk_to_s should be between 15 -- 75

    /// defining bulk relaxation time and additional transport coefficients
    Bulk_Relax_time    = 1./14.55/(1./3.-cs2)/(1./3.-cs2)/(grid_pt->epsilon + grid_pt->p)*bulk; // Bulk relaxation time from kinetic theory
    transport_coeff1   = 2.0/3.0*(Bulk_Relax_time);          /// from kinetic theory, small mass limit
    transport_coeff2   = 0.;                                 /// not known; put 0
    transport_coeff1_s = 8./5.*(1./3.-cs2)*Bulk_Relax_time;  /// from kinetic theory
    transport_coeff2_s = 0.;                                 /// not known;  put 0

    /// Computing Navier-Stokes term (-bulk viscosity * theta)
    NS_term = -bulk*(grid_pt->theta_u[rk_flag]);

    /// Computing relaxation term and nonlinear term: - Bulk - transport_coeff1*Bulk*theta
    tempf = -(grid_pt->pi_b[rk_flag]) - transport_coeff1*(grid_pt->theta_u[rk_flag])*(grid_pt->pi_b[rk_flag]);

    /// Computing nonlinear term: + transport_coeff2*Bulk*Bulk
    BB_term = transport_coeff2*(grid_pt->pi_b[rk_flag])*(grid_pt->pi_b[rk_flag]);


    /// Computing sigma^mu^nu
    int a, b;
    double sigma[4][4];

    for( a=0;a<=3;a++ ){
    for( b=0;b<=3;b++ ){

    sigma[a][b] = ( grid_pt->dUsup[rk_flag][a][b] + grid_pt->dUsup[rk_flag][b][a] )/2.

                - ( DATA->gmunu[a][b] + (grid_pt->u[rk_flag][a])*(grid_pt->u[rk_flag][b]) )*(grid_pt->theta_u[rk_flag])/3.

                + gamma/tau*(DATA->gmunu[a][3])*(DATA->gmunu[b][3])

                - ueta/tau/2.*( (DATA->gmunu[a][3])*(DATA->gmunu[b][0]) + (DATA->gmunu[b][3])*(DATA->gmunu[a][0]) )

                + ueta*gamma/tau/2.*( (DATA->gmunu[a][3])*(grid_pt->u[rk_flag][b]) +( DATA->gmunu[b][3])*(grid_pt->u[rk_flag][a]) )

                - ueta*ueta/tau/2.*( (DATA->gmunu[a][0])*(grid_pt->u[rk_flag][b]) + (DATA->gmunu[b][0])*(grid_pt->u[rk_flag][a]) )

                + ( grid_pt->u[rk_flag][a]*grid_pt->a[rk_flag][b] + grid_pt->u[rk_flag][b]*grid_pt->a[rk_flag][a] )/2. ;

    }
    }

    /// Computing terms that Couple with shear-stress tensor
    double Wsigma, WW, Shear_Sigma_term, Shear_Shear_term, Coupling_to_Shear;

    Wsigma =
  (
   grid_pt->Wmunu[rk_flag][0][0]*sigma[0][0]
  +grid_pt->Wmunu[rk_flag][1][1]*sigma[1][1]
  +grid_pt->Wmunu[rk_flag][2][2]*sigma[2][2]
  +grid_pt->Wmunu[rk_flag][3][3]*sigma[3][3]

  -2.*(
        grid_pt->Wmunu[rk_flag][0][1]*sigma[0][1]
       +grid_pt->Wmunu[rk_flag][0][2]*sigma[0][2]
       +grid_pt->Wmunu[rk_flag][0][3]*sigma[0][3]
       )

  +2.*(
        grid_pt->Wmunu[rk_flag][1][2]*sigma[1][2]
       +grid_pt->Wmunu[rk_flag][1][3]*sigma[1][3]
       +grid_pt->Wmunu[rk_flag][2][3]*sigma[2][3]
       )
   );

    WW =
  (
   grid_pt->Wmunu[rk_flag][0][0]*grid_pt->Wmunu[rk_flag][0][0]
  +grid_pt->Wmunu[rk_flag][1][1]*grid_pt->Wmunu[rk_flag][1][1]
  +grid_pt->Wmunu[rk_flag][2][2]*grid_pt->Wmunu[rk_flag][2][2]
  +grid_pt->Wmunu[rk_flag][3][3]*grid_pt->Wmunu[rk_flag][3][3]

  -2.*(
       grid_pt->Wmunu[rk_flag][0][1]*grid_pt->Wmunu[rk_flag][0][1]
       +grid_pt->Wmunu[rk_flag][0][2]*grid_pt->Wmunu[rk_flag][0][2]
       +grid_pt->Wmunu[rk_flag][0][3]*grid_pt->Wmunu[rk_flag][0][3]
       )

  +2.*(
       grid_pt->Wmunu[rk_flag][1][2]*grid_pt->Wmunu[rk_flag][1][2]
       +grid_pt->Wmunu[rk_flag][1][3]*grid_pt->Wmunu[rk_flag][1][3]
       +grid_pt->Wmunu[rk_flag][2][3]*grid_pt->Wmunu[rk_flag][2][3]
       )
   );

    /// multiply term by its respective transport coefficient
    Shear_Sigma_term = Wsigma*transport_coeff1_s;
    Shear_Shear_term = WW*transport_coeff2_s;

    /// full term that couples to shear is
    Coupling_to_Shear = -Shear_Sigma_term + Shear_Shear_term ;

    /// Final Answer
    Final_Answer = NS_term + tempf + BB_term + Coupling_to_Shear;

    return Final_Answer/(Bulk_Relax_time);
}/* Make_uPiSource */


