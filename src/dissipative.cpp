#include "util.h"
#include "grid.h"
#include "data.h"
#include "eos.h"
#include "dissipative.h"

using namespace std;

Diss::Diss(EOS *eosIn, InitData* DATA_in)
{
    eos = eosIn;
    minmod = new Minmod(DATA_in);
}

// destructor
Diss::~Diss()
{
    delete minmod;
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Dissipative parts */
/* Sangyong Nov 18 2014 */
/* change: alpha first which is the case
for everywhere else. also, this change is necessary
to use Wmunu[rk_flag][4][mu] as the dissipative baryon current*/
/* this is the only one that is being subtracted in the rhs */
double Diss::MakeWSource(double tau, int alpha, Grid *grid_pt, 
                         Grid *Lneighbor, Grid *Rneighbor, 
			       Grid *Lneighbor2, Grid *Rneighbor2, 
                         InitData *DATA,int rk_flag, int size, int rank)
{
    double delta[4], taufactor;
    double shear_on, bulk_on, diff_on;
    int i, nmax[4];

    if(DATA->turn_on_shear)
        shear_on = 1.0;
    else 
        shear_on = 0.0;
    
    if(DATA->turn_on_diff)
        diff_on = 1.0;
    else 
        diff_on = 0.0;

    if(DATA->turn_on_bulk)
        bulk_on = 1.0;
    else 
        bulk_on = 0.0;
    
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
    
    if(alpha == 4 && DATA->turn_on_diff == 0)
        return (0.0);

    /* Sangyong Nov 18 2014 */
    /* change: alpha first which is the case
               for everywhere else. also, this change is necessary
               to use Wmunu[rk_flag][4][mu] as the dissipative baryon current 
    */
    // dW/dtau
    double dWdtau;
    if(rk_flag==0)
    {
       // backward time derivative (first order is more stable) 
       dWdtau = (grid_pt->Wmunu[rk_flag][alpha][0] 
                       - grid_pt->prevWmunu[0][alpha][0])/DATA->delta_tau;
    }
    else if(rk_flag > 0)
    {
       /* first order since we don't know next values yet */
       dWdtau = (grid_pt->Wmunu[rk_flag][alpha][0] 
                       - grid_pt->Wmunu[0][alpha][0])/DATA->delta_tau;
    }
    else
    {
       fprintf(stderr,"rk_flag out of range. \n");
       exit(0);
    }

    /* bulk pressure term */
    double dPidtau;
    if(rk_flag==0)
    {
        /* first order since we don't know next values yet */
        dPidtau  = (grid_pt->Pimunu[rk_flag][alpha][0] 
                    - grid_pt->prevPimunu[rk_flag][alpha][0])/DATA->delta_tau;
    }
    else if(rk_flag > 0)
    {
        /* first order since we don't know next values yet */
        dPidtau  = (grid_pt->Pimunu[rk_flag][alpha][0] 
                    - grid_pt->Pimunu[0][alpha][0])/DATA->delta_tau;
    }

    double dWdx_perp = 0.0;
    double dPidx_perp = 0.0;
    for(i=1; i<=2; i++) // x and y
    {
        double sg = grid_pt->Wmunu[rk_flag][alpha][i];
        double bg = grid_pt->Pimunu[rk_flag][alpha][i];
        double sgp1, sgm1, bgp1, bgm1;
        if(grid_pt->position[i] == nmax[i])
        {
            sgp1 = sg;
            bgp1 = bg;
            sgm1 = grid_pt->nbr_m_1[i]->Wmunu[rk_flag][alpha][i];
            bgm1 = grid_pt->nbr_m_1[i]->Pimunu[rk_flag][alpha][i];
        } 
        else if(grid_pt->position[i] == 0)
        {
            sgp1 = grid_pt->nbr_p_1[i]->Wmunu[rk_flag][alpha][i];
            bgp1 = grid_pt->nbr_p_1[i]->Pimunu[rk_flag][alpha][i];
            sgm1 = sg;
            bgm1 = bg;
        }
        else
        {
            sgp1 = grid_pt->nbr_p_1[i]->Wmunu[rk_flag][alpha][i];
            bgp1 = grid_pt->nbr_p_1[i]->Pimunu[rk_flag][alpha][i];
            sgm1 = grid_pt->nbr_m_1[i]->Wmunu[rk_flag][alpha][i];
            bgm1 = grid_pt->nbr_m_1[i]->Pimunu[rk_flag][alpha][i];
        }
   
        dWdx_perp += minmod->minmod_dx(sgp1, sg, sgm1)/delta[i]; 
        dPidx_perp += minmod->minmod_dx(bgp1, bg, bgm1)/delta[i]; 
    }/* i */
    
    i=3;
    taufactor = tau;
    double dWdeta, dPideta;
  
    double sg = grid_pt->Wmunu[rk_flag][alpha][i];
    double bg = grid_pt->Pimunu[rk_flag][alpha][i];
    double sgp1, sgm1, bgp1, bgm1;
    if(grid_pt->position[i] == nmax[i])
    {
        // for the right most rank do boundary condition on the right
        if(rank == size-1) 
        {
            sgp1 = sg;
            bgp1 = bg;
            sgm1 = grid_pt->nbr_m_1[i]->Wmunu[rk_flag][alpha][i];
            bgm1 = grid_pt->nbr_m_1[i]->Pimunu[rk_flag][alpha][i];
        }
        else
        {
            sgp1 = Rneighbor->Wmunu[rk_flag][alpha][i];
            bgp1 = Rneighbor->Pimunu[rk_flag][alpha][i];
            sgm1 = grid_pt->nbr_m_1[i]->Wmunu[rk_flag][alpha][i];
            bgm1 = grid_pt->nbr_m_1[i]->Pimunu[rk_flag][alpha][i];
        }
    } 
    else if(grid_pt->position[i] == 0)
    {
        // for the left most rank do boundary condition on the left
        if(rank == 0) 
        {
            sgp1 = grid_pt->nbr_p_1[i]->Wmunu[rk_flag][alpha][i];
            bgp1 = grid_pt->nbr_p_1[i]->Pimunu[rk_flag][alpha][i];
            sgm1 = sg;
            bgm1 = bg;
        }
        else
        {
            sgp1 = grid_pt->nbr_p_1[i]->Wmunu[rk_flag][alpha][i];
            bgp1 = grid_pt->nbr_p_1[i]->Pimunu[rk_flag][alpha][i];
            sgm1 = Lneighbor->Wmunu[rk_flag][alpha][i];
            bgm1 = Lneighbor->Pimunu[rk_flag][alpha][i];
        }
    }
    else
    {
        sgp1 = grid_pt->nbr_p_1[i]->Wmunu[rk_flag][alpha][i];
        bgp1 = grid_pt->nbr_p_1[i]->Pimunu[rk_flag][alpha][i];
        sgm1 = grid_pt->nbr_m_1[i]->Wmunu[rk_flag][alpha][i];
        bgm1 = grid_pt->nbr_m_1[i]->Pimunu[rk_flag][alpha][i];
    }
    dWdeta = minmod->minmod_dx(sgp1, sg, sgm1)/delta[i]/taufactor; 
    dPideta = minmod->minmod_dx(bgp1, bg, bgm1)/delta[i]/taufactor; 

    /* partial_m (tau W^mn) = W^0n + tau partial_m W^mn */
    double sf = (tau*(dWdtau + dWdx_perp + dWdeta) 
                 + grid_pt->Wmunu[rk_flag][alpha][0]);
    double bf = (tau*(dPidtau + dPidx_perp + dPideta) 
                 + grid_pt->Pimunu[rk_flag][alpha][0]);

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

    // final result
    double result;
    if(alpha < 4)
        result = (sf*shear_on + bf*bulk_on);

    if(alpha == 4)
        result = sf*diff_on;

    if(isnan(result)) 
    {
        cout << "sf=" << sf << " bf=" << bf 
             << " Wmunu[" << rk_flag << "]=" 
             << grid_pt->Wmunu[rk_flag][alpha][0]
             << " Pimunu[" << rk_flag << "]=" 
             << grid_pt->Pimunu[rk_flag][alpha][0]
             << " prevWmunu=" << grid_pt->prevWmunu[0][alpha][0] << endl;
    }
    return result;
}/* MakeWSource */
/* MakeWSource is for Tmunu */


double Diss::Make_uWSource(double tau, Grid *grid_pt, int mu, int nu, 
                           InitData *DATA, int rk_flag)
{
    double tempf, tau_pi;
    double SW, shear, shear_to_s, T, epsilon, rhob, Ttr;
    int a, b;
    double sigma[4][4], gamma, ueta;
    double NS_term;
    
    if(DATA->turn_on_shear == 0)
        return 0.0;

    // Useful variables to define
    Ttr = 0.18/hbarc;  // phase transition temperature
    gamma = grid_pt->u[rk_flag][0];
    ueta  = grid_pt->u[rk_flag][3];
    epsilon = grid_pt->epsilon;
    rhob = grid_pt->rhob;
    T=eos->get_temperature(epsilon, rhob);

    if(DATA->T_dependent_shear_to_s == 1)
    {
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

    int include_WWterm = 1;
    int include_Vorticity_term = 0;
    int include_Wsigma_term = 1;

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///                 Defining transport coefficients                        ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
    //double s_den = eos->get_entropy(epsilon, rhob);
    shear = (shear_to_s)*(epsilon + grid_pt->p)/(T + 1e-15);
    tau_pi = 5.0*shear/(epsilon + grid_pt->p + 1e-15);
     
    //tau_pi = maxi(tau_pi, DATA->tau_pi);
    if(!isfinite(tau_pi))
    {
        tau_pi = DATA->delta_tau; 
        cout << "tau_pi was infinite ..." << endl;
    }
    if(tau_pi < DATA->delta_tau)
        tau_pi = DATA->delta_tau;

    // transport coefficient for nonlinear terms 
    // -- shear only terms -- 4Mar2013
    // transport coefficients of a massless gas of single component particles
    double transport_coefficient  = 9./70.*tau_pi/shear*(4./5.) ;
    double transport_coefficient2 = 4./3.*tau_pi;
    double transport_coefficient3 = 10./7.*tau_pi;
    double transport_coefficient4 = 2.*tau_pi;

    // transport coefficient for nonlinear terms 
    // -- coupling to bulk viscous pressure -- 4Mar2013
    // transport coefficients not yet known -- fixed to zero
    double transport_coefficient_b  = 6./5.*tau_pi ;
    double transport_coefficient2_b = 0.;


    /* This source has many terms */
    /* everything in the 1/(tau_pi) piece is here */
    /* third step in the split-operator time evol 
       use Wmunu[rk_flag] and u[rk_flag] with rk_flag = 0 */

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///            Wmunu + transport_coefficient2*Wmunu*theta                  ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///

    // full term is
    tempf = - (1.0 + transport_coefficient2*(grid_pt->theta_u[rk_flag]))
              *(grid_pt->Wmunu[rk_flag][mu][nu]);
    // FOR GUBSER ANALYTIC
    //tempf = -( transport_coefficient2*(grid_pt->theta_u[rk_flag]) )
    //         *(grid_pt->Wmunu[rk_flag][mu][nu]);

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///             Navier-Stokes Term -- -2.*shear*sigma^munu                 ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///

    for( a=0;a<4;a++ )
    {
        for( b=0;b<4;b++ )
        {
            sigma[a][b] = grid_pt->sigma[rk_flag][a][b];
        }
    }

    // full Navier-Stokes term is
    // sign changes according to metric sign convention
    NS_term = -2.*shear*sigma[mu][nu] ; 

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///                             Vorticity Term                             ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
    double omega[4][4];
    double term1_Vorticity;
    double Vorticity_term;

    // remember: dUsup[m][n] = partial^n u^m  ///
    // remember:  a[n]  =  u^m*partial_m u^n  ///
    if(include_Vorticity_term == 1)
    {
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

        // multiply term by its respective transport coefficient
        term1_Vorticity = transport_coefficient4*term1_Vorticity;

        // full term is
        Vorticity_term = term1_Vorticity ;
    }
    else
        Vorticity_term = 0.0;


/// //////////////////////////////////////////////////////////////////////////// ///
/// //////////////////////////////////////////////////////////////////////////// ///
///                   Add nonlinear term in shear-stress tensor                  ///
///   transport_coefficient3*Delta(mu nu)(alpha beta)*Wmu gamma sigma nu gamma   ///
/// //////////////////////////////////////////////////////////////////////////// ///
/// //////////////////////////////////////////////////////////////////////////// ///
    double Wsigma, Wsigma_term;
    double term1_Wsigma, term2_Wsigma;

    if(include_Wsigma_term == 1)
    {
        Wsigma = (
                    grid_pt->Wmunu[rk_flag][0][0]*sigma[0][0]
                  + grid_pt->Wmunu[rk_flag][1][1]*sigma[1][1]
                  + grid_pt->Wmunu[rk_flag][2][2]*sigma[2][2]
                  + grid_pt->Wmunu[rk_flag][3][3]*sigma[3][3]
                  -2.*(
                         grid_pt->Wmunu[rk_flag][0][1]*sigma[0][1]
                       + grid_pt->Wmunu[rk_flag][0][2]*sigma[0][2]
                       + grid_pt->Wmunu[rk_flag][0][3]*sigma[0][3]
                      )
                  +2.*(
                         grid_pt->Wmunu[rk_flag][1][2]*sigma[1][2]
                       + grid_pt->Wmunu[rk_flag][1][3]*sigma[1][3]
                       + grid_pt->Wmunu[rk_flag][2][3]*sigma[2][3]
                      )
                 );

        term1_Wsigma = (  
                  - grid_pt->Wmunu[rk_flag][mu][0]*sigma[nu][0] 
                  - grid_pt->Wmunu[rk_flag][nu][0]*sigma[mu][0]
                  + grid_pt->Wmunu[rk_flag][mu][1]*sigma[nu][1] 
                  + grid_pt->Wmunu[rk_flag][nu][1]*sigma[mu][1]
                  + grid_pt->Wmunu[rk_flag][mu][2]*sigma[nu][2] 
                  + grid_pt->Wmunu[rk_flag][nu][2]*sigma[mu][2]
                  + grid_pt->Wmunu[rk_flag][mu][3]*sigma[nu][3] 
                  + grid_pt->Wmunu[rk_flag][nu][3]*sigma[mu][3] )/2.;
        
        term2_Wsigma = -(1./3.)*(
            DATA->gmunu[mu][nu] 
            + grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][nu])*Wsigma; 

        // multiply term by its respective transport coefficient
        term1_Wsigma = transport_coefficient3*term1_Wsigma;
        term2_Wsigma = transport_coefficient3*term2_Wsigma;

        // full term is
        Wsigma_term = -term1_Wsigma - term2_Wsigma ;
    }
    else
        Wsigma_term = 0.0;

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///               Add nonlinear term in shear-stress tensor                ///
///   transport_coefficient*Delta(mu nu)(alpha beta)*Wmu gamma Wnu gamma   ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
    double Wsquare, WW_term;
    double term1_WW, term2_WW;
    if(include_WWterm == 1)
    {
        Wsquare = (
              grid_pt->Wmunu[rk_flag][0][0]*grid_pt->Wmunu[rk_flag][0][0]
            + grid_pt->Wmunu[rk_flag][1][1]*grid_pt->Wmunu[rk_flag][1][1]
            + grid_pt->Wmunu[rk_flag][2][2]*grid_pt->Wmunu[rk_flag][2][2]
            + grid_pt->Wmunu[rk_flag][3][3]*grid_pt->Wmunu[rk_flag][3][3]
            -2.*(
                   grid_pt->Wmunu[rk_flag][0][1]*grid_pt->Wmunu[rk_flag][0][1]
                 + grid_pt->Wmunu[rk_flag][0][2]*grid_pt->Wmunu[rk_flag][0][2]
                 + grid_pt->Wmunu[rk_flag][0][3]*grid_pt->Wmunu[rk_flag][0][3]
                )
            +2.*(
                   grid_pt->Wmunu[rk_flag][1][2]*grid_pt->Wmunu[rk_flag][1][2]
                 + grid_pt->Wmunu[rk_flag][1][3]*grid_pt->Wmunu[rk_flag][1][3]
                 + grid_pt->Wmunu[rk_flag][2][3]*grid_pt->Wmunu[rk_flag][2][3]
                ));

        term1_WW = (
            - grid_pt->Wmunu[rk_flag][mu][0]*grid_pt->Wmunu[rk_flag][nu][0]
            + grid_pt->Wmunu[rk_flag][mu][1]*grid_pt->Wmunu[rk_flag][nu][1]
            + grid_pt->Wmunu[rk_flag][mu][2]*grid_pt->Wmunu[rk_flag][nu][2]
            + grid_pt->Wmunu[rk_flag][mu][3]*grid_pt->Wmunu[rk_flag][nu][3] );

        term2_WW = -(1./3.)*( 
            DATA->gmunu[mu][nu] 
            + grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][nu] )*Wsquare;

        //multiply term by its respective transport coefficient
        term1_WW = term1_WW*transport_coefficient;
        term2_WW = term2_WW*transport_coefficient;

        // full term is
        WW_term = -term1_WW - term2_WW ; 
        // sign changes according to metric sign convention
    }
    else
        WW_term = 0.0;

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

    // multiply term by its respective transport coefficient
    Bulk_Sigma_term = Bulk_Sigma*transport_coefficient_b;
    Bulk_W_term     = Bulk_W*transport_coefficient2_b;

    // full term is
    Coupling_to_Bulk = -Bulk_Sigma_term + Bulk_W_term ;  
    // first term: sign changes according to metric sign convention


/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///


    // final answer is
    SW = (  NS_term + tempf + Vorticity_term + Wsigma_term 
          + WW_term + Coupling_to_Bulk)/(tau_pi);

    return SW;
}/* Make_uWSource */


int Diss::Make_uWRHS(double tau, Grid *grid_pt, 
                     Grid *Lneighbor, Grid *Rneighbor, 
			   Grid *Lneighbor2, Grid *Rneighbor2, 
                     double **w_rhs, InitData *DATA, int rk_flag, 
                     int size, int rank)
// Kurganov-Tadmor for Wmunu 
// implement 
// partial_tau (utau Wmn) + (1/tau)partial_eta (ueta Wmn) 
// + partial_x (ux Wmn) + partial_y (uy Wmn) + utau Wmn/tau = SW 
// or the right hand side of,
// partial_tau (utau Wmn) = 
// - (1/tau)partial_eta (ueta Wmn) - partial_x (ux Wmn) - partial_y (uy Wmn) 
// - utau Wmn/tau + SW 
//
// the local velocity is just u_x/u_tau, u_y/u_tau, u_eta/tau/u_tau 
// KT flux is given by 
// H_{j+1/2} = (fRph + fLph)/2 - ax(uRph - uLph) 
// Here fRph = ux WmnRph and ax uRph = |ux/utau|_max utau Wmn 
//
// This is the second step in the operator splitting. it uses
// rk_flag+1 as initial condition 
{
    int mu, nu, direc, nmax[4], ic;
    double f, fp1, fm1, fp2, fm2, delta[4];
    double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
    double uWphR, uWphL, uWmhR, uWmhL, WphR, WphL, WmhR, WmhL;
    double HWph, HWmh, taufactor, HW, ic_fac;
    double tempf, sum, shear_on;

    nmax[1] = DATA->nx;
    nmax[2] = DATA->ny;
    nmax[3] = DATA->neta-1;

    delta[1] = DATA->delta_x;
    delta[2] = DATA->delta_y;
    delta[3] = DATA->delta_eta;

    if(DATA->turn_on_shear)
        shear_on = 1.0;
    else 
        shear_on = 0.0;

    for(mu=0; mu<4; mu++)
    {
        for(nu=0; nu<4; nu++)
        {
            sum = 0.0;
            for(direc=1; direc<=3; direc++)
            {
                if(direc==3) 
                    taufactor = tau;
                else 
                    taufactor = 1.0;
                
                /* Get_uWmns */
                Get_uWmns(tau, grid_pt, Lneighbor, Rneighbor, 
                          Lneighbor2, Rneighbor2, mu, nu, direc, 
                          &g, &f, &gp1, &fp1, &gp2, &fp2, 
		              &gm1, &fm1, &gm2, &fm2, DATA, rk_flag, size, rank);
                
                /* MakeuWmnHalfs */
                // uWmn 
                uWphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f); 
                uWphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
                uWmhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
                uWmhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);
                
                /* just Wmn */
                WphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g); 
                WphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
                WmhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
                WmhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);

                a = fabs(grid_pt->u[rk_flag][direc])/grid_pt->u[rk_flag][0];

                if(direc<3) // x,y direction
	          {
                    if(grid_pt->position[direc] == 0)
                    {
                        am1 = a;
                    }
                    else
                    {
                        am1 = (fabs(grid_pt->nbr_m_1[direc]->u[rk_flag][direc])
                               /grid_pt->nbr_m_1[direc]->u[rk_flag][0]);
                    }
                    if(grid_pt->position[direc] == nmax[direc])
                    {
                        ap1 = a;
                    }
                    else
                    {
                        ap1 = (fabs(grid_pt->nbr_p_1[direc]->u[rk_flag][direc])
                               /grid_pt->nbr_p_1[direc]->u[rk_flag][0]);
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
                            am1 = (fabs(Lneighbor->u[rk_flag][direc])
                                    /Lneighbor->u[rk_flag][0]);
                        }
                    }
                    else
                    {
                        am1 = (fabs(grid_pt->nbr_m_1[direc]->u[rk_flag][direc])
                               /grid_pt->nbr_m_1[direc]->u[rk_flag][0]);
                    }
                    if(grid_pt->position[direc] == nmax[direc])
                    {
                        if (rank==size-1)
                        {
                            ap1 = a;
                        }
                        else
                        {
                            ap1 = (fabs(Rneighbor->u[rk_flag][direc])
                                   /Rneighbor->u[rk_flag][0]);
                        }
                    }
                    else
                    {
                        ap1 = (fabs(grid_pt->nbr_p_1[direc]->u[rk_flag][direc])
                               /grid_pt->nbr_p_1[direc]->u[rk_flag][0]);
                    }
                }

                ax = maxi(a, ap1);
                HWph = ((uWphR + uWphL) - ax*(WphR - WphL))*0.5;

                ax = maxi(a, am1); 
                HWmh = ((uWmhR + uWmhL) - ax*(WmhR - WmhL))*0.5;

                HW = (HWph - HWmh)/delta[direc]/taufactor;
                
                // make partial_i (u^i Wmn)
                sum += -HW;
            }/* direction */
            
            // add a source term -u^tau Wmn/tau due to 
            // the coordinate change to tau-eta
            sum += (
                - grid_pt->u[rk_flag][0]*grid_pt->Wmunu[rk_flag][mu][nu]/tau
                + grid_pt->theta_u[rk_flag]*grid_pt->Wmunu[rk_flag][mu][nu]  );
   
            // this is from udW = d(uW) - Wdu = RHS
            // or d(uW) = udW + Wdu
            // this term is being added to the rhs so that -4/3 + 1 = -1/3
            // other source terms due to the coordinate change to tau-eta
            tempf = 0.0;
            tempf = (( - DATA->gmunu[3][mu]*grid_pt->Wmunu[rk_flag][0][nu]
                       - DATA->gmunu[3][nu]*grid_pt->Wmunu[rk_flag][0][mu]
                       + DATA->gmunu[0][mu]*grid_pt->Wmunu[rk_flag][3][nu]
                       + DATA->gmunu[0][nu]*grid_pt->Wmunu[rk_flag][3][mu]
                       + grid_pt->Wmunu[rk_flag][3][nu]*grid_pt->u[rk_flag][mu]
                         *grid_pt->u[rk_flag][0]
                       + grid_pt->Wmunu[rk_flag][3][mu]*grid_pt->u[rk_flag][nu]
                         *grid_pt->u[rk_flag][0]
                       - grid_pt->Wmunu[rk_flag][0][nu]*grid_pt->u[rk_flag][mu]
                         *grid_pt->u[rk_flag][3]
                       - grid_pt->Wmunu[rk_flag][0][mu]*grid_pt->u[rk_flag][nu]
                         *grid_pt->u[rk_flag][3]
                     )*(grid_pt->u[rk_flag][3]/tau));
            
            for(ic=0; ic<4; ic++)
            {
                ic_fac = (ic==0 ? -1.0 : 1.0);
                
                tempf += (
                      grid_pt->Wmunu[rk_flag][ic][nu]*grid_pt->u[rk_flag][mu]
                      *grid_pt->a[rk_flag][ic]*ic_fac
                    + grid_pt->Wmunu[rk_flag][ic][mu]*grid_pt->u[rk_flag][nu]
                      *grid_pt->a[rk_flag][ic])*ic_fac;
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
	       
	       tgm2 = grid_pt->nbr_m_2[direc]->Wmunu[rk_flag][mu][nu];
	       tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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
	       
	       tgm2 = grid_pt->nbr_m_2[direc]->Wmunu[rk_flag][mu][nu];
	       tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
	     }
	   else if(grid_pt->position[direc] == 0)
	     {
	       tgp2 = grid_pt->nbr_p_2[direc]->Wmunu[rk_flag][mu][nu];
	       tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
	       
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
	       tgp2 = grid_pt->nbr_p_2[direc]->Wmunu[rk_flag][mu][nu];
	       tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
	       
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
	       tgp2 = grid_pt->nbr_p_2[direc]->Wmunu[rk_flag][mu][nu];
	       tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_2[direc]->Wmunu[rk_flag][mu][nu];
	       tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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
		   
		   tgm2 = grid_pt->nbr_m_2[direc]->Wmunu[rk_flag][mu][nu];
		   tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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
		   
		   tgm2 = grid_pt->nbr_m_2[direc]->Wmunu[rk_flag][mu][nu];
		   tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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
		   
		   tgm2 = grid_pt->nbr_m_2[direc]->Wmunu[rk_flag][mu][nu];
		   tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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
		   
		   tgm2 = grid_pt->nbr_m_2[direc]->Wmunu[rk_flag][mu][nu];
		   tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
		 
		 }
	     }
	   else if(grid_pt->position[direc] == 0)
	     {
	      if(rank == 0) // for the left most rank do boundary condition on the left
		{
		  tgp2 = grid_pt->nbr_p_2[direc]->Wmunu[rk_flag][mu][nu];
		  tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
		  tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
		  
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
		  tgp2 = grid_pt->nbr_p_2[direc]->Wmunu[rk_flag][mu][nu];
		  tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
		  tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
		  
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
		   tgp2 = grid_pt->nbr_p_2[direc]->Wmunu[rk_flag][mu][nu];
		   tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
		   
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
		   tgp2 = grid_pt->nbr_p_2[direc]->Wmunu[rk_flag][mu][nu];
		   tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
		   
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
	       tgp2 = grid_pt->nbr_p_2[direc]->Wmunu[rk_flag][mu][nu];
	       tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->Wmunu[rk_flag][mu][nu];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_2[direc]->Wmunu[rk_flag][mu][nu];
	       tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
	     }
	 }

       //cout << "tg=" << tg << ", tf=" << tf << ", tfp1=" << tfp1 
       //     << ", tfm1=" << tfm1 << endl;

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



int Diss::Make_uPRHS(
    double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
    Grid *Lneighbor2, Grid *Rneighbor2, double *p_rhs, InitData *DATA, 
    int rk_flag, int size, int rank)
{
    int direc, nmax[4];
    double f, fp1, fm1, fp2, fm2, delta[4];
    //double ux;
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

    if(DATA->turn_on_bulk)
        bulk_on = 1.0;
    else 
        bulk_on = 0.0;

     sum = 0.0;
     for(direc=1; direc<=3; direc++)
     {
         if(direc==3) 
            taufactor = tau;
         else 
            taufactor = 1.0;

         /* Get_uPis */
         Get_uPis(tau, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, 
                  direc, &g, &f, &gp1, &fp1, &gp2, &fp2, 
		      &gm1, &fm1, &gm2, &fm2, DATA, rk_flag, size, rank);

         /*  Make upi Halfs */
         /* uPi */
         uPiphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f); 
         uPiphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
         uPimhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
         uPimhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);

         /* just Pi */
         PiphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g); 
         PiphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
         PimhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
         PimhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);

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
         HPiph = ((uPiphR + uPiphL) - ax*(PiphR - PiphL))*0.5;
         
         ax = maxi(a, am1); 
         HPimh = ((uPimhR + uPimhL) - ax*(PimhR - PimhL))*0.5;
      
         HPi = (HPiph - HPimh)/delta[direc]/taufactor;

         /* make partial_i (u^i Pi) */
         sum += -HPi;

     }/* direction */
       
     /* add a source term due to the coordinate change to tau-eta */
     sum -= (grid_pt->pi_b[rk_flag])*(grid_pt->u[rk_flag][0])/tau;
     
     sum += (grid_pt->pi_b[rk_flag])*(grid_pt->theta_u[rk_flag]);

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
	       
	       tgm2 = grid_pt->nbr_m_2[direc]->pi_b[rk_flag];
	       tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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
	       
	       tgm2 = grid_pt->nbr_m_2[direc]->pi_b[rk_flag];
	       tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
	     }
	   else if(grid_pt->position[direc] == 0)
	     {
	       tgp2 = grid_pt->nbr_p_2[direc]->pi_b[rk_flag];
	       tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
	       
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
	       tgp2 = grid_pt->nbr_p_2[direc]->pi_b[rk_flag];
	       tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
	       
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
	       tgp2 = grid_pt->nbr_p_2[direc]->pi_b[rk_flag];
	       tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
	       
	       tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_2[direc]->pi_b[rk_flag];
	       tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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
		   
		   tgm2 = grid_pt->nbr_m_2[direc]->pi_b[rk_flag];
		   tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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
		   
		   tgm2 = grid_pt->nbr_m_2[direc]->pi_b[rk_flag];
		   tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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
		   
		   tgm2 = grid_pt->nbr_m_2[direc]->pi_b[rk_flag];
		   tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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
		   
		   tgm2 = grid_pt->nbr_m_2[direc]->pi_b[rk_flag];
		   tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
		   tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
		 }
	     }
	   else if(grid_pt->position[direc] == 0)
	     {
	       if(rank==0) // for left most rank use boundary condition
		 {
		   tgp2 = grid_pt->nbr_p_2[direc]->pi_b[rk_flag];
		   tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
		   
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
		   tgp2 = grid_pt->nbr_p_2[direc]->pi_b[rk_flag];
		   tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
		   
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
		   tgp2 = grid_pt->nbr_p_2[direc]->pi_b[rk_flag];
		   tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
		   
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
		   tgp2 = grid_pt->nbr_p_2[direc]->pi_b[rk_flag];
		   tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
		   tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];
		   
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
	       tgp2 = grid_pt->nbr_p_2[direc]->pi_b[rk_flag];
	       tfp2 = tgp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
	       tgp2 *=     grid_pt->nbr_p_2[direc]->u[rk_flag][0];

	       tgp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
	       tfp1 = tgp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
	       tgp1 *=     grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	       
	       tgm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
	       tfm1 = tgm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
	       tgm1 *=     grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	       
	       tgm2 = grid_pt->nbr_m_2[direc]->pi_b[rk_flag];
	       tfm2 = tgm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
	       tgm2 *=     grid_pt->nbr_m_2[direc]->u[rk_flag][0];
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


double Diss::Make_uPiSource(
    double tau, Grid *grid_pt, InitData *DATA, int rk_flag)
{
    double tempf;
    double shear, bulk;
    double Bulk_Relax_time;
    double transport_coeff1, transport_coeff2;
    double transport_coeff1_s, transport_coeff2_s;
    double NS_term, BB_term;
    double Final_Answer;

    // switch to include non-linear coupling terms in the bulk pi evolution
    int include_BBterm = 1;
    int include_coupling_to_shear = 1;
 
    // Useful variables to define
    double gamma, ueta, cs2;
    gamma = grid_pt->u[rk_flag][0];
    ueta  = grid_pt->u[rk_flag][3];

    if(DATA->turn_on_bulk == 0) return 0.0;

    // defining bulk viscosity coefficient

    // shear viscosity = constant * entropy density
    //s_den = eos->get_entropy(grid_pt->epsilon, grid_pt->rhob);
    //shear = (DATA->shear_to_s)*s_den;   

    // shear viscosity = constant * (e + P)/T
    double temperature = eos->get_temperature(grid_pt->epsilon, grid_pt->rhob);
    shear = (DATA->shear_to_s)*(grid_pt->epsilon + grid_pt->p)/temperature;  

    // cs2 is the velocity of sound squared
    cs2 = eos->get_cs2(grid_pt->epsilon, grid_pt->rhob);  

    // bulk viscosity = constant * shear viscosity * (1/3-cs2)**2
    // parameter DATA->bulk_to_s should be between 15 -- 75
    //bulk = (DATA->bulk_to_s)*shear*(1./3.-cs2)*(1./3.-cs2);  

    // T dependent bulk viscosity from Gabriel
    /////////////////////////////////////////////
    //           Parametrization 1             //
    /////////////////////////////////////////////
    double Ttr=0.18/0.1973;
    double dummy=temperature/Ttr;
    double A1=-13.77, A2=27.55, A3=13.45;
    double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;
 
    bulk = A1*dummy*dummy + A2*dummy - A3;
    if(temperature < 0.995*Ttr)
    {
        bulk = lambda3*exp((dummy-1)/sigma3)+ lambda4*exp((dummy-1)/sigma4)+0.03;
    }
    if(temperature > 1.05*Ttr)
    {
        bulk = lambda1*exp(-(dummy-1)/sigma1)+ lambda2*exp(-(dummy-1)/sigma2)+0.001;
    }
    bulk = bulk*(grid_pt->epsilon + grid_pt->p)/temperature;

    /////////////////////////////////////////////
    //           Parametrization 2             //
    /////////////////////////////////////////////
    //double Ttr=0.18/0.1973;
    //double dummy=temperature/Ttr;
    //double A1=-79.53, A2=159.067, A3=79.04;
    //double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    //double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;

    //bulk = A1*dummy*dummy + A2*dummy - A3;

    //if(temperature < 0.997*Ttr)
    //{
    //    bulk = lambda3*exp((dummy-1)/sigma3)+ lambda4*exp((dummy-1)/sigma4)+0.03;
    //}
    //if(temperature > 1.04*Ttr)
    //{
    //    bulk = lambda1*exp(-(dummy-1)/sigma1)+ lambda2*exp(-(dummy-1)/sigma2)+0.001;
    //}
    //bulk = bulk*s_den;

    ////////////////////////////////////////////
    //           Parametrization 3            //
    ////////////////////////////////////////////
    //double Ttr=0.18/0.1973;
    //double dummy=temperature/Ttr;
    //double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    //double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;
    
    //if(temperature<0.99945*Ttr)
    //{
    //    bulk = lambda3*exp((dummy-1)/sigma3)+ lambda4*exp((dummy-1)/sigma4)+0.03;
    //}
    //if(temperature>0.99945*Ttr)
    //{
    //    bulk = 0.901*exp(14.5*(1.0-dummy)) + 0.061/dummy/dummy;
    //}
    //bulk = bulk*s_den;

    // defining bulk relaxation time and additional transport coefficients
    // Bulk relaxation time from kinetic theory
    Bulk_Relax_time    = 1./14.55/(1./3.-cs2)/(1./3.-cs2)/(grid_pt->epsilon + grid_pt->p)*bulk; 

    transport_coeff1   = 2.0/3.0*(Bulk_Relax_time);          // from kinetic theory, small mass limit
    transport_coeff2   = 0.;                                 // not known; put 0
    transport_coeff1_s = 8./5.*(1./3.-cs2)*Bulk_Relax_time;  // from kinetic theory
    transport_coeff2_s = 0.;                                 // not known;  put 0

    // Computing Navier-Stokes term (-bulk viscosity * theta)
    NS_term = -bulk*(grid_pt->theta_u[rk_flag]);

    // Computing relaxation term and nonlinear term: - Bulk - transport_coeff1*Bulk*theta
    tempf = -(grid_pt->pi_b[rk_flag]) - transport_coeff1*(grid_pt->theta_u[rk_flag])*(grid_pt->pi_b[rk_flag]);

    // Computing nonlinear term: + transport_coeff2*Bulk*Bulk
    if (include_BBterm == 1)
        BB_term = transport_coeff2*(grid_pt->pi_b[rk_flag])*(grid_pt->pi_b[rk_flag]);
    else
        BB_term = 0.0;

    // Computing terms that Couple with shear-stress tensor
    double Wsigma, WW, Shear_Sigma_term, Shear_Shear_term, Coupling_to_Shear;

    if (include_coupling_to_shear == 1)
    {
        // Computing sigma^mu^nu
        double sigma[4][4];
        for(int a=0;a<4;a++ )
        {
            for(int b=0;b<4;b++ )
            {
                sigma[a][b] = grid_pt->sigma[rk_flag][a][b];

                //sigma[a][b] = ( grid_pt->dUsup[rk_flag][a][b] + grid_pt->dUsup[rk_flag][b][a] )/2.
                //            - ( DATA->gmunu[a][b] + (grid_pt->u[rk_flag][a])*(grid_pt->u[rk_flag][b]) )*(grid_pt->theta_u[rk_flag])/3.
                //            + gamma/tau*(DATA->gmunu[a][3])*(DATA->gmunu[b][3])
                //            - ueta/tau/2.*( (DATA->gmunu[a][3])*(DATA->gmunu[b][0]) + (DATA->gmunu[b][3])*(DATA->gmunu[a][0]) )
                //            + ueta*gamma/tau/2.*( (DATA->gmunu[a][3])*(grid_pt->u[rk_flag][b]) +( DATA->gmunu[b][3])*(grid_pt->u[rk_flag][a]) )
                //            - ueta*ueta/tau/2.*( (DATA->gmunu[a][0])*(grid_pt->u[rk_flag][b]) + (DATA->gmunu[b][0])*(grid_pt->u[rk_flag][a]) )
                //            + ( grid_pt->u[rk_flag][a]*grid_pt->a[rk_flag][b] + grid_pt->u[rk_flag][b]*grid_pt->a[rk_flag][a] )/2. ;

            }
        }

        Wsigma = (
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

        WW = (
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

        // multiply term by its respective transport coefficient
        Shear_Sigma_term = Wsigma*transport_coeff1_s;
        Shear_Shear_term = WW*transport_coeff2_s;

        // full term that couples to shear is
        Coupling_to_Shear = -Shear_Sigma_term + Shear_Shear_term ;
    }
    else
        Coupling_to_Shear = 0.0;
        
    // Final Answer
    Final_Answer = NS_term + tempf + BB_term + Coupling_to_Shear;

    return Final_Answer/(Bulk_Relax_time);
}/* Make_uPiSource */


/* Sangyong Nov 18 2014 */
/* baryon current parts */
/* this contains the source terms
   that is, all the terms that are not part of the current */
/* for the q part, we don't do tau*u*q we just do u*q 
this part contains 
-(1/tau_rho)(q[a] + kappa g[a][b]Dtildemu[b] + kappa u[a] u[b]g[b][c]Dtildemu[c])
+Delta[a][tau] u[eta] q[eta]/tau
-Delta[a][eta] u[eta] q[tau]/tau
-u[a]u[b]g[b][e] Dq[e]
*/
double Diss::Make_uqSource(double tau, Grid *grid_pt, int nu, InitData *DATA, int rk_flag)
{
    double tempf, tau_rho, tau_pi, shear, shear_to_s;
    double SW, kappa, T, epsilon, rhob, Ttr;
    int i;
    double q[4];
  
    if(DATA->turn_on_diff == 0) return 0.0;
 
    // Useful variables to define
    Ttr = 0.18/hbarc;  /// phase transition temperature
    epsilon = grid_pt->epsilon;
    rhob = grid_pt->rhob;

    T=eos->get_temperature(epsilon,rhob);
    if(DATA->T_dependent_shear_to_s == 1)
    {
      if(T < Ttr)
        shear_to_s=0.681-0.0594*T/Ttr-0.544*(T/Ttr)*(T/Ttr);
      else
        shear_to_s=-0.289+0.288*T/Ttr+0.0818*(T/Ttr)*(T/Ttr);
    }
    else
      shear_to_s = DATA->shear_to_s;

    //double s_den = eos->get_entropy(epsilon, rhob);
    shear = (shear_to_s)*(epsilon + grid_pt->p)/(T + 1e-15);
    tau_pi = 5.0*shear/(epsilon + grid_pt->p + 1e-15);
 
    if(!isfinite(tau_pi))
    {
        tau_pi = DATA->delta_tau; 
        cout << "tau_pi was infinite ..." << endl;
    }
    if(tau_pi < DATA->delta_tau)  // avoid tau_pi to be too small
        tau_pi = DATA->delta_tau;

    // Sangyong Nov 18 2014: From Gabriel
    // tau_rho = 27/20 * tau_shear
    // D = 9/64 * eta/T
    //tau_rho = (27.0/20.0)*tau_pi;
    //kappa = (9.0/64.0)*shear/T;
    tau_rho = 0.2/(T + 1e-15);
    double mub = eos->get_mu(epsilon, rhob);
    kappa = 0.2*rhob/(mub + 1e-15);

    // copy the value of \tilde{q^\mu}
    for(i=0; i<4; i++)
      q[i] = (grid_pt->Wmunu[rk_flag][4][i]);

    /*
      -(1/tau_rho)(q[a] + kappa g[a][b]Dtildemu[b] + kappa u[a] u[b]g[b][c]Dtildemu[c])
      + theta q[a] - q[a] u^\tau/tau
      +Delta[a][tau] u[eta] q[eta]/tau
      -Delta[a][eta] u[eta] q[tau]/tau
      -u[a] u[b]g[b][e] Dq[e] -> u[a] q[e] g[e][b] Du[b]
    */    
    SW = 0;  // record the final result
 
    // first: (1/tau_rho) part
    // recall that dUsup[4][i] = partial_i (muB/T) 
    // and dUsup[4][0] = -partial_tau (muB/T) = partial^tau (muB/T)
    // and a[4] = u^a partial_a (muB/T) = DmuB/T
    // -(1/tau_rho)(q[a] + kappa g[a][b]DmuB/T[b] + kappa u[a] u[b]g[b][c]DmuB/T[c])
    // a = nu 
    double NS = kappa*(grid_pt->dUsup[rk_flag][4][nu] 
                           + grid_pt->u[rk_flag][nu]*grid_pt->a[rk_flag][4]);
    if(isnan(NS))
    {
        cout << "Navier Stock term is nan! " << endl;
        cout << q[nu] << endl;
        cout << grid_pt->dUsup[rk_flag][4][nu] << endl; // derivative already upper index
        cout << grid_pt->a[rk_flag][4] << endl;
        cout << tau_rho << endl;
        cout << kappa << endl;
        cout << grid_pt->u[rk_flag][nu] << endl;
    }
  
    // add a new non-linear term (- q \theta)
    double transport_coeff = 1.0*tau_rho;   // from conformal kinetic theory
    double Nonlinear1 = -transport_coeff*(q[nu]*grid_pt->theta_u[rk_flag]);

    // add a new non-linear term (-q^\mu \sigma_\mu\nu)
    double transport_coeff_2 = 3./5.*tau_rho;   // from 14-momentum massless
    double temptemp = 0.0;
    for(int i = 0 ; i < 4; i++)
    {
        temptemp += q[i]*grid_pt->sigma[rk_flag][i][nu]*DATA->gmunu[i][i]; 
    }
    double Nonlinear2 = - transport_coeff_2*temptemp;

    SW = (-q[nu] - NS + Nonlinear1 + Nonlinear2)/(tau_rho + 1e-15);

    // all other geometric terms....

    // + theta q[a] - q[a] u^\tau/tau
    SW += (grid_pt->theta_u[rk_flag] - grid_pt->u[rk_flag][0]/tau)*q[nu];
 
    if(isnan(SW))
    {
        cout << "theta term is nan! " << endl;
    }

    // +Delta[a][tau] u[eta] q[eta]/tau 
    tempf = ((DATA->gmunu[nu][0] + grid_pt->u[rk_flag][nu]*grid_pt->u[rk_flag][0])
             *grid_pt->u[rk_flag][3]*q[3]/tau
             - (DATA->gmunu[nu][3] + grid_pt->u[rk_flag][nu]*grid_pt->u[rk_flag][3])
               *grid_pt->u[rk_flag][3]*q[0]/tau);
    SW += tempf;
 
    if(isnan(tempf))
    {
        cout << "Delta^{a \tau} and Delta^{a \eta} terms are nan!" << endl;
    }

    //-u[a] u[b]g[b][e] Dq[e] -> u[a] (q[e] g[e][b] Du[b])
    tempf = 0.0;
    for(i=0; i<4; i++)
    {
      tempf += q[i]*gmn(i)*(grid_pt->a[rk_flag][i]);
    }
    SW += (grid_pt->u[rk_flag][nu])*tempf;
    
    if(isnan(tempf))
    {
        cout << "u^a q_b Du^b term is nan! " << endl;
    }

    return SW;
}/* Make_uqSource */


int Diss::Make_uqRHS(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			 Grid *Lneighbor2, Grid *Rneighbor2, double **w_rhs, InitData *DATA, int rk_flag, int size, int rank)
{
  int mu, nu, direc, nmax[4];
  double f, fp1, fm1, fp2, fm2, delta[4];
  double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
  double uWphR, uWphL, uWmhR, uWmhL, WphR, WphL, WmhR, WmhL;
  double HWph, HWmh, taufactor, HW;
  double sum;

  /* Kurganov-Tadmor for q */
  /* implement 
    partial_tau (utau qmu) + (1/tau)partial_eta (ueta qmu) 
    + partial_x (ux qmu) + partial_y (uy qmu) + utau qmu/tau = SW 
  or the right hand side of,
    partial_tau (utau qmu) = 
    - (1/tau)partial_eta (ueta qmu) - partial_x (ux qmu) - partial_y (uy qmu) 
    - utau qmu/tau 
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


  // we use the Wmunu[4][nu] = q[nu] 
  mu = 4;
  for(nu=0; nu<4; nu++)
  {
    sum = 0.0;
    for(direc=1; direc<=3; direc++)
    {
      if(direc==3)
          taufactor = tau;
      else 
          taufactor = 1.0;

      /* Get_uWmns */
      Get_uWmns(tau, grid_pt, Lneighbor, Rneighbor, Lneighbor2, Rneighbor2, 
                mu, nu, direc, &g, &f, &gp1, &fp1, &gp2, &fp2, &gm1, &fm1, 
                &gm2, &fm2, DATA, rk_flag, size, rank);
 
      /*  MakeuWmnHalfs */
      /* uWmn */
      uWphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f); 
      uWphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
      uWmhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
      uWmhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);

      /* just Wmn */
      WphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g); 
      WphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
      WmhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
      WmhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);

      a = fabs(grid_pt->u[rk_flag][direc])/grid_pt->u[rk_flag][0];

      if(direc<3) // x,y direction
	{
	  if(grid_pt->position[direc] == 0)
	      am1 = a;
	  else
	  {
	      am1 = fabs(grid_pt->nbr_m_1[direc]->u[rk_flag][direc])
                  /grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	  }
	  if(grid_pt->position[direc] == nmax[direc])
	  {
	      ap1 = a;
	  }
	  else
	  {
	      ap1 = fabs(grid_pt->nbr_p_1[direc]->u[rk_flag][direc])
                  /grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	  }
      }
      else if (direc==3)
	{
	  if(grid_pt->position[direc] == 0)
	  {
	      if (rank==0)
	        am1 = a;
	      else
	      {
	        am1 = fabs(Lneighbor->u[rk_flag][direc])/Lneighbor->u[rk_flag][0];
	      }
	  }
	  else
	  {
	      am1 = fabs(grid_pt->nbr_m_1[direc]->u[rk_flag][direc])
                  /grid_pt->nbr_m_1[direc]->u[rk_flag][0];
	  }
	  
	  if(grid_pt->position[direc] == nmax[direc])
	  {
	      if (rank==size-1)
	        ap1 = a;
	      else
	      {
	        ap1 = fabs(Rneighbor->u[rk_flag][direc])/Rneighbor->u[rk_flag][0];
	      }		 
	  }
	  else
	  {
	      ap1 = fabs(grid_pt->nbr_p_1[direc]->u[rk_flag][direc])
                  /grid_pt->nbr_p_1[direc]->u[rk_flag][0];
	  }
	}
      ax = maxi(a, ap1);
      HWph = ((uWphR + uWphL) - ax*(WphR - WphL))*0.5;

      ax = maxi(a, am1); 
      HWmh = ((uWmhR + uWmhL) - ax*(WmhR - WmhL))*0.5;

      HW = (HWph - HWmh)/delta[direc]/taufactor;

      /* make partial_i (u^i Wmn) */
      sum += -HW;

    }/* direction */
    
    /* add a source term -u^tau Wmn/tau due to the coordinate change to tau-eta */
    /* Sangyong Nov 18 2014: don't need this. included in the uqSource. */
       /* this is from udW = d(uW) - Wdu = RHS */
       /* or d(uW) = udW + Wdu */
    /*    
        sum -= (grid_pt->u[rk_flag][0])*(grid_pt->Wmunu[rk_flag][mu][nu])/tau;
        sum += (grid_pt->theta_u[rk_flag])*(grid_pt->Wmunu[rk_flag][mu][nu]);
    */  
    w_rhs[mu][nu] = sum*(DATA->delta_tau);
  }/* nu */
  return 1; /* if successful */
}/* Make_uqRHS */

