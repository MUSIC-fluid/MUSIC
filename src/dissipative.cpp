// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <iomanip>
#include "./util.h"
#include "./cell.h"
#include "./grid.h"
#include "./data.h"
#include "./eos.h"
#include "./dissipative.h"

using namespace std;

Diss::Diss(EOS *eosIn, InitData* DATA_in) {
    eos = eosIn;
    minmod = new Minmod(DATA_in);
}

// destructor
Diss::~Diss() {
    delete minmod;
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* Dissipative parts */
/* Sangyong Nov 18 2014 */
/* change: alpha first which is the case
for everywhere else. also, this change is necessary
to use Wmunu[rk_flag][4][mu] as the dissipative baryon current*/
/* this is the only one that is being subtracted in the rhs */
double Diss::MakeWSource(double tau, int alpha, Grid &arena, int ix, int iy, int ieta,
                         InitData *DATA, int rk_flag) {
    const auto& grid_pt = arena(ix, iy, ieta);
    double delta[4], taufactor;
    double shear_on = DATA->turn_on_shear;
    double bulk_on  = DATA->turn_on_bulk;
    double diff_on  = DATA->turn_on_diff;
    int i;

    /* calculate d_m (tau W^{m,alpha}) + (geom source terms) */
    delta[1] = DATA->delta_x;
    delta[2] = DATA->delta_y;
    delta[3] = DATA->delta_eta;

    /* partial_tau W^tau alpha */
    /* this is partial_tau evaluated at tau */
    /* this is the first step. so rk_flag = 0 */
    if (alpha == 4 && DATA->turn_on_diff == 0)
        return (0.0);

    /* Sangyong Nov 18 2014 */
    /* change: alpha first which is the case
               for everywhere else. also, this change is necessary
               to use Wmunu[rk_flag][4][mu] as the dissipative baryon current*/
    // dW/dtau
    // backward time derivative (first order is more stable)
    int idx_1d_alpha0 = Util::map_2d_idx_to_1d(alpha, 0);
    double dWdtau;
    if (rk_flag == 0) {
        dWdtau = (grid_pt.Wmunu[rk_flag][idx_1d_alpha0] - grid_pt.prevWmunu[0][idx_1d_alpha0])/DATA->delta_tau;
    } else {
        dWdtau = (grid_pt.Wmunu[rk_flag][idx_1d_alpha0] - grid_pt.Wmunu    [0][idx_1d_alpha0])/DATA->delta_tau;
    }

    /* bulk pressure term */
    double dPidtau = 0.0;
    double Pi_alpha0 = 0.0;
    if (alpha < 4 && DATA->turn_on_bulk == 1) {
        double gfac = (alpha == 0 ? -1.0 : 0.0);
        Pi_alpha0 = (
            grid_pt.pi_b[rk_flag]*(gfac + grid_pt.u[rk_flag][alpha]
                                           *grid_pt.u[rk_flag][0]));
        if (rk_flag == 0) {
            dPidtau = ((Pi_alpha0 - grid_pt.prev_pi_b[0]
                                    *(gfac + grid_pt.prev_u[0][alpha]
                                             *grid_pt.prev_u[0][0]))
                       /DATA->delta_tau);
        } else {
            dPidtau = ((Pi_alpha0 - grid_pt.pi_b[0]
                                    *(gfac + grid_pt.u[0][alpha]
                                             *grid_pt.u[0][0]))
                       /DATA->delta_tau);
        }
    }

    // use central difference to preserve conservation law exactly
    int idx_1d;
    double dWdx_perp  = 0.0;
    double dPidx_perp = 0.0;

    // x-direction
    idx_1d          = Util::map_2d_idx_to_1d(alpha, 1);
    auto grid_pt_p1 = &arena(ix + 1, iy, ieta);
    auto grid_pt_m1 = &arena(ix - 1, iy, ieta);
    double sg       = grid_pt.Wmunu[rk_flag][idx_1d];
    double sgp1     = grid_pt_p1->Wmunu[rk_flag][idx_1d];
    double sgm1     = grid_pt_m1->Wmunu[rk_flag][idx_1d];
    dWdx_perp += minmod->minmod_dx(sgp1, sg, sgm1)/delta[1];
    //dWdx_perp += (sgp1 - sgm1)/(2.*delta[i]);
    if (alpha < 4 && DATA->turn_on_bulk == 1) {
        double gfac1 = (alpha == 1 ? 1.0 : 0.0);
        double bgp1  = grid_pt_p1->pi_b[rk_flag]*(gfac1 + grid_pt_p1->u[rk_flag][alpha]*grid_pt_p1->u[rk_flag][1]);
        double bg    = grid_pt    .pi_b[rk_flag]*(gfac1 + grid_pt    .u[rk_flag][alpha]*grid_pt    .u[rk_flag][1]);
        double bgm1  = grid_pt_m1->pi_b[rk_flag]*(gfac1 + grid_pt_m1->u[rk_flag][alpha]*grid_pt_m1->u[rk_flag][1]);
        dPidx_perp  += minmod->minmod_dx(bgp1, bg, bgm1)/delta[1];
        //dPidx_perp += (bgp1 - bgm1)/(2.*delta[i]);
    }

    // y-direction
    idx_1d     = Util::map_2d_idx_to_1d(alpha, 2);
    grid_pt_p1 = &arena(ix, iy + 1, ieta);
    grid_pt_m1 = &arena(ix, iy - 1, ieta);
    sg         = grid_pt.Wmunu[rk_flag][idx_1d];
    sgp1       = grid_pt_p1->Wmunu[rk_flag][idx_1d];
    sgm1       = grid_pt_m1->Wmunu[rk_flag][idx_1d];
    dWdx_perp += minmod->minmod_dx(sgp1, sg, sgm1)/delta[2];
    //dWdx_perp += (sgp1 - sgm1)/(2.*delta[i]);
    if (alpha < 4 && DATA->turn_on_bulk == 1) {
        double gfac1 = (alpha == 2 ? 1.0 : 0.0);
        double bgp1  = grid_pt_p1->pi_b[rk_flag]*(gfac1 + grid_pt_p1->u[rk_flag][alpha]*grid_pt_p1->u[rk_flag][2]);
        double bg    = grid_pt    .pi_b[rk_flag]*(gfac1 + grid_pt    .u[rk_flag][alpha]*grid_pt    .u[rk_flag][2]);
        double bgm1  = grid_pt_m1->pi_b[rk_flag]*(gfac1 + grid_pt_m1->u[rk_flag][alpha]*grid_pt_m1->u[rk_flag][2]);
        dPidx_perp  += minmod->minmod_dx(bgp1, bg, bgm1)/delta[2];
        //dPidx_perp += (bgp1 - bgm1)/(2.*delta[i]);
    }

    // eta-direction
    i              = 3;
    taufactor      = tau;
    double dWdeta  = 0.0;
    double dPideta = 0.0;
    idx_1d         = Util::map_2d_idx_to_1d(alpha, i);
    grid_pt_p1     = &arena(ix, iy, ieta + 1);
    grid_pt_m1     = &arena(ix, iy, ieta - 1);
    sg             = grid_pt.Wmunu[rk_flag][idx_1d];
    sgp1           = grid_pt_p1->Wmunu[rk_flag][idx_1d];
    sgm1           = grid_pt_m1->Wmunu[rk_flag][idx_1d];
    dWdeta         = minmod->minmod_dx(sgp1, sg, sgm1)/delta[i]/taufactor;
    //dWdeta = (sgp1 - sgm1)/(2.*delta[i]*taufactor);
    if (alpha < 4 && DATA->turn_on_bulk == 1) {
        double gfac3 = (alpha == i ? 1.0 : 0.0);
        double bgp1  = grid_pt_p1->pi_b[rk_flag]*(gfac3 + grid_pt_p1->u[rk_flag][alpha]*grid_pt_p1->u[rk_flag][i]);
        double bg    = grid_pt    .pi_b[rk_flag]*(gfac3 + grid_pt    .u[rk_flag][alpha]*grid_pt    .u[rk_flag][i]);
        double bgm1  = grid_pt_m1->pi_b[rk_flag]*(gfac3 + grid_pt_m1->u[rk_flag][alpha]*grid_pt_m1->u[rk_flag][i]);
        dPideta      = minmod->minmod_dx(bgp1, bg, bgm1)/delta[i]/taufactor;
        //dPideta    = (bgp1 - bgm1)/(2.*delta[i]*taufactor);
    }

    /* partial_m (tau W^mn) = W^0n + tau partial_m W^mn */
    double sf = (tau*(dWdtau + dWdx_perp + dWdeta)
                 + grid_pt.Wmunu[rk_flag][idx_1d_alpha0]);
    double bf = (tau*(dPidtau + dPidx_perp + dPideta)
                 + Pi_alpha0);

    /* sources due to coordinate transform this is added to partial_m W^mn */
    if (alpha == 0) {
        sf += grid_pt.Wmunu[rk_flag][9];
        bf += grid_pt.pi_b[rk_flag]*(1.0 + grid_pt.u[rk_flag][3]
                                            *grid_pt.u[rk_flag][3]);
    }
    if (alpha == 3) {
        sf += grid_pt.Wmunu[rk_flag][3];
        bf += grid_pt.pi_b[rk_flag]*(grid_pt.u[rk_flag][0]
                                      *grid_pt.u[rk_flag][3]);
    }

    // final result
    double result = 0.0;
    if (alpha < 4)
        result = (sf*shear_on + bf*bulk_on);
    else if (alpha == 4)
        result = sf*diff_on;

    if (isnan(result)) {
        cout << "[Error]Diss::MakeWSource: " << endl;
        cout << "rk_flag = " << rk_flag << endl;
        cout << "sf=" << sf << " bf=" << bf
             << " Wmunu[" << rk_flag << "]="
             << grid_pt.Wmunu[rk_flag][alpha]
             << " pi_b[" << rk_flag << "]="
             << grid_pt.pi_b[rk_flag]
             << " prev_pi_b=" << grid_pt.prev_pi_b[rk_flag] << endl;
    }
    return(result);
}/* MakeWSource */
/* MakeWSource is for Tmunu */

double Diss::Make_uWSource(double tau, Cell *grid_pt, int mu, int nu,
                           InitData *DATA, int rk_flag, double theta_local,
                           double *a_local, double *sigma_1d) {
    if (DATA->turn_on_shear == 0)
        return 0.0;

    double tempf, tau_pi;
    double SW, shear, shear_to_s, T, epsilon, rhob;
    int a, b;
    double gamma, ueta;
    double NS_term;
    double sigma[4][4];
    double Wmunu[4][4];
    for (int a = 0; a < 4; a++) {
        for (int b = a; b < 4; b++) {
            int idx_1d = Util::map_2d_idx_to_1d(a, b);
            Wmunu[a][b] = grid_pt->Wmunu[rk_flag][idx_1d];
            sigma[a][b] = sigma_1d[idx_1d];
        }
    }
    for (int a = 0; a < 4; a++) {
        for (int b = a+1; b < 4; b++) {
            Wmunu[b][a] = Wmunu[a][b];
            sigma[b][a] = sigma[a][b];
        }
    }

    // Useful variables to define
    gamma = grid_pt->u[rk_flag][0];
    ueta  = grid_pt->u[rk_flag][3];
    if (rk_flag == 0) {
        epsilon = grid_pt->epsilon;
        rhob = grid_pt->rhob;
    } else {
        epsilon = grid_pt->prev_epsilon;
        rhob = grid_pt->prev_rhob;
    }
    T = eos->get_temperature(epsilon, rhob);

    if (DATA->T_dependent_shear_to_s == 1) {
        shear_to_s = get_temperature_dependent_eta_s(DATA, T);
    } else {
        shear_to_s = DATA->shear_to_s;
    }

    int include_WWterm         = 0;
    int include_Vorticity_term = 0;
    int include_Wsigma_term    = 0;
    if (DATA->include_second_order_terms == 1) {
        include_WWterm      = 1;
        include_Wsigma_term = 1;
    }
    if (DATA->Initial_profile == 0) {
        include_WWterm         = 0;
        include_Wsigma_term    = 0;
        include_Vorticity_term = 0;
    }

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///                 Defining transport coefficients                        ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
    double pressure = eos->get_pressure(epsilon, rhob);
    shear = (shear_to_s)*(epsilon + pressure)/(T + 1e-15);
    tau_pi = 5.0*shear/(epsilon + pressure + 1e-15);

    // tau_pi = maxi(tau_pi, DATA->tau_pi);
    if (tau_pi > 1e20) {
        tau_pi = DATA->delta_tau;
        cout << "tau_pi was infinite ..." << endl;
    }
    if (tau_pi < DATA->delta_tau)
        tau_pi = DATA->delta_tau;

    // transport coefficient for nonlinear terms -- shear only terms -- 4Mar2013
    // transport coefficients of a massless gas of single component particles
    double transport_coefficient  = 9./70.*tau_pi/shear*(4./5.);
    double transport_coefficient2 = 4./3.*tau_pi;
    double transport_coefficient3 = 10./7.*tau_pi;
    double transport_coefficient4 = 2.*tau_pi;

    // transport coefficient for nonlinear terms
    // -- coupling to bulk viscous pressure -- 4Mar2013
    // transport coefficients not yet known -- fixed to zero
    double transport_coefficient_b  = 6./5.*tau_pi;
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
    tempf = (-(1.0 + transport_coefficient2*theta_local)*(Wmunu[mu][nu]));

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///             Navier-Stokes Term -- -2.*shear*sigma^munu                 ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///

    // full Navier-Stokes term is
    // sign changes according to metric sign convention
    NS_term = - 2.*shear*sigma[mu][nu];

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///                             Vorticity Term                             ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
    double term1_Vorticity;
    double Vorticity_term;

    // remember: dUsup[m][n] = partial^n u^m  ///
    // remember:  a[n]  =  u^m*partial_m u^n  ///
    if (include_Vorticity_term == 1) {
        double omega[4][4];
        for (a = 0; a < 4; a++) {
            for (b = 0; b <4; b++) {
                omega[a][b] = (
                    (grid_pt->dUsup[a][b]
                     - grid_pt->dUsup[b][a])/2.
                    + ueta/tau/2.*(  DATA->gmunu[a][0]*DATA->gmunu[b][3]
                                   - DATA->gmunu[b][0]*DATA->gmunu[a][3])
                    - ueta*gamma/tau/2.
                      *(  DATA->gmunu[a][3]*grid_pt->u[rk_flag][b]
                        - DATA->gmunu[b][3]*grid_pt->u[rk_flag][a])
                    + ueta*ueta/tau/2.
                      *(   DATA->gmunu[a][0]*grid_pt->u[rk_flag][b]
                         - DATA->gmunu[b][0]*grid_pt->u[rk_flag][a])
                    + (  grid_pt->u[rk_flag][a]*a_local[b]
                       - grid_pt->u[rk_flag][b]*a_local[a])/2.);
            }
        }
        term1_Vorticity = (- Wmunu[mu][0]*omega[nu][0]
                           - Wmunu[nu][0]*omega[mu][0]
                           + Wmunu[mu][1]*omega[nu][1]
                           + Wmunu[nu][1]*omega[mu][1]
                           + Wmunu[mu][2]*omega[nu][2]
                           + Wmunu[nu][2]*omega[mu][2]
                           + Wmunu[mu][3]*omega[nu][3]
                           + Wmunu[nu][3]*omega[mu][3])/2.;
        // multiply term by its respective transport coefficient
        term1_Vorticity = transport_coefficient4*term1_Vorticity;
        // full term is
        Vorticity_term = term1_Vorticity;
    } else {
        Vorticity_term = 0.0;
    }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                  Add nonlinear term in shear-stress tensor                ///
//  transport_coefficient3*Delta(mu nu)(alpha beta)*Wmu gamma sigma nu gamma ///
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
    double Wsigma, Wsigma_term;
    double term1_Wsigma, term2_Wsigma;
    if (include_Wsigma_term == 1) {
        Wsigma = (
               Wmunu[0][0]*sigma[0][0]
             + Wmunu[1][1]*sigma[1][1]
             + Wmunu[2][2]*sigma[2][2]
             + Wmunu[3][3]*sigma[3][3]
             - 2.*(  Wmunu[0][1]*sigma[0][1]
                   + Wmunu[0][2]*sigma[0][2]
                   + Wmunu[0][3]*sigma[0][3])
             +2.*(  Wmunu[1][2]*sigma[1][2]
                  + Wmunu[1][3]*sigma[1][3]
                  + Wmunu[2][3]*sigma[2][3]));
        term1_Wsigma = ( - Wmunu[mu][0]*sigma[nu][0]
                         - Wmunu[nu][0]*sigma[mu][0]
                         + Wmunu[mu][1]*sigma[nu][1]
                         + Wmunu[nu][1]*sigma[mu][1]
                         + Wmunu[mu][2]*sigma[nu][2]
                         + Wmunu[nu][2]*sigma[mu][2]
                         + Wmunu[mu][3]*sigma[nu][3]
                         + Wmunu[nu][3]*sigma[mu][3])/2.;

        term2_Wsigma = (-(1./3.)*(DATA->gmunu[mu][nu]
                                  + grid_pt->u[rk_flag][mu]
                                    *grid_pt->u[rk_flag][nu])*Wsigma);
        // multiply term by its respective transport coefficient
        term1_Wsigma = transport_coefficient3*term1_Wsigma;
        term2_Wsigma = transport_coefficient3*term2_Wsigma;

        // full term is
        Wsigma_term = -term1_Wsigma - term2_Wsigma;
    } else {
        Wsigma_term = 0.0;
    }

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///               Add nonlinear term in shear-stress tensor                ///
///   transport_coefficient*Delta(mu nu)(alpha beta)*Wmu gamma Wnu gamma   ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
    double Wsquare, WW_term;
    double term1_WW, term2_WW;
    if (include_WWterm == 1) {
        Wsquare = (  Wmunu[0][0]*Wmunu[0][0]
                   + Wmunu[1][1]*Wmunu[1][1]
                   + Wmunu[2][2]*Wmunu[2][2]
                   + Wmunu[3][3]*Wmunu[3][3]
            - 2.*(  Wmunu[0][1]*Wmunu[0][1]
                  + Wmunu[0][2]*Wmunu[0][2]
                  + Wmunu[0][3]*Wmunu[0][3])
            + 2.*(  Wmunu[1][2]*Wmunu[1][2]
                  + Wmunu[1][3]*Wmunu[1][3]
                  + Wmunu[2][3]*Wmunu[2][3]));
        term1_WW = ( - Wmunu[mu][0]*Wmunu[nu][0]
                     + Wmunu[mu][1]*Wmunu[nu][1]
                     + Wmunu[mu][2]*Wmunu[nu][2]
                     + Wmunu[mu][3]*Wmunu[nu][3]);

        term2_WW = (
            -(1./3.)*(DATA->gmunu[mu][nu]
                      + grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][nu])
            *Wsquare);

        // multiply term by its respective transport coefficient
        term1_WW = term1_WW*transport_coefficient;
        term2_WW = term2_WW*transport_coefficient;

        // full term is
        // sign changes according to metric sign convention
        WW_term = -term1_WW - term2_WW;
    } else {
        WW_term = 0.0;
    }

/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
///               Add coupling to bulk viscous pressure                    ///
///              transport_coefficient_b*Bulk*sigma^mu nu                  ///
///               transport_coefficient2_b*Bulk*W^mu nu                    ///
/// ////////////////////////////////////////////////////////////////////// ///
/// ////////////////////////////////////////////////////////////////////// ///
    double Coupling_to_Bulk;
    if (DATA->include_second_order_terms == 1) {
        double Bulk_Sigma = grid_pt->pi_b[rk_flag]*sigma[mu][nu];
        double Bulk_W = grid_pt->pi_b[rk_flag]*Wmunu[mu][nu];

        // multiply term by its respective transport coefficient
        double Bulk_Sigma_term = Bulk_Sigma*transport_coefficient_b;
        double Bulk_W_term = Bulk_W*transport_coefficient2_b;

        // full term is
        // first term: sign changes according to metric sign convention
        Coupling_to_Bulk = -Bulk_Sigma_term + Bulk_W_term;
    } else {
        Coupling_to_Bulk = 0.0;
    }

    // final answer is
    SW = (NS_term + tempf + Vorticity_term + Wsigma_term + WW_term
          + Coupling_to_Bulk)/(tau_pi);
    return(SW);
}/* Make_uWSource */


int Diss::Make_uWRHS(double tau, Grid &arena, int ix, int iy, int ieta,
                     double **w_rhs, InitData *DATA, int rk_flag,
                     double theta_local, double *a_local) {
    auto grid_pt = &(arena(ix, iy, ieta));
    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            w_rhs[a][b] = 0.0;
        }
    }
    if (DATA->turn_on_shear == 0)
        return(1);

    int mu, nu, direc, ic;
    double f, fp1, fm1, fp2, fm2;
    double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
    double uWphR, uWphL, uWmhR, uWmhL, WphR, WphL, WmhR, WmhL;
    double HWph, HWmh, taufactor, HW, ic_fac;
    double tempf;
    double Wmunu_local[4][4];
    for (int aa = 0; aa < 4; aa++) {
        for (int bb = aa; bb < 4; bb++) {
            int idx_1d = Util::map_2d_idx_to_1d(aa, bb);
            Wmunu_local[aa][bb] = grid_pt->Wmunu[rk_flag][idx_1d];
        }
    }
    for (int aa = 0; aa < 4; aa++) {
        for (int bb = aa+1; bb < 4; bb++) {
            Wmunu_local[bb][aa] = Wmunu_local[aa][bb];
        }
    }

    /* Kurganov-Tadmor for Wmunu */
    /* implement 
       partial_tau (utau Wmn) + (1/tau)partial_eta (ueta Wmn) 
       + partial_x (ux Wmn) + partial_y (uy Wmn) + utau Wmn/tau = SW 
       or the right hand side of,
       partial_tau (utau Wmn) = 
                        - (1/tau)partial_eta (ueta Wmn)
                        - partial_x (ux Wmn) - partial_y (uy Wmn) 
                        - utau Wmn/tau + SW*/

    /* the local velocity is just u_x/u_tau, u_y/u_tau, u_eta/tau/u_tau */
    /* KT flux is given by 
       H_{j+1/2} = (fRph + fLph)/2 - ax(uRph - uLph) 
       Here fRph = ux WmnRph and ax uRph = |ux/utau|_max utau Wmn */
    /* This is the second step in the operator splitting. it uses
       rk_flag+1 as initial condition */
    double delta[4];
    delta[1] = DATA->delta_x;
    delta[2] = DATA->delta_y;
    delta[3] = DATA->delta_eta;
    
    Cell *grid_pt_p1, *grid_pt_p2, *grid_pt_m1, *grid_pt_m2;
    double sum;
    // pi^\mu\nu is symmetric
    for (mu = 1; mu < 4; mu++) {
        for (nu = mu; nu < 4; nu++) {
            int idx_1d = Util::map_2d_idx_to_1d(mu, nu);
            sum = 0.0;
            for (direc = 1; direc <= 3; direc++) {
                if (direc == 1) {
                    grid_pt_p1 = &(arena(ix + 1, iy, ieta));
                    grid_pt_p2 = &(arena(ix + 2, iy, ieta));
                    grid_pt_m1 = &(arena(ix - 1, iy, ieta));
                    grid_pt_m2 = &(arena(ix - 2, iy, ieta));
                } else if (direc == 2) {
                    grid_pt_p1 = &(arena(ix, iy + 1, ieta));
                    grid_pt_p2 = &(arena(ix, iy + 2, ieta));
                    grid_pt_m1 = &(arena(ix, iy - 1, ieta));
                    grid_pt_m2 = &(arena(ix, iy - 2, ieta));
                } else {
                    grid_pt_p1 = &(arena(ix, iy, ieta + 1));
                    grid_pt_p2 = &(arena(ix, iy, ieta + 2));
                    grid_pt_m1 = &(arena(ix, iy, ieta - 1));
                    grid_pt_m2 = &(arena(ix, iy, ieta - 2));
                }

                if (direc == 3) {
                    taufactor = tau;
                } else {
                    taufactor = 1.0;
                }

                /* Get_uWmns */
                g = Wmunu_local[mu][nu];
                f = g*grid_pt->u[rk_flag][direc];
                g *=   grid_pt->u[rk_flag][0];
                   
                gp2 = grid_pt_p2->Wmunu[rk_flag][idx_1d];
                fp2 = gp2*grid_pt_p2->u[rk_flag][direc];
                gp2 *= grid_pt_p2->u[rk_flag][0];
                
                gp1 = grid_pt_p1->Wmunu[rk_flag][idx_1d];
                fp1 = gp1*grid_pt_p1->u[rk_flag][direc];
                gp1 *= grid_pt_p1->u[rk_flag][0];
                
                gm1 = grid_pt_m1->Wmunu[rk_flag][idx_1d];
                fm1 = gm1*grid_pt_m1->u[rk_flag][direc];
                gm1 *= grid_pt_m1->u[rk_flag][0];
                
                gm2 = grid_pt_m2->Wmunu[rk_flag][idx_1d];
                fm2 = gm2*grid_pt_m2->u[rk_flag][direc];
                gm2 *= grid_pt_m2->u[rk_flag][0];

                /* MakeuWmnHalfs */
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

                a = fabs(grid_pt->u[rk_flag][direc]);
                a /= grid_pt->u[rk_flag][0];
                am1 = (fabs(grid_pt_m1->u[rk_flag][direc])
                       /grid_pt_m1->u[rk_flag][0]);
                ap1 = (fabs(grid_pt_p1->u[rk_flag][direc])
                       /grid_pt_p1->u[rk_flag][0]);

                ax = maxi(a, ap1);
                HWph = ((uWphR + uWphL) - ax*(WphR - WphL))*0.5;

                ax = maxi(a, am1);
                HWmh = ((uWmhR + uWmhL) - ax*(WmhR - WmhL))*0.5;

                HW = (HWph - HWmh)/delta[direc]/taufactor;

                /* make partial_i (u^i Wmn) */
                sum += -HW;
            }/* direction */
            /* add a source term -u^tau Wmn/tau
               due to the coordinate change to tau-eta */
            sum += (- (grid_pt->u[rk_flag][0]*Wmunu_local[mu][nu])/tau
                    + (theta_local*Wmunu_local[mu][nu]));

            /* this is from udW = d(uW) - Wdu = RHS */
            /* or d(uW) = udW + Wdu */
            /* this term is being added to the rhs so that -4/3 + 1 = -1/3 */
            /* other source terms due to the coordinate change to tau-eta */
            tempf = 0.0;
            tempf = (
                - (DATA->gmunu[3][mu])*(Wmunu_local[0][nu])
                - (DATA->gmunu[3][nu])*(Wmunu_local[0][mu])
                + (DATA->gmunu[0][mu])*(Wmunu_local[3][nu])
                + (DATA->gmunu[0][nu])*(Wmunu_local[3][mu])
                + (Wmunu_local[3][nu])
                  *(grid_pt->u[rk_flag][mu])*(grid_pt->u[rk_flag][0])
                + (Wmunu_local[3][mu])
                  *(grid_pt->u[rk_flag][nu])*(grid_pt->u[rk_flag][0])
                - (Wmunu_local[0][nu])
                  *(grid_pt->u[rk_flag][mu])*(grid_pt->u[rk_flag][3])
                - (Wmunu_local[0][mu])
                  *(grid_pt->u[rk_flag][nu])*(grid_pt->u[rk_flag][3]))
                  *(grid_pt->u[rk_flag][3]/tau);
            for (ic = 0; ic < 4; ic++) {
                ic_fac = (ic == 0 ? -1.0 : 1.0);
                tempf += (
                      (Wmunu_local[ic][nu])*(grid_pt->u[rk_flag][mu])
                       *(a_local[ic])*ic_fac
                    + (Wmunu_local[ic][mu])*(grid_pt->u[rk_flag][nu])
                       *(a_local[ic])*ic_fac);
            }
            sum += tempf;
            w_rhs[mu][nu] = sum*(DATA->delta_tau);
        }  /* nu */
    }  /* mu */

    // pi^\mu\nu is symmetric
    for (int mu = 1; mu < 4; mu++) {
        for (int nu = mu+1; nu < 4; nu++) {
            w_rhs[nu][mu] = w_rhs[mu][nu];
        }
    }
    return(1);
}/* Make_uWRHS */



int Diss::Make_uPRHS(double tau, Grid &arena, int ix, int iy, int ieta,
                     double *p_rhs, InitData *DATA, 
                     int rk_flag, double theta_local) {
    auto grid_pt = &(arena(ix, iy, ieta));
    int direc;
    double f, fp1, fm1, fp2, fm2, delta[4];
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

    delta[1] = DATA->delta_x;
    delta[2] = DATA->delta_y;
    delta[3] = DATA->delta_eta;

    if (DATA->turn_on_bulk)
        bulk_on = 1.0;
    else 
        bulk_on = 0.0;

    sum = 0.0;
    Cell *grid_pt_p1, *grid_pt_p2, *grid_pt_m1, *grid_pt_m2;
    for (direc=1; direc<=3; direc++) {
        if (direc == 1) {
            grid_pt_p1 = &(arena(ix + 1, iy, ieta));
            grid_pt_p2 = &(arena(ix + 2, iy, ieta));
            grid_pt_m1 = &(arena(ix - 1, iy, ieta));
            grid_pt_m2 = &(arena(ix - 2, iy, ieta));
        } else if (direc == 2) {
            grid_pt_p1 = &(arena(ix, iy + 1, ieta));
            grid_pt_p2 = &(arena(ix, iy + 2, ieta));
            grid_pt_m1 = &(arena(ix, iy - 1, ieta));
            grid_pt_m2 = &(arena(ix, iy - 2, ieta));
        } else {
            grid_pt_p1 = &(arena(ix, iy, ieta + 1));
            grid_pt_p2 = &(arena(ix, iy, ieta + 2));
            grid_pt_m1 = &(arena(ix, iy, ieta - 1));
            grid_pt_m2 = &(arena(ix, iy, ieta - 2));
        }
        if (direc==3) 
            taufactor = tau;
        else 
            taufactor = 1.0;

        /* Get_uPis */
        g = grid_pt->pi_b[rk_flag];
        f = g*grid_pt->u[rk_flag][direc];
        g *= grid_pt->u[rk_flag][0];

        gp2 = grid_pt_p2->pi_b[rk_flag];
        fp2 = gp2*grid_pt_p2->u[rk_flag][direc];
        gp2 *= grid_pt_p2->u[rk_flag][0];
        
        gp1 = grid_pt_p1->pi_b[rk_flag];
        fp1 = gp1*grid_pt_p1->u[rk_flag][direc];
        gp1 *= grid_pt_p1->u[rk_flag][0];
        
        gm1 = grid_pt_m1->pi_b[rk_flag];
        fm1 = gm1*grid_pt_m1->u[rk_flag][direc];
        gm1 *= grid_pt_m1->u[rk_flag][0];
        
        gm2 = grid_pt_m2->pi_b[rk_flag];
        fm2 = gm2*grid_pt_m2->u[rk_flag][direc];
        gm2 *= grid_pt_m2->u[rk_flag][0];

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
  
        am1 = fabs(grid_pt_m1->u[rk_flag][direc]);
        am1 /= grid_pt_m1->u[rk_flag][0];

        ap1 = fabs(grid_pt_p1->u[rk_flag][direc]);
        ap1 /= grid_pt_p1->u[rk_flag][0];
        
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
     sum += (grid_pt->pi_b[rk_flag])*theta_local;
     *p_rhs = sum*(DATA->delta_tau)*bulk_on;

     return 1; /* if successful */
}/* Make_uPRHS */



double Diss::Make_uPiSource(double tau, Cell *grid_pt, InitData *DATA,
                        int rk_flag, double theta_local, double *sigma_1d) {
    if (DATA->turn_on_bulk == 0) return 0.0;

    double tempf;
    double bulk;
    double Bulk_Relax_time;
    double transport_coeff1, transport_coeff2;
    double transport_coeff1_s, transport_coeff2_s;
    double NS_term, BB_term;
    double Final_Answer;

    // switch to include non-linear coupling terms in the bulk pi evolution
    int include_BBterm = 0;
    int include_coupling_to_shear = 0;
    if (DATA->include_second_order_terms == 1) {
        include_BBterm = 1;
        include_coupling_to_shear = 1;
    }
 
    double epsilon, rhob;
    if (rk_flag == 0) {
        epsilon = grid_pt->epsilon;
        rhob = grid_pt->rhob;
    } else {
        epsilon = grid_pt->prev_epsilon;
        rhob = grid_pt->prev_rhob;
    }

    // defining bulk viscosity coefficient

    // shear viscosity = constant * entropy density
    //s_den = eos->get_entropy(epsilon, rhob);
    //shear = (DATA->shear_to_s)*s_den;   
    // shear viscosity = constant * (e + P)/T
    double temperature = eos->get_temperature(epsilon, rhob);

    // cs2 is the velocity of sound squared
    double cs2 = eos->get_cs2(epsilon, rhob);  
    double pressure = eos->get_pressure(epsilon, rhob);

    // T dependent bulk viscosity from Gabriel
    bulk = get_temperature_dependent_zeta_s(temperature);
    bulk = bulk*(epsilon + pressure)/temperature;

    // defining bulk relaxation time and additional transport coefficients
    // Bulk relaxation time from kinetic theory
    Bulk_Relax_time = (1./(14.55*(1./3. - cs2)*(1./3. - cs2))
                       /(epsilon + pressure)*bulk);

    // from kinetic theory, small mass limit
    transport_coeff1   = 2.0/3.0*(Bulk_Relax_time);
    transport_coeff2   = 0.;  // not known; put 0

    // from kinetic theory
    transport_coeff1_s = 8./5.*(1./3.-cs2)*Bulk_Relax_time;
    transport_coeff2_s = 0.;  // not known;  put 0

    // Computing Navier-Stokes term (-bulk viscosity * theta)
    NS_term = -bulk*theta_local;

    // Computing relaxation term and nonlinear term:
    // - Bulk - transport_coeff1*Bulk*theta
    tempf = (-(grid_pt->pi_b[rk_flag])
             - transport_coeff1*theta_local*(grid_pt->pi_b[rk_flag]));

    // Computing nonlinear term: + transport_coeff2*Bulk*Bulk
    if (include_BBterm == 1) {
        BB_term = (transport_coeff2*(grid_pt->pi_b[rk_flag])
                   *(grid_pt->pi_b[rk_flag]));
    } else {
        BB_term = 0.0;
    }

    // Computing terms that Couple with shear-stress tensor
    double Wsigma, WW, Shear_Sigma_term, Shear_Shear_term, Coupling_to_Shear;

    if (include_coupling_to_shear == 1) {
        // Computing sigma^mu^nu
        double sigma[4][4], Wmunu[4][4];
        for (int a = 0; a < 4 ; a++) {
            for (int b = a; b < 4; b++) {
                int idx_1d = Util::map_2d_idx_to_1d(a, b);
                sigma[a][b] = sigma_1d[idx_1d];
                Wmunu[a][b] = grid_pt->Wmunu[rk_flag][idx_1d];
            }
        }

        Wsigma = (  Wmunu[0][0]*sigma[0][0]
                  + Wmunu[1][1]*sigma[1][1]
                  + Wmunu[2][2]*sigma[2][2]
                  + Wmunu[3][3]*sigma[3][3]
                  - 2.*(  Wmunu[0][1]*sigma[0][1]
                        + Wmunu[0][2]*sigma[0][2]
                        + Wmunu[0][3]*sigma[0][3])
                  + 2.*(  Wmunu[1][2]*sigma[1][2]
                        + Wmunu[1][3]*sigma[1][3]
                        + Wmunu[2][3]*sigma[2][3]));

        WW = (   Wmunu[0][0]*Wmunu[0][0]
               + Wmunu[1][1]*Wmunu[1][1]
               + Wmunu[2][2]*Wmunu[2][2]
               + Wmunu[3][3]*Wmunu[3][3]
               - 2.*(  Wmunu[0][1]*Wmunu[0][1]
                     + Wmunu[0][2]*Wmunu[0][2]
                     + Wmunu[0][3]*Wmunu[0][3])
               + 2.*(  Wmunu[1][2]*Wmunu[1][2]
                     + Wmunu[1][3]*Wmunu[1][3]
                     + Wmunu[2][3]*Wmunu[2][3]));
        // multiply term by its respective transport coefficient
        Shear_Sigma_term = Wsigma*transport_coeff1_s;
        Shear_Shear_term = WW*transport_coeff2_s;

        // full term that couples to shear is
        Coupling_to_Shear = -Shear_Sigma_term + Shear_Shear_term ;
    } else {
        Coupling_to_Shear = 0.0;
    }
        
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
    -(1/tau_rho)(q[a] + kappa g[a][b]Dtildemu[b]
                 + kappa u[a] u[b]g[b][c]Dtildemu[c])
    +Delta[a][tau] u[eta] q[eta]/tau
    -Delta[a][eta] u[eta] q[tau]/tau
    -u[a]u[b]g[b][e] Dq[e]
*/
double Diss::Make_uqSource(double tau, Cell *grid_pt, int nu, InitData *DATA,
                           int rk_flag, double theta_local, double *a_local,
                           double *sigma_1d) {
    double q[4];
  
    if (DATA->turn_on_diff == 0) return 0.0;
 
    // Useful variables to define
    double epsilon, rhob;
    if (rk_flag == 0) {
        epsilon = grid_pt->epsilon;
        rhob = grid_pt->rhob;
    } else {
        epsilon = grid_pt->prev_epsilon;
        rhob = grid_pt->prev_rhob;
    }
    double pressure = eos->get_pressure(epsilon, rhob);
    double T        = eos->get_temperature(epsilon, rhob);

    double kappa_coefficient = DATA->kappa_coefficient;
    double tau_rho = kappa_coefficient/(T + 1e-15);
    double mub     = eos->get_mu(epsilon, rhob);
    double alpha   = mub/T;
    double kappa   = kappa_coefficient*(rhob/(3.*T*tanh(alpha) + 1e-15)
                                      - rhob*rhob/(epsilon + pressure));

    if (DATA->Initial_profile == 1) {
        // for 1+1D numerical test
        kappa = kappa_coefficient*(rhob/mub);
    }

    // copy the value of \tilde{q^\mu}
    for (int i = 0; i < 4; i++) {
        int idx_1d = Util::map_2d_idx_to_1d(4, i);
        q[i] = (grid_pt->Wmunu[rk_flag][idx_1d]);
    }

    /* -(1/tau_rho)(q[a] + kappa g[a][b]Dtildemu[b] 
     *              + kappa u[a] u[b]g[b][c]Dtildemu[c])
     * + theta q[a] - q[a] u^\tau/tau
     * + Delta[a][tau] u[eta] q[eta]/tau
     * - Delta[a][eta] u[eta] q[tau]/tau
     * - u[a] u[b]g[b][e] Dq[e] -> u[a] q[e] g[e][b] Du[b]
    */    
 
    // first: (1/tau_rho) part
    // recall that dUsup[4][i] = partial_i (muB/T) 
    // and dUsup[4][0] = -partial_tau (muB/T) = partial^tau (muB/T)
    // and a[4] = u^a partial_a (muB/T) = DmuB/T
    // -(1/tau_rho)(q[a] + kappa g[a][b]DmuB/T[b] 
    // + kappa u[a] u[b]g[b][c]DmuB/T[c])
    // a = nu 
    double NS = kappa*(grid_pt->dUsup[4][nu] 
                           + grid_pt->u[rk_flag][nu]*a_local[4]);
    if (isnan(NS)) {
        cout << "Navier Stock term is nan! " << endl;
        cout << q[nu] << endl;
        // derivative already upper index
        cout << grid_pt->dUsup[4][nu] << endl;
        cout << a_local[4] << endl;
        cout << tau_rho << endl;
        cout << kappa << endl;
        cout << grid_pt->u[rk_flag][nu] << endl;
    }
  
    // add a new non-linear term (- q \theta)
    double transport_coeff = 1.0*tau_rho;   // from conformal kinetic theory
    double Nonlinear1 = -transport_coeff*q[nu]*theta_local;

    // add a new non-linear term (-q^\mu \sigma_\mu\nu)
    double transport_coeff_2 = 3./5.*tau_rho;   // from 14-momentum massless
    double sigma[4][4];
    for (int ii = 0; ii < 4; ii++) {
        for (int jj = ii; jj < 4; jj++) {
            int idx_1d = Util::map_2d_idx_to_1d(ii, jj);
            sigma[ii][jj] = sigma_1d[idx_1d];
        }
    }
    for (int ii = 0; ii < 4; ii++) {
        for (int jj = ii+1; jj < 4; jj++) {
            sigma[jj][ii] = sigma[ii][jj];
        }
    }
    double temptemp = 0.0;
    for (int i = 0 ; i < 4; i++) {
        temptemp += q[i]*sigma[i][nu]*DATA->gmunu[i][i]; 
    }
    double Nonlinear2 = - transport_coeff_2*temptemp;

    double SW = (-q[nu] - NS + Nonlinear1 + Nonlinear2)/(tau_rho + 1e-15);
    if (DATA->Initial_profile == 1) {
        // for 1+1D numerical test
        SW = (-q[nu] - NS)/(tau_rho + 1e-15);
    }

    // all other geometric terms....
    // + theta q[a] - q[a] u^\tau/tau
    SW += (theta_local - grid_pt->u[rk_flag][0]/tau)*q[nu];
 
    if (isnan(SW)) {
        cout << "theta term is nan! " << endl;
    }

    // +Delta[a][tau] u[eta] q[eta]/tau 
    double tempf = ((DATA->gmunu[nu][0] 
                    + grid_pt->u[rk_flag][nu]*grid_pt->u[rk_flag][0])
                      *grid_pt->u[rk_flag][3]*q[3]/tau
                    - (DATA->gmunu[nu][3]
                       + grid_pt->u[rk_flag][nu]*grid_pt->u[rk_flag][3])
                      *grid_pt->u[rk_flag][3]*q[0]/tau);
    SW += tempf;
 
    if (isnan(tempf)) {
        cout << "Delta^{a \tau} and Delta^{a \eta} terms are nan!" << endl;
    }

    //-u[a] u[b]g[b][e] Dq[e] -> u[a] (q[e] g[e][b] Du[b])
    tempf = 0.0;
    for (int i = 0; i < 4; i++) {
        tempf += q[i]*gmn(i)*a_local[i];
    }
    SW += (grid_pt->u[rk_flag][nu])*tempf;
    
    if (isnan(tempf)) {
        cout << "u^a q_b Du^b term is nan! " << endl;
    }

    return SW;
}/* Make_uqSource */


int Diss::Make_uqRHS(double tau, Grid &arena, int ix, int iy, int ieta,
                     double **w_rhs, InitData *DATA, int rk_flag) {
    auto grid_pt = &(arena(ix, iy, ieta));
    int mu, nu, direc;
    double f, fp1, fm1, fp2, fm2, delta[4];
    double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
    double uWphR, uWphL, uWmhR, uWmhL, WphR, WphL, WmhR, WmhL;
    double HWph, HWmh, taufactor, HW;

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

    delta[1] = DATA->delta_x;
    delta[2] = DATA->delta_y;
    delta[3] = DATA->delta_eta;

    // we use the Wmunu[4][nu] = q[nu] 
    mu = 4;
    Cell *grid_pt_p1, *grid_pt_p2, *grid_pt_m1, *grid_pt_m2;
    for (nu = 1; nu < 4; nu++) {
        int idx_1d = Util::map_2d_idx_to_1d(mu, nu);
        double sum = 0.0;
        for (direc=1; direc<=3; direc++) {
            if (direc == 1) {
                grid_pt_p1 = &(arena(ix + 1, iy, ieta));
                grid_pt_p2 = &(arena(ix + 2, iy, ieta));
                grid_pt_m1 = &(arena(ix - 1, iy, ieta));
                grid_pt_m2 = &(arena(ix - 2, iy, ieta));
            } else if (direc == 2) {
                grid_pt_p1 = &(arena(ix, iy + 1, ieta));
                grid_pt_p2 = &(arena(ix, iy + 2, ieta));
                grid_pt_m1 = &(arena(ix, iy - 1, ieta));
                grid_pt_m2 = &(arena(ix, iy - 2, ieta));
            } else {
                grid_pt_p1 = &(arena(ix, iy, ieta + 1));
                grid_pt_p2 = &(arena(ix, iy, ieta + 2));
                grid_pt_m1 = &(arena(ix, iy, ieta - 1));
                grid_pt_m2 = &(arena(ix, iy, ieta - 2));
            }
            if (direc==3)
                taufactor = tau;
            else 
                taufactor = 1.0;

            /* Get_uWmns */
            g = grid_pt->Wmunu[rk_flag][idx_1d];
            f = g*grid_pt->u[rk_flag][direc];
            g *=   grid_pt->u[rk_flag][0];
               
            gp2 = grid_pt_p2->Wmunu[rk_flag][idx_1d];
            fp2 = gp2*grid_pt_p2->u[rk_flag][direc];
            gp2 *=     grid_pt_p2->u[rk_flag][0];
            
            gp1 = grid_pt_p1->Wmunu[rk_flag][idx_1d];
            fp1 = gp1*grid_pt_p1->u[rk_flag][direc];
            gp1 *=     grid_pt_p1->u[rk_flag][0];
            
            gm1 = grid_pt_m1->Wmunu[rk_flag][idx_1d];
            fm1 = gm1*grid_pt_m1->u[rk_flag][direc];
            gm1 *=     grid_pt_m1->u[rk_flag][0];
            
            gm2 = grid_pt_m2->Wmunu[rk_flag][idx_1d];
            fm2 = gm2*grid_pt_m2->u[rk_flag][direc];
            gm2 *=     grid_pt_m2->u[rk_flag][0];
 
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

            am1 = (fabs(grid_pt_m1->u[rk_flag][direc])
                   /grid_pt_m1->u[rk_flag][0]);
            ap1 = (fabs(grid_pt_p1->u[rk_flag][direc])
                   /grid_pt_p1->u[rk_flag][0]);
            ax = maxi(a, ap1);
            HWph = ((uWphR + uWphL) - ax*(WphR - WphL))*0.5;

            ax = maxi(a, am1); 
            HWmh = ((uWmhR + uWmhL) - ax*(WmhR - WmhL))*0.5;

            HW = (HWph - HWmh)/delta[direc]/taufactor;

            /* make partial_i (u^i Wmn) */
            sum += -HW;
        }/* direction */
    
        /* add a source term -u^tau Wmn/tau due to the coordinate 
         * change to tau-eta */
        /* Sangyong Nov 18 2014: don't need this. included in the uqSource. */
        /* this is from udW = d(uW) - Wdu = RHS */
        /* or d(uW) = udW + Wdu */
        /* 
         * sum -= (grid_pt->u[rk_flag][0])*(grid_pt->Wmunu[rk_flag][mu][nu])/tau;
         * sum += (grid_pt->theta_u[rk_flag])*(grid_pt->Wmunu[rk_flag][mu][nu]);
        */  
        w_rhs[mu][nu] = sum*(DATA->delta_tau);
    }/* nu */
    return 1; /* if successful */
} /* Make_uqRHS */

double Diss::get_temperature_dependent_eta_s(InitData *DATA, double T) {
    double Ttr = 0.18/hbarc;  // phase transition temperature
    double Tfrac = T/Ttr;
    double shear_to_s;
    if (T < Ttr) {
        shear_to_s = (DATA->shear_to_s + 0.0594*(1. - Tfrac)
                      + 0.544*(1. - Tfrac*Tfrac));
    } else {
        shear_to_s = (DATA->shear_to_s + 0.288*(Tfrac - 1.) 
                      + 0.0818*(Tfrac*Tfrac - 1.));
    }
    return(shear_to_s);
}

double Diss::get_temperature_dependent_zeta_s(double temperature) {
    // T dependent bulk viscosity from Gabriel
    /////////////////////////////////////////////
    //           Parametrization 1             //
    /////////////////////////////////////////////
    double Ttr=0.18/0.1973;
    double dummy=temperature/Ttr;
    double A1=-13.77, A2=27.55, A3=13.45;
    double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;
 
    double bulk = A1*dummy*dummy + A2*dummy - A3;
    if (temperature < 0.995*Ttr) {
        bulk = (lambda3*exp((dummy-1)/sigma3)
                + lambda4*exp((dummy-1)/sigma4) + 0.03);
    }
    if (temperature > 1.05*Ttr) {
        bulk = (lambda1*exp(-(dummy-1)/sigma1)
                + lambda2*exp(-(dummy-1)/sigma2) + 0.001);
    }
    
    /////////////////////////////////////////////
    //           Parametrization 2             //
    /////////////////////////////////////////////
    //double Ttr=0.18/0.1973;
    //double dummy=temperature/Ttr;
    //double A1=-79.53, A2=159.067, A3=79.04;
    //double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    //double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;

    //bulk = A1*dummy*dummy + A2*dummy - A3;

    //if (temperature < 0.997*Ttr) {
    //    bulk = (lambda3*exp((dummy-1)/sigma3)
    //            + lambda4*exp((dummy-1)/sigma4) + 0.03);
    //}
    //if (temperature > 1.04*Ttr) {
    //    bulk = (lambda1*exp(-(dummy-1)/sigma1)
    //            + lambda2*exp(-(dummy-1)/sigma2) + 0.001);
    //}

    ////////////////////////////////////////////
    //           Parametrization 3            //
    ////////////////////////////////////////////
    //double Ttr=0.18/0.1973;
    //double dummy=temperature/Ttr;
    //double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    //double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;
    
    //if (temperature<0.99945*Ttr) {
    //    bulk = (lambda3*exp((dummy-1)/sigma3)
    //            + lambda4*exp((dummy-1)/sigma4) + 0.03);
    //}
    //if (temperature>0.99945*Ttr) {
    //    bulk = 0.901*exp(14.5*(1.0-dummy)) + 0.061/dummy/dummy;
    //}

    return(bulk);
}


//! this function outputs the T and muB dependence of the baryon diffusion
//! coefficient, kappa
void Diss::output_kappa_T_and_muB_dependence(InitData *DATA) {
    cout << "output kappa_B(T, mu_B) ..." << endl;
    ofstream of("kappa_B_T_and_muB_dependence.dat");
    // write out the header of the file
    of << "# e (GeV/fm^3)  rhob (1/fm^3) T (GeV)  mu_B (GeV)  kappa (1/fm^2)"
       << endl;

    // define the grid
    double e_min = 1e-5;     // fm^-4
    double e_max = 100.0;    // fm^-4
    int ne = 1000;
    double de = (e_max - e_min)/(ne - 1.);
    double rhob_min = 0.0;   // fm^-3
    double rhob_max = sqrt(10.0);  // fm^-3
    int nrhob = 1000;
    double drhob = (rhob_max - rhob_min)/(nrhob - 1.);

    for (int i = 0; i < ne; i++) {
        double e_local = e_min + i*de;
        for (int j = 0; j < nrhob; j++) {
            double rhob_local = rhob_min + j*drhob;
            rhob_local *= rhob_local;
            double mu_B_local = eos->get_mu(e_local, rhob_local);
            if (mu_B_local*hbarc > 0.78)
                continue;  // discard points out of the table
            double p_local = eos->get_pressure(e_local, rhob_local);
            double T_local = eos->get_temperature(e_local, rhob_local);
            double alpha_local = mu_B_local/T_local;

            double kappa_local = (DATA->kappa_coefficient
                    *(rhob_local/(3.*T_local*tanh(alpha_local) + 1e-15)
                      - rhob_local*rhob_local/(e_local + p_local)));
            // output
            of << scientific << setw(18) << setprecision(8)
               << e_local*hbarc << "   " << rhob_local << "   "
               << T_local*hbarc << "   " << mu_B_local*hbarc << "   "
               << kappa_local << endl;
        }
    }
    of.close();  // close the file
}


//! this function outputs the T and muB dependence of the baryon diffusion
//! coefficient, kappa_B, along constant s/n_B trajectories
void Diss::output_kappa_along_const_sovernB(InitData *DATA) {
    cout << "output kappa_B(T, mu_B) along constant s/n_B trajectories..."
         << endl;
    
    double sovernB[] = {10.0, 20.0, 30.0, 51.0, 70.0, 94.0, 144.0, 420.0};
    int array_length = sizeof(sovernB)/sizeof(double);
    double s_0 = 0.00;         // 1/fm^3
    double s_max = 100.0;      // 1/fm^3
    double ds = 0.005;         // 1/fm^3
    int ns = static_cast<int>((s_max - s_0)/ds) + 1;
    for (int i = 0; i < array_length; i++) {
        ostringstream file_name;
        file_name << "kappa_B_sovernB_" << sovernB[i] << ".dat";
        ofstream of(file_name.str().c_str());
        // write out the header of the file
        of << "# e (GeV/fm^3)  rhob (1/fm^3) s (1/fm^3)  "
           << "T (GeV)  mu_B (GeV)  kappa (1/fm^2)" << endl;
        for (int j = 0; j < ns; j++) {
            double s_local = s_0 + j*ds;
            double nB_local = s_local/sovernB[i];
            double e_local = eos->get_s2e_finite_rhob(s_local, nB_local);
            double s_check = eos->get_entropy(e_local, nB_local);
            double p_local = eos->get_pressure(e_local, nB_local);
            double temperature = eos->get_temperature(e_local, nB_local);
            double mu_B = eos->get_mu(e_local, nB_local);
            if (mu_B*hbarc > 0.78)
                continue;  // discard points out of the table
            double alpha_local = mu_B/temperature;

            double kappa_local = (DATA->kappa_coefficient
                    *(nB_local/(3.*temperature*tanh(alpha_local) + 1e-15)
                      - nB_local*nB_local/(e_local + p_local)));
            // output
            of << scientific << setw(18) << setprecision(8)
               << e_local*hbarc << "   " << nB_local << "   "
               << s_check << "   "
               << temperature*hbarc << "   " << mu_B*hbarc << "   "
               << kappa_local << endl;
        }
        of.close();  // close the file
    }
}


//! this function outputs the T and muB dependence of the specific shear
//! viscosity eta/s
void Diss::output_eta_over_s_T_and_muB_dependence(InitData *DATA) {
    cout << "output eta/s(T, mu_B) ..." << endl;
    ofstream of("eta_over_s_T_and_muB_dependence.dat");
    // write out the header of the file
    of << "# e (GeV/fm^3)  rhob (1/fm^3) T (GeV)  mu_B (GeV)  eta/s"
       << endl;

    // define the grid
    double e_min    = 1e-5;     // fm^-4
    double e_max    = 100.0;    // fm^-4
    int ne          = 1000;
    double de       = (e_max - e_min)/(ne - 1.);
    double rhob_min = 0.0;   // fm^-3
    double rhob_max = sqrt(10.0);  // fm^-3
    int nrhob       = 1000;
    double drhob    = (rhob_max - rhob_min)/(nrhob - 1.);

    double etaT_over_enthropy = DATA->shear_to_s;
    for (int i = 0; i < ne; i++) {
        double e_local = e_min + i*de;
        for (int j = 0; j < nrhob; j++) {
            double rhob_local = rhob_min + j*drhob;
            rhob_local *= rhob_local;
            double mu_B_local = eos->get_mu(e_local, rhob_local);
            if (mu_B_local*hbarc > 0.78)
                continue;  // discard points out of the table
            double p_local = eos->get_pressure(e_local, rhob_local);
            double s_local = eos->get_entropy(e_local, rhob_local);
            double T_local = eos->get_temperature(e_local, rhob_local);

            double eta_over_s = (
                etaT_over_enthropy*(e_local + p_local)/(T_local*s_local));

            // output
            of << scientific << setw(18) << setprecision(8)
               << e_local*hbarc << "   " << rhob_local << "   "
               << T_local*hbarc << "   " << mu_B_local*hbarc << "   "
               << eta_over_s << endl;
        }
    }
    of.close();  // close the file
}


//! this function outputs the T and muB dependence of the specific shear
//! viscosity eta/s along constant s/n_B trajectories
void Diss::output_eta_over_s_along_const_sovernB(InitData *DATA) {
    cout << "output eta/s(T, mu_B) along constant s/n_B trajectories..."
         << endl;
    
    double etaT_over_enthropy = DATA->shear_to_s;
    double sovernB[] = {10.0, 20.0, 30.0, 51.0, 70.0, 94.0, 144.0, 420.0};
    int array_length = sizeof(sovernB)/sizeof(double);
    double s_0 = 0.00;         // 1/fm^3
    double s_max = 100.0;      // 1/fm^3
    double ds = 0.005;         // 1/fm^3
    int ns = static_cast<int>((s_max - s_0)/ds) + 1;
    for (int i = 0; i < array_length; i++) {
        ostringstream file_name;
        file_name << "eta_over_s_sovernB_" << sovernB[i] << ".dat";
        ofstream of(file_name.str().c_str());
        // write out the header of the file
        of << "# e (GeV/fm^3)  rhob (1/fm^3) s (1/fm^3)  "
           << "T (GeV)  mu_B (GeV)  eta/s" << endl;
        for (int j = 0; j < ns; j++) {
            double s_local = s_0 + j*ds;
            double nB_local = s_local/sovernB[i];
            double e_local = eos->get_s2e_finite_rhob(s_local, nB_local);
            double s_check = eos->get_entropy(e_local, nB_local);
            double p_local = eos->get_pressure(e_local, nB_local);
            double temperature = eos->get_temperature(e_local, nB_local);
            double mu_B = eos->get_mu(e_local, nB_local);
            if (mu_B*hbarc > 0.78)
                continue;  // discard points out of the table

            double eta_over_s = (
                etaT_over_enthropy*(e_local + p_local)/(temperature*s_local));

            // output
            of << scientific << setw(18) << setprecision(8)
               << e_local*hbarc << "   " << nB_local << "   "
               << s_check << "   "
               << temperature*hbarc << "   " << mu_B*hbarc << "   "
               << eta_over_s << endl;
        }
        of.close();  // close the file
    }
}
