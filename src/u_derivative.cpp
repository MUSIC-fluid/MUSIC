#include <omp.h>
#include "./util.h"
#include "./data.h"
#include "./cell.h"
#include "./grid.h"
#include "./minmod.h"
#include "./eos.h"
#include "./u_derivative.h"

using namespace std;

U_derivative::U_derivative(EOS *eosIn, InitData* DATA_in) :
    minmod(DATA_in)
{
   eos = eosIn;
   DATA_ptr = DATA_in;
   for (int i = 0; i < 5; i++) {
       for (int j = 0; j < 4; j++) {
           dUsup[i][j] = 0.0;
       }
   }
}

//! This function is a shell function to calculate parital^\nu u^\mu
int U_derivative::MakedU(double tau, InitData *DATA,
                         SCGrid &arena_prev, SCGrid &arena_current,
                         int ix, int iy, int ieta, int rk_flag) {
    // ideal hydro: no need to evaluate any flow derivatives
    if (DATA->viscosity_flag == 0) {
        return(1);
    }

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 4; j++) {
            dUsup[i][j] = 0.0;
        }
    }

	// this calculates du/dx, du/dy, (du/deta)/tau
    MakeDSpatial(tau, DATA, arena_current, ix, iy, ieta, rk_flag);
    // this calculates du/dtau
    MakeDTau(tau, DATA, &(arena_prev(ix, iy, ieta)), &(arena_current(ix,iy,ieta)), rk_flag); 

   return(1);
}


//! this function returns the expansion rate on the grid
double U_derivative::calculate_expansion_rate(
        double tau, SCGrid &arena, int ieta, int ix, int iy, int rk_flag) {
    double partial_mu_u_supmu = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        double gfac = (mu == 0 ? -1.0 : 1.0);
        // for expansion rate: theta
        partial_mu_u_supmu += dUsup[mu][mu]*gfac;
    }
    double theta = partial_mu_u_supmu + arena(ix,iy,ieta).u[0]/tau;
    return(theta);
}


//! this function returns Du^\mu
void U_derivative::calculate_Du_supmu(double tau, SCGrid &arena, int ieta,
                                      int ix, int iy, int rk_flag, DumuVec &a) {
    for (int mu = 0; mu <= 4; mu++) {
        double u_supnu_partial_nu_u_supmu = 0.0;
	    for (int nu = 0; nu < 4; nu++) {
            double tfac = (nu==0 ? -1.0 : 1.0);
            u_supnu_partial_nu_u_supmu += (
                tfac*arena(ix,iy,ieta).u[nu]
                *dUsup[mu][nu]);
        }
        a[mu] = u_supnu_partial_nu_u_supmu;
    }
}


//! This funciton returns the velocity shear tensor sigma^\mu\nu
void U_derivative::calculate_velocity_shear_tensor(double tau, SCGrid &arena, int ieta, int ix, int iy, int rk_flag, DumuVec &a_local, VelocityShearVec &sigma){
    FlowVec u_local = arena(ix, iy, ieta).u;
    double dUsup_local[4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            dUsup_local[i][j] = dUsup[i][j];
        }
    }
    double theta_u_local = calculate_expansion_rate(tau, arena, ieta, ix,
                                                    iy, rk_flag);
    double gfac = 0.0;
    double sigma_local[4][4];
    for (int a = 1; a < 4; a++) {
        for (int b = a; b < 4; b++) {
            if (b == a) {
                gfac = 1.0;
            } else {
                gfac = 0.0;
            }
            sigma_local[a][b] = ((dUsup_local[a][b] + dUsup_local[b][a])/2.
                - (gfac + u_local[a]*u_local[b])*theta_u_local/3.
                + u_local[0]/tau*DATA_ptr->gmunu[a][3]*DATA_ptr->gmunu[b][3]
                + u_local[3]*u_local[0]/tau/2.
                  *(DATA_ptr->gmunu[a][3]*u_local[b] 
                    + DATA_ptr->gmunu[b][3]*u_local[a])
                + (u_local[a]*a_local[b] + u_local[b]*a_local[a])/2.);
            sigma_local[b][a] = sigma_local[a][b];
        }
    }
    // make sigma[3][3] using traceless condition
    sigma_local[3][3] = (
        (  2.*(  u_local[1]*u_local[2]*sigma_local[1][2]
               + u_local[1]*u_local[3]*sigma_local[1][3]
               + u_local[2]*u_local[3]*sigma_local[2][3])
         - (u_local[0]*u_local[0] - u_local[1]*u_local[1])*sigma_local[1][1]
         - (u_local[0]*u_local[0] - u_local[2]*u_local[2])*sigma_local[2][2])
        /(u_local[0]*u_local[0] - u_local[3]*u_local[3]));
    // make sigma[0][i] using transversality
    for (int a = 1; a < 4; a++) {
        double temp = 0.0;
        for (int b = 1; b < 4; b++) {
            temp += sigma_local[a][b]*u_local[b];
        }
        sigma_local[0][a] = temp/u_local[0];
    }
    // make sigma[0][0]
    double temp = 0.0;
    for (int a = 1; a < 4; a++) {
        temp += sigma_local[0][a]*u_local[a];
    }
    sigma_local[0][0] = temp/u_local[0];

    sigma[0] = sigma_local[0][0];
    sigma[1] = sigma_local[0][1];
    sigma[2] = sigma_local[0][2];
    sigma[3] = sigma_local[0][3];
    sigma[4] = sigma_local[1][1];
    sigma[5] = sigma_local[1][2];
    sigma[6] = sigma_local[1][3];
    sigma[7] = sigma_local[2][2];
    sigma[8] = sigma_local[2][3];
    sigma[9] = sigma_local[3][3];
}


int U_derivative::MakeDSpatial(double tau, InitData *DATA, SCGrid &arena, int ix, int iy, int ieta,
                               int rk_flag) {
    const double delta[4] = {
      0.0,
      DATA->delta_x,
      DATA->delta_y,
      DATA->delta_eta*tau
    };  // taken care of the tau factor

    // calculate dUsup[m][n] = partial_n u_m
    Neighbourloop(arena, ix, iy, ieta, NLAMBDAS{
        for (int m = 1; m <= 3; m++) {
            const double f   = c.u[m];
            const double fp1 = p1.u[m];
            const double fm1 = m1.u[m];
            const double g   = minmod.minmod_dx(fp1, f, fm1) / delta[direction];
            dUsup[m][direction] = g;
        }
    });

    /* for u[0], use u[0]u[0] = 1 + u[i]u[i] */
    /* u[0]_m = u[i]_m (u[i]/u[0]) */
    /* for u[0] */

    for (int n = 1; n <= 3; n++) {
        double f = 0.0;
        for (int m = 1; m <= 3; m++) {
	        /* (partial_n u^m) u[m] */
	        f += dUsup[m][n]*(arena(ix,iy,ieta).u[m]);
        }
        f /= arena(ix,iy,ieta).u[0];
        dUsup[0][n] = f;
    }
    // Sangyong Nov 18 2014
    // Here we make derivatives of muB/T
    // dUsup[rk_flag][4][n] = partial_n (muB/T)
    // partial_x (muB/T) and partial_y (muB/T) first
    int m = 4; // means (muB/T)
    double rhob, eps, muB, T;
    rhob = arena(ix,iy,ieta).rhob;
    eps = arena(ix,iy,ieta).epsilon;
    muB = eos->get_mu(eps, rhob);
    T = eos->get_temperature(eps, rhob);
    double f = muB/T; 
    Neighbourloop(arena, ix, iy, ieta, NLAMBDAS{
        double fp1, fm1;
        fp1 = ( eos->get_mu(p1.epsilon, p1.rhob)
               /eos->get_temperature(p1.epsilon, p1.rhob));
        fm1 = ( eos->get_mu(m1.epsilon, m1.rhob)
               /eos->get_temperature(m1.epsilon, m1.rhob));
        double g = minmod.minmod_dx(fp1, f, fm1)/delta[direction];
        dUsup[m][direction] = g;
    });
    return 1;
}/* MakeDSpatial */

int U_derivative::MakeDTau(double tau, InitData *DATA,
                           Cell_small *grid_pt_prev, Cell_small *grid_pt,
                           int rk_flag) {
    int m;
    double f;
    // Sangyong Nov 18 2014: added these doubles 
    double tildemu, tildemu_prev, rhob, eps, muB, T;

    /* this makes dU[m][0] = partial^tau u^m */
    /* note the minus sign at the end because of g[0][0] = -1 */
    /* rk_flag is 0, 2, 4, ... */

    for (m=1; m<=3; m++) {
        /* first order is more stable */
        f = ((grid_pt->u[m] - grid_pt_prev->u[m])
             /DATA->delta_tau);
        dUsup[m][0] = -f; /* g00 = -1 */
    }/* m */

    /* I have now partial^tau u^i */
    /* I need to calculate (u^i partial^tau u^i) = u^0 partial^tau u^0 */
    /* u_0 d^0 u^0 + u_m d^0 u^m = 0 */
    /* -u^0 d^0 u^0 + u_m d^0 u^m = 0 */
    /* d^0 u^0 = u_m d^0 u^m/u^0 */

    f = 0.0;
    for (m=1; m<=3; m++) {
        /* (partial_0 u^m) u[m] */
        f += dUsup[m][0]*(grid_pt->u[m]);
    }
    f /= grid_pt->u[0];
    dUsup[0][0] = f;

    // Sangyong Nov 18 2014
    // Here we make the time derivative of (muB/T)
    m = 4;  
    // first order is more stable 
    // backward derivative
    // current values
    // f         = (grid_pt->rhob);
    rhob         = grid_pt->rhob;
    eps          = grid_pt->epsilon;
    muB          = eos->get_mu(eps, rhob);
    T            = eos->get_temperature(eps, rhob);
    tildemu      = muB/T;
    // f        -= (grid_pt->rhob_prev);
    rhob         = grid_pt_prev->rhob;
    eps          = grid_pt_prev->epsilon;
    muB          = eos->get_mu(eps, rhob);
    T            = eos->get_temperature(eps, rhob);
    tildemu_prev = muB/T;
    f            = (tildemu - tildemu_prev)/(DATA->delta_tau);
    dUsup[m][0] = -f; /* g00 = -1 */
    // Ends Sangyong's addition Nov 18 2014
    return 1;
}/* MakeDTau */

