#include "util.h"
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "minmod.h"
#include "eos.h"
#include "u_derivative.h"

U_derivative::U_derivative(const InitData &DATA_in, const EOS &eosIn) :
    DATA(DATA_in),
    eos(eosIn),
    minmod(DATA_in) {
    dUsup = {0.0};        // dUsup[m][n] = partial^n u^m
    dUoverTsup = {0.0};   // dUoverTsup[m][n] = partial^n (u^m/T)
}

//! This function is a shell function to calculate parital^\nu u^\mu
void U_derivative::MakedU(const double tau, SCGrid &arena_prev,
                          SCGrid &arena_current,
                          const int ix, const int iy, const int ieta) {
    dUsup = {0.0};
    dUoverTsup = {0.0};

    // this calculates du/dx, du/dy, (du/deta)/tau
    MakeDSpatial(tau, arena_current, ix, iy, ieta);
    // this calculates du/dtau
    MakeDTau(tau, &arena_prev(ix, iy, ieta), &arena_current(ix, iy, ieta));
}


//! this function returns the expansion rate on the grid
double U_derivative::calculate_expansion_rate(
        double tau, SCGrid &arena, int ieta, int ix, int iy) {
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
void U_derivative::calculate_Du_supmu(const double tau, SCGrid &arena,
                                      const int ieta, const int ix,
                                      const int iy, DumuVec &a) {
    for (int mu = 0; mu <= 4; mu++) {
        double u_supnu_partial_nu_u_supmu = 0.0;
        for (int nu = 0; nu < 4; nu++) {
            double tfac = (nu==0 ? -1.0 : 1.0);
            u_supnu_partial_nu_u_supmu += (
                tfac*arena(ix,iy,ieta).u[nu]*dUsup[mu][nu]);
        }
        a[mu] = u_supnu_partial_nu_u_supmu;
    }
}


void U_derivative::calculate_kinetic_vorticity(
            const double tau, SCGrid &arena,
            const int ieta, const int ix, const int iy,
            const DumuVec &a_local, VorticityVec &omega) {
    // this function computes the kinetic vorticity
    FlowVec u_local = arena(ix, iy, ieta).u;
    double dUsup_local[4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            dUsup_local[i][j] = dUsup[i][j];
        }
    }

    double omega_local[4][4];
    for (int mu = 0; mu < 4; mu++) {
        omega_local[mu][mu] = 0.0;
    }

    // compute the spatial components of omega^{\mu\nu}
    for (int mu = 1; mu < 4; mu++) {
        for (int nu = mu + 1; nu < 4; nu++) {
            omega_local[mu][nu] = 0.5*(
                  (dUsup_local[nu][mu] - dUsup_local[mu][nu])
                + (u_local[mu]*a_local[nu] - u_local[nu]*a_local[mu])
                + u_local[3]*u_local[0]/tau*(  u_local[mu]*DATA.gmunu[nu][3]
                                             - u_local[nu]*DATA.gmunu[mu][3])
                //- u_local[3]/tau*(- DATA.gmunu[mu][0]*DATA.gmunu[nu][3]
                //                  + DATA.gmunu[mu][3]*DATA.gmunu[nu][0])
                //+ u_local[3]*u_local[3]/tau*(- u_local[mu]*DATA.gmunu[nu][0]
                //                             + u_local[nu]*DATA.gmunu[mu][0])
            );
            omega_local[nu][mu] = -omega_local[mu][nu];
        }
    }

    // compute omega^{\tau x}, omega^{\tau y}, omega^{\tau\eta}
    // using u_\mu omega^{\mu\nu} = 0
    for (int mu = 1; mu < 4; mu++) {
        double temp = 0.0;
        for (int nu = 1; nu < 4; nu++) {
            temp += omega_local[mu][nu]*u_local[nu];
        }
        omega_local[0][mu] = temp/u_local[0];
    }

    omega[0] = omega_local[0][1];
    omega[1] = omega_local[0][2];
    omega[2] = omega_local[0][3];
    omega[3] = omega_local[1][2];
    omega[4] = omega_local[1][3];
    omega[5] = omega_local[2][3];
}


void U_derivative::calculate_thermal_vorticity(
            const double tau, SCGrid &arena, const int ieta,
            const int ix, const int iy, VorticityVec &omega) {
    // this function computes the thermal vorticity
    FlowVec u_local = arena(ix, iy, ieta).u;
    double T_local  = eos.get_temperature(arena(ix, iy, ieta).epsilon,
                                          arena(ix, iy, ieta).rhob);
    double dUsup_local[4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            dUsup_local[i][j] = dUoverTsup[i][j];
        }
    }

    double omega_thermal[4][4];
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = mu + 1; nu < 4; nu++) {
            omega_thermal[mu][nu] = -0.5*(
                (dUsup_local[nu][mu] - dUsup_local[mu][nu])
                - u_local[3]/(2.*tau*T_local)*(
                    - DATA.gmunu[mu][0]*DATA.gmunu[nu][3]
                    + DATA.gmunu[mu][3]*DATA.gmunu[nu][0])
            );
        }
    }

    omega[0] = omega_thermal[0][1];
    omega[1] = omega_thermal[0][2];
    omega[2] = omega_thermal[0][3];
    omega[3] = omega_thermal[1][2];
    omega[4] = omega_thermal[1][3];
    omega[5] = omega_thermal[2][3];
}


//! This funciton returns the velocity shear tensor sigma^\mu\nu
void U_derivative::calculate_velocity_shear_tensor(
        const double tau, SCGrid &arena, const int ieta, const int ix,
        const int iy, const DumuVec &a_local, VelocityShearVec &sigma) {
    FlowVec u_local = arena(ix, iy, ieta).u;
    double dUsup_local[4][4];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            dUsup_local[i][j] = dUsup[i][j];
        }
    }
    double theta_u_local = calculate_expansion_rate(tau, arena, ieta, ix, iy);
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
                + u_local[0]/tau*DATA.gmunu[a][3]*DATA.gmunu[b][3]
                + u_local[3]*u_local[0]/tau/2.
                  *(DATA.gmunu[a][3]*u_local[b] + DATA.gmunu[b][3]*u_local[a])
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

//! this function returns the vector D^\mu(\mu_B/T)
void U_derivative::get_DmuMuBoverTVec(DmuMuBoverTVec &vec) {
    for (int mu = 0; mu < 4; mu++)
        vec[mu] = dUsup[4][mu];
}


int U_derivative::MakeDSpatial(const double tau, SCGrid &arena,
                               const int ix, const int iy, const int ieta) {
    // taken care of the tau factor
    const double delta[4] = {0.0, DATA.delta_x, DATA.delta_y,
                             DATA.delta_eta*tau};

    // calculate dUsup[m][n] = partial_n u_m
    Neighbourloop(arena, ix, iy, ieta, NLAMBDAS{
        const double T   = eos.get_temperature(c.epsilon, c.rhob);
        const double Tp1 = eos.get_temperature(p1.epsilon, p1.rhob);
        const double Tm1 = eos.get_temperature(m1.epsilon, m1.rhob);
        for (int m = 0; m <= 3; m++) {
            const double f   = c.u[m];
            const double fp1 = p1.u[m];
            const double fm1 = m1.u[m];
            dUsup[m][direction] = (minmod.minmod_dx(fp1, f, fm1)
                                   /delta[direction]);
            dUoverTsup[m][direction] = (minmod.minmod_dx(fp1/Tp1, f/T, fm1/Tm1)
                                        /delta[direction]);
        }
    });

    /* for u[0], use u[0]u[0] = 1 + u[i]u[i] */
    /* u[0]_m = u[i]_m (u[i]/u[0]) */
    /* for u[0] */
    for (int n = 1; n <= 3; n++) {
        double f = 0.0;
        for (int m = 1; m <= 3; m++) {
            // (partial_n u^m) u[m]
            f += dUsup[m][n]*(arena(ix, iy, ieta).u[m]);
        }
        f /= arena(ix, iy, ieta).u[0];
        dUsup[0][n] = f;
    }

    // Sangyong Nov 18 2014
    // Here we make derivatives of muB/T
    // dUsup[rk_flag][4][n] = partial_n (muB/T)
    // partial_x (muB/T) and partial_y (muB/T) first
    const int m = 4;  // means (muB/T)
    const double rhob = arena(ix, iy, ieta).rhob;
    const double eps = arena(ix, iy, ieta).epsilon;
    const double muB = eos.get_muB(eps, rhob);
    const double T = eos.get_temperature(eps, rhob);
    const double f = muB/T;
    Neighbourloop(arena, ix, iy, ieta, NLAMBDAS{
        const double fp1 = (eos.get_muB(p1.epsilon, p1.rhob)
                            /eos.get_temperature(p1.epsilon, p1.rhob));
        const double fm1 = (eos.get_muB(m1.epsilon, m1.rhob)
                            /eos.get_temperature(m1.epsilon, m1.rhob));
        dUsup[m][direction] = minmod.minmod_dx(fp1, f, fm1)/delta[direction];
    });
    return 1;
}/* MakeDSpatial */

int U_derivative::MakeDTau(const double tau,
                           const Cell_small *grid_pt_prev,
                           const Cell_small *grid_pt) {
    /* this makes dU[m][0] = partial^tau u^m */
    /* note the minus sign at the end because of g[0][0] = -1 */

    const double eps  = grid_pt->epsilon;
    const double rhob = grid_pt->rhob;
    const double eps_prev  = grid_pt_prev->epsilon;
    const double rhob_prev = grid_pt_prev->rhob;
    const double T = eos.get_temperature(eps, rhob);
    const double T_prev = eos.get_temperature(eps_prev, rhob_prev);

    for (int m = 0; m < 4; m++) {
        // first order is more stable
        double f = (grid_pt->u[m] - grid_pt_prev->u[m])/DATA.delta_tau;
        dUsup[m][0] = -f;  // g^{00} = -1

        double duoverTdtau = ((grid_pt->u[m]/T - grid_pt_prev->u[m]/T_prev)
                              /DATA.delta_tau);
        dUoverTsup[m][0] = -duoverTdtau;   // g^{00} = -1
    }

    /* I have now partial^tau u^i */
    /* I need to calculate (u^i partial^tau u^i) = u^0 partial^tau u^0 */
    /* u_0 d^0 u^0 + u_m d^0 u^m = 0 */
    /* -u^0 d^0 u^0 + u_m d^0 u^m = 0 */
    /* d^0 u^0 = u_m d^0 u^m/u^0 */
    double f = 0.0;
    for (int m = 1; m < 4; m++) {
        /* (partial_0 u^m) u[m] */
        f += dUsup[m][0]*(grid_pt->u[m]);
    }
    dUsup[0][0] = f/grid_pt->u[0];

    // Sangyong Nov 18 2014
    // Here we make the time derivative of (muB/T)
    int m = 4;
    // first order is more stable backward derivative
    const double muB = eos.get_muB(eps, rhob);
    const double tildemu = muB/T;
    const double muB_prev = eos.get_muB(eps, rhob);
    const double tildemu_prev = muB_prev/T_prev;
    f = (tildemu - tildemu_prev)/(DATA.delta_tau);
    dUsup[m][0]  = -f;  // g^{00} = -1
    return 1;
}
