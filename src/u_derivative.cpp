#include "u_derivative.h"

#include "cell.h"
#include "data.h"
#include "minmod.h"
#include "util.h"

U_derivative::U_derivative(const InitData &DATA_in, const EOS &eosIn)
    : DATA(DATA_in), eos(eosIn), minmod(DATA_in) {
    dUsup = {{{0.0}}};       // dUsup[m][n] = partial^n u^m
    dUoverTsup = {{{0.0}}};  // dUoverTsup[m][n] = partial^n (u^m/T)
    dUTsup = {{{0.0}}};      // dUTsup[m][n] = partial^n (Tu^m)
}

//! This function is a shell function to calculate parital^\nu u^\mu
void U_derivative::MakedU(
    const double tau, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
    const int fieldIdx, const int ix, const int iy, const int ieta) {
    dUsup = {{{0.0}}};
    dUoverTsup = {{{0.0}}};
    dUTsup = {{{0.0}}};

    // this calculates du/dx, du/dy, (du/deta)/tau
    MakeDSpatial(tau, arenaFieldsCurr, fieldIdx, ix, iy, ieta);
    // this calculates du/dtau
    MakeDTau(tau, arenaFieldsPrev, arenaFieldsCurr, fieldIdx);
}

//! this function returns the expansion rate on the grid
double U_derivative::calculate_expansion_rate(
    const double tau, Fields &arena, const int fieldIdx) {
    double partial_mu_u_supmu = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        double gfac = (mu == 0 ? -1.0 : 1.0);
        // for expansion rate: theta
        partial_mu_u_supmu += dUsup[mu][mu] * gfac;
    }
    double theta = partial_mu_u_supmu + arena.u_[0][fieldIdx] / tau;
    return (theta);
}

//! this function returns Du^\mu
void U_derivative::calculate_Du_supmu(
    const double tau, Fields &arena, const int fieldIdx, DumuVec &a) {
    for (int mu = 0; mu <= 4; mu++) {
        double u_supnu_partial_nu_u_supmu = 0.0;
        for (int nu = 0; nu < 4; nu++) {
            double tfac = (nu == 0 ? -1.0 : 1.0);
            u_supnu_partial_nu_u_supmu +=
                (tfac * arena.u_[nu][fieldIdx] * dUsup[mu][nu]);
        }
        a[mu] = u_supnu_partial_nu_u_supmu;
    }
}

// This is a shell function to compute all 4 kinds of vorticity tensors
void U_derivative::compute_vorticity_shell(
    const double tau, Fields &arena_prev, Fields &arena_curr, const int ieta,
    const int ix, const int iy, const int fieldIdx, const double eta,
    VorticityVec &omega_local_kSP, VorticityVec &omega_local_knoSP,
    VorticityVec &omega_local_th, VorticityVec &omega_local_T,
    VelocityShearVec &sigma_th_local, DmuMuBoverTVec &DbetaMu) {
    MakedU(tau, arena_prev, arena_curr, fieldIdx, ix, iy, ieta);
    DumuVec a_local;
    calculate_Du_supmu(tau, arena_curr, fieldIdx, a_local);

    VorticityVec omega_local;
    calculate_kinetic_vorticity_with_spatial_projector(
        tau, arena_curr, fieldIdx, a_local, omega_local);
    omega_local_kSP = transform_vorticity_to_tz(omega_local, eta);
    calculate_kinetic_vorticity_no_spatial_projection(
        tau, arena_curr, fieldIdx, omega_local);
    omega_local_knoSP = transform_vorticity_to_tz(omega_local, eta);
    calculate_thermal_vorticity(tau, arena_curr, fieldIdx, omega_local);
    omega_local_th = transform_vorticity_to_tz(omega_local, eta);
    calculate_T_vorticity(tau, arena_curr, fieldIdx, omega_local);
    omega_local_T = transform_vorticity_to_tz(omega_local, eta);

    VelocityShearVec sigma_th_Mline;
    calculate_thermal_shear_tensor(tau, arena_curr, fieldIdx, sigma_th_Mline);
    sigma_th_local = transform_SigmaMuNu_to_tz(sigma_th_Mline, eta);
    DmuMuBoverTVec DbetaMu_local;
    get_DmuMuBoverTVec(DbetaMu_local);
    DbetaMu = transform_vector_to_tz(DbetaMu_local, eta);
}

// This function transforms the vorticity tensor from tau-eta to tz
// It outputs (sigma^tt, sigma^tx, sigma^ty, sigma^tz,
//             sigma^xx, sigma^xy, sigma^xz,
//             sigma^yy, sigma^yz,
//             sigma^zz)
VelocityShearVec U_derivative::transform_SigmaMuNu_to_tz(
    const VelocityShearVec sigma_Mline, const double eta) {
    VelocityShearVec sigma_Cart;
    const double cosh_eta = cosh(eta);
    const double sinh_eta = sinh(eta);
    sigma_Cart[0] =
        (sigma_Mline[0] * cosh_eta * cosh_eta
         + 2. * sigma_Mline[3] * cosh_eta * sinh_eta
         + sigma_Mline[9] * sinh_eta * sinh_eta);
    sigma_Cart[1] = sigma_Mline[1] * cosh_eta + sigma_Mline[6] * sinh_eta;
    sigma_Cart[2] = sigma_Mline[2] * cosh_eta + sigma_Mline[8] * sinh_eta;
    sigma_Cart[3] =
        (sigma_Mline[0] * cosh_eta * sinh_eta
         + sigma_Mline[3] * (cosh_eta * cosh_eta + sinh_eta * sinh_eta)
         + sigma_Mline[9] * sinh_eta * cosh_eta);
    sigma_Cart[4] = sigma_Mline[4];
    sigma_Cart[5] = sigma_Mline[5];
    sigma_Cart[6] = sigma_Mline[1] * sinh_eta + sigma_Mline[6] * cosh_eta;
    sigma_Cart[7] = sigma_Mline[7];
    sigma_Cart[8] = sigma_Mline[2] * sinh_eta + sigma_Mline[8] * cosh_eta;
    sigma_Cart[9] =
        (sigma_Mline[0] * sinh_eta * sinh_eta
         + 2. * sigma_Mline[3] * cosh_eta * sinh_eta
         + sigma_Mline[9] * cosh_eta * cosh_eta);
    return (sigma_Cart);
}

// This function transforms the vorticity tensor from tau-eta to tz
// It outputs (omega^tx, omega^ty, omega^tz, omega^xy, omega^xz, omega^yz)
VorticityVec U_derivative::transform_vorticity_to_tz(
    const VorticityVec omega_Mline, const double eta) {
    VorticityVec omega_Cart;
    const double cosh_eta = cosh(eta);
    const double sinh_eta = sinh(eta);
    omega_Cart[0] = omega_Mline[0] * cosh_eta - omega_Mline[4] * sinh_eta;
    omega_Cart[1] = omega_Mline[1] * cosh_eta - omega_Mline[5] * sinh_eta;
    omega_Cart[2] = omega_Mline[2];
    omega_Cart[3] = omega_Mline[3];
    omega_Cart[4] = omega_Mline[4] * cosh_eta - omega_Mline[0] * sinh_eta;
    omega_Cart[5] = omega_Mline[5] * cosh_eta - omega_Mline[1] * sinh_eta;
    return (omega_Cart);
}

// This function transforms the vorticity tensor from tau-eta to tz
// It outputs (vec^t, vec^x, vec^y, vec^z)
DmuMuBoverTVec U_derivative::transform_vector_to_tz(
    const DmuMuBoverTVec vec_Mline, const double eta) {
    DmuMuBoverTVec vec_Cart;
    const double cosh_eta = cosh(eta);
    const double sinh_eta = sinh(eta);
    vec_Cart[0] = vec_Mline[0] * cosh_eta + vec_Mline[3] * sinh_eta;
    vec_Cart[1] = vec_Mline[1];
    vec_Cart[2] = vec_Mline[2];
    vec_Cart[3] = vec_Mline[0] * sinh_eta + vec_Mline[3] * cosh_eta;
    return (vec_Cart);
}

//! this function computes the kinetic vorticity with the spatial projectors
//! the output omega^{\mu\nu} can be directly used in the EoM for the
//! shear viscous tensor and diffusion current
//! Please note that it is defined with the opposite sign compared to the
//! conventional kinetic vorticity in the literature
//! Because MUSIC use the metric g = (-1, 1, 1, 1), the output omega^{\mu\nu}
//! differs from the ones with g = (1, -1, -1, -1) by a minus sign
void U_derivative::calculate_kinetic_vorticity_with_spatial_projector(
    const double tau, Fields &arena, const int fieldIdx, const DumuVec &a_local,
    VorticityVec &omega) {
    FlowVec u_local;
    double dUsup_local[4][4];
    for (int i = 0; i < 4; i++) {
        u_local[i] = arena.u_[i][fieldIdx];
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
            omega_local[mu][nu] =
                0.5
                * ((dUsup_local[nu][mu] - dUsup_local[mu][nu])
                   + (u_local[mu] * a_local[nu] - u_local[nu] * a_local[mu])
                   + u_local[3] * u_local[0] / tau
                         * (u_local[mu] * DATA.gmunu[nu][3]
                            - u_local[nu] * DATA.gmunu[mu][3])
                   //- u_local[3]/tau*(- DATA.gmunu[mu][0]*DATA.gmunu[nu][3]
                   //                  + DATA.gmunu[mu][3]*DATA.gmunu[nu][0])
                   //+ u_local[3]*u_local[3]/tau*(-
                   // u_local[mu]*DATA.gmunu[nu][0]
                   //                             +
                   //                             u_local[nu]*DATA.gmunu[mu][0])
                );
            omega_local[nu][mu] = -omega_local[mu][nu];
        }
    }

    // compute omega^{\tau x}, omega^{\tau y}, omega^{\tau\eta}
    // using u_\mu omega^{\mu\nu} = 0
    for (int mu = 1; mu < 4; mu++) {
        double temp = 0.0;
        for (int nu = 1; nu < 4; nu++) {
            temp += omega_local[mu][nu] * u_local[nu];
        }
        omega_local[0][mu] = temp / u_local[0];
    }

    omega[0] = omega_local[0][1];
    omega[1] = omega_local[0][2];
    omega[2] = omega_local[0][3];
    omega[3] = omega_local[1][2];
    omega[4] = omega_local[1][3];
    omega[5] = omega_local[2][3];
}

//! this function computes the thermal vorticity
//! it outputs omega^{\mu\nu} in the metric g = (-1, 1, 1, 1) which differs
//! from the ones with g = (1, -1, -1, -1) by a minus sign
void U_derivative::calculate_thermal_vorticity(
    const double tau, Fields &arena, const int fieldIdx, VorticityVec &omega) {
    // this function computes the thermal vorticity
    double T_local =
        eos.get_temperature(arena.e_[fieldIdx], arena.rhob_[fieldIdx]);
    if (T_local > T_tol) {
        FlowVec u_local;
        double dUsup_local[4][4];
        for (int i = 0; i < 4; i++) {
            u_local[i] = arena.u_[i][fieldIdx];
            for (int j = 0; j < 4; j++) {
                dUsup_local[i][j] = dUoverTsup[i][j];
            }
        }

        double omega_thermal[4][4];
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
                omega_thermal[mu][nu] =
                    (-0.5 * (dUsup_local[nu][mu] - dUsup_local[mu][nu])
                     - u_local[3] / (2. * tau * T_local)
                           * (DATA.gmunu[mu][0] * DATA.gmunu[nu][3]
                              - DATA.gmunu[mu][3] * DATA.gmunu[nu][0]));
            }
        }

        omega[0] = omega_thermal[0][1];
        omega[1] = omega_thermal[0][2];
        omega[2] = omega_thermal[0][3];
        omega[3] = omega_thermal[1][2];
        omega[4] = omega_thermal[1][3];
        omega[5] = omega_thermal[2][3];
    } else {
        omega[0] = 0.;
        omega[1] = 0.;
        omega[2] = 0.;
        omega[3] = 0.;
        omega[4] = 0.;
        omega[5] = 0.;
    }
}

//! this function computes the thermal shear tensor
//! it outputs sigma_th^{\mu\nu} in the metric g = (-1, 1, 1, 1) which differs
//! from the ones with g = (1, -1, -1, -1) by a minus sign
void U_derivative::calculate_thermal_shear_tensor(
    const double tau, Fields &arena, const int fieldIdx,
    VelocityShearVec &sigma_th) {
    // this function computes the thermal shear tensor
    double T_local =
        eos.get_temperature(arena.e_[fieldIdx], arena.rhob_[fieldIdx]);
    if (T_local > T_tol) {
        FlowVec u_local;
        double dUsup_local[4][4];
        for (int i = 0; i < 4; i++) {
            u_local[i] = arena.u_[i][fieldIdx];
            for (int j = 0; j < 4; j++) {
                dUsup_local[i][j] = dUoverTsup[i][j];
            }
        }

        double sigma_thermal[4][4];
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu; nu < 4; nu++) {
                sigma_thermal[mu][nu] =
                    (0.5 * (dUsup_local[nu][mu] + dUsup_local[mu][nu])
                     - u_local[3] / (2. * tau * T_local)
                           * (DATA.gmunu[mu][0] * DATA.gmunu[nu][3]
                              + DATA.gmunu[mu][3] * DATA.gmunu[nu][0])
                     + u_local[0] / (tau * T_local)
                           * (DATA.gmunu[mu][3] * DATA.gmunu[nu][3]));
            }
        }

        sigma_th[0] = sigma_thermal[0][0];
        sigma_th[1] = sigma_thermal[0][1];
        sigma_th[2] = sigma_thermal[0][2];
        sigma_th[3] = sigma_thermal[0][3];
        sigma_th[4] = sigma_thermal[1][1];
        sigma_th[5] = sigma_thermal[1][2];
        sigma_th[6] = sigma_thermal[1][3];
        sigma_th[7] = sigma_thermal[2][2];
        sigma_th[8] = sigma_thermal[2][3];
        sigma_th[9] = sigma_thermal[3][3];
    } else {
        for (auto &sigma_i : sigma_th) sigma_i = 0.;
    }
}

//! this function computes the temperature- (T-)vorticity
//! it outputs omega^{\mu\nu} in the metric g = (-1, 1, 1, 1) which differs
//! from the ones with g = (1, -1, -1, -1) by a minus sign
void U_derivative::calculate_T_vorticity(
    const double tau, Fields &arena, const int fieldIdx, VorticityVec &omega) {
    // this function computes the T-vorticity
    FlowVec u_local;
    double T_local =
        eos.get_temperature(arena.e_[fieldIdx], arena.rhob_[fieldIdx]);
    double dUsup_local[4][4];
    for (int i = 0; i < 4; i++) {
        u_local[i] = arena.u_[i][fieldIdx];
        for (int j = 0; j < 4; j++) {
            dUsup_local[i][j] = dUTsup[i][j];
        }
    }

    double omega_thermal[4][4];
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = mu + 1; nu < 4; nu++) {
            omega_thermal[mu][nu] =
                (-0.5 * (dUsup_local[nu][mu] - dUsup_local[mu][nu])
                 - u_local[3] * T_local / (2. * tau)
                       * (DATA.gmunu[mu][0] * DATA.gmunu[nu][3]
                          - DATA.gmunu[mu][3] * DATA.gmunu[nu][0]));
        }
    }

    omega[0] = omega_thermal[0][1];
    omega[1] = omega_thermal[0][2];
    omega[2] = omega_thermal[0][3];
    omega[3] = omega_thermal[1][2];
    omega[4] = omega_thermal[1][3];
    omega[5] = omega_thermal[2][3];
}

//! this function computes the conventional kinetic vorticity
//! it outputs omega^{\mu\nu} in the metric g = (-1, 1, 1, 1) which differs
//! from the ones with g = (1, -1, -1, -1) by a minus sign
void U_derivative::calculate_kinetic_vorticity_no_spatial_projection(
    const double tau, Fields &arena, const int fieldIdx, VorticityVec &omega) {
    // this function computes the full kinetic vorticity without the spatial
    // projection
    // omega^{\mu\nu} = \partial^\mu u^\nu - \partial^\nu u^\mu
    FlowVec u_local;
    double dUsup_local[4][4];
    for (int i = 0; i < 4; i++) {
        u_local[i] = arena.u_[i][fieldIdx];
        for (int j = 0; j < 4; j++) {
            dUsup_local[i][j] = dUsup[i][j];
        }
    }

    double omega_thermal[4][4];
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = mu + 1; nu < 4; nu++) {
            omega_thermal[mu][nu] =
                (-0.5 * (dUsup_local[nu][mu] - dUsup_local[mu][nu])
                 - u_local[3] / (2. * tau)
                       * (+DATA.gmunu[mu][0] * DATA.gmunu[nu][3]
                          - DATA.gmunu[mu][3] * DATA.gmunu[nu][0]));
        }
    }

    omega[0] = omega_thermal[0][1];
    omega[1] = omega_thermal[0][2];
    omega[2] = omega_thermal[0][3];
    omega[3] = omega_thermal[1][2];
    omega[4] = omega_thermal[1][3];
    omega[5] = omega_thermal[2][3];
}

//! This funciton returns the velocity shear tensor sigma^{\mu\nu}
//! it outputs sigma^{\mu\nu} in the metric g = (-1, 1, 1, 1)
//! Please note that this output differs from the sigma^{\mu\nu} in the metric
//! g = (1, -1, -1, -1) by a minus sign
void U_derivative::calculate_velocity_shear_tensor(
    const double tau, Fields &arena, const int fieldIdx,
    const double theta_u_local, const DumuVec &a_local,
    VelocityShearVec &sigma) {
    FlowVec u_local;
    double dUsup_local[4][4];
    for (int i = 0; i < 4; i++) {
        u_local[i] = arena.u_[i][fieldIdx];
        for (int j = 0; j < 4; j++) {
            dUsup_local[i][j] = dUsup[i][j];
        }
    }
    double gfac = 0.0;
    double sigma_local[4][4];
    for (int a = 1; a < 4; a++) {
        for (int b = a; b < 4; b++) {
            if (b == a) {
                gfac = 1.0;
            } else {
                gfac = 0.0;
            }
            sigma_local[a][b] =
                ((dUsup_local[a][b] + dUsup_local[b][a]) / 2.
                 - (gfac + u_local[a] * u_local[b]) * theta_u_local / 3.
                 + u_local[0] / tau * DATA.gmunu[a][3] * DATA.gmunu[b][3]
                 + u_local[3] * u_local[0] / tau / 2.
                       * (DATA.gmunu[a][3] * u_local[b]
                          + DATA.gmunu[b][3] * u_local[a])
                 + (u_local[a] * a_local[b] + u_local[b] * a_local[a]) / 2.);
            sigma_local[b][a] = sigma_local[a][b];
        }
    }
    // make sigma[3][3] using traceless condition
    sigma_local[3][3] =
        ((2.
              * (u_local[1] * u_local[2] * sigma_local[1][2]
                 + u_local[1] * u_local[3] * sigma_local[1][3]
                 + u_local[2] * u_local[3] * sigma_local[2][3])
          - (u_local[0] * u_local[0] - u_local[1] * u_local[1])
                * sigma_local[1][1]
          - (u_local[0] * u_local[0] - u_local[2] * u_local[2])
                * sigma_local[2][2])
         / (u_local[0] * u_local[0] - u_local[3] * u_local[3]));
    // make sigma[0][i] using transversality
    for (int a = 1; a < 4; a++) {
        double temp = 0.0;
        for (int b = 1; b < 4; b++) {
            temp += sigma_local[a][b] * u_local[b];
        }
        sigma_local[0][a] = temp / u_local[0];
    }
    // make sigma[0][0]
    double temp = 0.0;
    for (int a = 1; a < 4; a++) {
        temp += sigma_local[0][a] * u_local[a];
    }
    sigma_local[0][0] = temp / u_local[0];

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
    for (int mu = 0; mu < 4; mu++) vec[mu] = dUsup[4][mu];
}

int U_derivative::MakeDSpatial(
    const double tau, Fields &arena, const int fieldIdx, const int ix,
    const int iy, const int ieta) {
    // taken care of the tau factor
    const double delta[4] = {
        0.0, DATA.delta_x, DATA.delta_y, DATA.delta_eta * tau};

    // calculate dUsup[m][n] = partial^n u^m
    FieldNeighbourLoopIdeal1(
        arena, ix, iy, ieta, FNLILAMBDAS1 {
            for (int m = 1; m < 4; m++) {
                const double f = arena.u_[m][Ic];
                const double fp1 = arena.u_[m][Ip1];
                const double fm1 = arena.u_[m][Im1];
                dUsup[m][direction] =
                    (minmod.minmod_dx(fp1, f, fm1) / delta[direction]);
            }
            if (DATA.output_vorticity == 1) {
                for (int m = 0; m < 4; m++) {
                    const double f = arena.u_[m][Ic];
                    const double fp1 = arena.u_[m][Ip1];
                    const double fm1 = arena.u_[m][Im1];
                    double T =
                        eos.get_temperature(arena.e_[Ic], arena.rhob_[Ic]);
                    double Tp1 =
                        eos.get_temperature(arena.e_[Ip1], arena.rhob_[Ip1]);
                    double Tm1 =
                        eos.get_temperature(arena.e_[Im1], arena.rhob_[Im1]);
                    if (T > T_tol && Tp1 > T_tol && Tm1 > T_tol) {
                        dUoverTsup[m][direction] =
                            (minmod.minmod_dx(fp1 / Tp1, f / T, fm1 / Tm1)
                             / delta[direction]);
                    } else {
                        dUoverTsup[m][direction] = 0.;
                    }
                    dUTsup[m][direction] =
                        (minmod.minmod_dx(fp1 * Tp1, f * T, fm1 * Tm1)
                         / delta[direction]);
                }
            }

            // Sangyong Nov 18 2014
            // Here we make derivatives of muB/T
            // dUsup[rk_flag][4][n] = partial_n (muB/T)
            // partial_x (muB/T) and partial_y (muB/T) first
            if (DATA.turn_on_diff == 1 || DATA.output_vorticity == 1) {
                double f =
                    (eos.get_muB(arena.e_[Ic], arena.rhob_[Ic])
                     / eos.get_temperature(arena.e_[Ic], arena.rhob_[Ic]));
                double fp1 =
                    (eos.get_muB(arena.e_[Ip1], arena.rhob_[Ip1])
                     / eos.get_temperature(arena.e_[Ip1], arena.rhob_[Ip1]));
                double fm1 =
                    (eos.get_muB(arena.e_[Im1], arena.rhob_[Im1])
                     / eos.get_temperature(arena.e_[Im1], arena.rhob_[Im1]));
                dUsup[4][direction] =
                    (minmod.minmod_dx(fp1, f, fm1) / delta[direction]);
            }
        });

    /* for u[0], use u[0]u[0] = 1 + u[i]u[i] */
    /* u[0]_m = u[i]_m (u[i]/u[0]) */
    /* for u[0] */
    for (int n = 1; n < 4; n++) {
        double f = 0.0;
        for (int m = 1; m < 4; m++) {
            // (partial_n u^m) u[m]
            f += dUsup[m][n] * (arena.u_[m][fieldIdx]);
        }
        f /= arena.u_[0][fieldIdx];
        dUsup[0][n] = f;
    }
    return 1;
} /* MakeDSpatial */

int U_derivative::MakeDTau(
    const double tau, const Fields &arenaFieldsPrev,
    const Fields &arenaFieldsCurr, const int fieldIdx) {
    /* this makes dU[m][0] = partial^tau u^m */
    /* note the minus sign at the end because of g[0][0] = -1 */

    const double eps = arenaFieldsCurr.e_[fieldIdx];
    const double rhob = arenaFieldsCurr.rhob_[fieldIdx];
    const double eps_prev = arenaFieldsPrev.e_[fieldIdx];
    const double rhob_prev = arenaFieldsPrev.rhob_[fieldIdx];
    const double T = eos.get_temperature(eps, rhob);
    const double T_prev = eos.get_temperature(eps_prev, rhob_prev);

    for (int m = 0; m < 4; m++) {
        // first order is more stable
        double uCurr = arenaFieldsCurr.u_[m][fieldIdx];
        double uPrev = arenaFieldsPrev.u_[m][fieldIdx];
        double f = (uCurr - uPrev) / DATA.delta_tau;
        dUsup[m][0] = -f;  // g^{00} = -1

        if (DATA.output_vorticity == 1) {
            if (T > T_tol && T_prev > T_tol) {
                double duoverTdtau =
                    (uCurr / T - uPrev / T_prev) / DATA.delta_tau;
                dUoverTsup[m][0] = -duoverTdtau;  // g^{00} = -1
            } else {
                dUoverTsup[m][0] = 0.;
            }
            double duTdtau = (uCurr * T - uPrev * T_prev) / DATA.delta_tau;
            dUTsup[m][0] = -duTdtau;  // g^{00} = -1
        }
    }

    /* I have now partial^tau u^i */
    /* I need to calculate (u^i partial^tau u^i) = u^0 partial^tau u^0 */
    /* u_0 d^0 u^0 + u_m d^0 u^m = 0 */
    /* -u^0 d^0 u^0 + u_m d^0 u^m = 0 */
    /* d^0 u^0 = u_m d^0 u^m/u^0 */
    double f = 0.0;
    for (int m = 1; m < 4; m++) {
        /* (partial_0 u^m) u[m] */
        f += dUsup[m][0] * arenaFieldsCurr.u_[m][fieldIdx];
    }
    dUsup[0][0] = f / arenaFieldsCurr.u_[0][fieldIdx];

    // Sangyong Nov 18 2014
    // Here we make the time derivative of (muB/T)
    int m = 4;
    // first order is more stable backward derivative
    const double muB = eos.get_muB(eps, rhob);
    const double tildemu = muB / T;
    const double muB_prev = eos.get_muB(eps_prev, rhob_prev);
    const double tildemu_prev = muB_prev / T_prev;
    f = (tildemu - tildemu_prev) / (DATA.delta_tau);
    dUsup[m][0] = -f;  // g^{00} = -1
    return 1;
}
