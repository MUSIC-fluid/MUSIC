// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <iostream>
#include <algorithm>
#include <cmath>
#include "cell.h"
#include "eos.h"
#include "reconst.h"

Reconst::Reconst(const EOS &eosIn, const int echo_level_in, int beastMode) :
    eos(eosIn),
    max_iter(100),
    v_critical(0.563624),
    echo_level(echo_level_in) {
    if (beastMode != 0) {
        abs_err = 1e-8;
        rel_err = 1e-6;
    } else {
        abs_err = 1e-16;
        rel_err = 1e-14;
    }
}


ReconstCell Reconst::ReconstIt_shell(double tau, const TJbVec &tauq_vec,
                                     const ReconstCell &grid_pt) {
    ReconstCell grid_p1;

    TJbVec q_vec;
    for (int i = 0; i < 7; i++) {
        q_vec[i] = tauq_vec[i]/tau;
    }

    int flag = ReconstIt_velocity_Newton(grid_p1, tau, q_vec, grid_pt);

    if (flag == -1) {
        revert_grid(grid_p1, grid_pt);
    } else if (flag == -2) {
        regulate_grid(grid_p1, q_vec[0]);
    }

    return grid_p1;
}


//! This function reverts the grid information back its values
//! at the previous time step
void Reconst::revert_grid(ReconstCell &grid_current,
                          const ReconstCell &grid_prev) const {
    grid_current.e    = grid_prev.e;
    grid_current.rhob = grid_prev.rhob;
    grid_current.rhoq = grid_prev.rhoq;
    grid_current.rhos = grid_prev.rhos;
    grid_current.u    = grid_prev.u;
}


//! reconstruct TJb from q[0] - q[4]
//! reconstruct velocity first for finite mu_B case
//! use Newton's method to solve v and u0
int Reconst::ReconstIt_velocity_Newton(ReconstCell &grid_p, double tau,
                                       const TJbVec &q,
                                       const ReconstCell &grid_pt) {

    double K00 = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    double M   = sqrt(K00);
    double T00 = q[0];

    double J0B = q[4];
    double J0Q = q[5];
    double J0S = q[6];

    if ((T00 < abs_err)) {
        // T^{0\mu} is too small, directly set it to
        // e = abs_err, u^\mu = (1, 0, 0, 0)
        return(-2);
    }

    if (T00 < M) {
        if (echo_level > 9) {
            music_message.warning(
                            "Reconst:: can not find solution! Revert back~");
            music_message << "T00 = " << T00 << ", M = " << M;
            music_message.flush("warning");
        }
        return(-1);
    }

    double u[4], epsilon, pressure, rhob, rhoq, rhos;

    double v_solution = 0.0;
    double v_guess = sqrt(1. - 1./(grid_pt.u[0]*grid_pt.u[0] + abs_err));
    //if (v_guess != v_guess) {
    //    v_guess = 0.0;
    //}
    int v_status = solve_v_Hybrid(v_guess, T00, M, J0B, J0Q, J0S, v_solution);

    if (v_status == 0) {
        return(-1);
    }

    u[0] = 1./(sqrt(1. - v_solution*v_solution) + v_solution*abs_err);
    epsilon = T00 - v_solution*sqrt(K00);
    rhob = J0B/u[0];
    rhoq = J0Q/u[0];
    rhos = J0S/u[0];

    if (std::isnan(epsilon)) {
        std::cout << "reconst: e is nan! "
                  << T00 << ", " << v_solution << ", " << K00 << std::endl;
        exit(1);
    }

    if (v_solution > v_critical && epsilon > 1e-8) {
        // for large velocity, solve u0
        double u0_solution = u[0];
        int u0_status = solve_u0_Hybrid(u[0], T00, K00, M, J0B, J0Q, J0S,
                                        u0_solution);

        if (u0_status == 1) {
            u[0] = u0_solution;
            epsilon = T00 - sqrt((1. - 1./(u0_solution*u0_solution))*K00);
            rhob = J0B/u0_solution;
            rhoq = J0Q/u0_solution;
            rhos = J0S/u0_solution;
        }
    }

    double check_u0_var = std::abs(u[0] - grid_pt.u[0])/grid_pt.u[0];
    if (check_u0_var > 100.) {
        if (grid_pt.e > 1e-6 && echo_level > 2) {
            music_message << "Reconst velocity Newton:: "
                          << "u0 varies more than 100 times compared to "
                          << "its value at previous time step";
            music_message.flush("warning");
            music_message << "e = " << grid_pt.e
                          << ", u[0] = " << u[0]
                          << ", prev_u[0] = " << grid_pt.u[0];
            music_message.flush("warning");
        }
        return(-1);
    }

    grid_p.e = epsilon;
    grid_p.rhob = rhob;
    grid_p.rhoq = rhoq;
    grid_p.rhos = rhos;

    pressure = eos.get_pressure(epsilon, rhob, rhoq, rhos);

    // individual components of velocity
    double velocity_inverse_factor = u[0]/(T00 + pressure);

    u[1] = q[1]*velocity_inverse_factor;
    u[2] = q[2]*velocity_inverse_factor;
    u[3] = q[3]*velocity_inverse_factor;

    // Correcting normalization of 4-velocity
    double u_mag_sq = u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
    if (std::abs(u[0]*u[0] - u_mag_sq - 1.0) > abs_err) {
        double scalef = sqrt((u[0]*u[0] - 1.)/(u_mag_sq + abs_err));
        u[1] *= scalef;
        u[2] *= scalef;
        u[3] *= scalef;
    }

    for (int mu = 0; mu < 4; mu++) {
        grid_p.u[mu] = u[mu];
    }

    return(1);
}


//! This function regulate the grid information
void Reconst::regulate_grid(ReconstCell &grid_cell, double elocal) const {
    grid_cell.e = std::max(abs_err, elocal);
    grid_cell.rhob = 0.0;
    grid_cell.rhoq = 0.0;
    grid_cell.rhos = 0.0;

    grid_cell.u[0] = 1.0;
    grid_cell.u[1] = 0.0;
    grid_cell.u[2] = 0.0;
    grid_cell.u[3] = 0.0;
}


int Reconst::solve_velocity_Newton(const double v_guess, const double T00,
                                   const double M, const double J0B,
                                   const double J0Q, const double J0S,
                                   double &v_solution) {
    int v_status       = 1;
    int iter           = 0;
    double rel_error_v = 10.0;
    double abs_error_v = 10.0;
    double v_next      = v_guess;
    double v_prev      = v_guess;
    double fv, dfdv;
    do {
        iter++;
        reconst_velocity_fdf(v_prev, T00, M, J0B, J0Q, J0S, fv, dfdv);
        v_next = v_prev - (fv/dfdv);
        v_next = std::max(0.0, std::min(1.0, v_next));
        abs_error_v = fv;
        rel_error_v = 2.*abs_error_v/(v_next + v_prev + abs_err);
        v_prev = v_next;
        if (iter > max_iter) {
            v_status = 0;
            break;
        }
    } while (std::abs(abs_error_v) > abs_err
             && std::abs(rel_error_v) > rel_err);

    v_solution = v_next;
    if (v_status == 0 && echo_level > 5) {
        music_message.warning(
                "Reconst velocity Newton:: can not find solution!");
        music_message.warning("output the results at the last iteration:");
        music_message.warning("iter  [lower, upper]  root  err(est)");
        music_message << iter << "   [" << v_prev << ",  " << v_next
                      << "]  " << abs_error_v << "  " << rel_error_v;
        music_message.flush("warning");
    }
    return(v_status);
}


int Reconst::solve_v_Hybrid(const double v_guess, const double T00,
                            const double M, const double J0B,
                            const double J0Q, const double J0S,
                            double &v_solution) {
    int v_status = 1;
    double v_l = std::max(0., v_guess - 0.05);
    double v_h = std::min(1., v_guess + 0.05);
    double fv_l, fv_h;
    reconst_velocity_f(v_l, T00, M, J0B, J0Q, J0S, fv_l);
    reconst_velocity_f(v_h, T00, M, J0B, J0Q, J0S,fv_h);
    if (fv_l*fv_h > 0.) {
        v_l = 0;
        v_h = 1;
        reconst_velocity_f(v_l, T00, M, J0B, J0Q, J0S, fv_l);
        reconst_velocity_f(v_h, T00, M, J0B, J0Q, J0S, fv_h);
        if (fv_l*fv_h > 0.) {
            v_status = 0;
            music_message.error(
                    "Reconst velocity Hybrid:: can not find solution!");
            return(v_status);
        }
    }

    if (std::abs(fv_l) < abs_err) {
        v_solution = v_l;
        return(1);
    }
    if (std::abs(fv_h) < abs_err) {
        v_solution = v_h;
        return(1);
    }

    double dv_prev = v_h - v_l;
    double dv_curr = dv_prev;
    double v_root = (v_h + v_l)/2.;
    double fv = 0;
    double dfdv = 0;
    if (dv_prev > 0.1*v_root) {
        reconst_velocity_fdf(v_root, T00, M, J0B, J0Q, J0S, fv, dfdv);
    } else {
        reconst_velocity_f(v_root, T00, M, J0B, J0Q, J0S, fv);
        dfdv = (fv_h - fv_l)/(v_h - v_l);
    }
    double abs_error_v = 10.0;
    double rel_error_v = 10.0;
    int iter_v = 0;
    do {
        iter_v++;
        double v_prev = v_root;
        double fv_prev = fv;
        if (((v_root - v_h)*dfdv - fv)*((v_root - v_l)*dfdv - fv) > 0.
            || (std::abs(2.*fv) > std::abs(dv_prev*dfdv))) {
            dv_prev = dv_curr;
            dv_curr = (v_h - v_l)/2.;
            v_root  = v_l + dv_curr;
        } else {
            dv_prev = dv_curr;
            dv_curr = fv/dfdv;
            v_root  = v_root - dv_curr;
        }
        abs_error_v = std::abs(dv_curr);
        rel_error_v = abs_error_v/(v_root + abs_err);
        if (std::abs(v_root - v_prev) > 0.1*v_root) {
            reconst_velocity_fdf(v_root, T00, M, J0B, J0Q, J0S, fv, dfdv);
        } else {
            reconst_velocity_f(v_root, T00, M, J0B, J0Q, J0S, fv);
            dfdv = (fv - fv_prev)/(v_root - v_prev + abs_err);
        }
        if (fv*fv_l < 0.) {
            v_h  = v_root;
            fv_h = fv;
        } else {
            v_l  = v_root;
            fv_l = fv;
        }
        if (iter_v > max_iter) {
            v_status = 0;
            break;
        }
        //if (iter_v > 50) {
        //    std::cout << "iter_v = " << iter_v
        //              << ", v = " << v_root
        //              << ", abs_err = " << abs_error_v
        //              << ", rel_err = " << rel_error_v << std::endl;
        //}
    } while (abs_error_v > abs_err && rel_error_v > rel_err);
    v_solution = v_root;

    if (v_status == 0 && echo_level > 5) {
        music_message.warning(
                "Reconst velocity Hybrid:: can not find solution!");
        music_message.warning(
                "output the results at the last iteration:");
        music_message.warning("iter  [lower, upper]  root  err(est)");
        music_message << iter_v << "   [" << v_l << ",  "
                      << v_h << "]  " << abs_error_v << "  "
                      << rel_error_v;
        music_message.flush("warning");
    }
    return(v_status);
}


int Reconst::solve_u0_Newton(const double u0_guess, const double T00,
                             const double K00, const double M, const double J0B,
                             const double J0Q, const double J0S,
                             double &u0_solution) {
    int u0_status = 1;
    double u0_prev = u0_guess;
    double u0_next = u0_prev;
    double abs_error_u0 = 10.0;
    double rel_error_u0 = 10.0;
    double fu0, dfdu0;
    int iter_u0 = 0;
    do {
        iter_u0++;
        reconst_u0_fdf(u0_prev, T00, K00, M, J0B, J0Q, J0S, fu0, dfdu0);
        u0_next = u0_prev - fu0/dfdu0;
        u0_next = std::max(1.0, u0_next);
        abs_error_u0 = fu0;
        rel_error_u0 = 2.*abs_error_u0/(u0_next + u0_prev + abs_err);
        u0_prev = u0_next;
        if (iter_u0 > max_iter) {
            u0_status = 0;
            break;
        }
    } while (std::abs(abs_error_u0) > abs_err
             && std::abs(rel_error_u0) > rel_err);

    u0_solution = u0_next;
    if (u0_status == 0 && echo_level > 5) {
        music_message.warning(
                "Reconst velocity Newton:: can not find solution!");
        music_message.warning(
                "output the results at the last iteration:");
        music_message.warning("iter  [lower, upper]  root  err(est)");
        music_message << iter_u0 << "   [" << u0_prev << ",  "
                      << u0_next << "]  " << abs_error_u0 << "  "
                      << rel_error_u0;
        music_message.flush("warning");
    }
    return(u0_status);
}


int Reconst::solve_u0_Hybrid(const double u0_guess, const double T00,
                             const double K00, const double M,
                             const double J0B,
                             const double J0Q, const double J0S,
                             double &u0_solution) {
    int u0_status = 1;
    double u0_l = std::max(1., 0.99*u0_guess);
    double u0_h = 1.01*u0_guess;
    double fu0_l, fu0_h;
    reconst_u0_f(u0_l, T00, K00, M, J0B, J0Q, J0S, fu0_l);
    reconst_u0_f(u0_h, T00, K00, M, J0B, J0Q, J0S, fu0_h);
    if (fu0_l*fu0_h > 0.) {
        u0_l = 1.0;
        u0_h = 1e5;
        reconst_u0_f(u0_l, T00, K00, M, J0B, J0Q, J0S, fu0_l);
        reconst_u0_f(u0_h, T00, K00, M, J0B, J0Q, J0S, fu0_h);
        if (fu0_l*fu0_h > 0.) {
            if (echo_level > 5) {
                music_message << "Reconst u0 Hybrid:: can not find solution!";
                music_message.flush("error");
                music_message << "u0_guess = " << u0_guess
                              << ", T00 = " << T00 << ", M = " << M
                              << ", J0B = " << J0B << ", J0Q = " << J0Q
                              << ", J0S = " << J0S;
                music_message.flush("error");
                music_message << "u0_l = " << u0_l << ", fu0_l = " << fu0_l
                              << ", u0_h = " << u0_h << ", fu0_h = " << fu0_h;
                music_message.flush("error");
            }
            return(0);
        }
    }

    if (std::abs(fu0_l) < abs_err) {
        u0_solution = u0_l;
        return(1);
    }

    if (std::abs(fu0_h) < abs_err) {
        u0_solution = u0_h;
        return(1);
    }

    double du0_prev = u0_h - u0_l;
    double du0_curr = du0_prev;
    double u0_root = (u0_h + u0_l)/2.;
    double fu0 = 0;
    reconst_u0_f(u0_root, T00, K00, M, J0B, J0Q, J0S, fu0);
    double dfdu0 = (fu0_h - fu0_l)/(u0_h - u0_l);
    double abs_error_u0 = 10.0;
    double rel_error_u0 = 10.0;
    int iter_u0 = 0;
    do {
        iter_u0++;
        double u0_prev = u0_root;
        double fu0_prev = fu0;
        if (((u0_root - u0_h)*dfdu0 - fu0)*((u0_root - u0_l)*dfdu0 - fu0) > 0.
            || (std::abs(2.*fu0) > std::abs(du0_prev*dfdu0))) {
            du0_prev = du0_curr;
            du0_curr = (u0_h - u0_l)/2.;
            u0_root  = u0_l + du0_curr;
        } else {
            du0_prev = du0_curr;
            du0_curr = fu0/dfdu0;
            u0_root  = u0_root - du0_curr;
        }
        if (u0_root != u0_root) {
            u0_status = 0;
            break;
        }
        reconst_u0_f(u0_root, T00, K00, M, J0B, J0Q, J0S, fu0);
        dfdu0 = (fu0 - fu0_prev)/(u0_root - u0_prev);
        abs_error_u0 = du0_curr;
        rel_error_u0 = du0_curr/u0_root;
        if (fu0*fu0_l < 0.) {
            u0_h  = u0_root;
            fu0_h = fu0;
        } else {
            u0_l  = u0_root;
            fu0_l = fu0;
        }
        if (iter_u0 > max_iter) {
            u0_status = 0;
            break;
        }
    } while (   std::abs(abs_error_u0) > abs_err
             && std::abs(rel_error_u0) > rel_err);
    u0_solution = u0_root;
    //std::cout << "iter_u0 = " << iter_u0 << std::endl;

    if (u0_status == 0 && echo_level > 5) {
        music_message.warning(
                "Reconst velocity Hybrid:: can not find solution!");
        music_message.warning(
                "output the results at the last iteration:");
        music_message.warning("iter  [lower, upper]  root  err(est)");
        music_message << iter_u0 << "   [" << u0_l << ",  "
                      << u0_h << "]  " << abs_error_u0 << "  "
                      << rel_error_u0;
        music_message.flush("warning");
    }
    return(u0_status);
}


void Reconst::reconst_velocity_f(const double v, const double T00,
                                 const double M, const double J0B,
                                 const double J0Q, const double J0S,
                                 double &fv) const {
    const double epsilon = T00 - v*M;
    const double temp = sqrt(1. - v*v);
    const double rhob = J0B*temp;
    const double rhoq = J0Q*temp;
    const double rhos = J0S*temp;

    double pressure = eos.get_pressure(epsilon, rhob, rhoq, rhos);

    //fv = v*(T00 + pressure) - M;
    fv = v - M/(T00 + pressure);
}


//void Reconst::reconst_velocity_fdf(const double v, const double T00,
//                                   const double M, const double J0B,
//                                   const double J0Q, const double J0S,
//                                   double &fv, double &dfdv) const {
//    const double epsilon = T00 - v*M;
//    const double temp = sqrt(1. - v*v);
//    const double rhob = J0B*temp;
//
//    const double rhoq = J0Q*temp;
//    const double rhos = J0S*temp;
//
//    double pressure = 0;
//    double dPde = 0;
//    double dPdrhob = 0;
//    double dPdrhoq = 0;
//    double dPdrhos = 0;
//    eos.get_pressure_with_gradients(epsilon, rhob, rhoq, rhos, pressure,
//                                    dPde, dPdrhob, dPdrhoq, dPdrhos);
//
//    const double temp1 = T00 + pressure;
//    const double temp2 = v/std::max(abs_err, temp);
//
//    fv   = v - M/temp1;
//    dfdv = 1. - M/(temp1*temp1)*(M*dPde + J0B*temp2*dPdrhob
//                                 + J0Q*temp2*dPdrhoq + J0S*temp2*dPdrhos);
//}


void Reconst::reconst_velocity_fdf(const double v, const double T00,
                                   const double M, const double J0B,
                                   const double J0Q, const double J0S,
                                   double &fv, double &dfdv) const {
    reconst_velocity_f(v, T00, M, J0B, J0Q, J0S, fv);
    double fv1;
    double v1 = v*0.999;
    reconst_velocity_f(v1, T00, M, J0B, J0Q, J0S, fv1);
    dfdv = (fv - fv1)/(v - v1);
}


void Reconst::reconst_u0_fdf(const double u0, const double T00,
                             const double K00, const double M, const double J0B,
                             const double J0Q, const double J0S,
                             double &fu0, double &dfdu0) const {
    reconst_u0_f(u0, T00, K00, M, J0B, J0Q, J0S, fu0);

    double u02 = u0*1.001;
    double fu02;
    reconst_u0_f(u02, T00, K00, M, J0B, J0Q, J0S, fu02);
    dfdu0 = (fu02 - fu0)/(u02 - u0);
}


void Reconst::reconst_u0_f(const double u0, const double T00,
                           const double K00, const double M, const double J0B,
                           const double J0Q, const double J0S,
                           double &fu0) const {
    const double u0_inv = 1./u0;
    const double v = sqrt(1. - u0_inv*u0_inv);

    const double epsilon = T00 - v*M;
    const double rhob = J0B*u0_inv;
    const double rhoq = J0Q*u0_inv;
    const double rhos = J0S*u0_inv;

    double pressure = eos.get_pressure(epsilon, rhob, rhoq, rhos);
    const double temp1 = 1. - K00/((T00 + pressure)*(T00 + pressure));
    fu0 = u0 - 1./sqrt(temp1);
}


//void Reconst::reconst_u0_fdf(const double u0, const double T00,
//                           const double K00, const double M, const double J0B,
//                           const double J0Q, const double J0S,
//                           double &fu0, double &dfdu0) const {
//    const double u0_inv = 1./u0;
//    const double v = sqrt(1. - u0_inv*u0_inv);
//
//    const double epsilon = T00 - v*M;
//    const double rhob = J0B*u0_inv;
//    const double rhoq = J0Q*u0_inv;
//    const double rhos = J0S*u0_inv;
//
//    const double dedu0 = - M/std::max(abs_err, u0*u0*u0*v);
//    const double drhobdu0 = - J0B/(u0*u0);
//    const double drhoqdu0 = - J0Q/(u0*u0);
//    const double drhosdu0 = - J0S/(u0*u0);
//
//    double pressure = 0;
//    double dPde = 0;
//    double dPdrhob = 0;
//    double dPdrhoq = 0;
//    double dPdrhos = 0;
//    eos.get_pressure_with_gradients(epsilon, rhob, rhoq, rhos, pressure,
//                                    dPde, dPdrhob, dPdrhoq, dPdrhos);
//
//    const double temp1 = std::max(abs_err*abs_err,
//                                  (T00 + pressure)*(T00 + pressure) - K00);
//    const double denorm1 = sqrt(temp1);
//    const double temp = (T00 + pressure)/denorm1;
//
//    fu0    = u0 - temp;
//    dfdu0  = 1. + (dedu0*dPde + drhobdu0*dPdrhob + drhoqdu0*dPdrhoq
//                   + drhosdu0*dPdrhos)*K00/(temp1*denorm1);
//}

