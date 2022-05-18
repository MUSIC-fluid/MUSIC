// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <iostream>
#include <algorithm>
#include <cmath>
#include "cell.h"
#include "grid.h"
#include "eos.h"
#include "reconst.h"

Reconst::Reconst(const EOS &eosIn, const int echo_level_in) :
    eos(eosIn),
    max_iter(100),
    rel_err(1e-16),
    abs_err(1e-16),
    LARGE(1e20),
    v_critical(0.563624),
    echo_level(echo_level_in) {}



ReconstCell Reconst::ReconstIt_shell(double tau, const TJbVec &tauq_vec,
                                     const Cell_small &grid_pt) {
    ReconstCell grid_p1;

    TJbVec q_vec;
    for (int i = 0; i < 5; i++) {
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
                          const Cell_small &grid_prev) const {
    grid_current.e    = grid_prev.epsilon;
    grid_current.rhob = grid_prev.rhob;
    grid_current.u    = grid_prev.u;
}

//! reconstruct TJb from q[0] - q[4]
//! reconstruct velocity first for finite mu_B case
//! use Newton's method to solve v and u0
int Reconst::ReconstIt_velocity_Newton(ReconstCell &grid_p, double tau,
                                       const TJbVec &q,
                                       const Cell_small &grid_pt) {
    double K00 = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    double M   = sqrt(K00);
    double T00 = q[0];
    double J0  = q[4];
    
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

    double u[4], epsilon, pressure, rhob;
    
    double v_guess = sqrt(1. - 1./(grid_pt.u[0]*grid_pt.u[0] + abs_err));
    if (v_guess != v_guess) {
        v_guess = 0.0;
    }
    double v_solution = 0.0;
    //int v_status = solve_velocity_Newton(v_guess, T00, M, J0, v_solution);
    int v_status = solve_v_Hybrid(v_guess, T00, M, J0, v_solution);
    if (v_status == 0) {
        return(-1);
    }

    u[0] = 1./(sqrt(1. - v_solution*v_solution) + v_solution*abs_err);
    epsilon = T00 - v_solution*sqrt(K00);
    rhob = J0/u[0];
    if (v_solution > v_critical) {
        // for large velocity, solve u0
        double u0_guess    = u[0];
        double u0_solution = u0_guess;
        //int u0_status = solve_u0_Newton(u0_guess, T00, K00, M, J0, u0_solution);
        int u0_status = solve_u0_Hybrid(u0_guess, T00, K00, M, J0, u0_solution);
        if (u0_status == 1) {
            u[0] = u0_solution;
            epsilon = T00 - sqrt((1. - 1./(u0_solution*u0_solution))*K00);
            rhob = J0/u0_solution;
        }
    }

    double check_u0_var = std::abs(u[0] - grid_pt.u[0])/grid_pt.u[0];
    if (check_u0_var > 100.) {
        if (grid_pt.epsilon > 1e-6 && echo_level > 2) {
            music_message << "Reconst velocity Newton:: "
                          << "u0 varies more than 100 times compared to "
                          << "its value at previous time step";
            music_message.flush("warning");
            music_message << "e = " << grid_pt.epsilon
                          << ", u[0] = " << u[0]
                          << ", prev_u[0] = " << grid_pt.u[0];
            music_message.flush("warning");
        }
        return(-1);
    }

    grid_p.e = epsilon;
    grid_p.rhob = rhob;
    pressure = eos.get_pressure(epsilon, rhob);

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
    grid_cell.u[0] = 1.0;
    grid_cell.u[1] = 0.0;
    grid_cell.u[2] = 0.0;
    grid_cell.u[3] = 0.0;
}


int Reconst::solve_velocity_Newton(const double v_guess, const double T00,
                                   const double M, const double J0,
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
        reconst_velocity_fdf(v_prev, T00, M, J0, fv, dfdv);
        v_next = v_prev - (fv/dfdv);
        v_next = std::max(0.0, std::min(1.0, v_next));
        abs_error_v = fv;
        rel_error_v = 2.*abs_error_v/(v_next + v_prev + abs_err);
        v_prev = v_next;
        if (iter > max_iter) {
            v_status = 0;
            break;
        }
    } while (std::abs(abs_error_v) > abs_err && std::abs(rel_error_v) > rel_err);

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
                            const double M, const double J0,
                            double &v_solution) {
    int v_status = 1;
    double v_l = 0.0;
    double v_h = 1.0;
    double fv_l, fv_h;
    double dfdv_l, dfdv_h;
    reconst_velocity_fdf(v_l, T00, M, J0, fv_l, dfdv_l);
    reconst_velocity_fdf(v_h, T00, M, J0, fv_h, dfdv_h);
    if (std::abs(fv_l) < abs_err) {
        v_solution = v_l;
        return(1);
    }
    if (std::abs(fv_h) < abs_err) {
        v_solution = v_h;
        return(1);
    }

    if (fv_l*fv_h > 0.) {
        v_status = 0;
        music_message.error(
                "Reconst velocity Hybrid:: can not find solution!");
        return(v_status);
    }

    double dv_prev = v_h - v_l;
    double dv_curr = dv_prev;
    double v_root = (v_h + v_l)/2.;
    double fv, dfdv;
    reconst_velocity_fdf(v_root, T00, M, J0, fv, dfdv);
    double abs_error_v = 10.0;
    double rel_error_v = 10.0;
    int iter_v = 0;
    do {
        iter_v++;
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
        abs_error_v = dv_curr;
        rel_error_v = abs_error_v/(v_root + abs_err);
        reconst_velocity_fdf(v_root, T00, M, J0, fv, dfdv);
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
    } while (   std::abs(abs_error_v) > abs_err
             && std::abs(rel_error_v) > rel_err);
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
                             const double K00, const double M, const double J0,
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
        reconst_u0_fdf(u0_prev, T00, K00, M, J0, fu0, dfdu0);
        u0_next = u0_prev - fu0/dfdu0;
        u0_next = std::max(1.0, u0_next);
        abs_error_u0 = fu0;
        rel_error_u0 = 2.*abs_error_u0/(u0_next + u0_prev + abs_err);
        u0_prev = u0_next;
        if (iter_u0 > max_iter) {
            u0_status = 0;
            break;
        }
    } while (std::abs(abs_error_u0) > abs_err && std::abs(rel_error_u0) > rel_err);

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
                             const double K00, const double M, const double J0,
                             double &u0_solution) {
    int u0_status = 1;
    double u0_l = 1.0;
    double u0_h = 1e5;
    if (u0_guess >= 1.0) {
        u0_l = std::max(u0_l, 0.5*u0_guess);
        u0_h = std::min(u0_h, 1.5*u0_guess);
    }
    double fu0_l, fu0_h;
    double dfdu0_l, dfdu0_h;
    reconst_u0_fdf(u0_l, T00, K00, M, J0, fu0_l, dfdu0_l);
    reconst_u0_fdf(u0_h, T00, K00, M, J0, fu0_h, dfdu0_h);
    if (std::abs(fu0_l) < abs_err) {
        u0_solution = u0_l;
        return(1);
    }
    if (std::abs(fu0_h) < abs_err) {
        u0_solution = u0_h;
        return(1);
    }

    if (fu0_l*fu0_h > 0.) {
        u0_status = 0;
        music_message << "Reconst u0 Hybrid:: can not find solution!";
        music_message.flush("error");
        music_message << "u0_guess = " << u0_guess << ", T00 = " << T00
                      << ", M = " << M << ", J0 = " << J0;
        music_message.flush("error");
        return(u0_status);
    }

    double du0_prev = u0_h - u0_l;
    double du0_curr = du0_prev;
    double u0_root = (u0_h + u0_l)/2.;
    double fu0, dfdu0;
    reconst_u0_fdf(u0_root, T00, K00, M, J0, fu0, dfdu0);
    double abs_error_u0 = 10.0;
    double rel_error_u0 = 10.0;
    int iter_u0 = 0;
    do {
        iter_u0++;
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
        reconst_u0_fdf(u0_root, T00, K00, M, J0, fu0, dfdu0);
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



void Reconst::reconst_velocity_fdf(const double v, const double T00,
                                   const double M, const double J0,
                                   double &fv, double &dfdv) const {
    const double epsilon = T00 - v*M;
    const double temp    = sqrt(1. - v*v);
    const double rho     = J0*temp;

    const double pressure = eos.get_pressure(epsilon, rho);
    const double temp1    = T00 + pressure;
    const double temp2    = v/temp;
    const double dPde     = eos.get_dpde(epsilon, rho);
    const double dPdrho   = eos.get_dpdrhob(epsilon, rho);

    fv   = v - M/temp1;
    dfdv = 1. - M/(temp1*temp1)*(M*dPde + J0*temp2*dPdrho);
}

void Reconst::reconst_u0_fdf(const double u0, const double T00,
                             const double K00, const double M, const double J0,
                             double &fu0, double &dfdu0) const {
    const double v       = sqrt(1. - 1./(u0*u0));
    const double epsilon = T00 - v*M;
    const double rho     = J0/u0;

    const double dedu0   = - M/(u0*u0*u0*v + abs_err);
    const double drhodu0 = - J0/(u0*u0);

    const double pressure = eos.get_pressure(epsilon, rho);
    const double dPde     = eos.get_dpde(epsilon, rho);
    const double dPdrho   = eos.get_dpdrhob(epsilon, rho);

    const double temp1 = (T00 + pressure)*(T00 + pressure) - K00;
    const double denorm1 = sqrt(temp1);
    const double temp = (T00 + pressure)/denorm1;

    fu0    = u0 - temp;
    dfdu0  = 1. + (dedu0*dPde + drhodu0*dPdrho)*K00/(temp1*denorm1);
}
