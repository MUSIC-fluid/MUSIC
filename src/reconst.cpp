// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <iostream>
#include <algorithm>
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "eos.h"
#include "reconst.h"

Reconst::Reconst(const EOS &eosIn, const InitData &DATA_in) :
    eos(eosIn),
    DATA(DATA_in),
    max_iter(100),
    rel_err(1e-9),
    abs_err(1e-10),
    LARGE(1e20),
    v_critical(0.563624) {
    eos_eps_max = eos.get_eps_max();
    echo_level = DATA.echo_level;
}

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
    grid_current.e = grid_prev.epsilon;
    grid_current.rhob = grid_prev.rhob;
    grid_current.u = grid_prev.u;
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

    if ((T00 < abs_err) || ((T00 - K00/T00) < 0.0)) {
        // can't make Tmunu with this. restore the previous value
        // remember that uq are eigher halfway cells or the final q_next
        // at this point, the original values in grid_pt->TJb are not touched.
        if (echo_level > 9) {
            music_message.warning(
                    "Reconst velocity Newton:: can not find solution!");
            music_message << "T00 = " << T00 << ", K00 = " << K00;
            music_message.flush("warning");
        }
        return(-2);
    }

    double u[4], epsilon, pressure, rhob;
    
    double v_guess = sqrt(1. - 1./(grid_pt.u[0]*grid_pt.u[0] + 1e-15));
    if (v_guess != v_guess) {
        v_guess = 0.0;
    }
    double v_solution = 0.0;
    int v_status = solve_velocity_Newton(v_guess, T00, M, J0, v_solution);
    if (v_status == 0) {
        return(-1);
    }
    
    if (v_solution < v_critical) {
        u[0] = 1./(sqrt(1. - v_solution*v_solution) + v_solution*abs_err);
        epsilon = T00 - v_solution*sqrt(K00);
        rhob = J0/u[0];
    } else {  // for large velocity, solve u0
        double u0_guess = 1./sqrt(1. - v_solution*v_solution);
        double u0_solution = u0_guess;
        int u0_status = solve_u0_Newton(u0_guess, T00, K00, M, J0, u0_solution);
        if (u0_status == 0) {
            return(-1);
        }
        u[0] = u0_solution;
        epsilon = T00 - sqrt((1. - 1./(u0_solution*u0_solution))*K00);
        rhob = J0/u0_solution;
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
    if (fabs(u[0]*u[0] - u_mag_sq - 1.0) > abs_err) {
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
    grid_cell.e = std::max(1e-12, elocal);
    grid_cell.rhob = 0.0;
    grid_cell.u[0] = 1.0;
    grid_cell.u[1] = 0.0;
    grid_cell.u[2] = 0.0;
    grid_cell.u[3] = 0.0;
}


int Reconst::solve_velocity_Newton(const double v_guess, const double T00,
                                   const double M, const double J0,
                                   double &v_solution) {
    int v_status      = 1;
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
        rel_error_v = 2.*abs_error_v/(v_next + v_prev + 1e-15);
        v_prev = v_next;
        if (iter > max_iter) {
            v_status = 0;
            break;
        }
    } while (fabs(abs_error_v) > abs_err && fabs(rel_error_v) > rel_err);

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
        rel_error_u0 = 2.*abs_error_u0/(u0_next + u0_prev + 1e-15);
        u0_prev = u0_next;
        if (iter_u0 > max_iter) {
            u0_status = 0;
            break;
        }
    } while (fabs(abs_error_u0) > abs_err && fabs(rel_error_u0) > rel_err);

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



void Reconst::reconst_velocity_fdf(const double v, const double T00,
                                   const double M, const double J0,
                                   double &fv, double &dfdv) const {
    const double epsilon = T00 - v*M;
    const double temp    = sqrt(1. - v*v);
    const double rho     = J0*temp;

    const double pressure = eos.get_pressure(epsilon, rho);
    const double temp1    = T00 + pressure;
    const double temp2    = v/temp;
    const double dPde     = eos.p_e_func(epsilon, rho);
    const double dPdrho   = eos.p_rho_func(epsilon, rho);

    fv   = v - M/temp1;
    dfdv = 1. - M/(temp1*temp1)*(M*dPde + J0*temp2*dPdrho);
}

void Reconst::reconst_u0_fdf(const double u0, const double T00,
                             const double K00, const double M, const double J0,
                             double &fu0, double &dfdu0) const {
    const double v       = sqrt(1. - 1./(u0*u0));
    const double epsilon = T00 - v*M;
    const double rho     = J0/u0;

    const double dedu0   = - M/(u0*u0*u0*v + 1e-15);
    const double drhodu0 = - J0/(u0*u0);

    const double pressure = eos.get_pressure(epsilon, rho);
    const double dPde     = eos.p_e_func(epsilon, rho);
    const double dPdrho   = eos.p_rho_func(epsilon, rho);

    const double temp1 = (T00 + pressure)*(T00 + pressure) - K00;
    const double denorm1 = sqrt(temp1);
    const double temp = (T00 + pressure)/denorm1;

    fu0    = u0 - temp;
    dfdu0  = 1. + (dedu0*dPde + drhodu0*dPdrho)*K00/(temp1*denorm1);
}
