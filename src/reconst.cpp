// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <iostream>
#include "util.h"
#include "data.h"
#include "grid.h"
#include "eos.h"
#include "reconst.h"

using namespace std;

Reconst::Reconst(EOS *eosIn) {
    eos = eosIn;
    util = new Util;

    // initialize gsl root finding solver
    gsl_rootfinding_max_iter = 100;
    gsl_rootfinding_relerr = 1e-9;
    gsl_rootfinding_abserr = 1e-10;
    gsl_solverType = gsl_root_fsolver_brent;   // Brent-Dekker method
    gsl_rootfinding_solver = gsl_root_fsolver_alloc (gsl_solverType);
}

// destructor
Reconst::~Reconst() {
    delete util;
    gsl_root_fsolver_free (gsl_rootfinding_solver);
}

int Reconst::ReconstIt(Grid *grid_p, int direc, double tau, double **uq,
                       Grid *grid_pt, double eps_init, double rhob_init,
                       InitData *DATA, int rk_flag) {
    /* reconstruct TJb from q[0] - q[4] */
    double K00, T00, J0, u[4], epsilon, p, h, rhob;
    double epsilon_prev, rhob_prev, p_prev, p_guess, temperr;
    double epsilon_next, rhob_next, p_next, err, tempf, temph, cs2;
    double eps_guess, scalef;
    int iter, mu, nu, alpha;
    double q[5];
    const double RECONST_PRECISION = 1e-8;

    /* prepare for the iteration */
    /* uq = qiphL, qiphR, etc 
       qiphL[alpha][direc] means, for instance, TJ[alpha][0] 
       in the cell at x+dx/2 calculated from the left 
       */

    /* uq are the conserved charges. That is, the ones appearing in
       d_tau (Ttautau/tau) + d_eta(Ttaueta/tau) + d_perp(Tperptau) = -Tetaeta
       d_tau (Ttaueta) + d_eta(Tetaeta) + d_v(tau Tveta) = -Ttaueta/tau 
       d_tau (Ttauv) + d_eta(Tetav) + d_w(tau Twv) = 0
       d_tau (Jtau) + d_eta Jeta + d_perp(tau Jperp) = 0
       */

    /* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta
       q[4] = Jtau */
    /* uq = qiphL, qiphR, qimhL, qimhR, qirk */


    for (alpha=0; alpha<5; alpha++) {
        q[alpha] = uq[alpha][direc];
    }

    K00 = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    K00 /= (tau*tau);
 
    T00 = q[0]/tau;
    J0 = q[4]/tau;
 
    if ( (T00 < 0.0) || ((T00 - K00/T00) < 0.0) || (T00 < (SMALL)) ) {
        // can't make Tmunu with this. restore the previous value
        // remember that uq are eigher halfway cells or the final q_next
        // at this point, the original values in grid_pt->TJb are not touched.
        grid_p->epsilon = grid_pt->epsilon;
        grid_p->rhob = grid_pt->rhob;
     
        grid_p->p = grid_pt->p;
        for (mu=0; mu<4; mu++) {
            grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
            grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
 
            for (nu=0; nu<4; nu++) {
                grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
            }/* nu */
        }/* mu */
        return -1;
    } /* if t00-k00/t00 < 0.0 */

    /* Iteration scheme */
 
    cs2 = eos->p_e_func(eps_init, rhob_init);
    eps_guess = GuessEps(T00, K00, cs2);
    epsilon_next = eps_guess;
     
    if(isnan(epsilon_next))
        cout << "problem " << eps_guess << " T00=" << T00 << " K00=" << K00
             << " cs2=" << cs2 << " q[0]=" << q[0] << " uq[0][" 
             << direc << "]=" << uq[0][direc] 
			 << " q[1]=" << q[1] << " q[2]=" << q[2] << endl;
    p_guess = eos->get_pressure(epsilon_next, rhob_init);
    p_next = p_guess;

    /* rhob = J0*sqrt( (eps + p)/(T00 + p) ) */

    if (J0 == 0.0) {
        rhob_next = 0.0;
    } else {
        rhob_next = J0*sqrt((eps_guess + p_guess)/(T00 + p_guess));
        if (!isfinite(rhob_next) || rhob_next < 0) {
            rhob_next = 0.0;
        }
    }
    err = 1.0;
    for (iter=0; iter<100; iter++) {
        if (err < (RECONST_PRECISION)*0.01) {
            if (isnan(epsilon_next))
                cout << "problem2" << endl;
            p_next = eos->get_pressure(epsilon_next, rhob_next);
            break;
        } else {
            epsilon_prev = epsilon_next;
            rhob_prev = rhob_next;
        }
   
        if (isnan(epsilon_prev))
            cout << "problem3" << endl;

        p_prev = eos->get_pressure(epsilon_prev, rhob_prev);
        epsilon_next = T00 - K00/(T00 + p_prev);
        err = 0.0;
   
        if (DATA->turn_on_rhob == 1) {
            rhob_next = J0*sqrt((epsilon_prev + p_prev)/(T00 + p_prev));
            temperr = fabs((rhob_next-rhob_prev)/(rhob_prev+SMALL));
            if (isfinite(temperr))
                err += temperr;
            else if (isinf(temperr))
                err += 1000.0; /* big enough */
            else if (isnan(temperr))
                err += 1000.0;
        } else {
            rhob_next = 0.0;
        }
   
        temperr = fabs((epsilon_next-epsilon_prev)/(epsilon_prev+SMALL));
        if (isfinite(temperr))
            err += temperr;
        else if (isinf(temperr))
            err += 1000.0; /* big enough */
        else if (isnan(temperr))
            err += 1000.0;
    }/* iter */ 

    if (iter == 100) {
        fprintf(stderr, "Reconst didn't converge.\n");
        cout << grid_p->epsilon << endl;
        fprintf(stderr, "Reverting to the previous TJb...\n"); 
        for (mu=0; mu<4; mu++) {
            grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
            grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
        
            for (nu=0; nu<4; nu++) {
                grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
            }/* nu */
        }/* mu */
     
        grid_p->epsilon = grid_pt->epsilon;
        grid_p->p = grid_pt->p;
        grid_p->rhob = grid_pt->rhob;
        return -2;
    } /* if iteration is unsuccessful, revert */

    /* update */
 
    epsilon = grid_p->epsilon = epsilon_next;
    p = grid_p->p = p_next;
    rhob = grid_p->rhob = rhob_next;
    h = p+epsilon;

    /* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta,
       q[4] = Jtau */

    u[0] = sqrt((q[0]/tau + p)/h);
    //remove if for speed
    if (!isfinite(u[0])) {
        u[0] = 1.0;
        u[1] = 0.0;
        u[2] = 0.0;
        u[3] = 0.0;
    } else {
        u[1] = q[1]/tau/h/u[0]; 
        u[2] = q[2]/tau/h/u[0]; 
        u[3] = q[3]/tau/h/u[0]; 
    }

    if (u[0] > cosh(DATA->local_y_max)) {
        fprintf(stderr, "Reconst: u[0] = %e is too large.\n", u[0]);
        if(grid_pt->epsilon > 0.3) {
	        fprintf(stderr, "Reconst: u[0] = %e is too large.\n", u[0]);
	        fprintf(stderr, "epsilon = %e\n", grid_pt->epsilon);
	        fprintf(stderr, "Reverting to the previous TJb...\n"); 
	    }
        for (mu=0; mu<4; mu++) {
	        grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
	        grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
	  
	        for(nu=0; nu<4; nu++) {
	            grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
	        }/* nu */
	    }/* mu */
        grid_p->epsilon = grid_pt->epsilon;
        grid_p->p = grid_pt->p;
        grid_p->rhob = grid_pt->rhob;
      
        return -2;
    }/* if u[0] is too large, revert */

    /* Correcting normalization of 4-velocity */
    temph = u[0]*u[0] - u[1]*u[1] - u[2]*u[2] - u[3]*u[3];
    // Correct velocity when unitarity is not satisfied to numerical accuracy
    // (constant "SMALL")
    if (fabs(temph - 1.0) > SMALL) {
        // If the deviation is too large, exit MUSIC
        if (fabs(temph - 1.0) > 0.1) {
            fprintf(stderr, "In Reconst, reconstructed : u2 = %e\n", temph);
            fprintf(stderr, "Can't happen.\n");
            exit(0);
        } else if(fabs(temph - 1.0) > sqrt(SMALL)) {
            // Warn only when the deviation from 1 is relatively large
            fprintf(stderr, "In Reconst, reconstructed : u2 = %e\n", temph);
            fprintf(stderr, "with u[0] = %e\n", u[0]);
            fprintf(stderr, "Correcting it...\n");
        }   

        // Rescaling spatial components of velocity so that unitarity 
        // is exactly satisfied (u[0] is not modified)
        scalef = (u[0]-1.0);
        scalef *= (u[0]+1.0);
        scalef /= (u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
        scalef = sqrt(scalef);
        u[1] *= scalef;
        u[2] *= scalef;
        u[3] *= scalef;
    }/* if u^mu u_\mu != 1 */
    /* End: Correcting normalization of 4-velocity */

    for (mu=0; mu<4; mu++) {
        tempf = grid_p->TJb[0][4][mu] = rhob*u[mu];
        tempf = grid_p->u[0][mu] = u[mu];
   
        for (nu=0; nu<4; nu++) {
            tempf = grid_p->TJb[0][nu][mu] 
                  = ((epsilon + p)*u[nu]*u[mu] + p*(DATA->gmunu)[nu][mu]);
            if (!isfinite(tempf)) {
                fprintf(stderr, "Update: TJb[0][%d][%d] is %e.\n",
                        nu, mu, grid_p->TJb[0][nu][mu]);
                fprintf(stderr, "Update: epsilon is %e.\n", epsilon);
                exit(0);
            }
        }/* nu */
    }/* mu */
    return 1; /* on successful execution */
}/* Reconst */

int Reconst::ReconstIt_velocity(
    Grid *grid_p, int direc, double tau, double **uq, Grid *grid_pt,
    double eps_init, double rhob_init, InitData *DATA, int rk_flag) {
    /* reconstruct TJb from q[0] - q[4] */
    /* reconstruct velocity first for finite mu_B case
     * (add by C. Shen Nov. 2014) */
    double K00, T00, J0, u[4], epsilon, pressure, rhob;
    int iter, mu, nu, alpha;
    double q[5];

    double v_critical = 0.563624;
    int echo_level = DATA->echo_level;

    /* prepare for the iteration */
    /* uq = qiphL, qiphR, etc 
       qiphL[alpha][direc] means, for instance, TJ[alpha][0] 
       in the cell at x+dx/2 calculated from the left */
    
    /* uq are the conserved charges. That is, the ones appearing in
       d_tau (Ttautau/tau) + d_eta(Ttaueta/tau) + d_perp(Tperptau) = -Tetaeta
       d_tau (Ttaueta) + d_eta(Tetaeta) + d_v(tau Tveta) = -Ttaueta/tau 
       d_tau (Ttauv) + d_eta(Tetav) + d_w(tau Twv) = 0
       d_tau (Jtau) + d_eta Jeta + d_perp(tau Jperp) = 0 */
    
    /* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta
       q[4] = Jtau */
    /* uq = qiphL, qiphR, qimhL, qimhR, qirk */


    for (alpha=0; alpha<5; alpha++)
        q[alpha] = uq[alpha][direc]/tau;

    K00 = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    T00 = q[0];
    J0 = q[4];
 
    if ( (T00 < SMALL) || ((T00 - K00/T00) < 0.0) ) {
        // can't make Tmunu with this. restore the previous value 
        // remember that uq are eigher halfway cells or the final q_next 
        // at this point, the original values in grid_pt->TJb are not touched. 
        grid_p->epsilon = grid_pt->epsilon;
        grid_p->rhob = grid_pt->rhob;
        grid_p->p = grid_pt->p;
        for (mu=0; mu<4; mu++) {
            grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
            grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
 
            for (nu=0; nu<4; nu++) {
                grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
            }/* nu */
        }/* mu */
        return -1;
    }/* if t00-k00/t00 < 0.0 */

    // solving velocity with gsl
    reconst_v_params params;
    params.T00 = T00;
    params.K00 = K00;
    params.J0 = J0;

    CCallbackHolder *Callback_params = new CCallbackHolder;
    Callback_params->cls = this;
    Callback_params->params = &params;

    // solve velocity first
    gslFunc.function = this->CCallback_reconst_v;
    gslFunc.params = Callback_params;
    gsl_root_fsolver_set (gsl_rootfinding_solver, &gslFunc, 0.0, 1.0);
 
    int status;
    iter = 0;
    do {
       iter++;
       status = gsl_root_fsolver_iterate (gsl_rootfinding_solver);
       double x_lo = gsl_root_fsolver_x_lower (gsl_rootfinding_solver);
       double x_hi = gsl_root_fsolver_x_upper (gsl_rootfinding_solver);
       status = gsl_root_test_interval(x_lo, x_hi, gsl_rootfinding_abserr,
                                       gsl_rootfinding_relerr);
    } while (status == GSL_CONTINUE && iter < gsl_rootfinding_max_iter);

    double v_solution;
    if (status == GSL_SUCCESS) {
        v_solution = gsl_root_fsolver_root(gsl_rootfinding_solver);
    } else {
        if (echo_level > 5) {
           double x_lo = gsl_root_fsolver_x_lower (gsl_rootfinding_solver);
           double x_hi = gsl_root_fsolver_x_upper (gsl_rootfinding_solver);
           double result = gsl_root_fsolver_root (gsl_rootfinding_solver);
           fprintf(
                stderr, 
                "***Warning: Reconst velocity:: can not find solution!!!\n");
           fprintf(stderr, "***output the results at the last iteration: \n");
           fprintf(stderr, "%5s [%9s, %9s] %9s %10s \n",
                   "iter", "lower", "upper", "root", "err(est)");
           fprintf(stderr, "%5d [%.7f, %.7f] %.7f %.7f\n", 
                   iter, x_lo, x_hi, result, x_hi - x_lo);
           fprintf(stderr, "Reverting to the previous TJb...\n"); 
        }
        for (mu=0; mu<4; mu++) {
            grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
            grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
            for (nu=0; nu<4; nu++)
                grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
        }
        grid_p->epsilon = grid_pt->epsilon;
        grid_p->p = grid_pt->p;
        grid_p->rhob = grid_pt->rhob;
        return -2;
    }/* if iteration is unsuccessful, revert */
   
    // for large velocity, solve u0
    double u0_solution = 1.0;
    if (v_solution > v_critical) {
        double u0_max_guess = max(5e3, 2*T00/(T00*T00 - K00));
        gslFunc.function = this->CCallback_reconst_u0;
        gslFunc.params = Callback_params;
        gsl_root_fsolver_set (gsl_rootfinding_solver, &gslFunc, 1.0,
                              u0_max_guess);
 
        int status;
        iter = 0;
        do {
            iter++;
            status = gsl_root_fsolver_iterate (gsl_rootfinding_solver);
            double x_lo = gsl_root_fsolver_x_lower (gsl_rootfinding_solver);
            double x_hi = gsl_root_fsolver_x_upper (gsl_rootfinding_solver);
            status = gsl_root_test_interval(x_lo, x_hi, gsl_rootfinding_abserr,
                                            gsl_rootfinding_relerr);
        } while (status == GSL_CONTINUE && iter < gsl_rootfinding_max_iter);

        if (status == GSL_SUCCESS) {
            u0_solution = gsl_root_fsolver_root(gsl_rootfinding_solver);
        } else {
            if (echo_level > 5) {
                double x_lo = gsl_root_fsolver_x_lower(gsl_rootfinding_solver);
                double x_hi = gsl_root_fsolver_x_upper(gsl_rootfinding_solver);
                double result = gsl_root_fsolver_root(gsl_rootfinding_solver);
                fprintf(
                    stderr,
                    "***Warning: Reconst velocity:: can not find solution!\n");
                fprintf(stderr,
                        "***output the results at the last iteration:\n");
                fprintf(stderr, "%5s [%9s, %9s] %9s %10s \n",
                        "iter", "lower", "upper", "root", "err(est)");
                fprintf(stderr, "%5d [%.7f, %.7f] %.7f %.7f\n", 
                        iter, x_lo, x_hi, result, x_hi - x_lo);
                fprintf(stderr, "Reverting to the previous TJb...\n"); 
            }
            for (mu=0; mu<4; mu++) {
                grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
                grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
                for (nu=0; nu<4; nu++)
                    grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
            }
            grid_p->epsilon = grid_pt->epsilon;
            grid_p->p = grid_pt->p;
            grid_p->rhob = grid_pt->rhob;
            return -2;
        }/* if iteration is unsuccessful, revert */
        v_solution = sqrt(1. - 1./(u0_solution*u0_solution));
    }

    // successfully found velocity, now update everything else
    if (v_solution < v_critical) {
        u[0] = 1./(sqrt(1. - v_solution*v_solution) + v_solution*SMALL);
        epsilon = T00 - v_solution*sqrt(K00);
        rhob = J0/u[0];
    } else {
        u[0] = u0_solution;
        epsilon = T00 - sqrt((1. - 1./(u0_solution*u0_solution))*K00);
        rhob = J0/u0_solution;
    }
    grid_p->epsilon = epsilon;
    grid_p->rhob = rhob;

    pressure = eos->get_pressure(epsilon, rhob);
    grid_p->p = pressure;

    // individual components of velocity
    double velocity_inverse_factor = u[0]/(T00 + pressure);

    double u_max = 242582597.70489514; // cosh(20)
    //remove if for speed
    if (!isfinite(u[0])) {
        u[0] = 1.0;
        u[1] = 0.0;
        u[2] = 0.0;
        u[3] = 0.0;
    } else if(u[0] > u_max) {
        // check whether velocity is too large
        if (echo_level > 5) {
            fprintf(stderr, "Reconst velocity: u[0] = %e is too large.\n",
                    u[0]);
            if (grid_pt->epsilon > 0.3) {
	            fprintf(stderr, "Reconst velocity: u[0] = %e is too large.\n",
                        u[0]);
	            fprintf(stderr, "epsilon = %e\n", grid_pt->epsilon);
	            fprintf(stderr, "Reverting to the previous TJb...\n"); 
	        }
        }
        for (mu=0; mu<4; mu++) {
	        grid_p->TJb[0][4][mu] = grid_pt->TJb[rk_flag][4][mu];
	        grid_p->u[0][mu] = grid_pt->u[rk_flag][mu];
	  
	        for (nu=0; nu<4; nu++) {
	            grid_p->TJb[0][nu][mu] = grid_pt->TJb[rk_flag][nu][mu];
	        }/* nu */
	    }/* mu */
        grid_p->epsilon = grid_pt->epsilon;
        grid_p->p = grid_pt->p;
        grid_p->rhob = grid_pt->rhob;

        return -2;
    } else {
        u[1] = q[1]*velocity_inverse_factor; 
        u[2] = q[2]*velocity_inverse_factor; 
        u[3] = q[3]*velocity_inverse_factor; 
    }

    // Correcting normalization of 4-velocity
    double temp_usq = u[0]*u[0] - u[1]*u[1] - u[2]*u[2] - u[3]*u[3];
    // Correct velocity when unitarity is not satisfied to numerical accuracy
    if (fabs(temp_usq - 1.0) > SMALL) {
        // If the deviation is too large, exit MUSIC
        if (fabs(temp_usq - 1.0) > 0.1*u[0]) {
            fprintf(stderr,
                    "In Reconst velocity, reconstructed: u^2 - 1= %e\n",
                    temp_usq - 1.0);
            fprintf(stderr, "Can't happen.\n");
            fprintf(stderr, "u[0]=%.6e, u[1]=%.6e, u[2]=%.6e, u[3]=%.6e\n",
                    u[0], u[1], u[2], u[3]);
            fprintf(stderr, "e=%.6e, rhob=%.6e, p=%.6e\n",
                    epsilon, rhob, pressure);
            exit(0);
        } else if (fabs(temp_usq - 1.0) > sqrt(SMALL)*u[0] && echo_level > 5) {
            // Warn only when the deviation from 1 is relatively large
            fprintf(stderr,
                    "In Reconst velocity, reconstructed: u^2 - 1 = %.8e \n",
                    temp_usq - 1.0);
            double f_res;
            if (v_solution < v_critical)
                f_res = fabs(reconst_velocity_function(v_solution, &params));
            else
                f_res = fabs(reconst_u0_function(u0_solution, &params));
            fprintf(stderr, "with v = %.8e, u[0] = %.8e, res = %.8e \n",
                    v_solution, u[0], f_res);
            fprintf(stderr, "with u[1] = %e\n", u[1]);
            fprintf(stderr, "with u[2] = %e\n", u[2]);
            fprintf(stderr, "with u[3] = %e\n", u[3]);
            fprintf(stderr, "with T00 = %e, K = %e \n", T00, K00);
            fprintf(stderr, "with q1 = %e, q2 = %e, q3 = %e \n",
                    q[1], q[2], q[3]);
            fprintf(stderr, "Correcting it...\n");
        }
        // Rescaling spatial components of velocity so that unitarity 
        // is exactly satisfied (u[0] is not modified)
        double scalef = sqrt(
                (u[0]*u[0] - 1.0)/(u[1]*u[1] + u[2]*u[2] + u[3]*u[3] + SMALL));
        u[1] *= scalef;
        u[2] *= scalef;
        u[3] *= scalef;
    }// if u^mu u_\mu != 1 
    // End: Correcting normalization of 4-velocity
   
    for (mu=0; mu<4; mu++) {
        double tempf;
        grid_p->TJb[0][4][mu] = rhob*u[mu];
        grid_p->u[0][mu] = u[mu];
        for (nu=0; nu<4; nu++) {
            tempf = grid_p->TJb[0][nu][mu]
                  = ((epsilon + pressure)*u[nu]*u[mu]
                          + pressure*(DATA->gmunu)[nu][mu]);
            if (!isfinite(tempf)) {
                fprintf(stderr, "Update: TJb[0][%d][%d] is %e.\n",
                        nu, mu, grid_p->TJb[0][nu][mu]);
                fprintf(stderr, "Update: epsilon is %e.\n", epsilon);
                exit(0);
            }
            grid_p->TJb[0][nu][mu] = tempf;
        }/* nu */
    }/* mu */

    // clean up
    delete Callback_params;

    return 1;  /* on successful execution */
}/* Reconst */

void Reconst::ReconstError(const char *str, int i, int rk_flag, double *qi,
                           double **qi2, Grid *grid_pt) {
    int alpha;
    fprintf(stderr, "Reconst %s in the direction = %d reports an error.\n", 
            str, i); 
    fprintf(stderr, "grid_pt position = (%d, %d, %d).\n", 
            grid_pt->position[1], grid_pt->position[2], grid_pt->position[3]);
    fprintf(stderr, "rk_flag = %d\n", rk_flag); 
 
    for (alpha=0; alpha<5; alpha++) {
        fprintf(stderr, "qi[%d] = %e\n", alpha, qi[alpha]);
    }
    for (alpha=0; alpha<5; alpha++) {
        fprintf(stderr, "qi2[%d][%d] = %e\n", alpha, i, qi2[alpha][i]);
    }
    return;
}/* ReconstErr */

double Reconst::GuessEps(double T00, double K00, double cs2) {
    double f;
 
    if (cs2 < SMALL) {
        f = ((-K00 + util->Power(T00,2))
             *(cs2*K00*util->Power(T00,2) + util->Power(T00,4) 
               + util->Power(cs2,2)*K00*(2*K00 - util->Power(T00,2))))
            /util->Power(T00,5);
    } else {
        f = ((-1.0 + cs2)*T00 
            + sqrt(-4.0*cs2*K00 + util->Power(T00,2)
                   + 2.0*cs2*util->Power(T00,2)
                   + util->Power(cs2,2)*util->Power(T00,2)))/(2.0*cs2);
    }
    return f;
}/*  GuessEps */

double Reconst::reconst_velocity_function(double v, void *params) {
    reconst_v_params *params_ptr = (reconst_v_params *) params;

    double T00 = params_ptr->T00;
    double K00 = params_ptr->K00;
    double J0 = params_ptr->J0;

    double M = sqrt(K00);

    double epsilon = T00 - v*M;
    double rho = J0*sqrt(1 - v*v);
   
    double pressure = eos->get_pressure(epsilon, rho);
    //double f = v*(T00 + pressure) - K00;
    double f = v - M/(T00 + pressure);

    return(f);
}

double Reconst::reconst_u0_function(double u0, void *params) {
    reconst_v_params *params_ptr = (reconst_v_params *) params;

    double T00 = params_ptr->T00;
    double K00 = params_ptr->K00;
    double J0 = params_ptr->J0;

    double M = sqrt(K00);

    double epsilon = T00 - sqrt(1. - 1./u0/u0)*M;
    double rho = J0/u0;
    
    double pressure = eos->get_pressure(epsilon, rho);
    double f = (
        u0 - (T00 + pressure)/sqrt((T00 + pressure)*(T00 + pressure) - K00));

    return(f);
}
