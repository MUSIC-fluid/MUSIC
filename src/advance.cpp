// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include <omp.h>
#include "./util.h"
#include "./data.h"
#include "./cell.h"
#include "./grid.h"
#include "./reconst.h"
#include "./eos.h"
#include "./evolve.h"
#include "./advance.h"

using namespace std;

Advance::Advance(EOS *eosIn, InitData *DATA_in,
                 hydro_source *hydro_source_in) {
  DATA_ptr         = DATA_in;
  eos              = eosIn;
  util             = new Util;
  reconst_ptr      = new Reconst(eos, DATA_in);
  diss             = new Diss(eosIn, DATA_in);
  minmod           = new Minmod(DATA_in);
  u_derivative_ptr = new U_derivative(eosIn, DATA_in);
  if (DATA_in->Initial_profile == 12 || DATA_in->Initial_profile == 13
      || DATA_in->Initial_profile == 30) {
    flag_add_hydro_source = true;
    hydro_source_ptr = hydro_source_in;
  } else {
        flag_add_hydro_source = false;
        hydro_source_ptr = NULL;
  }
  
  grid_nx   = DATA_in->nx;
  grid_ny   = DATA_in->ny;
  grid_neta = DATA_in->neta;
  rk_order  = DATA_in->rk_order;
}

// destructor
Advance::~Advance() {
  delete util;
  delete diss;
  delete reconst_ptr;
  delete minmod;
  delete u_derivative_ptr;
}


//! this function evolves one Runge-Kutta step in tau
int Advance::AdvanceIt(double tau, InitData *DATA, Grid &arena,
                       int rk_flag) {
  int ieta, ix, iy;
#pragma omp parallel private(ieta, ix, iy)
  {
#pragma omp for collapse(3)
    for (ieta = 0; ieta < grid_neta; ieta++) {
      for (ix = 0; ix <= grid_nx; ix++) {
	for (iy = 0; iy <= grid_ny; iy++) {
	  double eta_s_local = (- DATA_ptr->eta_size/2.
				+ ieta*DATA_ptr->delta_eta);
	  double x_local = (- DATA_ptr->x_size/2.
			    + ix*DATA_ptr->delta_x);
	  double y_local = (- DATA_ptr->y_size/2.
			    + iy*DATA_ptr->delta_y);
	  
	  //		    FirstRKStepT(tau, x_local, y_local, eta_s_local,
	  //             DATA, &(arena(ix,iy,ieta)), rk_flag);
	  FirstRKStepT(tau, x_local, y_local, eta_s_local, DATA, arena, ix, iy, ieta, rk_flag);
          
	  if (DATA->viscosity_flag == 1) {
	    //double tau_rk = tau;
	    //if (rk_flag == 1) {
	    //    tau_rk = tau + DATA_ptr->delta_tau;
	    //}
	    //double theta_local = (
	    //        u_derivative_ptr->calculate_expansion_rate(
	    //            tau_rk, arena, ieta, ix, iy, rk_flag));
	    double theta_local = (
				  u_derivative_ptr->calculate_expansion_rate(
									     tau, arena, ieta, ix, iy, rk_flag));
	    double *a_local = new double[5];
	    double *sigma_local = new double[10];
	    //u_derivative_ptr->calculate_Du_supmu(
	    //        tau_rk, arena, ieta, ix, iy, rk_flag, a_local);
	    //u_derivative_ptr->calculate_velocity_shear_tensor(
	    //        tau_rk, arena, ieta, ix, iy, rk_flag, a_local,
	    //        sigma_local);
	    u_derivative_ptr->calculate_Du_supmu(
						 tau, arena, ieta, ix, iy, rk_flag, a_local);
	    u_derivative_ptr->calculate_velocity_shear_tensor(
							      tau, arena, ieta, ix, iy, rk_flag, a_local,
							      sigma_local);
	    
	    FirstRKStepW(tau, DATA, &(arena(ix,iy,ieta)),
			 rk_flag, theta_local, a_local,
			 sigma_local, ieta, ix, iy);
	    
	    delete[] a_local;
	    delete[] sigma_local;
	  }
	}
      }
    }
#pragma omp barrier
  }
  return(1);
}


/* %%%%%%%%%%%%%%%%%%%%%% First steps begins here %%%%%%%%%%%%%%%%%% */
int Advance::FirstRKStepT(double tau, double x_local, double y_local,
			  double eta_s_local, InitData *DATA, Grid &arena, int ix, int iy, int ieta, int rk_flag) {
  // this advances the ideal part
  double tau_now  = tau;
  double tau_next = tau + (DATA_ptr->delta_tau);
  double tau_rk;
  if (rk_flag == 0) {
    tau_rk = tau_now;
  } else if (rk_flag == 1) {
    tau_rk = tau_next;
  } else {
    fprintf(stderr,"rk_flag = %d out of range.\n", rk_flag);
    exit(0);
  }
  
  // Solve partial_a T^{a mu} = -partial_a W^{a mu}
  // Update T^{mu nu}
  
  // MakeDelatQI gets
  //   qi = q0 if rk_flag = 0 or
  //   qi = q0 + k1 if rk_flag = 1
  // rhs[alpha] is what MakeDeltaQI outputs. 
  // It is the spatial derivative part of partial_a T^{a mu}
  // (including geometric terms)
  double *qi = new double[5];
  MakeDeltaQI(tau_rk, arena, ix, iy, ieta, qi, rk_flag);
  
  double *j_mu = new double[4];
  for (int ii = 0; ii < 4; ii++) {
    j_mu[ii] = 0.0;
  }
  double rhob_source = 0.0;
  if (flag_add_hydro_source) {
    double *u_local = new double[4];
    for (int ii = 0; ii < 4; ii++) {
      //            u_local[ii] = grid_pt->u[rk_flag][ii];
      u_local[ii] = arena(ix,iy,ieta).u[rk_flag][ii];
    }
    hydro_source_ptr->get_hydro_energy_source(
					      tau_rk, x_local, y_local, eta_s_local, u_local, j_mu);
    for (int ii = 0; ii < 4; ii++) {
      j_mu[ii] *= tau_rk;
    }
    if (DATA->turn_on_rhob == 1) {
      rhob_source = tau_rk*hydro_source_ptr->get_hydro_rhob_source(
								   tau_rk, x_local, y_local, eta_s_local, u_local);
    }
    delete[] u_local;
  }
  
  for (int alpha = 0; alpha < 5; alpha++) {
    // now MakeWSource returns partial_a W^{a mu}
    // (including geometric terms) 
    double dwmn = diss->MakeWSource(tau_rk, alpha, &arena(ix,iy,ieta), DATA, rk_flag);
    /* dwmn is the only one with the minus sign */
    qi[alpha] -= dwmn*(DATA->delta_tau);
    
    if (flag_add_hydro_source) {
      // adding hydro_source terms
      if (alpha < 4) {
	qi[alpha] += j_mu[alpha]*DATA->delta_tau;
      } else {
	qi[alpha] += rhob_source*DATA->delta_tau;
      }
    }
    
    // set baryon density back to zero if viscous correction made it
    // non-zero remove/modify if rho_b!=0 
    // - this is only to remove the viscous correction that 
    // can make rho_b negative which we do not want.
    if (DATA->turn_on_rhob == 0) {
      if (alpha == 4 && fabs(qi[alpha]) > 1e-12)
	qi[alpha] = 0.;
    }
    
    /* if rk_flag > 0, we now have q0 + k1 + k2. 
     * So add q0 and multiply by 1/2 */
    if (rk_flag > 0) {
      qi[alpha] += get_TJb(&arena(ix,iy,ieta), 0, alpha, 0)*tau_now;
      qi[alpha] *= 0.5;
    }
  }
  delete[] j_mu;
  
  int flag = 0;
  Cell grid_rk_t;
  //grid_rk_t.u = util->mtx_malloc(1, 4);
  flag = reconst_ptr->ReconstIt_shell(&grid_rk_t, tau_next, qi, &arena(ix,iy,ieta),
				      rk_flag); 
  
  delete[] qi;
  
  UpdateTJbRK(&grid_rk_t, &arena(ix,iy,ieta), rk_flag); 
  //    util->mtx_free(grid_rk_t.u, 1, 4);
  return(flag);
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/*
  Done with T 
  Start W
*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

int Advance::FirstRKStepW(double tau, InitData *DATA, Cell *grid_pt,
                          int rk_flag, double theta_local, double* a_local,
                          double *sigma_local, int ieta, int ix, int iy) {
  double tau_now = tau;
  double tau_next = tau + (DATA->delta_tau);
  
  int trk_flag = rk_flag + 1;
  if (rk_flag == 1) {
    trk_flag = 0;
  }
  
  double **w_rhs;
  w_rhs = new double* [5];
  for (int i = 0; i < 5; i++) {
    w_rhs[i] = new double[4];
    for (int j = 0; j < 4; j++) {
      w_rhs[i][j] = 0.;
    }
  }
  
  // Sangyong Nov 18 2014 implemented mu_max
  int mu_max;
  if (DATA->turn_on_rhob == 1)
    mu_max = 4;
  else 
    mu_max = 3;
  
  // Solve partial_a (u^a W^{mu nu}) = 0
  // Update W^{mu nu}
  // mu = 4 is the baryon current qmu
  
  // calculate delta uWmunu
  // need to use u[0][mu], remember rk_flag = 0 here
  // with the KT flux 
  // solve partial_tau (u^0 W^{kl}) = -partial_i (u^i W^{kl}
  
  /* Advance uWmunu */
  double tempf, temps;
  if (rk_flag == 0) {
    diss->Make_uWRHS(tau_now, grid_pt, w_rhs, DATA, rk_flag,
		     theta_local, a_local);
    for (int mu = 1; mu < 4; mu++) {
      for (int nu = mu; nu < 4; nu++) {
	int idx_1d = util->map_2d_idx_to_1d(mu, nu);
	tempf = ((grid_pt->Wmunu[rk_flag][idx_1d])
		 *(grid_pt->u[rk_flag][0]));
	temps = diss->Make_uWSource(tau_now, grid_pt, mu, nu, DATA,
				    rk_flag, theta_local, a_local,
				    sigma_local); 
	tempf += temps*(DATA->delta_tau);
	tempf += w_rhs[mu][nu];
	grid_pt->Wmunu[trk_flag][idx_1d] = (
					    tempf/(grid_pt->u[trk_flag][0]));
      }
    }
  } else if (rk_flag > 0) {
    diss->Make_uWRHS(tau_next, grid_pt, w_rhs, DATA, rk_flag,
		     theta_local, a_local);
    for (int mu = 1; mu < 4; mu++) {
      for (int nu = mu; nu < 4; nu++) {
	int idx_1d = util->map_2d_idx_to_1d(mu, nu);
	tempf = (grid_pt->Wmunu[0][idx_1d])*(grid_pt->prev_u[0][0]);
	temps = diss->Make_uWSource(tau_next, grid_pt, mu, nu, DATA,
				    rk_flag, theta_local, a_local,
				    sigma_local); 
	tempf += temps*(DATA->delta_tau);
	tempf += w_rhs[mu][nu];
	
	tempf += ((grid_pt->Wmunu[rk_flag][idx_1d])
		  *(grid_pt->u[rk_flag][0]));
	tempf *= 0.5;
	
	grid_pt->Wmunu[trk_flag][idx_1d] = (
                                            tempf/(grid_pt->u[trk_flag][0]));
      }
    }
  } /* rk_flag > 0 */
  
  if (DATA->turn_on_bulk == 1) {
    /* calculate delta u pi */
    double p_rhs;
    if (rk_flag == 0) {
      /* calculate delta u^0 pi */
      diss->Make_uPRHS(tau_now, grid_pt, &p_rhs, DATA, rk_flag,
		       theta_local);
      
      tempf = (grid_pt->pi_b[rk_flag])*(grid_pt->u[rk_flag][0]);
      temps = diss->Make_uPiSource(tau_now, grid_pt, DATA, rk_flag,
				   theta_local, sigma_local);
      tempf += temps*(DATA->delta_tau);
      tempf += p_rhs;
      
      grid_pt->pi_b[trk_flag] = tempf/(grid_pt->u[trk_flag][0]);
    } else if (rk_flag > 0) {
      /* calculate delta u^0 pi */
      diss->Make_uPRHS(tau_next, grid_pt, &p_rhs, DATA, rk_flag,
		       theta_local);
      
      tempf = (grid_pt->pi_b[0])*(grid_pt->prev_u[0][0]);
      temps = diss->Make_uPiSource(tau_next, grid_pt, DATA, rk_flag,
				   theta_local, sigma_local);
      tempf += temps*(DATA->delta_tau);
      tempf += p_rhs;
      
      tempf += (grid_pt->pi_b[rk_flag])*(grid_pt->u[rk_flag][0]);
      tempf *= 0.5;
      
      grid_pt->pi_b[trk_flag] = tempf/(grid_pt->u[trk_flag][0]);
    }
  } else {
    grid_pt->pi_b[trk_flag] = 0.0;
  }
  
  // CShen: add source term for baryon diffusion
  if (DATA->turn_on_diff == 1) {
    if (rk_flag == 0) {
      diss->Make_uqRHS(tau_now, grid_pt, w_rhs, DATA, rk_flag);
      int mu = 4;
      for (int nu = 1; nu < 4; nu++) {
	int idx_1d = util->map_2d_idx_to_1d(mu, nu);
	tempf = ((grid_pt->Wmunu[rk_flag][idx_1d])
		 *(grid_pt->u[rk_flag][0]));
	temps = diss->Make_uqSource(tau_now, grid_pt, nu, DATA,
				    rk_flag, theta_local, a_local,
				    sigma_local); 
	tempf += temps*(DATA->delta_tau);
	tempf += w_rhs[mu][nu];
	
	grid_pt->Wmunu[trk_flag][idx_1d] = (
                                            tempf/(grid_pt->u[trk_flag][0]));
      }
    } else if (rk_flag > 0) {
      diss->Make_uqRHS(tau_next, grid_pt, w_rhs, DATA, rk_flag);
      int mu = 4;
      for (int nu = 1; nu < 4; nu++) {
	int idx_1d = util->map_2d_idx_to_1d(mu, nu);
	tempf = (grid_pt->Wmunu[0][idx_1d])*(grid_pt->prev_u[0][0]);
	temps = diss->Make_uqSource(tau_next, grid_pt, nu, DATA,
				    rk_flag, theta_local, a_local,
				    sigma_local); 
	tempf += temps*(DATA->delta_tau);
	tempf += w_rhs[mu][nu];
	
	tempf += ((grid_pt->Wmunu[rk_flag][idx_1d])
		  *(grid_pt->u[rk_flag][0]));
	tempf *= 0.5;
	
	grid_pt->Wmunu[trk_flag][idx_1d] = (
					    tempf/(grid_pt->u[trk_flag][0]));
      }
    } /* rk_flag > 0 */
  } else {
    for (int nu = 0; nu < 4; nu++) {
      int idx_1d = util->map_2d_idx_to_1d(4, nu);
      grid_pt->Wmunu[trk_flag][idx_1d] = 0.0;
    }
  }
  
  // re-make Wmunu[3][3] so that Wmunu[mu][nu] is traceless
  grid_pt->Wmunu[trk_flag][9] = (
				 (2.*(grid_pt->u[trk_flag][1]*grid_pt->u[trk_flag][2]
				      *grid_pt->Wmunu[trk_flag][5]
				      + grid_pt->u[trk_flag][1]*grid_pt->u[trk_flag][3]
				      *grid_pt->Wmunu[trk_flag][6]
				      + grid_pt->u[trk_flag][2]*grid_pt->u[trk_flag][3]
				      *grid_pt->Wmunu[trk_flag][8])
				  - (grid_pt->u[trk_flag][0]*grid_pt->u[trk_flag][0] 
				     - grid_pt->u[trk_flag][1]*grid_pt->u[trk_flag][1])
				  *grid_pt->Wmunu[trk_flag][4] 
				  - (grid_pt->u[trk_flag][0]*grid_pt->u[trk_flag][0] 
				     - grid_pt->u[trk_flag][2]*grid_pt->u[trk_flag][2])
				  *grid_pt->Wmunu[trk_flag][7])
				 /(grid_pt->u[trk_flag][0]*grid_pt->u[trk_flag][0] 
				   - grid_pt->u[trk_flag][3]*grid_pt->u[trk_flag][3]));
  
  // make Wmunu[i][0] using the transversality
  for (int mu = 1; mu < 4; mu++) {
    tempf = 0.0;
    for (int nu = 1; nu < 4; nu++) {
      int idx_1d = util->map_2d_idx_to_1d(mu, nu);
      tempf += (
                grid_pt->Wmunu[trk_flag][idx_1d]*grid_pt->u[trk_flag][nu]);
    }
    grid_pt->Wmunu[trk_flag][mu] = tempf/(grid_pt->u[trk_flag][0]);
  }
  
  // make Wmunu[0][0]
  tempf = 0.0;
  for (int nu=1; nu<4; nu++)
    tempf += grid_pt->Wmunu[trk_flag][nu]*grid_pt->u[trk_flag][nu]; 
  grid_pt->Wmunu[trk_flag][0] = tempf/(grid_pt->u[trk_flag][0]);
  
  if (DATA->turn_on_diff == 1) {
    // make qmu[0] using transversality
    for (int mu = 4; mu < mu_max + 1; mu++) {
      tempf = 0.0;
      for (int nu = 1; nu < 4; nu++) {
	int idx_1d = util->map_2d_idx_to_1d(mu, nu);
	tempf += (grid_pt->Wmunu[trk_flag][idx_1d]
		  *grid_pt->u[trk_flag][nu]);
      }
      grid_pt->Wmunu[trk_flag][10] = (
				      tempf/(grid_pt->u[trk_flag][0]));
    }
  } else {
    grid_pt->Wmunu[trk_flag][10] = 0.0;
  }
  
  // If the energy density of the fluid element is smaller than 0.01GeV
  // reduce Wmunu using the QuestRevert algorithm
  int revert_flag = 0;
  int revert_q_flag = 0;
  if (DATA->Initial_profile != 0) {
    revert_flag = QuestRevert(tau, grid_pt, rk_flag, DATA, ieta, ix, iy);
    if (DATA->turn_on_diff == 1) {
      revert_q_flag = QuestRevert_qmu(tau, grid_pt, rk_flag, DATA,
				      ieta, ix, iy);
    }
  }
  
  for (int i = 0; i < 5; i++) {
    delete[] w_rhs[i];
  }
  delete[] w_rhs;
  
  if (revert_flag == 1 || revert_q_flag == 1) {
    return(-1);
  } else {
    return(1);
  }
}/* FirstRKStepW */

// update results after RK evolution to grid_pt
void Advance::UpdateTJbRK(Cell *grid_rk, Cell *grid_pt, int rk_flag) {
  int trk_flag = rk_flag + 1;
  if (rk_flag == 1) {
    trk_flag = 0;
  }
  
  if (rk_flag == 0) {
    grid_pt->epsilon_t = grid_rk->epsilon;
    grid_pt->rhob_t = grid_rk->rhob;
  } else {
    grid_pt->epsilon = grid_rk->epsilon;
    grid_pt->rhob = grid_rk->rhob;
  }
  
  // reconstructed grid_rk uses rk_flag 0 only
  for (int mu=0; mu<4; mu++) {
    grid_pt->u[trk_flag][mu] = grid_rk->u[0][mu];
  }
}/* UpdateTJbRK */


//! this function reduce the size of shear stress tensor and bulk pressure
//! in the dilute region to stablize numerical simulations
int Advance::QuestRevert(double tau, Cell *grid_pt, int rk_flag,
                         InitData *DATA, int ieta, int ix, int iy) {
  int revert_flag = 0;
  double eps_scale = 0.5;   // 1/fm^4
  double e_local = grid_pt->epsilon;
  double rhob = grid_pt->rhob;
  
  int trk_flag = rk_flag + 1;
  if (rk_flag == 1) {
    trk_flag = 0;
  }
  
  // regulation factor in the default MUSIC
  //double factor = 300.*tanh(grid_pt->epsilon/eps_scale);
  double xi = 0.05;
  double factor = 100.*(1./(exp(-(e_local - eps_scale)/xi) + 1.)
			- 1./(exp(eps_scale/xi) + 1.));
  double factor_bulk = factor;
  
  // regulation factor from Bjoern
  //double eps_scale = 0.1;  // GeV/fm^3
  //double factor = max(
  //    1e-15, 4.*(1. + tanh(3.*(grid_pt->epsilon - eps_scale/hbarc - 0.5))));
  //double factor_bulk = 2.*factor;
  
  double pi_00 = grid_pt->Wmunu[trk_flag][0];
  double pi_01 = grid_pt->Wmunu[trk_flag][1];
  double pi_02 = grid_pt->Wmunu[trk_flag][2];
  double pi_03 = grid_pt->Wmunu[trk_flag][3];
  double pi_11 = grid_pt->Wmunu[trk_flag][4];
  double pi_12 = grid_pt->Wmunu[trk_flag][5];
  double pi_13 = grid_pt->Wmunu[trk_flag][6];
  double pi_22 = grid_pt->Wmunu[trk_flag][7];
  double pi_23 = grid_pt->Wmunu[trk_flag][8];
  double pi_33 = grid_pt->Wmunu[trk_flag][9];
  
  double pisize = (pi_00*pi_00 + pi_11*pi_11 + pi_22*pi_22 + pi_33*pi_33
		   - 2.*(pi_01*pi_01 + pi_02*pi_02 + pi_03*pi_03)
		   + 2.*(pi_12*pi_12 + pi_13*pi_13 + pi_23*pi_23));
  
  double pi_local = grid_pt->pi_b[trk_flag];
  double bulksize = 3.*pi_local*pi_local;
  
  double p_local = eos->get_pressure(e_local, rhob);
  double eq_size = e_local*e_local + 3.*p_local*p_local;
  
  // In default MUSIC
  double rho_shear = sqrt(pisize/eq_size)/factor; 
  double rho_bulk  = sqrt(bulksize/eq_size)/factor_bulk;
  
  // Bjoern's quest_revert
  //double rho_shear = pisize/eq_size/factor; 
  //double rho_bulk  = bulksize/eq_size/factor_bulk;
  
  // Reducing the shear stress tensor 
  double rho_shear_max = 0.1;
  if (rho_shear > rho_shear_max) {
    if (e_local > eps_scale && DATA_ptr->echo_level > 5) {
      music_message << "ieta = " << ieta << ", ix = " << ix
		    << ", iy = " << iy
		    << ", energy density = " << e_local*hbarc
		    << " GeV/fm^3, shear |pi/(epsilon+3*P)| = "
		    << rho_shear;
      music_message.flush("warning");
    }
    for (int mu = 0; mu < 10; mu++) {
      grid_pt->Wmunu[trk_flag][mu] = (
				      (rho_shear_max/rho_shear)*grid_pt->Wmunu[trk_flag][mu]);
    }
    revert_flag = 1;
  }
  
  // Reducing bulk viscous pressure 
  double rho_bulk_max = 0.1;
  if (rho_bulk > rho_bulk_max) {
    if (e_local > eps_scale && DATA_ptr->echo_level > 5) {
      music_message << "ieta = " << ieta << ", ix = " << ix
		    << ", iy = " << iy
		    << ", energy density = " << e_local*hbarc
		    << " GeV/fm^3, bulk |Pi/(epsilon+3*P)| = "
		    << rho_bulk;
      music_message.flush("warning");
    }
    grid_pt->pi_b[trk_flag] = (
			       (rho_bulk_max/rho_bulk)*grid_pt->pi_b[trk_flag]);
    revert_flag = 1;
  }
  return(revert_flag);
}


//! this function reduce the size of net baryon diffusion current
//! in the dilute region to stablize numerical simulations
int Advance::QuestRevert_qmu(double tau, Cell *grid_pt, int rk_flag,
                             InitData *DATA, int ieta, int ix, int iy) {
  
  int trk_flag = rk_flag + 1;
  if (rk_flag == 1) {
    trk_flag = 0;
  }
  int revert_flag = 0;
  double eps_scale = 0.5;   // in 1/fm^4
  
  double xi = 0.05;
  double factor = 100.*(1./(exp(-(grid_pt->epsilon - eps_scale)/xi) + 1.)
			- 1./(exp(eps_scale/xi) + 1.));
  
  double q_mu_local[4];
  for (int i = 0; i < 4; i++) {
    // copy the value from the grid
    q_mu_local[i] = grid_pt->Wmunu[trk_flag][10+i];
  }
  
  // calculate the size of q^\mu
  double q_size = 0.0;
  for (int i = 0; i < 4; i++) {
    double gfac = (i == 0 ? -1.0 : 1.0);
    q_size += gfac*q_mu_local[i]*q_mu_local[i];
  }
  
  // first check the positivity of q^mu q_mu 
  // (in the conversion of gmn = diag(-+++))
  if (q_size < 0.0) {
    music_message << "Advance::QuestRevert_qmu: q^mu q_mu = " << q_size
		  << " < 0!";
    music_message.flush("warning");
    music_message << "Reset it to zero!!!!";
    music_message.flush("warning");
    for (int i = 0; i < 4; i++) {
      int idx_1d = util->map_2d_idx_to_1d(4, i);
      grid_pt->Wmunu[trk_flag][idx_1d] = 0.0;
    }
    revert_flag = 1;
  }
  
  // reduce the size of q^mu according to rhoB
  double e_local = grid_pt->epsilon;
  double rhob_local = grid_pt->rhob;
  double rho_q = sqrt(q_size/(rhob_local*rhob_local))/factor;
  double rho_q_max = 0.1;
  if (rho_q > rho_q_max) {
    if (e_local > eps_scale && DATA_ptr->echo_level > 5) {
      music_message << "ieta = " << ieta << ", ix = " << ix
		    << ", iy = " << iy
		    << ", energy density = " << e_local*hbarc
		    << "GeV/fm^3"
		    << ", rhob = " << rhob_local << "1/fm^3"
		    << "-- diffusion |q/rhob| = " << rho_q;
      music_message.flush("warning");
    }
    for (int i = 0; i < 4; i++) {
      grid_pt->Wmunu[trk_flag][10+i] = (
					(rho_q_max/rho_q)*q_mu_local[i]);
    }
    revert_flag = 1;
  }
  return(revert_flag);
}


//! This function computes the rhs array. It computes the spatial
//! derivatives of T^\mu\nu using the KT algorithm
void Advance::MakeDeltaQI(double tau, Grid &arena, int ix, int iy, int ieta, double *qi, int rk_flag) {
  double delta[4];
  delta[1] = DATA_ptr->delta_x;
  delta[2] = DATA_ptr->delta_y;
  delta[3] = DATA_ptr->delta_eta;
  
  double rhs[5];
  for (int alpha = 0; alpha < 5; alpha++) 
    {
      qi[alpha] = get_TJb(&arena(ix,iy,ieta), rk_flag, alpha, 0)*tau;
      rhs[alpha] = 0.0; }/* get qi first */
  
  double *qiphL = new double[5];
  double *qiphR = new double[5];
  double *qimhL = new double[5];
  double *qimhR = new double[5];
  
  Cell grid_phL, grid_phR, grid_mhL, grid_mhR;
  for (int direc = 1; direc < 4; direc++)
    {
      double tau_fac = tau;
      int ixp, iyp, ietap, ixp2, iyp2, ietap2;
      int ixm, iym, ietam, ixm2, iym2, ietam2;
      // this could be made more elegant
      if (direc == 1) {
	ixp = ix+1; iyp = iy; ietap = ieta;
	ixm = ix-1; iym = iy; ietam = ieta;
	ixp2 = ix+2; iyp2 = iy; ietap2 = ieta;
	ixm2 = ix-2; iym2 = iy; ietam2 = ieta;
      }
      if (direc == 2) {
	ixp = ix; iyp = iy+1; ietap = ieta;
	ixm = ix; iym = iy-1; ietam = ieta;
	ixp2 = ix; iyp2 = iy+2; ietap2 = ieta;
	ixm2 = ix; iym2 = iy-2; ietam2 = ieta;
      }
      if (direc == 3) {
	tau_fac = 1.0;
	ixp = ix; iyp = iy; ietap = ieta+1;
	ixm = ix; iym = iy; ietam = ieta-1;
	ixp2 = ix; iyp2 = iy; ietap2 = ieta+2;
	ixm2 = ix; iym2 = iy; ietam2 = ieta-2;
     }
      for (int alpha = 0; alpha < 5; alpha++) {
	double gphL = qi[alpha];
	double gphR = (
		       tau*get_TJb(&arena(ixp,iyp,ietap), rk_flag, alpha, 0));
	double gmhL = (
		       tau*get_TJb(&arena(ixm,iym,ietam), rk_flag, alpha, 0));
	double gmhR = qi[alpha];
	double fphL = 0.5*minmod->minmod_dx(gphR, qi[alpha], gmhL);
	double fphR = -0.5*minmod->minmod_dx(
					     tau*get_TJb(&arena(ixp2,iyp2,ietap2), rk_flag, alpha, 0),
					     gphR, qi[alpha]);
	double fmhL = 0.5*minmod->minmod_dx(qi[alpha], gmhL,
					    tau*get_TJb(&arena(ixm2,iym2,ietam2), rk_flag, alpha, 0));
	double fmhR = -0.5*minmod->minmod_dx(gphR, qi[alpha], gmhL);
	qiphL[alpha] = gphL + fphL;
	qiphR[alpha] = gphR + fphR;
	qimhL[alpha] = gmhL + fmhL;
	qimhR[alpha] = gmhR + fmhR;
      }
      // for each direction, reconstruct half-way cells
      // reconstruct e, rhob, and u[4] for half way cells
      int flag = reconst_ptr->ReconstIt_shell(
					      &grid_phL, tau, qiphL, &arena(ix,iy,ieta), 0);
      flag *= reconst_ptr->ReconstIt_shell(
					   &grid_phR, tau, qiphR, &arena(ix,iy,ieta), 0); 
      flag *= reconst_ptr->ReconstIt_shell(
					   &grid_mhL, tau, qimhL, &arena(ix,iy,ieta), 0);
      flag *= reconst_ptr->ReconstIt_shell(
					   &grid_mhR, tau, qimhR, &arena(ix,iy,ieta), 0);
      double aiphL = MaxSpeed(tau, direc, &grid_phL);
      double aiphR = MaxSpeed(tau, direc, &grid_phR);
      double aimhL = MaxSpeed(tau, direc, &grid_mhL);
      double aimhR = MaxSpeed(tau, direc, &grid_mhR);
      
      double aiph = maxi(aiphL, aiphR);
      double aimh = maxi(aimhL, aimhR);
      
      for (int alpha = 0; alpha < 5; alpha++) {
	double FiphL = get_TJb(&grid_phL, 0, alpha, direc)*tau_fac;
	double FiphR = get_TJb(&grid_phR, 0, alpha, direc)*tau_fac;
	double FimhL = get_TJb(&grid_mhL, 0, alpha, direc)*tau_fac;
	double FimhR = get_TJb(&grid_mhR, 0, alpha, direc)*tau_fac;
        
	// KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2})/2
	//                  - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2
	double Fiph = 0.5*((FiphL + FiphR)
			   - aiph*(qiphR[alpha] - qiphL[alpha]));
	double Fimh = 0.5*((FimhL + FimhR)
			   - aimh*(qimhR[alpha] - qimhL[alpha]));
	double DFmmp = (Fimh - Fiph)/delta[direc];
        
	rhs[alpha] += DFmmp*(DATA_ptr->delta_tau);
      }
    }
  
  // geometric terms
  rhs[0] -= get_TJb(&arena(ix,iy,ieta), rk_flag, 3, 3)*DATA_ptr->delta_tau;
  rhs[3] -= get_TJb(&arena(ix,iy,ieta), rk_flag, 3, 0)*DATA_ptr->delta_tau;
  
  for (int i = 0; i < 5; i++) {
    qi[i] += rhs[i];
  }
  
  delete[] qiphL;
  delete[] qiphR;
  delete[] qimhL;
  delete[] qimhR;
  
  // util->mtx_free(grid_phL.u, 1, 4);
  // util->mtx_free(grid_phR.u, 1, 4);
  // util->mtx_free(grid_mhL.u, 1, 4);
  // util->mtx_free(grid_mhR.u, 1, 4);
}/* MakeDeltaQI */


/* Calculate the right-hand-side */
/*
 du/dt = (1/Delta x)(F_imh - F_iph) + D
 F_iph(a) = (1/2)(f_a(qiphR) + f_a(qiphL))- (1/2)aiph(qiphR - qiphL)
 F_imh(a) = (1/2)(f_a(qimhR) + f_a(qimhL))- (1/2)aimh(qimhR - qimhL)
 RK 2nd order Heun's rules:

 k1 = f(t, u);
 u1 = u + k1*h
 k2 = f(t+h, u1);
 u2 = u1 + k2*h = u + (k1+k2)*h;
 u_next = u + (1/2)(k1+k2)*h
        = (1/2)(u + u2);
*/

// determine the maximum signal propagation speed at the given direction
double Advance::MaxSpeed(double tau, int direc, Cell *grid_p) {
    //grid_p = grid_p_h_L, grid_p_h_R, grid_m_h_L, grid_m_h_R
    //these are reconstructed by Reconst which only uses u[0] and TJb[0]
    double utau = (grid_p->u[0][0]);
    double utau2 = utau*utau;
    double ux = fabs((grid_p->u[0][direc]));
    double ux2 = ux*ux;
    double ut2mux2 = utau2 - ux2;
  
    double eps = grid_p->epsilon;
    double rhob = grid_p->rhob;
  
    double vs2 = eos->get_cs2(eps, rhob);

    double den = utau2*(1. - vs2) + vs2;
    double num_temp_sqrt = (ut2mux2 - (ut2mux2 - 1.)*vs2)*vs2;
    double num;
    if (num_temp_sqrt >= 0) {
        num = utau*ux*(1. - vs2) + sqrt(num_temp_sqrt);
    } else {
        double dpde = eos->p_e_func(eps, rhob);
        double p = eos->get_pressure(eps, rhob);
        double h = p+eps;
        if (dpde < 0.001) {
            num = (sqrt(-(h*dpde*h*(dpde*(-1.0 + ut2mux2) - ut2mux2))) 
                   - h*(-1.0 + dpde)*utau*ux);
        } else {
            fprintf(stderr,"WARNING: in MaxSpeed. \n");
            fprintf(stderr, "Expression under sqrt in num=%lf. \n",
                    num_temp_sqrt);
            fprintf(stderr,"at value e=%lf. \n",eps);
            fprintf(stderr,"at value p=%lf. \n",p);
            fprintf(stderr,"at value h=%lf. \n",h);
            fprintf(stderr,"at value rhob=%lf. \n",rhob);
            fprintf(stderr,"at value utau=%lf. \n", utau);
            fprintf(stderr,"at value uk=%lf. \n", ux);
            fprintf(stderr,"at value vs^2=%lf. \n", vs2);
            fprintf(stderr,"at value dpde=%lf. \n", eos->p_e_func(eps, rhob));
            fprintf(stderr,"at value dpdrhob=%lf. \n",
                    eos->p_rho_func(eps, rhob));
            fprintf(stderr, "MaxSpeed: exiting.\n");
            exit(1);
        }
    }
    
    double f = num/(den + 1e-15);
    if (f < 0.0) {
        fprintf(stderr, "SpeedMax = %e\n is negative.\n", f);
        fprintf(stderr, "Can't happen.\n");
        exit(0);
    } else if (f <  ux/utau) {
        if (num != 0.0) {
            if (fabs(f-ux/utau)<0.0001) {
                f = ux/utau;
            } else {
	            fprintf(stderr, "SpeedMax-v = %lf\n", f-ux/utau);
	            fprintf(stderr, "SpeedMax = %e\n is smaller than v = %e.\n",
                        f, ux/utau);
	            fprintf(stderr, "Can't happen.\n");
	            exit(0);
            }
        }
    } else if (f >  1.0) {
        fprintf(stderr, "SpeedMax = %e\n is bigger than 1.\n", f);
        fprintf(stderr, "Can't happen.\n");
        fprintf(stderr, "SpeedMax = num/den, num = %e, den = %e \n", num, den);
        fprintf(stderr, "cs2 = %e \n", vs2);
        f =1.;
        exit(1);
    }
    if (direc == 3)
        f /= tau;

    return f;
}/* MaxSpeed */

double Advance::get_TJb(Cell *grid_p, int rk_flag, int mu, int nu) {
    double rhob = grid_p->rhob;
    if (rk_flag == 1) {
        rhob = grid_p->rhob_t;
    }
    double u_nu = grid_p->u[rk_flag][nu];
    if (mu == 4) {
        double J_nu = rhob*u_nu;
        return(J_nu);
    } else if (mu < 4) {
        double e = grid_p->epsilon;
        if (rk_flag == 1) {
            e = grid_p->epsilon_t;
        }
        double gfac = 0.0;
        double u_mu = 0.0;
        if (mu == nu) {
            u_mu = u_nu;
            if (mu == 0) {
                gfac = -1.0;
            } else {
                gfac = 1.0;
            }
        } else {
            u_mu = grid_p->u[rk_flag][mu];
        }
        double pressure = eos->get_pressure(e, rhob);
        double T_munu = (e + pressure)*u_mu*u_nu + pressure*gfac;
        return(T_munu);
    } else {
        return(0.0);
    }
}
