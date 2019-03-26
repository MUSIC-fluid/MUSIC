// Copyright @ Chun Shen

#include "critical_modes.h"
#include "eos.h"
#include "data.h"
#include "data_struct.h"
#include "minmod.h"
#include <fstream>
#include <cmath>
#include <iomanip>
#include <iostream>

CriticalSlowModes::CriticalSlowModes(
        const EOS &eos_in, const InitData &DATA_in) :
    DATA(DATA_in), eos(eos_in), minmod(DATA_in), n_renorm(4) {
}

CriticalSlowModes::~CriticalSlowModes() {
    Qvec.clear();
}


void CriticalSlowModes::InitializeFields_Gubser(SCGrid &arena_current) {
    Qvec.clear();
    const int nQ = 3;
    dQ = 1.0;
    Qvec.resize(nQ);
    Qvec[0] = 1.0;
    Qvec[1] = 2.0;
    Qvec[2] = 5.0;

    std::string input_filename = (
                        "tests/Gubser_flow/Gubser_ideal_phiQ_init_tau_1.dat");
    std::ifstream profile(input_filename.c_str());
    if (!profile.good()) {
        music_message << "CriticalSlowModes::InitializeFields_Gubser: "
                      << "Can not open the initial file: " << input_filename;
        music_message.flush("error");
        exit(1);
    }
    // read the information line
    std::string dummy;
    std::getline(profile, dummy);

    for (int ix = 0; ix < arena_current.nX(); ix++)
    for (int iy = 0; iy < arena_current.nY(); iy++)
    for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
        arena_current(ix, iy, ieta).phi_Q.resize(nQ);
        double x_local, y_local;
        profile >> x_local >> y_local;
        for (int iQ = 0; iQ < nQ; iQ++) {
            double phiQ_tmp = 0.0;
            profile >> phiQ_tmp;
            if (phiQ_tmp < 1e-16) {
                music_message << "ix = " << ix << ", iy = " << iy
                              << "phiQ_init = " << phiQ_tmp << " is too small";
                music_message.flush("warning");
            }
            arena_current(ix, iy, ieta).phi_Q[iQ] = std::max(1e-16, phiQ_tmp);
        }
    }
    profile.close();
    music_message.info(
        "The critical slow modes phiQ are initialized with Gubser solution.");
}


void CriticalSlowModes::InitializeFields(const int nQ, SCGrid &arena_current) {
    Qvec.clear();
    const double Q_min = 0.1;
    const double Q_max = 10.0;
    dQ = (Q_max - Q_min)/(nQ - 1);
    Qvec.resize(nQ);
    for (int i = 0; i < nQ; i++) {
        Qvec[i] = Q_min + i*dQ;
    }
    for (int ix = 0; ix < arena_current.nX(); ix++)
    for (int iy = 0; iy < arena_current.nY(); iy++)
    for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
        arena_current(ix, iy, ieta).phi_Q.resize(nQ);
        const double eps  = arena_current(ix, iy, ieta).epsilon;
        const double rhob = arena_current(ix, iy, ieta).rhob;
        const double xi   = eos.get_correlation_length(eps, rhob);
        for (int iQ = 0; iQ < nQ; iQ++) {
            arena_current(ix, iy, ieta).phi_Q[iQ] = 1.0*(
                compute_phiQ_equilibrium(Qvec[iQ]*xi, eps, rhob));
        }
    }
    music_message.info("The critical slow modes phiQ are initialized.");
}


double CriticalSlowModes::phiQbar_f2(const double x) const {
    return(1./(1. + x*x));
}


double CriticalSlowModes::phiQbar_0(const double e, const double rho_b) const {
    const double c_p = eos.get_dedT(e, rho_b);
    return(c_p/(rho_b*rho_b + 1e-16));
}


//! This function gets the local correlation length
double CriticalSlowModes::get_xi(const double e, const double rho_b) const {
    //const double xi = 2.0;
    const double xi = eos.get_correlation_length(e, rho_b);
    return(xi);
}


//! This function computes the equilibrium value of phi_Q
double CriticalSlowModes::compute_phiQ_equilibrium(
        const double Qxi, const double e, const double rho_b) const {
    const double phi_0 = phiQbar_0(e, rho_b);
    const double f2    = phiQbar_f2(Qxi);
    return(phi_0*f2);
}
    

//! This function computes the relaxation rate for the phiQ fields
double CriticalSlowModes::get_GammaQ(const double Q, const double xi,
                                     const double T, const double eta) const {
    const double Gamma_xi = T/(6.*M_PI*eta*xi*xi*xi);  // [1/fm]
    const double Qxi      = Q*xi;
    const double Kawasaki = 0.75*(
            1. + Qxi*Qxi + (Qxi*Qxi*Qxi - 1./(Qxi + 1e-16))*atan(Qxi + 1e-16));
    const double GammaQ = 2.*Gamma_xi*Kawasaki;  // [1/fm]
    return(GammaQ);
}


//! This function evolves the phi_Q field at ix, iy, ieta by one RK step in time
void CriticalSlowModes::evolve_phiQfields(
        const double tau, SCGrid &arena_prev, SCGrid &arena_current,
        SCGrid &arena_future, const double theta_local,
        const int ix, const int iy, const int ieta,
        const int rk_flag) {
    
    auto grid_pt_prev = &(arena_prev   (ix, iy, ieta));
    auto grid_pt_c    = &(arena_current(ix, iy, ieta));
    auto grid_pt_f    = &(arena_future (ix, iy, ieta));

    const double tau_rk = tau + rk_flag*(DATA.delta_tau);
    DeltaXVec delta = {DATA.delta_tau, DATA.delta_x, DATA.delta_y,
                       DATA.delta_eta*tau_rk};

    for (unsigned int iQ = 0; iQ < Qvec.size(); iQ++) {
        double flux_term = 0.0;
        compute_KTflux(tau_rk, arena_current, theta_local, ix, iy, ieta, iQ,
                       delta, flux_term);
        double tempf = (
                   (1. - rk_flag)*(grid_pt_c->phi_Q[iQ]*grid_pt_c->u[0])
                 + rk_flag*(grid_pt_prev->phi_Q[iQ]*grid_pt_prev->u[0]));
        double source_term = compute_relaxation_source_term(
                            tau_rk, grid_pt_c, grid_pt_prev, iQ, rk_flag);
        tempf += (flux_term + source_term)*(DATA.delta_tau);
        tempf += rk_flag*((grid_pt_c->phi_Q[iQ])*(grid_pt_c->u[0]));
        tempf *= 1./(1. + rk_flag);
        grid_pt_f->phi_Q[iQ] = tempf/(grid_pt_f->u[0]);
    }
}


void CriticalSlowModes::compute_KTflux(
        const double tau, SCGrid &arena, const double theta_local,
        const int ix, const int iy, const int ieta, const int iQ,
        const DeltaXVec delta, double &flux_term) const {
    auto grid_pt = &(arena(ix, iy, ieta));
    flux_term = 0.0;
    
    /* Kurganov-Tadmor for phi_Q */
    /* implement 
      u^mu partial_mu (phi_Q) = (u^mu phi_Q)_{,mu} - u^mu_{,mu} phi_Q
                              = partial_mu (u^\mu phi_Q) + utau phi_Q/tau
                                - phi_Q theta
      partial_mu (u^mu phi_Q) = partial_tau (utau phi_Q)
                                + (1/tau)partial_eta (ueta phi_Q) 
                                + partial_x (ux phi_Q)
                                + partial_y (uy phi_Q)
      So the KT flux is the right hand side of the following equation
      partial_tau (utau phi_Q) = - (1/tau)partial_eta (ueta phi_Q)
                                 - partial_x (ux phi_Q)
                                 - partial_y (uy phi_Q)
                                 - utau phi_Q/tau
                                 + phi_Q theta
      the local velocity is just u_x/u_tau, u_y/u_tau, u_eta/tau/u_tau
      KT flux is given by 
      H_{j+1/2} = (fRph + fLph)/2 - ax(uRph - uLph) 
      Here fRph = ux phiQRph and ax uRph = |ux/utau|_max utau Phi_Q
    */
    Neighbourloop(arena, ix, iy, ieta, NLAMBDAS{
        double g = c.phi_Q[iQ];
        double f = g*c.u[direction];
        g *= c.u[0];

        double gp2 = p2.phi_Q[iQ];
        double fp2 = gp2*p2.u[direction];
        gp2 *= p2.u[0];

        double gp1 = p1.phi_Q[iQ];
        double fp1 = gp1*p1.u[direction];
        gp1 *= p1.u[0];

        double gm1 = m1.phi_Q[iQ];
        double fm1 = gm1*m1.u[direction];
        gm1 *= m1.u[0];

        double gm2 = m2.phi_Q[iQ];
        double fm2 = gm2*m2.u[direction];
        gm2 *= m2.u[0];

        double temp;
        /*  Make u*phi_Q Halfs */
        const double uPhiQphR = fp1 - 0.5*minmod.minmod_dx(fp2, fp1, f);
        temp                  = 0.5*minmod.minmod_dx(fp1, f, fm1);
        const double uPhiQphL = f + temp;
        const double uPhiQmhR = f - temp;
        const double uPhiQmhL = fm1 + 0.5*minmod.minmod_dx(f, fm1, fm2);

        /* just Phi_Q */
        const double PhiQphR = gp1 - 0.5*minmod.minmod_dx(gp2, gp1, g);
        temp                 = 0.5*minmod.minmod_dx(gp1, g, gm1);
        const double PhiQphL = g + temp;
        const double PhiQmhR = g - temp;
        const double PhiQmhL = gm1 + 0.5*minmod.minmod_dx(g, gm1, gm2);

        /* compute flux following Kurganov-Tadmor */
        const double a   = std::abs(c.u[direction])/c.u[0];
        const double am1 = std::abs(m1.u[direction])/m1.u[0];
        const double ap1 = std::abs(p1.u[direction])/p1.u[0];

        double ax = std::max(a, ap1);
        const double HPhiQph = ((uPhiQphR + uPhiQphL)
                                - ax*(PhiQphR - PhiQphL))*0.5;
        ax = std::max(a, am1);
        const double HPhiQmh = ((uPhiQmhR + uPhiQmhL)
                                - ax*(PhiQmhR - PhiQmhL))*0.5;

        const double HPhiQ = (HPhiQph - HPhiQmh)/delta[direction];
        /* make partial_i (u^i Phi_Q) */
        flux_term += -HPhiQ;
    });

    /* add a source term due to the coordinate change to tau-eta */
    flux_term += (grid_pt->phi_Q[iQ])*(-(grid_pt->u[0])/tau + theta_local);
}


double CriticalSlowModes::compute_relaxation_source_term(
        const double tau, Cell_small *grid_pt, Cell_small *grid_pt_prev,
        const int iQ, const int rk_flag) const {

    double epsilon, rhob;
    if (rk_flag == 0) {
        epsilon = grid_pt->epsilon;
        rhob    = grid_pt->rhob;
    } else {
        epsilon = grid_pt_prev->epsilon;
        rhob    = grid_pt_prev->rhob;
    }
    const double temperature = eos.get_temperature(epsilon, rhob);
    const double s_local     = eos.get_entropy(epsilon, rhob);
    const double shear_eta   = std::max(DATA.shear_to_s, 0.08)*s_local;
    const double xi          = get_xi(epsilon, rhob);

    const double phiQ_relax_rate = std::min(
            1./(DATA.delta_tau), get_GammaQ(Qvec[iQ], xi,
                                            temperature, shear_eta));
    const double phiQ_eq = compute_phiQ_equilibrium(Qvec[iQ]*xi, epsilon, rhob);
    double source_term   = - phiQ_relax_rate*(
            (phiQ_eq/(grid_pt->phi_Q[iQ])*(grid_pt->phi_Q[iQ] - phiQ_eq)));
    return(source_term);
}


//! This function computes the renormalizations for EoS and
//! transport coefficients (future)
void CriticalSlowModes::compute_renormalizations(
        const double tau, SCGrid &arena_current,
        const int ix, const int iy, const int ieta) const {
    auto grid_pt_c = &(arena_current(ix, iy, ieta));
    
    const double eps    = grid_pt_c->epsilon;
    const double rhob   = grid_pt_c->rhob;
    const double P0     = eos.get_pressure(eps, rhob);
    const double beta0  = 1./eos.get_temperature(eps, rhob);
    const double xi     = eos.get_correlation_length(eps, rhob);

    const double phase_factor = 1./(2.*2.*M_PI*M_PI)*dQ;
    double Delta_s     = 0.0;
    double Delta_beta  = 0.0;
    double Delta_alpha = 0.0;
    for (unsigned int iQ = 0; iQ < Qvec.size(); iQ++) {
        const double Qxi       = Qvec[iQ]*xi;
        const double phiQ_eq   = compute_phiQ_equilibrium(Qxi, eps, rhob);
        const double ratio     = grid_pt_c->phi_Q[iQ]/phiQ_eq;
        const double Q2        = Qvec[iQ]*Qvec[iQ];

        const double delta_s_i     = entropy_intergrand(ratio);
        const double delta_alpha_i = alpha_intergrand(eps, rhob, Qxi, ratio);
        const double delta_beta_i  = beta_intergrand(eps, rhob, Qxi, ratio);

        Delta_s     += Q2*delta_s_i;
        Delta_alpha += Q2*delta_alpha_i;
        Delta_beta  += Q2*delta_beta_i;
    }
    Delta_s     *= phase_factor;
    Delta_alpha *= phase_factor;
    Delta_beta  *= phase_factor;
    
    double Delta_P = ((Delta_s - (eps + P0)*Delta_beta + rhob*Delta_alpha)
                      /(beta0 + Delta_beta));
}


double CriticalSlowModes::entropy_intergrand(const double x) const {
    double integrand = log(x) + 1. - x;
    return(integrand);
}


double CriticalSlowModes::alpha_intergrand(
                        const double e, const double rhob,
                        const double Qxi, const double phi_ratio) const {
    double dCpde = 0.0;
    double dxide = 0.0;
    double Cp = eos.get_dedT(e, rhob);
    double factor = dCpde/(Cp + 1e-16) - 2.*Qxi/(1. + Qxi*Qxi)*dxide;
    double integrand = factor*(phi_ratio - 1.);
    return(integrand);
}


double CriticalSlowModes::beta_intergrand(
                        const double e, const double rhob,
                        const double Qxi, const double phi_ratio) const {
    double dCpdrhob = 0.0;
    double dxidrhob = 0.0;
    double Cp = eos.get_dedT(e, rhob);
    double factor = (dCpdrhob/(Cp + 1e-16) - 2./(rhob + 1e-16)
                     - 2.*Qxi/(1. + Qxi*Qxi)*dxidrhob);
    double integrand = factor*(phi_ratio - 1.);
    return(integrand);
}
