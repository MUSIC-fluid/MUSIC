// Copyright @ Chun Shen

#include "critical_modes.h"
#include "eos.h"
#include "data.h"
#include "data_struct.h"
#include "minmod.h"
#include <cmath>
#include <iomanip>
#include <iostream>

CriticalSlowModes::CriticalSlowModes(
        const EOS &eos_in, const InitData &DATA_in) :
    DATA(DATA_in), eos(eos_in), minmod(DATA_in) {
}

CriticalSlowModes::~CriticalSlowModes() {
    Qvec.clear();
}


void CriticalSlowModes::InitializeFields(const int nQ, SCGrid &arena_current) {
    Qvec.clear();
    const double Q_min = 0.1;
    const double Q_max = 10.0;
    const double dQ    = (Q_max - Q_min)/(nQ - 1);
    Qvec.resize(nQ);
    for (int i = 0; i < nQ; i++) {
        Qvec[i] = Q_min + i*dQ;
        for (int ix = 0; ix < arena_current.nX(); ix++)
        for (int iy = 0; iy < arena_current.nY(); iy++)
        for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
            arena_current(ix, iy, ieta).phi_Q.resize(nQ);
            const double e_local    = arena_current(ix, iy, ieta).epsilon;
            const double rhob_local = arena_current(ix, iy, ieta).rhob;
            arena_current(ix, iy, ieta).phi_Q[i] = 0.1*(
                compute_phiQ_equilibrium(Qvec[i], e_local, rhob_local));
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
        const double Q, const double e, const double rho_b) const {
    const double phi_0 = phiQbar_0(e, rho_b);
    const double Qxi   = Q*get_xi(e, rho_b);
    const double f2    = phiQbar_f2(Qxi);
    return(phi_0*f2);
}
    

//! This function computes the relaxation rate for the phiQ fields
double CriticalSlowModes::get_GammaQ(const double Q, const double xi,
                                     const double T, const double eta) const {
    const double Gamma_xi = T/(6.*M_PI*eta*xi*xi*xi);  // [1/fm]
    const double Qxi      = Q*xi;
    const double Kawasaki = 3./4.*(
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
    
    auto grid_pt_prev = &(arena_prev(ix, iy, ieta));
    auto grid_pt_c    = &(arena_current(ix, iy, ieta));
    auto grid_pt_f    = &(arena_future(ix, iy, ieta));

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
        tempf += source_term*(DATA.delta_tau);
        tempf += flux_term;
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
    flux_term -= (grid_pt->phi_Q[iQ])*(grid_pt->u[0])/tau;
    flux_term += (grid_pt->phi_Q[iQ])*theta_local;

    flux_term *= delta[0];
}


double CriticalSlowModes::compute_relaxation_source_term(
        const double tau, Cell_small *grid_pt, Cell_small *grid_pt_prev,
        const int iQ, const int rk_flag) const {

    double epsilon, rhob;
    if (rk_flag == 0) {
        epsilon = grid_pt->epsilon;
        rhob = grid_pt->rhob;
    } else {
        epsilon = grid_pt_prev->epsilon;
        rhob = grid_pt_prev->rhob;
    }
    const double temperature = eos.get_temperature(epsilon, rhob);
    const double s_local     = eos.get_entropy(epsilon, rhob);
    const double shear_eta   = std::max(DATA.shear_to_s, 0.08)*s_local;
    const double xi          = get_xi(epsilon, rhob);

    const double phiQ_relax_time = std::max(3.*DATA.delta_tau,
                                            1./(get_GammaQ(Qvec[iQ], xi,
                                                temperature, shear_eta)));
    const double phiQ_eq = compute_phiQ_equilibrium(Qvec[iQ], epsilon, rhob);
    double source_term   = - (
            (phiQ_eq/(grid_pt->phi_Q[iQ])*(grid_pt->phi_Q[iQ] - phiQ_eq))
            /phiQ_relax_time);
    return(source_term);
}
