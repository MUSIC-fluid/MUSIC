// Copyright @ Chun Shen

#include "critical_modes.h"
#include "eos.h"
#include "data.h"
#include "data_struct.h"
#include <cmath>

CriticalSlowModes::CriticalSlowModes(const EOS &eos_in, InitData &DATA_in) :
    eos(eos_in), DATA(DATA_in) {
}

CriticalSlowModes::~CriticalSlowModes() {
    Qvec.clear();
}


void CriticalSlowModes::InitializeFields(const int nQ, SCGrid &arena_current) {
    Qvec.clear();
    const double Q_min = 0.1;
    const double Q_max = 1.0;
    const double dQ = (Q_max - Q_min)/(nQ - 1);
    Qvec.resize(nQ);
    for (int i = 0; i < nQ; i++) {
        Qvec[i] = Q_min + i*dQ;
        for (int ix = 0; ix < arena_current.nX(); ix++) {
            for (int iy = 0; iy < arena_current.nY(); iy++) {
                for (int ieta = 0; ieta < arena_current.nEta(); ieta++) {
                    arena_current(ix, iy, ieta).phi_Q.resize(nQ);
                    const double e_local = arena_current(ix, iy, ieta).epsilon;
                    const double rhob_local = arena_current(ix, iy, ieta).rhob;
                    arena_current(ix, iy, ieta).phi_Q[i] = (
                        compute_phiQ_equilibrium(Qvec[i], e_local, rhob_local));
                }
            }
        }
    }
    music_message.info("The critical slow modes phiQ are initialized.");
}


double CriticalSlowModes::phiQbar_f2(const double x) const {
    return(1./(1. + x*x));
}


double CriticalSlowModes::phiQbar_0(const double e, const double rho_b) const {
    const double c_p = 1.0;
    //return(c_p/(rho_b*rho_b + 1e-16));
    return(1.0);
}


//! This function gets the local correlation length
double CriticalSlowModes::get_xi(const double e, const double rho_b) const {
    double xi = 1.0;
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
    const double Gamma_xi = T/(6.*M_PI*eta*xi*xi*xi);
    const double Qxi      = Q*xi;
    const double Kawasaki = 3./4.*(
            1. + Qxi*Qxi + (Qxi*Qxi*Qxi - 1./(Qxi + 1e-16))*atan(Qxi + 1e-16));
    const double GammaQ = 2.*Gamma_xi*Kawasaki;
    return(GammaQ);
}


//! This function evolves the phi_Q field at ix, iy, ieta by one RK step in time
void CriticalSlowModes::evolve_phiQfields(
        const double tau, const int ix, const int iy, const int ieta,
        const int rk_flag, const double T_local, FlowVec umu) {
    const double tau_rk = tau + rk_flag*(DATA.delta_tau);
}
