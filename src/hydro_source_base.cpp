// Copyright 2019 Chun Shen

#include "hydro_source_base.h"
#include "data_struct.h"

void HydroSourceBase::get_hydro_energy_source_before_tau(
    const double tau, const double x, const double y, const double eta_s,
    EnergyFlowVec &j_mu) const {

    FlowVec u                   = {0};
    u[0]                        = 1.0;
    EnergyFlowVec j_mu_one_step = {0};

    double tau0 = 0.0;
    double dtau = 0.005;
    int n_tau_steps = static_cast<int>((tau - tau0)/dtau);
    for (int i = 0; i < n_tau_steps; i++) {
        j_mu_one_step = {0};
        const double tau_local = tau0 + (i + 0.5)*dtau;
        get_hydro_energy_source(tau_local, x, y, eta_s, u, j_mu_one_step);
        for (int j = 0; j < 4; j++) {
            j_mu[j] += tau_local*j_mu_one_step[j]*dtau;
        }
    }
    for (int j = 0; j < 4; j++) {
        j_mu[j] /= tau;
    }
}

double HydroSourceBase::get_hydro_rhob_source_before_tau(
                                const double tau, const double x,
                                const double y, const double eta_s) const {
    FlowVec u = {0};
    u[0] = 1.0;

    double res  = 0.;
    double tau0 = 0.0;
    double dtau = 0.005;

    int n_tau_steps = static_cast<int>((tau - tau0)/dtau);
    for (int i = 0; i < n_tau_steps; i++) {
        const double tau_local = tau0 + (i + 0.5)*dtau;
        const double res_local = get_hydro_rhob_source(
                                                tau_local, x, y, eta_s, u);
        res += tau_local*res_local*dtau;
    }

    return(res/tau);
}
