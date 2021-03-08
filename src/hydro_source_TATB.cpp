// Copyright 2019 Chun Shen

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <memory>

#include "hydro_source_TATB.h"
#include "util.h"

using std::string;


HydroSourceTATB::HydroSourceTATB(const InitData &DATA_in) :
    DATA_(DATA_in) {
    double tau_overlap = 2.*7./(sinh(DATA_.beam_rapidity));
    tau_source = std::max(DATA_.tau0, tau_overlap);
    set_source_tau_min(tau_source);
    set_source_tau_max(tau_source);
    yL_frac_ = DATA_.yL_frac;
    TA_ = 0.;
    TB_ = 0.;
    read_in_TATB();
}


HydroSourceTATB::~HydroSourceTATB() {
    profile_TA.clear();
    profile_TB.clear();
}


//! This function reads in the spatal information of the nuclear thickness
//! functions
void HydroSourceTATB::read_in_TATB() {
    music_message << "read in TA and TB from " << DATA_.initName_TA
                  << " and " << DATA_.initName_TB;
    music_message.flush("info");

    string text_string;

    std::ifstream TAfile(DATA_.initName_TA.c_str());
    if (!TAfile) {
        music_message << "hydro_source::read_in_TATB: "
                      << "can not open TA file: " << DATA_.initName_TA;
        music_message.flush("error");
        exit(1);
    }
    std::ifstream TBfile(DATA_.initName_TB.c_str());
    if (!TBfile) {
        music_message << "hydro_source::read_in_TATB: "
                      << "can not open TB file: " << DATA_.initName_TB;
        music_message.flush("error");
        exit(1);
    }

    const int nx = DATA_.nx;
    const int ny = DATA_.ny;
    for (int i = 0; i < nx; i++) {
        std::vector<double> TA_temp;
        std::vector<double> TB_temp;
        for (int j = 0; j < ny; j++) {
            double TA, TB;
            TAfile >> TA;
            TBfile >> TB;
            TA_temp.push_back(TA);
            TB_temp.push_back(TB);
            TA_ += TA;
            TB_ += TB;
        }
        profile_TA.push_back(TA_temp);
        profile_TB.push_back(TB_temp);
    }
    TAfile.close();
    TBfile.close();
    TA_ *= DATA_.delta_x*DATA_.delta_y;
    TB_ *= DATA_.delta_x*DATA_.delta_y;
    double N_B = TA_ + TB_;
    double total_energy = DATA_.ecm/2.*N_B;
    music_message << "sqrt{s} = " << DATA_.ecm << " GeV, "
                  << "beam rapidity = " << DATA_.beam_rapidity << ", "
                  << "total energy = " << total_energy << " GeV, "
                  << "N_B = " << N_B;
    music_message.flush("info");

    music_message << "HydroSourceTATB: tau_min = " << get_source_tau_min()
                  << " fm/c.";
    music_message.flush("info");
    music_message << "HydroSourceTATB: tau_max = " << get_source_tau_max()
                  << " fm/c.";
    music_message.flush("info");
    music_message << "Longitudinal velocity fraction yL_frac = " << yL_frac_;
    music_message.flush("info");
}


void HydroSourceTATB::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s,
    const FlowVec &u_mu, EnergyFlowVec &j_mu) const {
    const double dtau = DATA_.delta_tau;
    j_mu = {0};
    if (std::abs((tau - tau_source)) > 1./2.*dtau) return;

    const int ix = static_cast<int>((x + DATA_.x_size/2.)/DATA_.delta_x + 0.1);
    const int iy = static_cast<int>((y + DATA_.y_size/2.)/DATA_.delta_y + 0.1);

    double eta_flat = DATA_.eta_flat;
    if (DATA_.ecm > 2000.) {
        // at LHC energies, we vary eta_flat according to TA + TB
        double total_num_nucleons = 208.*2.;  // for PbPb runs
        double slope = -0.4;
        eta_flat = eta_flat + ((TA_ + TB_)/total_num_nucleons - 0.5)*slope;
    }

    double y_CM = atanh((profile_TA[ix][iy] - profile_TB[ix][iy])
                        /(profile_TA[ix][iy] + profile_TB[ix][iy]
                          + Util::small_eps)
                        *tanh(DATA_.beam_rapidity));
    double y_L = yL_frac_*y_CM;
    double M_inv = ((profile_TA[ix][iy] + profile_TB[ix][iy])
                    *Util::m_N*cosh(DATA_.beam_rapidity)/Util::hbarc);  // [1/fm^3]
    double eta0 = std::min(eta_flat/2.0,
                           std::abs(DATA_.beam_rapidity - (y_CM - y_L)));
    double eta_envelop = eta_profile_plateau(eta_s - (y_CM - y_L), eta0,
                                             DATA_.eta_fall_off);
    double E_norm = tau_source*energy_eta_profile_normalisation(
                                    y_CM, eta0, DATA_.eta_fall_off);
    //double eta_envelop = eta_profile_plateau_frag(eta_s - (y_CM - y_L), eta0,
    //                                              DATA_.eta_fall_off);
    //double E_norm = tau_source*energy_eta_profile_normalisation_numerical(
    //                                y_CM, eta0, DATA_.eta_fall_off);
    //std::cout << "check E_norm1 = " << E_norm1 << ", E_norm = " << E_norm << std::endl;
    double epsilon = M_inv*eta_envelop/E_norm/dtau;  // [1/fm^5]
    j_mu[0] = epsilon*cosh(y_L);  // [1/fm^5]
    j_mu[3] = epsilon*sinh(y_L);  // [1/fm^5]
}


double HydroSourceTATB::get_hydro_rhob_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu) const {
    const double dtau = DATA_.delta_tau;
    double res = 0.;
    if (std::abs((tau - tau_source)) > 1./2.*dtau) return(res);

    const int ix = static_cast<int>((x + DATA_.x_size/2.)/DATA_.delta_x + 0.1);
    const int iy = static_cast<int>((y + DATA_.y_size/2.)/DATA_.delta_y + 0.1);
    double eta_rhob_left  = eta_rhob_left_factor(eta_s);
    double eta_rhob_right = eta_rhob_right_factor(eta_s);
    res = profile_TA[ix][iy]*eta_rhob_right + profile_TB[ix][iy]*eta_rhob_left;  // [1/fm^3]
    res /= dtau;  // [1/fm^4]
    return(res);
}


double HydroSourceTATB::eta_rhob_left_factor(const double eta) const {
    double eta_0       = -std::abs(DATA_.eta_rhob_0);
    double delta_eta_1 = DATA_.eta_rhob_width_1;
    double delta_eta_2 = DATA_.eta_rhob_width_2;
    double norm        = 2./(sqrt(M_PI)*tau_source*(delta_eta_1 + delta_eta_2));
    double exp_arg     = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_1;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_2;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}


double HydroSourceTATB::eta_rhob_right_factor(const double eta) const {
    double eta_0       = std::abs(DATA_.eta_rhob_0);
    double delta_eta_1 = DATA_.eta_rhob_width_1;
    double delta_eta_2 = DATA_.eta_rhob_width_2;
    double norm        = 2./(sqrt(M_PI)*tau_source*(delta_eta_1 + delta_eta_2));
    double exp_arg     = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_2;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_1;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}


double HydroSourceTATB::eta_profile_plateau(
        const double eta, const double eta_0, const double sigma_eta) const {
    // this function return the eta envelope profile for energy density
    // Hirano's plateau + Gaussian fall-off
    double res;
    double exparg1 = (std::abs(eta) - eta_0)/sigma_eta;
    double exparg = exparg1*exparg1/2.0;
    res = exp(-exparg*Util::theta(exparg1));
    return res;
}


double HydroSourceTATB::energy_eta_profile_normalisation(
        const double y_CM, const double eta_0, const double sigma_eta) const {
    // this function returns the normalization of the eta envelope profile
    // for energy density
    //  \int deta eta_profile_plateau(eta - y_CM, eta_0, sigma_eta)*cosh(eta)
    double f1 = (exp(eta_0)*erfc(-sqrt(0.5)*sigma_eta)
                 + exp(-eta_0)*erfc(sqrt(0.5*sigma_eta)));
    double f2 = sqrt(M_PI/2.)*sigma_eta*exp(sigma_eta*sigma_eta/2.);
    double f3 = sinh(eta_0 + y_CM) - sinh(-eta_0 + y_CM);
    double norm = cosh(y_CM)*f2*f1 + f3;
    return(norm);
}


double HydroSourceTATB::eta_profile_plateau_frag(
        const double eta, const double eta_0, const double sigma_eta) const {
    // this function return the eta envelope profile for energy density
    // Hirano's plateau + Gaussian fall-off
    double res;
    double exparg1 = (std::abs(eta) - eta_0)/sigma_eta;
    double exparg = exparg1*exparg1/2.0;
    double eta_fragmentation = std::max(eta_0, DATA_.beam_rapidity - 2.0);
    double exparg_frag1 = (std::abs(eta) - eta_fragmentation)/0.5;
    double exparg_frag = exparg_frag1*exparg_frag1/2.0;
    res = (exp(-exparg*Util::theta(exparg1))
           *exp(-exparg_frag*Util::theta(exparg_frag1)));
    return res;
}


double HydroSourceTATB::energy_eta_profile_normalisation_numerical(
        const double y_CM, const double eta_0, const double sigma_eta) const {
    // this function returns the normalization of the eta envelope profile
    // for energy density by performing numerical integral
    //  \int deta eta_profile_plateau(eta - y_CM, eta_0, sigma_eta)*cosh(eta)
    const int npoints = 200;
    const double eta_max = DATA_.beam_rapidity + 3.0;
    const double deta = 2.*eta_max/(npoints - 1);
    double f_eta = 0.;
    for (int i = 0; i < npoints; i++) {
        double eta_i = - eta_max + i*deta + y_CM;
        f_eta += eta_profile_plateau_frag(eta_i - y_CM, eta_0, sigma_eta)*cosh(eta_i);
    }
    f_eta *= deta;
    return(f_eta);
}
