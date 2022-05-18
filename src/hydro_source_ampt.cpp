// Copyright 2019 Chun Shen

#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <memory>

#include "hydro_source_ampt.h"
#include "util.h"

using std::string;

HydroSourceAMPT::HydroSourceAMPT(const InitData &DATA_in) :
    DATA(DATA_in) {
    set_source_tau_min(100.0);
    set_source_tau_max(0.0);
    set_sigma_tau(0.1);
    set_sigma_x  (0.5);
    set_sigma_eta(0.2);
    parton_quench_factor = 1.0;    // no diffusion current from the source
    read_in_AMPT_partons();
}


HydroSourceAMPT::~HydroSourceAMPT() {
    parton_list.clear();
}


//! This function reads in the partons information from the AMPT model
void HydroSourceAMPT::read_in_AMPT_partons() {
    parton_list.clear();
    string AMPT_filename = DATA.initName_AMPT;
    music_message << "hydro_source: "
                  << "read in AMPT parton list from " << AMPT_filename;
    music_message.flush("info");

    string text_string;
    std::ifstream AMPT_file(AMPT_filename.c_str());
    if (!AMPT_file) {
        music_message << "hydro_source::read_in_AMPT_partons: "
                      << "can not open the AMPT file: " << AMPT_filename;
        music_message.flush("error");
        exit(1);
    }

    int n_partons = 0;
    int event_id  = 0;
    int dummy;
    getline(AMPT_file, text_string);
    std::stringstream text_stream1(text_string);
    text_stream1 >> event_id >> dummy >> n_partons;

    // now we read in data
    for (int ipart = 0; ipart < n_partons; ipart++) {
        getline(AMPT_file, text_string);
        std::stringstream text_stream(text_string);
        std::shared_ptr<parton> new_parton(new parton);
        double t_local, z_local, pz_local;
        int pid;
        text_stream >> pid >> new_parton->px >> new_parton->py >> pz_local
                    >> new_parton->mass
                    >> new_parton->x >> new_parton->y >> z_local >> t_local;
        if (t_local < z_local) continue;
        if (std::abs(pid) > 3) continue;

        // Now the parton is inside the light cone
        new_parton->E = sqrt(  new_parton->mass*new_parton->mass
                            + new_parton->px*new_parton->px
                            + new_parton->py*new_parton->py
                            + pz_local*pz_local);
        new_parton->tau      = sqrt(t_local*t_local - z_local*z_local);
        new_parton->eta_s    = 0.5*log( (t_local + z_local)
                                      /(t_local - z_local + Util::small_eps));
        new_parton->rapidity = 0.5*log( (new_parton->E + pz_local)
                                      /(new_parton->E - pz_local));
        double u_perp = (sqrt(  new_parton->px*new_parton->px
                              + new_parton->py*new_parton->py)
                         /new_parton->mass);
        new_parton->rapidity_perp = asinh(u_perp);
        if (pid == 1) {
            // d quark
            new_parton->baryon_number   =  1./3.;
            new_parton->strangness      =  0.0  ;
            new_parton->electric_charge = -1./3.;
        } else if (pid == -1) {
            // anti-d quark
            new_parton->baryon_number   = -1./3.;
            new_parton->strangness      =  0.0  ;
            new_parton->electric_charge =  1./3.;
        } else if (pid == 2) {
            // u quark
            new_parton->baryon_number   =  1./3.;
            new_parton->strangness      =  0.0  ;
            new_parton->electric_charge =  2./3.;
        } else if (pid == -2) {
            // anti-u quark
            new_parton->baryon_number   = -1./3.;
            new_parton->strangness      =  0.0  ;
            new_parton->electric_charge = -2./3.;
        } else if (pid == 3) {
            // s quark
            new_parton->baryon_number   =  1./3.;
            new_parton->strangness      = -1.0  ;
            new_parton->electric_charge = -1./3.;
        } else if (pid == -3) {
            // anti-s quark
            new_parton->baryon_number   = -1./3.;
            new_parton->strangness      =  1.0  ;
            new_parton->electric_charge =  1./3.;
        } else {
            music_message << "pid = " << pid;
            music_message.flush("info");
        }
        parton_list.push_back(new_parton);
        if (get_source_tau_max() < new_parton->tau) {
            set_source_tau_max(new_parton->tau);
        }
        if (get_source_tau_min() > new_parton->tau) {
            set_source_tau_min(new_parton->tau);
        }
    }
    AMPT_file.close();
    music_message << "HydroSourceAMPT:: read in " << parton_list.size() << "/"
                  << n_partons << " partons.";
    music_message.flush("info");
    music_message << "HydroSourceAMPT:: tau_min = " << get_source_tau_min()
                  << " fm.";
    music_message.flush("info");
    music_message << "HydroSourceAMPT:: tau_max = " << get_source_tau_max()
                  << " fm.";
    music_message.flush("info");
}


void HydroSourceAMPT::prepare_list_for_current_tau_frame(
                                                const double tau_local) {
    double dtau = DATA.delta_tau;
    parton_list_current_tau.clear();
    for (auto &it: parton_list) {
        double tau_dis = it->tau - tau_local;
        if (tau_dis > 0. && tau_dis < dtau) {
            parton_list_current_tau.push_back(it);
        }
    }
    music_message << "hydro_source: tau = " << tau_local
                  << " number of source: "
                  << parton_list_current_tau.size();
    music_message.flush("info");
}


void HydroSourceAMPT::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s, 
    const FlowVec &u_mu, EnergyFlowVec &j_mu) const {
    j_mu = {0};
    if (parton_list_current_tau.size() == 0) return;
    
    const double sigma_tau = get_sigma_tau();
    const double sigma_x   = get_sigma_x();
    const double sigma_eta = get_sigma_eta();

    // flow velocity
    const double dtau = DATA.delta_tau;

    const double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
    const double prefactor_tau  = 1./dtau;
    const double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
    const double n_sigma_skip   = 5.;
    const double skip_dis_x     = n_sigma_skip*sigma_x;
    const double skip_dis_eta   = n_sigma_skip*sigma_eta;
    const double exp_tau = 1./tau;

    // AMPT parton sources
    double tau_dis_max = tau - get_source_tau_max();
    if (tau_dis_max < n_sigma_skip*sigma_tau) {
        for (auto &it: parton_list_current_tau) {
            double x_dis = x - it->x;
            if (std::abs(x_dis) > skip_dis_x) continue;

            double y_dis = y - it->y;
            if (std::abs(y_dis) > skip_dis_x) continue;

            double eta_s_dis = eta_s - it->eta_s;
            if (std::abs(eta_s_dis) > skip_dis_eta) continue;

            double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                    /(sigma_x*sigma_x));
            double exp_eta_s = (
                    exp(-eta_s_dis*eta_s_dis/(sigma_eta*sigma_eta)));

            double f_smear = exp_tau*exp_xperp*exp_eta_s;
            double p_perp_sq = it->px*it->px + it->py*it->py;
            double m_perp = sqrt(it->mass*it->mass + p_perp_sq);
            j_mu[0] += m_perp*cosh(it->rapidity - eta_s)*f_smear;
            j_mu[1] += it->px*f_smear;
            j_mu[2] += it->py*f_smear;
            j_mu[3] += m_perp*sinh(it->rapidity - eta_s)*f_smear;
        }
        double norm = DATA.sFactor/Util::hbarc;     // 1/fm^4
        double prefactor = norm*prefactor_tau*prefactor_prep*prefactor_etas;
        j_mu[0] *= prefactor;
        j_mu[1] *= prefactor;
        j_mu[2] *= prefactor;
        j_mu[3] *= prefactor;
    }
}


double HydroSourceAMPT::get_hydro_rhob_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu) const {
    double res = 0.;
    if (parton_list_current_tau.size() == 0) return(res);
    
    const double sigma_tau = get_sigma_tau();
    const double sigma_x   = get_sigma_x();
    const double sigma_eta = get_sigma_eta();

    // flow velocity
    const double gamma_perp_flow  = sqrt(1. + u_mu[1]*u_mu[1] + u_mu[2]*u_mu[2]);
    const double y_perp_flow      = acosh(gamma_perp_flow);
    const double y_long_flow      = asinh(u_mu[3]/gamma_perp_flow) + eta_s;
    const double sinh_y_perp_flow = sinh(y_perp_flow);
    const double dtau             = DATA.delta_tau;

    const double exp_tau        = 1.0/tau;
    const double n_sigma_skip   = 5.;
    const double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
    const double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
    const double prefactor_tau  = 1./dtau;
    const double skip_dis_x     = n_sigma_skip*sigma_x;
    const double skip_dis_eta   = n_sigma_skip*sigma_eta;

    double tau_dis_max = tau - get_source_tau_max();
    if (tau_dis_max < n_sigma_skip*sigma_tau) {
        for (auto &it: parton_list_current_tau) {
            // skip the evaluation if the strings is too far away in the
            // space-time grid
            double x_dis = x - it->x;
            if (std::abs(x_dis) > skip_dis_x) continue;

            double y_dis = y - it->y;
            if (std::abs(y_dis) > skip_dis_x) continue;

            double eta_s_dis = eta_s - it->eta_s;
            if (std::abs(eta_s_dis) > skip_dis_eta) continue;

            double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                    /(sigma_x*sigma_x));
            double exp_eta_s = (
                    exp(-eta_s_dis*eta_s_dis/(sigma_eta*sigma_eta)));
            double f_smear = exp_tau*exp_xperp*exp_eta_s;
            double y_dump = ((1. - parton_quench_factor)*it->rapidity
                             + parton_quench_factor*y_long_flow);
            double y_dump_perp = ((1. - parton_quench_factor)*it->rapidity_perp
                                  + parton_quench_factor*y_perp_flow);
            double p_dot_u = (u_mu[0]
                - tanh(y_dump_perp)*sinh_y_perp_flow/cosh(y_dump - eta_s)
                - tanh(y_dump - eta_s)*u_mu[3]);
            res += p_dot_u*f_smear;
        }
        res *= prefactor_tau*prefactor_prep*prefactor_etas;
    }
    return(res);
}
