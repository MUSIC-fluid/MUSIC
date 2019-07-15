// Copyright 2019 Chun Shen

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <memory>

#include "hydro_source_strings.h"
#include "util.h"

using std::string;

HydroSourceStrings::HydroSourceStrings(const InitData &DATA_in) :
    DATA(DATA_in) {
    set_source_tau_min(100.0);
    set_source_tau_max(0.0);
    set_sigma_tau(0.1);
    set_sigma_x  (0.5);
    set_sigma_eta(0.5);
    string_dump_mode     = DATA.string_dump_mode;
    string_quench_factor = DATA.string_quench_factor;
    parton_quench_factor = 1.0;    // no diffusion current from the source
    read_in_QCD_strings_and_partons();
}


HydroSourceStrings::~HydroSourceStrings() {
    QCD_strings_list.clear();
}


//! This function reads in the spatal information of the strings and partons
//! which are produced from the MC-Glauber-LEXUS model
void HydroSourceStrings::read_in_QCD_strings_and_partons() {
    string QCD_strings_filename = DATA.initName;
    music_message << "read in QCD strings list from " << QCD_strings_filename;
    music_message.flush("info");
    string text_string;

    std::ifstream QCD_strings_file(QCD_strings_filename.c_str());
    if (!QCD_strings_file) {
        music_message << "hydro_source::read_in_QCD_strings_and_partons: "
                      << "can not open QCD strings file: "
                      << QCD_strings_filename;
        music_message.flush("error");
        exit(1);
    }

    for (int i = 0; i < 2; i++) {
        // read the header
        getline(QCD_strings_file, text_string);
    }

    // now we read in data
    getline(QCD_strings_file, text_string);
    while (!QCD_strings_file.eof()) {
        std::stringstream text_stream(text_string);
        std::shared_ptr<QCD_string> new_string(new QCD_string);
        text_stream >> new_string->norm >> new_string->m_over_sigma
                    >> new_string->tau_form
                    >> new_string->tau_0 >> new_string->eta_s_0
                    >> new_string->x_perp >> new_string->y_perp
                    >> new_string->eta_s_left >> new_string->eta_s_right
                    >> new_string->y_l >> new_string->y_r
                    >> new_string->frac_l >> new_string->frac_r
                    >> new_string->y_l_i >> new_string->y_r_i
                    >> new_string->eta_s_baryon_left
                    >> new_string->eta_s_baryon_right
                    >> new_string->y_l_baryon
                    >> new_string->y_r_baryon
                    >> new_string->baryon_frac_l;
        if (!text_stream.eof()) {
            // read in the last element
            text_stream >> new_string->baryon_frac_r;
        } else {
            // the string is too short
            music_message << "read_in_QCD_strings_and_partons: "
                          << "the format of file"
                          << QCD_strings_filename << "is wrong~";
            music_message.flush("error");
            exit(1);
        }
        if (!text_stream.eof()) {
            std::string str_temp;
            text_stream >> str_temp;
            if (str_temp.find_first_not_of(' ') != std::string::npos) {
                // the string is too long
                music_message << "read_in_QCD_strings_and_partons: "
                              << "the format of file"
                              << QCD_strings_filename << "is wrong~";
                music_message.flush("error");
                exit(1);
            }
        }

        if (DATA.Initial_profile == 131) {
            new_string->tau_0 = 0.;
            new_string->eta_s_0 = 0.;
        }

        // compute the string end tau
        double temp_factor1 = (new_string->tau_0*new_string->tau_0
                               - new_string->tau_form*new_string->tau_form);
        double temp_factor2 = (new_string->tau_0
                        *cosh(new_string->eta_s_left - new_string->eta_s_0));
        double temp_factor3 = (new_string->tau_0
                    *cosh(new_string->eta_s_right - new_string->eta_s_0));
        double tau_end_left_local = (
            temp_factor2 + sqrt(temp_factor2*temp_factor2 - temp_factor1));
        double tau_end_right_local = (
            temp_factor3 + sqrt(temp_factor3*temp_factor3 - temp_factor1));
        new_string->tau_end_left = tau_end_left_local;
        new_string->tau_end_right = tau_end_right_local;

        // compute the baryon number tau
        temp_factor2 = (new_string->tau_0*cosh(
                        new_string->eta_s_baryon_left - new_string->eta_s_0));
        temp_factor3 = (new_string->tau_0*cosh(
                        new_string->eta_s_baryon_right - new_string->eta_s_0));
        new_string->tau_baryon_left = (
            temp_factor2 + sqrt(temp_factor2*temp_factor2 - temp_factor1));
        new_string->tau_baryon_right = (
            temp_factor3 + sqrt(temp_factor2*temp_factor2 - temp_factor1));
        
        // determine the tau_start and eta_s_start of the string
        if (new_string->eta_s_left > new_string->eta_s_0) {
            new_string->tau_start = tau_end_left_local;
            new_string->eta_s_start = new_string->eta_s_left;
        } else if (new_string->eta_s_right < new_string->eta_s_0) {
            new_string->tau_start = tau_end_right_local;
            new_string->eta_s_start = new_string->eta_s_right;
        } else {
            new_string->tau_start = new_string->tau_0 + new_string->tau_form;
            new_string->eta_s_start = new_string->eta_s_0;
        }

        // read in one string properly
        QCD_strings_list.push_back(new_string);

        // record the proper time of the first and last string sources
        double source_tau = new_string->tau_form;
        if (new_string->tau_end_left > new_string->tau_end_right) {
            source_tau = new_string->tau_end_left;
        } else {
            source_tau = new_string->tau_end_right;
        }

        if (get_source_tau_max() < source_tau) {
            set_source_tau_max(source_tau);
        }

        if (get_source_tau_min() > new_string->tau_start) {
            set_source_tau_min(new_string->tau_start);
        }

        getline(QCD_strings_file, text_string);
    }
    QCD_strings_file.close();
    music_message << "HydroSourceStrings: tau_min = " << get_source_tau_min()
                  << " fm/c.";
    music_message.flush("info");
    music_message << "HydroSourceStrings: tau_max = " << get_source_tau_max()
                  << " fm/c.";
    music_message.flush("info");
    
    double total_baryon_number = 0;
    for (auto const& it : QCD_strings_list) {
        total_baryon_number += it->baryon_frac_l + it->baryon_frac_r;
    }
    music_message << "total baryon number = " << total_baryon_number;
    music_message.flush("info");
    compute_norm_for_strings();
}


void HydroSourceStrings::compute_norm_for_strings() {
    const int neta              = 500;
    const double eta_range      = 12.;
    const double deta           = 2.*eta_range/(neta - 1);

    const double sigma_eta = get_sigma_eta();
    const double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);

    double E_string_total = 0.0;
    double E_baryon_total = 0.0;
    for (auto &it: QCD_strings_list) {
        double E_string_norm   = 0.;
        double E_baryon_L_norm = 0.;
        double E_baryon_R_norm = 0.;
        for (int ieta = 0; ieta < neta; ieta++) {
            double eta_local = - eta_range + ieta*deta;
            double f_eta = (it->frac_l
                + (it->frac_l - it->frac_r)/(it->eta_s_left - it->eta_s_right)
                  *(eta_local - it->eta_s_left));
            double y_eta = (it->y_l
                + (it->y_l - it->y_r)/(it->eta_s_left - it->eta_s_right)
                  *(eta_local - it->eta_s_left));

            double expon_left  = (it->eta_s_left  - eta_local)/sigma_eta;
            double expon_right = (it->eta_s_right - eta_local)/sigma_eta;
            double e_eta = 0.5*(- erf(expon_left) + erf(expon_right));
            E_string_norm += f_eta*e_eta*cosh(y_eta);

            double e_baryon_L = exp(-expon_left*expon_left);
            double e_baryon_R = exp(-expon_right*expon_right);
            E_baryon_L_norm += e_baryon_L*cosh(eta_local);
            E_baryon_R_norm += e_baryon_R*cosh(eta_local);
        }
        E_string_norm  *= prefactor_etas*deta;
        double E_string = (  it->frac_l*cosh(it->y_l_i)
                           + it->frac_r*cosh(it->y_r_i)
                           - it->frac_l*cosh(it->y_l)
                           - it->frac_r*cosh(it->y_r));
        it->norm = E_string/(E_string_norm + 1e-16);
        E_string_total += E_string;

        // here the E_norm is for the energy of remnants at the string ends
        // frac_l and frac_r should be used here
        E_baryon_L_norm *= it->frac_l*prefactor_etas*deta;
        E_baryon_R_norm *= it->frac_r*prefactor_etas*deta;
        double E_baryon_L   = it->frac_l*cosh(it->y_l);
        double E_baryon_R   = it->frac_r*cosh(it->y_r);
        it->E_baryon_norm_L = E_baryon_L/(E_baryon_L_norm + 1e-16);
        it->E_baryon_norm_R = E_baryon_R/(E_baryon_R_norm + 1e-16);
        E_baryon_total += E_baryon_L + E_baryon_R;
    }
    music_message << "E_total = "
                  << (E_string_total + E_baryon_total)*DATA.sFactor << " GeV. "
                  << "E_string_total = " << E_string_total*DATA.sFactor
                  << " GeV" << ", E_baryon_total = "
                  << E_baryon_total*DATA.sFactor << " GeV.";
    music_message.flush("info");
}


void HydroSourceStrings::prepare_list_for_current_tau_frame(
                                                    const double tau_local) {
    double dtau = DATA.delta_tau;
    QCD_strings_list_current_tau.clear();
    QCD_strings_baryon_list_current_tau.clear();
    QCD_strings_remnant_list_current_tau.clear();
    for (auto &it: QCD_strings_list) {
        if ((   it->tau_baryon_left >= (tau_local - 1./2.*dtau)
             && it->tau_baryon_left <  (tau_local + 3./2.*dtau)
             && it->baryon_frac_l > 0.)
            || (   it->tau_baryon_right >= (tau_local - 1./2.*dtau)
                && it->tau_baryon_right <  (tau_local + 3./2.*dtau)
                && it->baryon_frac_r > 0.)
            ) {
            QCD_strings_baryon_list_current_tau.push_back(it);
        }
        if ((   it->tau_end_left >= (tau_local - 1./2.*dtau)
             && it->tau_end_left <  (tau_local + 3./2.*dtau))
            || (   it->tau_end_right >= (tau_local - 1./2.*dtau)
                && it->tau_end_right <  (tau_local + 3./2.*dtau))) {
            QCD_strings_remnant_list_current_tau.push_back(it);
        }
        if (   it->tau_start <= tau_local + 3./2.*dtau
            && std::max(it->tau_end_left, it->tau_end_right) >= tau_local - dtau/2.) {
            QCD_strings_list_current_tau.push_back(it);
        }
    }
    music_message << "hydro_source: tau = " << tau_local << " fm."
                  << " number of strings for energy density: "
                  << QCD_strings_list_current_tau.size()
                  << " number of strings for net baryon density: "
                  << QCD_strings_baryon_list_current_tau.size();
    music_message.flush("info");
}


void HydroSourceStrings::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s, 
    const FlowVec &u_mu, EnergyFlowVec &j_mu) const {
    j_mu = {0};
    if (   QCD_strings_list_current_tau.size() == 0
        && QCD_strings_remnant_list_current_tau.size() == 0) return;

    const double sigma_x   = get_sigma_x();
    const double sigma_eta = get_sigma_eta();

    // flow velocity
    const double gamma_perp_flow = sqrt(1. + u_mu[1]*u_mu[1] + u_mu[2]*u_mu[2]);
    const double y_perp_flow     = acosh(gamma_perp_flow);
    const double y_long_flow     = asinh(u_mu[3]/gamma_perp_flow) + eta_s;
    const double sin_phi_flow    = u_mu[1]/gamma_perp_flow;
    const double cos_phi_flow    = u_mu[2]/gamma_perp_flow;
    const double dtau            = DATA.delta_tau;

    const double n_sigma_skip   = 5.;
    const double skip_dis_x     = n_sigma_skip*sigma_x;
    const double skip_dis_eta   = n_sigma_skip*sigma_eta;
    const double sfactor        = DATA.sFactor/Util::hbarc;
    const double exp_tau = 1./tau;
    for (auto const&it: QCD_strings_list_current_tau) {
        // energy source from strings
        const double tau_0     = it.lock()->tau_0;
        const double delta_tau = it.lock()->tau_form;
        
        double x_dis = x - it.lock()->x_perp;
        if (std::abs(x_dis) > skip_dis_x) continue;

        double y_dis = y - it.lock()->y_perp;
        if (std::abs(y_dis) > skip_dis_x) continue;

        // calculate the crossed string segments in the eta direction
        // normally, there will be two segments
        // [eta_L_next, eta_L] and [eta_R, eta_R_next]
        // the envelop profile for a segment [eta_L, eta_R] is
        // f(eta) = 0.5*(- Erf((eta_L - eta)/sigma)
        //               + Erf((eta_R - eta)/sigma))
        double eta_s_shift = 0.0;
        double tau_L = tau - dtau/2.;
        if (tau_L > tau_0 + delta_tau) {
            eta_s_shift = acosh((tau_L*tau_L + tau_0*tau_0
                                    - delta_tau*delta_tau)
                                   /(2.*tau_L*tau_0 + 1e-16));
        }
        double eta_s_L = std::min(it.lock()->eta_s_right,
                                  it.lock()->eta_s_0 - eta_s_shift);
        double eta_s_R = std::max(it.lock()->eta_s_left,
                                  it.lock()->eta_s_0 + eta_s_shift);

        double eta_s_next_shift = 0.0;
        double tau_next = tau + dtau/2.;
        if (tau_next > tau_0 + delta_tau) {
            eta_s_next_shift = acosh((tau_next*tau_next + tau_0*tau_0
                                      - delta_tau*delta_tau)
                                     /(2.*tau_next*tau_0 + 1e-16));
        }
        double eta_s_L_next = std::max(it.lock()->eta_s_left,
                                       it.lock()->eta_s_0 - eta_s_next_shift);
        double eta_s_R_next = std::min(it.lock()->eta_s_right,
                                       it.lock()->eta_s_0 + eta_s_next_shift);

        bool flag_left = true;  // the left string segment is valid
        if (eta_s_L_next > eta_s_L) flag_left = false;
        
        bool flag_right = true;  // the right string segment is valid
        if (eta_s_R_next < eta_s_R) flag_right = false;

        double exp_eta_s = 0.;
        if (flag_left) {
            if (   eta_s > eta_s_L_next - skip_dis_eta 
                && eta_s < eta_s_L + skip_dis_eta) {
                exp_eta_s += 0.5*(- erf((eta_s_L_next - eta_s)/sigma_eta)
                                  + erf((eta_s_L - eta_s)/sigma_eta));
            }
        }
        if (flag_right) {
            if (   eta_s > eta_s_R - skip_dis_eta 
                && eta_s < eta_s_R_next + skip_dis_eta) {
                exp_eta_s += 0.5*(- erf((eta_s_R - eta_s)/sigma_eta)
                                  + erf((eta_s_R_next - eta_s)/sigma_eta));
            }
        }

        double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                /(sigma_x*sigma_x));

        double e_frac = 1.0;
        if (eta_s < it.lock()->eta_s_left) {
            e_frac = it.lock()->frac_l;
        } else if (eta_s < it.lock()->eta_s_right) {
            e_frac = (it.lock()->frac_l
                      + (it.lock()->frac_r - it.lock()->frac_l)
                        /(it.lock()->eta_s_right - it.lock()->eta_s_left)
                        *(eta_s - it.lock()->eta_s_left));
        } else {
            e_frac = it.lock()->frac_r;
        }
        double e_local = e_frac*exp_tau*exp_xperp*exp_eta_s;
        e_local *= it.lock()->norm*sfactor;  // 1/fm^4
        double y_string = (
                it.lock()->y_l + (it.lock()->y_r - it.lock()->y_l)
                            /(it.lock()->eta_s_right - it.lock()->eta_s_left)
                            *(eta_s - it.lock()->eta_s_left));
        double y_dump = ((1. - string_quench_factor)*y_string
                         + string_quench_factor*y_long_flow);
        double y_dump_perp = string_quench_factor*y_perp_flow;
        double cosh_long = cosh(y_dump - eta_s);
        double sinh_long = sinh(y_dump - eta_s);
        double cosh_perp = 1.0;
        double sinh_perp = 0.0;
        if (std::abs(y_dump_perp) > 1e-6) {
            cosh_perp = cosh(y_dump_perp);
            sinh_perp = sinh(y_dump_perp);
        }
        j_mu[0] += e_local*cosh_long*cosh_perp;
        j_mu[1] += e_local*sinh_perp*cos_phi_flow;
        j_mu[2] += e_local*sinh_perp*sin_phi_flow;
        j_mu[3] += e_local*sinh_long*cosh_perp;
    }

    for (auto const&it: QCD_strings_remnant_list_current_tau) {
        // add remnant energy at the string ends
        bool flag_left = false;
        if (   it.lock()->tau_end_left >= tau - dtau/2.
            && it.lock()->tau_end_left <  tau + dtau/2.) {
            flag_left = true;
        }

        bool flag_right = false;
        if (   it.lock()->tau_end_right >= tau - dtau/2.
            && it.lock()->tau_end_right <  tau + dtau/2.) {
            flag_right = true;
        }
        
        double x_dis = x - it.lock()->x_perp;
        if (std::abs(x_dis) > skip_dis_x) continue;
        
        double y_dis = y - it.lock()->y_perp;
        if (std::abs(y_dis) > skip_dis_x) continue;

        double exp_eta_s_left = 0.0;
        if (flag_left) {
            double eta_dis_left = std::abs(eta_s - it.lock()->eta_s_left);
            if (eta_dis_left < skip_dis_eta) {
                exp_eta_s_left = (exp(-eta_dis_left*eta_dis_left
                                      /(sigma_eta*sigma_eta)));
            }
        }
        double exp_eta_s_right = 0.0;
        if (flag_right) {
            double eta_dis_right = std::abs(eta_s - it.lock()->eta_s_right);
            if (eta_dis_right < skip_dis_eta) {
                exp_eta_s_right = (exp(-eta_dis_right*eta_dis_right
                                       /(sigma_eta*sigma_eta)));
            }
        }
        double exp_factors = exp_tau*(
              exp_eta_s_left*(it.lock()->frac_l)*(it.lock()->E_baryon_norm_L)
            + exp_eta_s_right*(it.lock()->frac_r)*(it.lock()->E_baryon_norm_R)
        );
        double e_baryon_local = 0.0;
        if (exp_factors > 0) {
            double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                    /(sigma_x*sigma_x));
            e_baryon_local = exp_xperp*exp_factors;
        }
        j_mu[0] += e_baryon_local*sfactor;
    }
    const double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
    const double prefactor_tau  = 1./dtau;
    const double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
    const double prefactors = prefactor_tau*prefactor_prep*prefactor_etas;
    j_mu[0] *= prefactors;
    j_mu[1] *= prefactors;
    j_mu[2] *= prefactors;
    j_mu[3] *= prefactors;
}


double HydroSourceStrings::get_hydro_rhob_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu) const {
    double res = 0.;
    if (QCD_strings_baryon_list_current_tau.size() == 0) return(res);
    
    const double sigma_x   = get_sigma_x();
    const double sigma_eta = get_sigma_eta();

    // flow velocity
    const double gamma_perp_flow  = sqrt(1. + u_mu[1]*u_mu[1] + u_mu[2]*u_mu[2]);
    const double y_perp_flow      = acosh(gamma_perp_flow);
    const double y_long_flow      = asinh(u_mu[3]/gamma_perp_flow) + eta_s;
    const double sin_phi_flow     = u_mu[1]/gamma_perp_flow;
    const double cos_phi_flow     = u_mu[2]/gamma_perp_flow;
    const double dtau             = DATA.delta_tau;

    const double exp_tau        = 1.0/tau;
    const double n_sigma_skip   = 5.;
    const double skip_dis_x     = n_sigma_skip*sigma_x;
    const double skip_dis_eta   = n_sigma_skip*sigma_eta;
    for (auto &it: QCD_strings_baryon_list_current_tau) {
        // skip the evaluation if the strings is too far away in the
        // space-time grid
        // dumping energy into the medium from the active strings
        //double tau_dis_left = fabs(tau - it->tau_end_left);
        //double tau_dis_right = fabs(tau - it->tau_end_right);
        int flag_left = 0;
        if (   it.lock()->tau_baryon_left >= tau - dtau/2.
            && it.lock()->tau_baryon_left <  tau + dtau/2.
            && it.lock()->baryon_frac_l > 0.) {
            flag_left = 1;
        }

        int flag_right = 0;
        if (   it.lock()->tau_baryon_right >= tau - dtau/2.
            && it.lock()->tau_baryon_right <  tau + dtau/2.
            && it.lock()->baryon_frac_r > 0.) {
            flag_right = 1;
        }

        if (flag_left == 0 && flag_right == 0) continue;

        double x_dis = x - it.lock()->x_perp;
        if (std::abs(x_dis) > skip_dis_x) continue;
        
        double y_dis = y - it.lock()->y_perp;
        if (std::abs(y_dis) > skip_dis_x) continue;

        double exp_eta_s_left = 0.0;
        if (flag_left == 1) {
            double eta_dis_left = std::abs(eta_s
                                           - it.lock()->eta_s_baryon_left);
            if (eta_dis_left < skip_dis_eta) {
                exp_eta_s_left = (exp(-eta_dis_left*eta_dis_left
                                      /(sigma_eta*sigma_eta)));
            }
        }

        double exp_eta_s_right = 0.0;
        if (flag_right == 1) {
            double eta_dis_right = std::abs(eta_s
                                            - it.lock()->eta_s_baryon_right);
            if (eta_dis_right < skip_dis_eta) {
                exp_eta_s_right = (exp(-eta_dis_right*eta_dis_right
                                       /(sigma_eta*sigma_eta)));
            }
        }
        
        double exp_factors = exp_tau*(
                exp_eta_s_left*it.lock()->baryon_frac_l
                + exp_eta_s_right*it.lock()->baryon_frac_r);
        if (exp_factors > 0) {
            double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                    /(sigma_x*sigma_x));
            double fsmear = exp_xperp*exp_factors;
            double rapidity_local = (
                (  exp_eta_s_left*(it.lock()->baryon_frac_l)*(it.lock()->y_l_baryon)
                 + exp_eta_s_right*(it.lock()->baryon_frac_r)*(it.lock()->y_r_baryon))
                /(  exp_eta_s_left*(it.lock()->baryon_frac_l)
                  + exp_eta_s_right*(it.lock()->baryon_frac_r) + 1e-16));
            double y_dump = ((1. - parton_quench_factor)*rapidity_local
                             + parton_quench_factor*y_long_flow);
            double y_dump_perp = parton_quench_factor*y_perp_flow;
            double p_dot_u = 1.;
            if (parton_quench_factor < 1.) {
                p_dot_u = (  u_mu[0]*cosh(y_dump)*cosh(y_dump_perp)
                           - u_mu[1]*sinh(y_dump_perp)*cos_phi_flow
                           - u_mu[2]*sinh(y_dump_perp)*sin_phi_flow
                           - u_mu[3]*sinh(y_dump)*cosh(y_dump_perp));
            }
            res += p_dot_u*fsmear;
        }
    }
    const double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
    const double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
    const double prefactor_tau  = 1./dtau;
    res *= prefactor_tau*prefactor_prep*prefactor_etas;
    return(res);
}
