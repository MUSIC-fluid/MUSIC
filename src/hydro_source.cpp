// Copyright 2017 Chun Shen

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "./hydro_source.h"
#include "./util.h"

using namespace std;

hydro_source::hydro_source(const InitData &DATA_in) :
    DATA(DATA_in) {
    source_tau_max = 0.0;
    source_tau_min = 100.0;
    if (DATA.Initial_profile == 13) {  // MC-Glauber-LEXUS
        sigma_tau            = 0.1;
        sigma_x              = 0.5;
        sigma_eta            = 0.5;
        volume               = DATA.delta_x*DATA.delta_y*DATA.delta_eta;
        string_dump_mode     = DATA.string_dump_mode;
        string_quench_factor = DATA.string_quench_factor;
        parton_quench_factor = 1.0;
        read_in_QCD_strings_and_partons();
    }
    if (DATA.Initial_profile == 30) {  // AMPT
        sigma_tau = 0.1;
        sigma_x   = 0.5;
        sigma_eta = 0.2;
        volume = DATA.delta_x*DATA.delta_y*DATA.delta_eta;
        parton_quench_factor = 1.0;
        read_in_AMPT_partons();
    }
}

hydro_source::~hydro_source() {
    if (DATA.Initial_profile == 13) {
        QCD_strings_list.clear();
    }
    if (DATA.Initial_profile == 30) {
        parton_list.clear();
    }
}

//! This function reads in the spatal information of the strings and partons
//! which are produced from the MC-Glauber-LEXUS model
void hydro_source::read_in_QCD_strings_and_partons() {
    string QCD_strings_filename = DATA.initName;
    string partons_filename = DATA.initName_rhob;
    music_message << "read in QCD strings list from " << QCD_strings_filename
                  << " and partons list from " << partons_filename;
    music_message.flush("info");
    string text_string;

    ifstream QCD_strings_file(QCD_strings_filename.c_str());
    if (!QCD_strings_file) {
        music_message << "hydro_source::read_in_QCD_strings_and_partons: "
                      << "can not open QCD strings file: "
                      << QCD_strings_filename;
        music_message.flush("error");
        exit(1);
    }
    getline(QCD_strings_file, text_string);  // read the header
    // now we read in data
    getline(QCD_strings_file, text_string);
    while (!QCD_strings_file.eof()) {
        stringstream text_stream(text_string);
        std::shared_ptr<QCD_string> new_string(new QCD_string);
        text_stream >> new_string->norm >> new_string->delta_E
                    >> new_string->tau_form
                    >> new_string->tau_0 >> new_string->eta_s_0
                    >> new_string->x_perp >> new_string->y_perp
                    >> new_string->eta_s_left >> new_string->eta_s_right
                    >> new_string->y_l >> new_string->y_r
                    >> new_string->frac_l >> new_string->frac_r
                    >> new_string->y_l_i;
        if (!text_stream.eof()) {
            // read in the last element
            text_stream >> new_string->y_r_i;
        } else {
            // the string is too short
            music_message << "read_in_QCD_strings_and_partons: "
                          << "the format of file"
                          << QCD_strings_filename << "is wrong~";
            music_message.flush("error");
            exit(1);
        }
        if (!text_stream.eof()) {
            // the string is too long
            music_message << "read_in_QCD_strings_and_partons: "
                          << "the format of file"
                          << QCD_strings_filename << "is wrong~";
            music_message.flush("error");
            exit(1);
        }

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

        if (source_tau_max < source_tau) {
            source_tau_max = source_tau;
        }

        if (source_tau_min > (new_string->tau_0 + new_string->tau_form)) {
            source_tau_min = new_string->tau_0 + new_string->tau_form;
        }
        getline(QCD_strings_file, text_string);
    }
    QCD_strings_file.close();
    music_message << "hydro_source: tau_min = " << source_tau_min << " fm/c.";
    music_message.flush("info");
    music_message << "hydro_source: tau_max = " << source_tau_max << " fm/c.";
    music_message.flush("info");
    
    double total_baryon_number = 0;
    for (auto const& it : QCD_strings_list)
        total_baryon_number += it->frac_l + it->frac_r;
    music_message << "total baryon number = " << total_baryon_number;
    music_message.flush("info");
}


//! This function reads in the partons information from the AMPT model
void hydro_source::read_in_AMPT_partons() {
    parton_list.clear();
    string AMPT_filename = DATA.initName_AMPT;
    music_message << "hydro_source: "
                  << "read in AMPT parton list from " << AMPT_filename;
    music_message.flush("info");

    string text_string;
    ifstream AMPT_file(AMPT_filename.c_str());
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
    stringstream text_stream1(text_string);
    text_stream1 >> event_id >> dummy >> n_partons;

    // now we read in data
    for (int ipart = 0; ipart < n_partons; ipart++) {
        getline(AMPT_file, text_string);
        stringstream text_stream(text_string);
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
                                      /(t_local - z_local + 1e-15));
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
            cout << "pid = " << pid << endl;
        }
        parton_list.push_back(new_parton);
        if (source_tau_max < new_parton->tau) {
            source_tau_max = new_parton->tau;
        }
        if (source_tau_min > new_parton->tau) {
            source_tau_min = new_parton->tau;
        }
    }
    AMPT_file.close();
    music_message << "hydro_source:: read in " << parton_list.size() << "/"
                  << n_partons << " partons.";
    music_message.flush("info");
    music_message << "hydro_source:: tau_min = " << source_tau_min << " fm.";
    music_message.flush("info");
    music_message << "hydro_source:: tau_max = " << source_tau_max << " fm.";
    music_message.flush("info");
}


void hydro_source::prepare_list_for_current_tau_frame(double tau_local) {
    double dtau = DATA.delta_tau;
    QCD_strings_list_current_tau.clear();
    QCD_strings_baryon_list_current_tau.clear();
    parton_list_current_tau.clear();
    if (DATA.Initial_profile == 13) {
        for (auto &it: QCD_strings_list) {
            if ((   it->tau_end_left >= (tau_local - 1./2.*dtau)
                 && it->tau_end_left <  (tau_local + 3./2.*dtau))
                || (   it->tau_end_right >= (tau_local - 1./2.*dtau)
                    && it->tau_end_right <  (tau_local + 3./2.*dtau))) {
                QCD_strings_baryon_list_current_tau.push_back(it);
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
    } else if (DATA.Initial_profile == 30) {
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
}

void hydro_source::get_hydro_energy_source(
    double tau, double x, double y, double eta_s, 
    FlowVec &u_mu, EnergyFlowVec &j_mu) {
    j_mu = {0};
    // flow velocity
    const double gamma_perp_flow = sqrt(1. + u_mu[1]*u_mu[1] + u_mu[2]*u_mu[2]);
    const double y_perp_flow     = acosh(gamma_perp_flow);
    const double y_long_flow     = asinh(u_mu[3]/gamma_perp_flow) + eta_s;
    const double sin_phi_flow    = u_mu[1]/gamma_perp_flow;
    const double cos_phi_flow    = u_mu[2]/gamma_perp_flow;
    const double dtau            = DATA.delta_tau;

    const double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
    const double prefactor_tau = 1./dtau;
    const double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
    const double n_sigma_skip   = 5.;
    const double skip_dis_x     = n_sigma_skip*sigma_x;
    const double skip_dis_eta   = n_sigma_skip*sigma_eta;
    const double sfactor        = DATA.sFactor/hbarc;
    if (DATA.Initial_profile == 13) {
        // energy source from strings
        // double prefactor_tau = 1./(sqrt(M_PI)*sigma_tau);
        for (auto &it: QCD_strings_list_current_tau) {
            const double tau_0     = it->tau_0;
            const double delta_tau = it->tau_form;
            
            double x_dis = x - it->x_perp;
            if (std::abs(x_dis) > skip_dis_x) continue;

            double y_dis = y - it->y_perp;
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
                                       /(2.*tau_L*tau_0 + 1e-10));
            }
            double eta_s_L = std::min(it->eta_s_right, it->eta_s_0 - eta_s_shift);
            double eta_s_R = std::max(it->eta_s_left,  it->eta_s_0 + eta_s_shift);

            double eta_s_next_shift = 0.0;
            double tau_next = tau + dtau/2.;
            if (tau_next > tau_0 + delta_tau) {
                eta_s_next_shift = acosh((tau_next*tau_next + tau_0*tau_0
                                          - delta_tau*delta_tau)
                                         /(2.*tau_next*tau_0 + 1e-10));
            }
            double eta_s_L_next = std::max(it->eta_s_left,  it->eta_s_0 - eta_s_next_shift);
            double eta_s_R_next = std::min(it->eta_s_right, it->eta_s_0 + eta_s_next_shift);

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

            double exp_tau = 1./tau;
            double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                    /(sigma_x*sigma_x));

            double e_frac = 1.0;
            if (eta_s < it->eta_s_left) {
                e_frac = it->frac_l;
            } else if (eta_s < it->eta_s_right) {
                e_frac = (it->frac_l
                          + (it->frac_r - it->frac_l)
                            /(it->eta_s_right - it->eta_s_left)
                            *(eta_s - it->eta_s_left));
            } else {
                e_frac = it->frac_r;
            }
            double e_local = e_frac*exp_tau*exp_xperp*exp_eta_s;
            e_local *= sfactor;  // 1/fm^4
            double y_string = (
                    it->y_l + (it->y_r - it->y_l)
                                /(it->eta_s_right - it->eta_s_left)
                                *(eta_s - it->eta_s_left));
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
        double prefactors = prefactor_tau*prefactor_prep*prefactor_etas;
        j_mu[0] *= prefactors;
        j_mu[1] *= prefactors;
        j_mu[2] *= prefactors;
        j_mu[3] *= prefactors;
    } else if (DATA.Initial_profile == 30) {
        // AMPT parton sources
        double n_sigma_skip = 5.;
        double tau_dis_max = tau - source_tau_max;
        if (tau_dis_max < n_sigma_skip*sigma_tau) {
            for (auto &it: parton_list_current_tau) {
                double x_dis = x - it->x;
                if (fabs(x_dis) > skip_dis_x) {
                    continue;
                }
                double y_dis = y - it->y;
                if (fabs(y_dis) > skip_dis_x) {
                    continue;
                }
                double eta_s_dis = eta_s - it->eta_s;
                if (fabs(eta_s_dis) > skip_dis_eta) {
                    continue;
                }
                double exp_tau = 1./tau;
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
            double norm = DATA.sFactor/hbarc;     // 1/fm^4
            double prefactor = norm*prefactor_tau*prefactor_prep*prefactor_etas;
            j_mu[0] *= prefactor;
            j_mu[1] *= prefactor;
            j_mu[2] *= prefactor;
            j_mu[3] *= prefactor;
        }
    }
}

double hydro_source::get_hydro_rhob_source(double tau, double x, double y,
                                           double eta_s, FlowVec &u_mu) {
    double res = 0.;

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
    if (DATA.Initial_profile == 13) {
        for (auto &it: QCD_strings_baryon_list_current_tau) {
            // skip the evaluation if the strings is too far away in the
            // space-time grid
            // dumping energy into the medium from the active strings
            //double tau_dis_left = fabs(tau - it->tau_end_left);
            //double tau_dis_right = fabs(tau - it->tau_end_right);
            int flag_left = 0;
            if (   it->tau_end_left >= tau - dtau/2.
                && it->tau_end_left <  tau + dtau/2.) {
                flag_left = 1;
            }

            int flag_right = 0;
            if (   it->tau_end_right >= tau - dtau/2.
                && it->tau_end_right <  tau + dtau/2.) {
                flag_right = 1;
            }

            if (flag_left == 0 && flag_right == 0) continue;

            double x_dis = x - it->x_perp;
            if (std::abs(x_dis) > skip_dis_x) continue;
            
            double y_dis = y - it->y_perp;
            if (std::abs(y_dis) > skip_dis_x) continue;

            double exp_eta_s_left = 0.0;
            if (flag_left == 1) {
                double eta_dis_left = std::abs(eta_s - it->eta_s_left);
                if (eta_dis_left < skip_dis_eta) {
                    exp_eta_s_left = (exp(-eta_dis_left*eta_dis_left
                                          /(sigma_eta*sigma_eta)));
                }
            }

            double exp_eta_s_right = 0.0;
            if (flag_right == 1) {
                double eta_dis_right = std::abs(eta_s - it->eta_s_right);
                if (eta_dis_right < skip_dis_eta) {
                    exp_eta_s_right = (exp(-eta_dis_right*eta_dis_right
                                           /(sigma_eta*sigma_eta)));
                }
            }
            
            double exp_factors = exp_tau*(
                    exp_eta_s_left*it->frac_l + exp_eta_s_right*it->frac_r);
            if (exp_factors > 0) {
                double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(sigma_x*sigma_x));
                double fsmear = exp_xperp*exp_factors;
                double rapidity_local = (
                    (  exp_eta_s_left*it->frac_l*it->y_l
                     + exp_eta_s_right*it->frac_r*it->y_r)
                    /(  exp_eta_s_left*it->frac_l
                      + exp_eta_s_right*it->frac_r));
                double y_dump = ((1. - parton_quench_factor)*rapidity_local
                                 + parton_quench_factor*y_long_flow);
                double y_dump_perp = parton_quench_factor*y_perp_flow;
                double p_dot_u = (u_mu[0]
                    - tanh(y_dump_perp)*sinh_y_perp_flow/cosh(y_dump - eta_s)
                    - tanh(y_dump - eta_s)*u_mu[3]);
                res += p_dot_u*fsmear;
            }
        }
        res *= prefactor_tau*prefactor_prep*prefactor_etas;
    } else if (DATA.Initial_profile == 30) {
        double tau_dis_max = tau - source_tau_max;
        if (tau_dis_max < n_sigma_skip*sigma_tau) {
            for (auto &it: parton_list_current_tau) {
                // skip the evaluation if the strings is too far away in the
                // space-time grid
                double x_dis = x - it->x;
                if (fabs(x_dis) > skip_dis_x) {
                    continue;
                }
                double y_dis = y - it->y;
                if (fabs(y_dis) > skip_dis_x) {
                    continue;
                }
                double eta_s_dis = eta_s - it->eta_s;
                if (fabs(eta_s_dis) > skip_dis_eta) {
                    continue;
                }
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
    }
    return(res);
}

void hydro_source::get_hydro_energy_source_before_tau(
    double tau, double x, double y, double eta_s, double *j_mu) {
    FlowVec u                   = {0};
    u[0]                        = 1.0;
    EnergyFlowVec j_mu_one_step = {0};

    double tau0 = 0.0;
    double dtau = DATA.delta_tau;
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

double hydro_source::get_hydro_rhob_source_before_tau(
        double tau, double x, double y, double eta_s) {
    FlowVec u = {0};
    u[0] = 1.0;

    double res  = 0.;
    double tau0 = 0.0;
    double dtau = DATA.delta_tau;

    int n_tau_steps = static_cast<int>((tau - tau0)/dtau);
    for (int i = 0; i < n_tau_steps; i++) {
        const double tau_local = tau0 + (i + 0.5)*dtau;
        const double res_local = get_hydro_rhob_source(tau_local, x, y, eta_s, u);
        res += tau_local*res_local*dtau;
    }

    return res/tau;
}

