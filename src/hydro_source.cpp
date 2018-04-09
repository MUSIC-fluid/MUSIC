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
        sigma_eta = 0.5;
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
        QCD_string new_string;
        text_stream >> new_string.norm >> new_string.delta_E
                    >> new_string.tau_form
                    >> new_string.tau_0 >> new_string.eta_s_0
                    >> new_string.x_perp >> new_string.y_perp
                    >> new_string.eta_s_left >> new_string.eta_s_right
                    >> new_string.y_l >> new_string.y_r
                    >> new_string.frac_l >> new_string.frac_r
                    >> new_string.y_l_i;
        if (!text_stream.eof()) {
            // read in the last element
            text_stream >> new_string.y_r_i;
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

        new_string.status = 0;
        double temp_factor1 = (new_string.tau_0*new_string.tau_0
                               - new_string.tau_form*new_string.tau_form);
        double temp_factor2 = (new_string.tau_0
                        *cosh(new_string.eta_s_left - new_string.eta_s_0));
        double temp_factor3 = (new_string.tau_0
                    *cosh(new_string.eta_s_right - new_string.eta_s_0));
        double tau_end_left_local = (
            temp_factor2 + sqrt(temp_factor2*temp_factor2 - temp_factor1));
        double tau_end_right_local = (
            temp_factor3 + sqrt(temp_factor3*temp_factor3 - temp_factor1));
        new_string.tau_end_left = tau_end_left_local;
        new_string.tau_end_right = tau_end_right_local;

        // read in one string properly
        QCD_strings_list.push_back(new_string);

        // record the proper time of the first and last string sources
        double source_tau = new_string.tau_form;
        if (new_string.tau_end_left > new_string.tau_end_right) {
            source_tau = new_string.tau_end_left;
        } else {
            source_tau = new_string.tau_end_right;
        }

        if (source_tau_max < source_tau) {
            source_tau_max = source_tau;
        }

        if (source_tau_min > (new_string.tau_0 + new_string.tau_form)) {
            source_tau_min = new_string.tau_0 + new_string.tau_form;
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
        total_baryon_number += it.frac_l + it.frac_r;
    music_message << "total baryon number = " << total_baryon_number;
    music_message.flush("info");
}


//! This function reads in the partons information from the AMPT model
void hydro_source::read_in_AMPT_partons() {
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
    AMPT_file >> n_partons;

    // now we read in data
    getline(AMPT_file, text_string);
    while (!AMPT_file.eof()) {
        stringstream text_stream(text_string);
        parton new_parton;
        double t_local, z_local, pz_local;
        text_stream >> t_local >> new_parton.x >> new_parton.y >> z_local
                    >> new_parton.E >> new_parton.px >> new_parton.py
                    >> pz_local;
        if (t_local > z_local) {
            // the parton is inside the light cone
            double p_perp_sq = (new_parton.px*new_parton.px
                                + new_parton.py*new_parton.py);
            double mass_sq = (new_parton.E*new_parton.E
                              - p_perp_sq - pz_local*pz_local);
            if (mass_sq > 0.) {
                new_parton.mass = sqrt(mass_sq);
                new_parton.y_perp = asinh(sqrt(p_perp_sq)/new_parton.mass);
                new_parton.tau = sqrt(t_local*t_local - z_local*z_local);
                new_parton.eta_s = 0.5*log((t_local + z_local)
                                           /(t_local - z_local + 1e-15));
                new_parton.rapidity = 0.5*log((new_parton.E + pz_local)
                                              /(new_parton.E - pz_local));
                parton_list.push_back(new_parton);
                if (source_tau_max < new_parton.tau) {
                    source_tau_max = new_parton.tau;
                }
            }
        }
        getline(AMPT_file, text_string);
    }
    AMPT_file.close();
    music_message << "hydro_source:: read in " << parton_list.size() << "/"
                  << n_partons << " partons.";
    music_message.flush("info");
    music_message << "hydro_source:: tau_max = " << source_tau_max << " fm.";
    music_message.flush("info");
}

void hydro_source::get_hydro_energy_source(
    double tau, double x, double y, double eta_s, 
    FlowVec &u_mu, EnergyFlowVec &j_mu) {
    j_mu = {0};
    // flow velocity
    double gamma_perp_flow = sqrt(1. + u_mu[1]*u_mu[1] + u_mu[2]*u_mu[2]);
    double y_perp_flow     = acosh(gamma_perp_flow);
    double y_long_flow     = asinh(u_mu[3]/gamma_perp_flow) + eta_s;
    double sin_phi_flow    = u_mu[1]/gamma_perp_flow;
    double cos_phi_flow    = u_mu[2]/gamma_perp_flow;

    if (DATA.Initial_profile == 13) {
        // energy source from strings
        double n_sigma_skip = 5.;
        double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
        double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
        double dtau = DATA.delta_tau;
        // double prefactor_tau = 1./(sqrt(M_PI)*sigma_tau);
        for (auto &it: QCD_strings_list) {
            double tau_0 = it.tau_0;
            double delta_tau = it.tau_form;
            double tau_start = tau_0 + delta_tau;
            if (it.eta_s_left > it.eta_s_0) {
                tau_start = it.tau_end_left;
            } else if (it.eta_s_right < it.eta_s_0) {
                tau_start = it.tau_end_right;
            }
            if (tau > tau_start && it.status == 0) {
                // activiate the string when the constant tau hypersurface
                // starts to cross it
                it.status = 1;
            }
            if (it.status == 1) {
                // dumping energy into the medium from the active strings
                if (tau > it.tau_end_left && tau > it.tau_end_right) {
                    it.status = 2;
                    continue;
                }
                double x_dis = x - it.x_perp;
                if (fabs(x_dis) > n_sigma_skip*sigma_x) {
                    continue;
                }
                double y_dis = y - it.y_perp;
                if (fabs(y_dis) > n_sigma_skip*sigma_x) {
                    continue;
                }
                int flag_left = 1;
                int flag_right = 1;

                if (it.eta_s_0 < it.eta_s_left) {
                    flag_left = 0;
                }
                if (it.eta_s_0 > it.eta_s_right) {
                    flag_right = 0;
                }

                double eta_s_shift = acosh((tau*tau + tau_0*tau_0
                                            - delta_tau*delta_tau)
                                           /(2.*tau*tau_0));
                double eta_s_left = it.eta_s_0 - eta_s_shift;
                double eta_s_right = it.eta_s_0 + eta_s_shift;

                double eta_s_prev_shift = 0.0;
                double tau_prev = tau - dtau;
                if (tau_prev > tau_0 + delta_tau) {
                    eta_s_prev_shift = acosh((tau_prev*tau_prev + tau_0*tau_0
                                              - delta_tau*delta_tau)
                                             /(2.*tau_prev*tau_0));
                }
                double eta_s_left_prev = it.eta_s_0 - eta_s_prev_shift;
                double eta_s_right_prev = it.eta_s_0 + eta_s_prev_shift;

                if (eta_s_left_prev < it.eta_s_left) {
                    flag_left = 0;
                } else if (eta_s_left < it.eta_s_left) {
                    eta_s_left = it.eta_s_left;
                }
                
                if (eta_s_right_prev > it.eta_s_right) {
                    flag_right = 0;
                } else if (eta_s_right > it.eta_s_right) {
                    eta_s_right = it.eta_s_right;
                }

                double eta_s_left_dis = (
                        fabs(eta_s - (eta_s_left + eta_s_left_prev)/2.));
                double eta_s_right_dis = (
                        fabs(eta_s - (eta_s_right + eta_s_right_prev)/2.));

                if (eta_s_left_dis > n_sigma_skip*sigma_eta
                     && eta_s_right_dis > n_sigma_skip*sigma_eta) {
                    continue;
                }
                double exp_eta_s = 0.;
                if (flag_left == 1
                        && eta_s_left_dis < n_sigma_skip*sigma_eta) {
                    double deta_local = fabs(eta_s_left - eta_s_left_prev);
                    exp_eta_s += (exp(-(eta_s_left_dis*eta_s_left_dis)
                                      /(sigma_eta*sigma_eta))*deta_local);
                }
                if (flag_right == 1
                        && eta_s_right_dis < n_sigma_skip*sigma_eta) {
                    double deta_local = fabs(eta_s_right - eta_s_right_prev);
                    exp_eta_s += (exp(-(eta_s_right_dis*eta_s_right_dis)
                                      /(sigma_eta*sigma_eta))*deta_local);
                }
                double exp_tau = 1./tau;
                double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(sigma_x*sigma_x));

                double e_frac = 1.0;
                if (eta_s < it.eta_s_left) {
                    e_frac = it.frac_l;
                } else if (eta_s < it.eta_s_right) {
                    e_frac = (it.frac_l
                              + (it.frac_r - it.frac_l)
                                /(it.eta_s_right - it.eta_s_left)
                                *(eta_s - it.eta_s_left));
                } else {
                    e_frac = it.frac_r;
                }
                double e_local = e_frac*exp_tau*exp_xperp*exp_eta_s;
                e_local *= DATA.sFactor/hbarc;  // 1/fm^4
                double y_string = (
                        it.y_l + (it.y_r - it.y_l)
                                    /(it.eta_s_right - it.eta_s_left)
                                    *(eta_s - it.eta_s_left));
                double y_dump = ((1. - string_quench_factor)*y_string
                                 + string_quench_factor*y_long_flow);
                double y_dump_perp = string_quench_factor*y_perp_flow;
                double cosh_long = cosh(y_dump - eta_s);
                double sinh_long = sinh(y_dump - eta_s);
                double cosh_perp = cosh(y_dump_perp);
                double sinh_perp = sinh(y_dump_perp);
                j_mu[0] += e_local*cosh_long*cosh_perp;
                j_mu[1] += e_local*sinh_perp*cos_phi_flow;
                j_mu[2] += e_local*sinh_perp*sin_phi_flow;
                j_mu[3] += e_local*sinh_long*cosh_perp;
            }
        }
        double prefactors = prefactor_prep*prefactor_etas;
        j_mu[0] *= prefactors;
        j_mu[1] *= prefactors;
        j_mu[2] *= prefactors;
        j_mu[3] *= prefactors;
    } else if (DATA.Initial_profile == 30) {
        // AMPT parton sources
        double n_sigma_skip = 5.;
        double tau_dis_max = tau - source_tau_max;
        if (tau_dis_max < n_sigma_skip*sigma_tau) {
            double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
            double prefactor_tau = 1./(sqrt(M_PI)*sigma_tau);
            double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
            for (vector<parton>::iterator it = parton_list.begin();
                 it != parton_list.end(); it++) {
                double tau_dis = tau - (*it).tau;
                if (fabs(tau_dis) > n_sigma_skip*sigma_tau) {
                    continue;
                }
                double x_dis = x - (*it).x;
                if (fabs(x_dis) > n_sigma_skip*sigma_x) {
                    continue;
                }
                double y_dis = y - (*it).y;
                if (fabs(y_dis) > n_sigma_skip*sigma_x) {
                    continue;
                }
                double eta_s_dis = eta_s - (*it).eta_s;
                if (fabs(eta_s_dis) > n_sigma_skip*sigma_eta) {
                    continue;
                }
                double exp_tau = (
                    1./((*it).tau)
                    *exp(-tau_dis*tau_dis/(sigma_tau*sigma_tau)));
                double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(sigma_x*sigma_x));
                double exp_eta_s = (
                        exp(-eta_s_dis*eta_s_dis/(sigma_eta*sigma_eta)));

                double f_smear = exp_tau*exp_xperp*exp_eta_s;
                double p_perp_sq = (*it).px*(*it).px + (*it).py*(*it).py;
                double m_perp = sqrt((*it).mass*(*it).mass + p_perp_sq);
                j_mu[0] += m_perp*cosh((*it).rapidity - eta_s)*f_smear;
                j_mu[1] += (*it).px*f_smear;
                j_mu[2] += (*it).py*f_smear;
                j_mu[3] += m_perp*sinh((*it).rapidity - eta_s)*f_smear;
            }
            double norm = DATA.sFactor/hbarc;     // 1/fm^4
            j_mu[0] *= norm*prefactor_tau*prefactor_prep*prefactor_etas;
            j_mu[1] *= norm*prefactor_tau*prefactor_prep*prefactor_etas;
            j_mu[2] *= norm*prefactor_tau*prefactor_prep*prefactor_etas;
            j_mu[3] *= norm*prefactor_tau*prefactor_prep*prefactor_etas;
        }
    }
}

double hydro_source::get_hydro_rhob_source(double tau, double x, double y,
                                           double eta_s, FlowVec &u_mu) {
    double res = 0.;

    // flow velocity
    double gamma_perp_flow  = sqrt(1. + u_mu[1]*u_mu[1] + u_mu[2]*u_mu[2]);
    double y_perp_flow      = acosh(gamma_perp_flow);
    double y_long_flow      = asinh(u_mu[3]/gamma_perp_flow) + eta_s;
    double sinh_y_perp_flow = sinh(y_perp_flow);
    double dtau             = DATA.delta_tau;

    if (DATA.Initial_profile == 13) {
        double n_sigma_skip = 5.;
        double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
        double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
        double prefactor_tau = 1./dtau;
        for (auto &it: QCD_strings_list) {
            // skip the evaluation if the strings is too far away in the
            // space-time grid
            // dumping energy into the medium from the active strings
            //double tau_dis_left = fabs(tau - it.tau_end_left);
            //double tau_dis_right = fabs(tau - it.tau_end_right);
            //if (tau_dis_left > n_sigma_skip*sigma_tau
            //        && tau_dis_right > n_sigma_skip*sigma_tau) {
            //    continue;
            //}
            int flag_left = 0;
            int flag_right = 0;
            if (tau > it.tau_end_left
                    && tau < it.tau_end_left + dtau) {
                flag_left = 1;
            }
            if (tau > it.tau_end_right
                    && tau < it.tau_end_right + dtau) {
                flag_right = 1;
            }

            if (flag_left == 0 && flag_right == 0) {
                continue;
            }

            double x_dis = x - it.x_perp;
            if (fabs(x_dis) > n_sigma_skip*sigma_x) {
                continue;
            }
            double y_dis = y - it.y_perp;
            if (fabs(y_dis) > n_sigma_skip*sigma_x) {
                continue;
            }

            double exp_tau_left = 1.0/tau;
            double exp_eta_s_left = 0.0;
            if (flag_left == 1) {
                double eta_dis_left = fabs(eta_s - it.eta_s_left);
                if (eta_dis_left < n_sigma_skip*sigma_eta) {
                    exp_eta_s_left = (exp(-eta_dis_left*eta_dis_left
                                          /(sigma_eta*sigma_eta)));
                }
            }

            double exp_tau_right = 1.0/tau;
            double exp_eta_s_right = 0.0;
            if (flag_right == 1) {
                double eta_dis_right = fabs(eta_s - it.eta_s_right);
                if (eta_dis_right < n_sigma_skip*sigma_eta) {
                    exp_eta_s_right = (exp(-eta_dis_right*eta_dis_right
                                           /(sigma_eta*sigma_eta)));
                }
            }
            
            double exp_factors = (exp_tau_left*exp_eta_s_left*it.frac_l
                                + exp_tau_right*exp_eta_s_right*it.frac_r);
            if (exp_factors > 0) {
                double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(sigma_x*sigma_x));
                double fsmear = exp_xperp*exp_factors;
                double rapidity_local = (
                    (exp_tau_left*exp_eta_s_left*it.frac_l*it.y_l
                     + exp_tau_right*exp_eta_s_right*it.frac_r*it.y_r)
                    /(exp_tau_left*exp_eta_s_left*it.frac_l
                      + exp_tau_right*exp_eta_s_right*it.frac_r));
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
        double n_sigma_skip = 5.;
        double tau_dis_max = tau - source_tau_max;
        if (tau_dis_max < n_sigma_skip*sigma_tau) {
            double prefactor_prep = 1./(M_PI*sigma_x*sigma_x);
            double prefactor_etas = 1./(sqrt(M_PI)*sigma_eta);
            double prefactor_tau = 1./(sqrt(M_PI)*sigma_tau);
            for (vector<parton>::iterator it = parton_list.begin();
                 it != parton_list.end(); it++) {
                // skip the evaluation if the strings is too far away in the
                // space-time grid
                double tau_dis = tau - (*it).tau;
                if (fabs(tau_dis) > n_sigma_skip*sigma_tau) {
                    continue;
                }
                double x_dis = x - (*it).x;
                if (fabs(x_dis) > n_sigma_skip*sigma_x) {
                    continue;
                }
                double y_dis = y - (*it).y;
                if (fabs(y_dis) > n_sigma_skip*sigma_x) {
                    continue;
                }
                double eta_s_dis = eta_s - (*it).eta_s;
                if (fabs(eta_s_dis) > n_sigma_skip*sigma_eta) {
                    continue;
                }
                double exp_tau = (
                    1./((*it).tau)
                    *exp(-tau_dis*tau_dis/(sigma_tau*sigma_tau)));
                double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(sigma_x*sigma_x));
                double exp_eta_s = (
                        exp(-eta_s_dis*eta_s_dis/(sigma_eta*sigma_eta)));
                double f_smear = exp_tau*exp_xperp*exp_eta_s;
                double y_dump = ((1. - parton_quench_factor)*(*it).rapidity
                                 + parton_quench_factor*y_long_flow);
                double y_dump_perp = ((1. - parton_quench_factor)*(*it).y_perp
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
