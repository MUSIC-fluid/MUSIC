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

HydroSourceStrings::HydroSourceStrings(InitData &DATA_in) :
    DATA(DATA_in) {
    set_source_tau_min(100.0);
    set_source_tau_max(0.0);
    set_source_tauStart_max(0.0);
    set_sigma_tau(0.1);
    set_sigma_x  (DATA.stringSourceSigmaX);
    set_sigma_eta(DATA.stringSourceSigmaEta);
    string_dump_mode     = DATA.string_dump_mode;
    string_quench_factor = DATA.string_quench_factor;
    parton_quench_factor = 1.0;    // no diffusion current from the source
    stringTransverseShiftFrac_ = DATA.stringTransverseShiftFrac;
    preEqFlowFactor_ = DATA.preEqFlowFactor;
    read_in_QCD_strings_and_partons();
}


HydroSourceStrings::HydroSourceStrings(
        InitData &DATA_in, std::vector< std::vector<double> > QCDStringList) :
    DATA(DATA_in) {
    set_source_tau_min(100.0);
    set_source_tau_max(0.0);
    set_source_tauStart_max(0.0);
    set_sigma_tau(0.1);
    set_sigma_x  (DATA.stringSourceSigmaX);
    set_sigma_eta(DATA.stringSourceSigmaEta);
    string_dump_mode     = DATA.string_dump_mode;
    string_quench_factor = DATA.string_quench_factor;
    parton_quench_factor = 1.0;    // no diffusion current from the source
    stringTransverseShiftFrac_ = DATA.stringTransverseShiftFrac;
    preEqFlowFactor_ = DATA.preEqFlowFactor;
    read_in_QCD_strings_and_partons(QCDStringList);
}


HydroSourceStrings::~HydroSourceStrings() {
    QCD_strings_list.clear();
}


double HydroSourceStrings::getStringEndTau(
        const double tau_0, const double tau_form,
        const double eta_s_0, const double eta_s) const {
    double temp_factor1 = tau_0*tau_0 - tau_form*tau_form;
    double temp_factor2 = tau_0*cosh(eta_s - eta_s_0);
    double tau_end_local = (temp_factor2
                            + sqrt(temp_factor2*temp_factor2 - temp_factor1));
    return(tau_end_local);
}


//! This function reads in the spatal information of the strings and partons
//! which are produced from the MC-Glauber-LEXUS model
void HydroSourceStrings::read_in_QCD_strings_and_partons(
        std::vector< std::vector<double> > QCDStringList_in) {
    music_message.info("read in QCD strings list from vector");

    for (auto string_i: QCDStringList_in) {
        std::shared_ptr<QCD_string> new_string(new QCD_string);
        new_string->mass = string_i[0];
        new_string->m_over_sigma = string_i[1];
        new_string->tau_form = string_i[2];
        new_string->tau_0 = string_i[3];
        new_string->eta_s_0 = string_i[4];
        new_string->x_perp = string_i[5];
        new_string->y_perp = string_i[6];
        new_string->x_pl = string_i[7];
        new_string->y_pl = string_i[8];
        new_string->x_pr = string_i[9];
        new_string->y_pr = string_i[10];
        new_string->eta_s_left = string_i[11];
        new_string->eta_s_right = string_i[12];
        new_string->y_l = string_i[13];
        new_string->y_r = string_i[14];
        new_string->remnant_l = string_i[15];
        new_string->remnant_r = string_i[16];
        new_string->y_l_i = string_i[17];
        new_string->y_r_i = string_i[18];
        new_string->eta_s_baryon_left = string_i[19];
        new_string->eta_s_baryon_right = string_i[20];
        new_string->y_l_baryon = string_i[21];
        new_string->y_r_baryon = string_i[22];
        new_string->baryon_frac_l = string_i[23];
        new_string->baryon_frac_r = string_i[24];
        new_string->px_i = string_i[25];
        new_string->py_i = string_i[26];

        if (DATA.Initial_profile == 131) {
            new_string->tau_0 = 0.;
            new_string->eta_s_0 = 0.;
            //new_string->tau_form = 0.5;
        }
        new_string->sigma_x = get_sigma_x();
        new_string->sigma_eta = get_sigma_eta();

        // compute the string end tau
        new_string->tau_end_left = getStringEndTau(
                new_string->tau_0, new_string->tau_form,
                new_string->eta_s_0, new_string->eta_s_left);
        new_string->tau_end_right = getStringEndTau(
                new_string->tau_0, new_string->tau_form,
                new_string->eta_s_0, new_string->eta_s_right);

        // compute the baryon number tau
        new_string->tau_baryon_left = getStringEndTau(
                new_string->tau_0, new_string->tau_form,
                new_string->eta_s_0, new_string->eta_s_baryon_left);
        new_string->tau_baryon_right = getStringEndTau(
                new_string->tau_0, new_string->tau_form,
                new_string->eta_s_0, new_string->eta_s_baryon_right);

        // determine the tau_start and eta_s_start of the string
        if (new_string->eta_s_left > new_string->eta_s_0) {
            new_string->tau_start = new_string->tau_end_left;
            new_string->eta_s_start = new_string->eta_s_left;
        } else if (new_string->eta_s_right < new_string->eta_s_0) {
            new_string->tau_start = new_string->tau_end_right;
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

        if (new_string->eta_s_left < 1 && new_string->eta_s_right > -1) {
            // the string in the mid-rapidity
            if (new_string->tau_start > get_source_tauStart_max()) {
                set_source_tauStart_max(new_string->tau_start);
            }
        }

        // set the maximum tau = 10 fm/c for source terms
        source_tau = std::min(10., source_tau);

        if (get_source_tau_max() < source_tau) {
            set_source_tau_max(source_tau);
        }

        if (get_source_tau_min() > new_string->tau_start) {
            set_source_tau_min(new_string->tau_start);
        }
    }
    music_message << "HydroSourceStrings: tau_min = " << get_source_tau_min()
                  << " fm/c.";
    music_message.flush("info");
    music_message << "HydroSourceStrings: tau_max = " << get_source_tau_max()
                  << " fm/c.";
    music_message.flush("info");

    double total_baryon_number = 0;
    double total_energy = 0.;
    for (auto const& it : QCD_strings_list) {
        total_baryon_number += it->baryon_frac_l + it->baryon_frac_r;
        total_energy += it->mass*(cosh(it->y_l_i) + cosh(it->y_r_i));
    }
    music_message << "total baryon number = " << total_baryon_number;
    music_message.flush("info");
    music_message << "total energy = " << total_energy << " GeV";
    music_message.flush("info");
    compute_norm_for_strings();

    double xMax = 0.;
    double yMax = 0.;
    for (auto const& it : QCD_strings_list) {
        xMax = std::max(xMax, std::abs(it->x_perp));
        xMax = std::max(xMax, std::abs(it->x_pl));
        xMax = std::max(xMax, std::abs(it->x_pr));
        yMax = std::max(yMax, std::abs(it->y_perp));
        yMax = std::max(yMax, std::abs(it->y_pl));
        yMax = std::max(yMax, std::abs(it->y_pr));
    }

    // adjust transverse grid size
    double energyAddRadius = sqrt(total_energy
                                  /std::max(1., total_baryon_number)/2760.);
    double npartAddRadius = total_baryon_number/500.;
    double reRunAddRadius = DATA.reRunCount*2.;
    double gridOffset = 3 + energyAddRadius + npartAddRadius + reRunAddRadius;
    gridOffset = std::max(gridOffset, 5.*DATA.stringSourceSigmaX);
    DATA.x_size = 2.*(xMax + gridOffset);
    DATA.y_size = 2.*(yMax + gridOffset);
    DATA.delta_x = DATA.x_size/(DATA.nx - 1);
    DATA.delta_y = DATA.y_size/(DATA.ny - 1);
    if (DATA.resetDtau) {
        // make sure delta_tau is not too large for delta_x and delta_y
        DATA.delta_tau = std::min(DATA.delta_tau,
                                  std::min(DATA.delta_x*DATA.dtaudxRatio,
                                           DATA.delta_y*DATA.dtaudxRatio));
        if (DATA.delta_tau > 0.001)
            DATA.delta_tau = (static_cast<int>(DATA.delta_tau*1000))/1000.;
        DATA.nt = static_cast<int>(DATA.tau_size/DATA.delta_tau + 0.5);
    }
    music_message << "[HydroSource] Grid info: x_size = "
                  << DATA.x_size << ", y_size = " << DATA.y_size
                  << ", dx = " << DATA.delta_x << " fm, dy = "
                  << DATA.delta_y << " fm, dtau = " << DATA.delta_tau << " fm";
    music_message.flush("info");
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

    // read in collision information
    getline(QCD_strings_file, text_string);
    std::stringstream coll_info(text_string);
    std::string dummy;
    double total_energy = 0;
    for (int ii = 0; ii < 16; ii++) {
        coll_info >> dummy;
    }
    coll_info >> total_energy;

    // read the header
    getline(QCD_strings_file, text_string);

    // now we read in data
    getline(QCD_strings_file, text_string);
    while (!QCD_strings_file.eof()) {
        std::stringstream text_stream(text_string);
        std::shared_ptr<QCD_string> new_string(new QCD_string);
        text_stream >> new_string->mass >> new_string->m_over_sigma
                    >> new_string->tau_form
                    >> new_string->tau_0 >> new_string->eta_s_0
                    >> new_string->x_perp >> new_string->y_perp
                    >> new_string->x_pl >> new_string->y_pl
                    >> new_string->x_pr >> new_string->y_pr
                    >> new_string->eta_s_left >> new_string->eta_s_right
                    >> new_string->y_l >> new_string->y_r
                    >> new_string->remnant_l >> new_string->remnant_r
                    >> new_string->y_l_i >> new_string->y_r_i
                    >> new_string->eta_s_baryon_left
                    >> new_string->eta_s_baryon_right
                    >> new_string->y_l_baryon
                    >> new_string->y_r_baryon
                    >> new_string->baryon_frac_l;
        if (!text_stream.eof()) {
            // read in the last three elements
            text_stream >> new_string->baryon_frac_r
                        >> new_string->px_i
                        >> new_string->py_i;
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
            //new_string->tau_form = 0.5;
        }
        new_string->sigma_x = get_sigma_x();
        new_string->sigma_eta = get_sigma_eta();

        // compute the string end tau
        new_string->tau_end_left = getStringEndTau(
                new_string->tau_0, new_string->tau_form,
                new_string->eta_s_0, new_string->eta_s_left);
        new_string->tau_end_right = getStringEndTau(
                new_string->tau_0, new_string->tau_form,
                new_string->eta_s_0, new_string->eta_s_right);

        // compute the baryon number tau
        new_string->tau_baryon_left = getStringEndTau(
                new_string->tau_0, new_string->tau_form,
                new_string->eta_s_0, new_string->eta_s_baryon_left);
        new_string->tau_baryon_right = getStringEndTau(
                new_string->tau_0, new_string->tau_form,
                new_string->eta_s_0, new_string->eta_s_baryon_right);

        // determine the tau_start and eta_s_start of the string
        if (new_string->eta_s_left > new_string->eta_s_0) {
            new_string->tau_start = new_string->tau_end_left;
            new_string->eta_s_start = new_string->eta_s_left;
        } else if (new_string->eta_s_right < new_string->eta_s_0) {
            new_string->tau_start = new_string->tau_end_right;
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

        if (new_string->eta_s_left < 1 && new_string->eta_s_right > -1) {
            // the string in the mid-rapidity
            if (new_string->tau_start > get_source_tauStart_max()) {
                set_source_tauStart_max(new_string->tau_start);
            }
        }

        // set the maximum tau = 10 fm/c for source terms
        source_tau = std::min(10., source_tau);

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
    music_message << "total energy = " << total_energy << " GeV";
    music_message.flush("info");
    compute_norm_for_strings();

    double xMax = 0.;
    double yMax = 0.;
    for (auto const& it : QCD_strings_list) {
        xMax = std::max(xMax, std::abs(it->x_perp));
        xMax = std::max(xMax, std::abs(it->x_pl));
        xMax = std::max(xMax, std::abs(it->x_pr));
        yMax = std::max(yMax, std::abs(it->y_perp));
        yMax = std::max(yMax, std::abs(it->y_pl));
        yMax = std::max(yMax, std::abs(it->y_pr));
    }

    // adjust transverse grid size
    double energyAddRadius = sqrt(total_energy
                                  /std::max(1., total_baryon_number)/2760.);
    double npartAddRadius = total_baryon_number/500.;
    double reRunAddRadius = DATA.reRunCount*2.;
    double gridOffset = 3 + energyAddRadius + npartAddRadius + reRunAddRadius;
    gridOffset = std::max(gridOffset, 5.*DATA.stringSourceSigmaX);
    DATA.x_size = 2.*(xMax + gridOffset);
    DATA.y_size = 2.*(yMax + gridOffset);
    DATA.delta_x = DATA.x_size/(DATA.nx - 1);
    DATA.delta_y = DATA.y_size/(DATA.ny - 1);
    if (DATA.resetDtau) {
        // make sure delta_tau is not too large for delta_x and delta_y
        DATA.delta_tau = std::min(DATA.delta_tau,
                                  std::min(DATA.delta_x*DATA.dtaudxRatio,
                                           DATA.delta_y*DATA.dtaudxRatio));
        if (DATA.delta_tau > 0.001)
            DATA.delta_tau = (static_cast<int>(DATA.delta_tau*1000))/1000.;
        DATA.nt = static_cast<int>(DATA.tau_size/DATA.delta_tau + 0.5);
    }
    music_message << "[HydroSource] Grid info: x_size = "
                  << DATA.x_size << ", y_size = " << DATA.y_size
                  << ", dx = " << DATA.delta_x << " fm, dy = "
                  << DATA.delta_y << " fm, dtau = " << DATA.delta_tau << " fm";
    music_message.flush("info");
}


void HydroSourceStrings::compute_norm_for_strings() {
    const int neta = 500;
    const double eta_range = 12.;
    const double deta = 2.*eta_range/(neta - 1);

    double E_string_total  = 0.0;
    double E_remnant_total = 0.0;
    double Pz_string_total  = 0.0;
    double Pz_remnant_total = 0.0;
    for (auto &it: QCD_strings_list) {
        const double sigma_eta = it->sigma_eta;
        const double prefactor_etas = 1./(sqrt(2.*M_PI)*sigma_eta);

        double E_string_norm    = 0.;
        double E_remnant_L_norm = 0.;
        double E_remnant_R_norm = 0.;
        for (int ieta = 0; ieta < neta; ieta++) {
            double eta_local = - eta_range + ieta*deta;
            double Delta_eta = it->eta_s_right - it->eta_s_left;
            double denorm_safe = std::copysign(
                std::max(Util::small_eps, std::abs(Delta_eta)), Delta_eta);
            double y_eta = (it->y_l + (it->y_r - it->y_l)/denorm_safe
                                      *(eta_local - it->eta_s_left));
            double expon_left  = (
                        (it->eta_s_left  - eta_local)/(sqrt(2.)*sigma_eta));
            double expon_right = (
                        (it->eta_s_right - eta_local)/(sqrt(2.)*sigma_eta));
            double e_eta = 0.5*(- erf(expon_left) + erf(expon_right));
            E_string_norm += e_eta*cosh(y_eta);

            double e_remnant_L = exp(-expon_left*expon_left);
            double e_remnant_R = exp(-expon_right*expon_right);
            E_remnant_L_norm += e_remnant_L*cosh(it->y_l);
            E_remnant_R_norm += e_remnant_R*cosh(it->y_r);
        }
        E_string_norm *= prefactor_etas*deta;
        double E_string = (it->mass*cosh(it->y_l_i) + it->mass*cosh(it->y_r_i)
                           - it->mass*cosh(it->y_l) - it->mass*cosh(it->y_r));
        it->norm = E_string/std::max(E_string_norm, Util::small_eps);
        E_string_total += E_string;
        double Pz_string = (  it->mass*sinh(it->y_l_i) - it->mass*sinh(it->y_l)
                            + it->mass*sinh(it->y_r_i) - it->mass*sinh(it->y_r));
        Pz_string_total += Pz_string;

        // here the E_norm is for the energy of remnants at the string ends
        E_remnant_L_norm *= prefactor_etas*deta;
        E_remnant_R_norm *= prefactor_etas*deta;
        double E_remnant_L   = it->remnant_l*it->mass*cosh(it->y_l);
        double E_remnant_R   = it->remnant_r*it->mass*cosh(it->y_r);
        it->E_remnant_norm_L = (E_remnant_L
                                /std::max(E_remnant_L_norm, Util::small_eps));
        it->E_remnant_norm_R = (E_remnant_R
                                /std::max(E_remnant_R_norm, Util::small_eps));
        E_remnant_total += E_remnant_L + E_remnant_R;
        double Pz_remnant_L = it->remnant_l*it->mass*sinh(it->y_l);
        double Pz_remnant_R = it->remnant_r*it->mass*sinh(it->y_r);
        Pz_remnant_total += Pz_remnant_L + Pz_remnant_R;
    }
    music_message << "E_total = "
                  << E_string_total + E_remnant_total << " GeV. "
                  << "E_string_total = " << E_string_total << " GeV, "
                  << "E_remnant_total = " << E_remnant_total << " GeV.";
    music_message.flush("info");
    music_message << "Pz_total = "
                  << Pz_string_total + Pz_remnant_total << " GeV. "
                  << "Pz_string_total = " << Pz_string_total << " GeV, "
                  << "Pz_remnant_total = " << Pz_remnant_total << " GeV.";
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

    const double dtau = DATA.delta_tau;
    const double n_sigma_skip = 8.;
    const double exp_tau = 1./tau;
    for (auto const&it: QCD_strings_list_current_tau) {
        const double sigma_x = it->sigma_x;
        const double sigma_x_sq = sigma_x*sigma_x;
        const double sigma_eta = it->sigma_eta;
        const double alpsig = preEqFlowFactor_*sigma_x;
        const double prefactor_prep = (
            1./(2.*M_PI*(sigma_x_sq
                         + exp(alpsig*alpsig/2.)*sqrt(M_PI/2)
                           *alpsig*sigma_x_sq*erf(alpsig/sqrt(2.)))
               )
        );
        const double prefactor_etas = 1./(sqrt(2.*M_PI)*sigma_eta);
        const double skip_dis_x = n_sigma_skip*sigma_x;
        const double skip_dis_eta = n_sigma_skip*sigma_eta;
        // energy source from strings
        const double tau_0     = it->tau_0;
        const double delta_tau = it->tau_form;

        if (eta_s < it->eta_s_left - skip_dis_eta
                || eta_s > it->eta_s_right + skip_dis_eta) continue;
        double eta_frac = ((eta_s - it->eta_s_left)
                           /std::max(Util::small_eps,
                                     it->eta_s_right - it->eta_s_left));
        eta_frac = std::max(0., std::min(1., eta_frac));

        const double x_perp = getStringTransverseCoord(it->x_pl, it->x_pr,
                                                       eta_frac);
        double x_dis = x - x_perp;
        if (std::abs(x_dis) > skip_dis_x) continue;

        const double y_perp = getStringTransverseCoord(it->y_pl, it->y_pr,
                                                       eta_frac);
        double y_dis = y - y_perp;
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
                                   /std::max(Util::small_eps, 2.*tau_L*tau_0));
        }
        double eta_s_L = std::min(it->eta_s_right,
                                  it->eta_s_0 - eta_s_shift);
        double eta_s_R = std::max(it->eta_s_left,
                                  it->eta_s_0 + eta_s_shift);

        double eta_s_next_shift = 0.0;
        double tau_next = tau + dtau/2.;
        if (tau_next > tau_0 + delta_tau) {
            eta_s_next_shift = acosh(
                (tau_next*tau_next + tau_0*tau_0 - delta_tau*delta_tau)
                /std::max(Util::small_eps, 2.*tau_next*tau_0));
        }
        double eta_s_L_next = std::max(it->eta_s_left,
                                       it->eta_s_0 - eta_s_next_shift);
        double eta_s_R_next = std::min(it->eta_s_right,
                                       it->eta_s_0 + eta_s_next_shift);

        bool flag_left = true;  // the left string segment is valid
        if (eta_s_L_next > eta_s_L) flag_left = false;

        bool flag_right = true;  // the right string segment is valid
        if (eta_s_R_next < eta_s_R) flag_right = false;

        double exp_eta_s = 0.;
        if (flag_left) {
            if (   eta_s > eta_s_L_next - skip_dis_eta
                && eta_s < eta_s_L + skip_dis_eta) {
                exp_eta_s += 0.5*(
                        - erf((eta_s_L_next - eta_s)/(sqrt(2.)*sigma_eta))
                        + erf((eta_s_L - eta_s)/(sqrt(2.)*sigma_eta)));
            }
        }
        if (flag_right) {
            if (   eta_s > eta_s_R - skip_dis_eta
                && eta_s < eta_s_R_next + skip_dis_eta) {
                exp_eta_s += 0.5*(
                        - erf((eta_s_R - eta_s)/(sqrt(2.)*sigma_eta))
                        + erf((eta_s_R_next - eta_s)/(sqrt(2.)*sigma_eta)));
            }
        }

        double exp_xperp = exp(-(x_dis*x_dis + y_dis*y_dis)
                                /(2.*sigma_x*sigma_x));
        double cosh_perp = (
            cosh(preEqFlowFactor_*sqrt(x_dis*x_dis + y_dis*y_dis)));
        double sinh_perp = (
            sinh(preEqFlowFactor_*sqrt(x_dis*x_dis + y_dis*y_dis)));
        double phi_perp = atan2(y_dis, x_dis);

        double e_local = exp_tau*exp_xperp*exp_eta_s*it->norm;
        double Delta_eta = it->eta_s_right - it->eta_s_left;
        double denorm_safe = std::copysign(
                std::max(Util::small_eps, std::abs(Delta_eta)), Delta_eta);
        double y_string = (it->y_l + (it->y_r - it->y_l)/denorm_safe
                                     *(eta_s - it->eta_s_left));
        double cosh_long = cosh(y_string - eta_s);
        double sinh_long = sinh(y_string - eta_s);
        //double cosh_perp = 1.0;
        double local_eperp = prefactor_etas*prefactor_prep*e_local;
        j_mu[0] += local_eperp*cosh_long*cosh_perp;
        if (std::isnan(j_mu[0])) {
            std::cout << local_eperp << "  " << cosh_long << std::endl;
            std::cout << prefactor_etas << "  " << prefactor_prep << "  "
                      << e_local << "  " << cosh_perp
                      << std::endl;
            std::cout << exp_tau << "  " << exp_xperp << "  "
                      << exp_eta_s << "  " << it->norm << std::endl;
            std::cout << x_dis << "  " << y_dis << "  " << sigma_x << std::endl;
            std::cout << it->x_pl << "  " << it->x_pr << "  " << eta_frac
                      << "  " << it->eta_s_right << "  "
                      << it->eta_s_left << std::endl;
        }
        j_mu[1] += local_eperp*sinh_perp*cos(phi_perp);
        j_mu[2] += local_eperp*sinh_perp*sin(phi_perp);
        j_mu[3] += local_eperp*sinh_long*cosh_perp;
    }

    for (auto const&it: QCD_strings_remnant_list_current_tau) {
        const double sigma_x = it->sigma_x;
        const double sigma_x_sq = sigma_x*sigma_x;
        const double sigma_eta = it->sigma_eta;
        const double alpsig = preEqFlowFactor_*sigma_x;
        const double prefactor_prep = (
            1./(2.*M_PI*(sigma_x_sq
                         + exp(alpsig*alpsig/2.)*sqrt(M_PI/2)
                           *alpsig*sigma_x_sq*erf(alpsig/sqrt(2.)))
               )
        );
        const double prefactor_etas = 1./(sqrt(2.*M_PI)*sigma_eta);
        const double prefactors = prefactor_prep*prefactor_etas;
        const double skip_dis_x = n_sigma_skip*sigma_x;
        const double skip_dis_eta = n_sigma_skip*sigma_eta;
        // add remnant energy at the string ends
        bool flag_left = false;
        if (   it->tau_end_left >= tau - dtau/2.
            && it->tau_end_left <  tau + dtau/2.
            && it->remnant_l > 0.) {
            flag_left = true;
        }

        bool flag_right = false;
        if (   it->tau_end_right >= tau - dtau/2.
            && it->tau_end_right <  tau + dtau/2.
            && it->remnant_r > 0.) {
            flag_right = true;
        }

        double x_dis = x - it->x_perp;
        if (std::abs(x_dis) > skip_dis_x) continue;

        double y_dis = y - it->y_perp;
        if (std::abs(y_dis) > skip_dis_x) continue;

        double exp_eta_s_left = 0.0;
        if (flag_left) {
            double eta_dis_left = std::abs(eta_s - it->eta_s_left);
            if (eta_dis_left < skip_dis_eta) {
                exp_eta_s_left = (exp(-eta_dis_left*eta_dis_left
                                      /(2.*sigma_eta*sigma_eta)));
            }
        }
        double exp_eta_s_right = 0.0;
        if (flag_right) {
            double eta_dis_right = std::abs(eta_s - it->eta_s_right);
            if (eta_dis_right < skip_dis_eta) {
                exp_eta_s_right = (exp(-eta_dis_right*eta_dis_right
                                       /(2.*sigma_eta*sigma_eta)));
            }
        }
        double cosh_long = exp_tau*(
              exp_eta_s_left*(it->remnant_l)*(it->E_remnant_norm_L)*cosh(it->y_l - eta_s)
            + exp_eta_s_right*(it->remnant_r)*(it->E_remnant_norm_R)*cosh(it->y_r - eta_s)
        );
        double sinh_long = exp_tau*(
              exp_eta_s_left*(it->remnant_l)*(it->E_remnant_norm_L)*sinh(it->y_l - eta_s)
            + exp_eta_s_right*(it->remnant_r)*(it->E_remnant_norm_R)*sinh(it->y_r - eta_s)
        );

        double cosh_perp = 1.;
        double sinh_perp = 0.;
        double phi_perp = 0.;
        double exp_xperp = 0.;
        if (cosh_long > 0) {
            exp_xperp = (prefactors*exp(-(x_dis*x_dis + y_dis*y_dis)
                                        /(2.*sigma_x*sigma_x)));
            cosh_perp = (
                cosh(preEqFlowFactor_*sqrt(x_dis*x_dis + y_dis*y_dis)));
            sinh_perp = (
                sinh(preEqFlowFactor_*sqrt(x_dis*x_dis + y_dis*y_dis)));
            phi_perp = atan2(y_dis, x_dis);
        }
        j_mu[0] += exp_xperp*cosh_long*cosh_perp;
        j_mu[1] += exp_xperp*sinh_perp*cos(phi_perp);
        j_mu[2] += exp_xperp*sinh_perp*sin(phi_perp);
        j_mu[3] += exp_xperp*sinh_long*cosh_perp;
    }
    const double prefactor_tau = 1./dtau;
    const double unit_convert = 1.0/Util::hbarc;
    const double prefactors = prefactor_tau*unit_convert;
    j_mu[0] *= prefactors;  // 1/fm^4
    j_mu[1] *= prefactors;  // 1/fm^4
    j_mu[2] *= prefactors;  // 1/fm^4
    j_mu[3] *= prefactors;  // 1/fm^4
}


double HydroSourceStrings::get_hydro_rhob_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu) const {
    double res = 0.;
    if (QCD_strings_baryon_list_current_tau.size() == 0) return(res);

    // flow velocity
    const double gamma_perp_flow = sqrt(1. + u_mu[1]*u_mu[1] + u_mu[2]*u_mu[2]);
    const double y_perp_flow     = acosh(gamma_perp_flow);
    const double y_long_flow     = asinh(u_mu[3]/gamma_perp_flow) + eta_s;
    const double sin_phi_flow    = u_mu[1]/gamma_perp_flow;
    const double cos_phi_flow    = u_mu[2]/gamma_perp_flow;
    const double dtau            = DATA.delta_tau;

    const double exp_tau      = 1.0/tau;
    const double n_sigma_skip = 8.;
    for (auto &it: QCD_strings_baryon_list_current_tau) {
        const double sigma_x = it->sigma_x;
        const double sigma_eta = it->sigma_eta;
        const double prefactor_prep = 1./(2.*M_PI*sigma_x*sigma_x);
        const double prefactor_etas = 1./(sqrt(2.*M_PI)*sigma_eta);
        const double skip_dis_x = n_sigma_skip*sigma_x;
        const double skip_dis_eta = n_sigma_skip*sigma_eta;

        // skip the evaluation if the strings is too far away in the
        // space-time grid
        // dumping energy into the medium from the active strings
        //double tau_dis_left = fabs(tau - it->tau_end_left);
        //double tau_dis_right = fabs(tau - it->tau_end_right);
        int flag_left = 0;
        if (   it->tau_baryon_left >= tau - dtau/2.
            && it->tau_baryon_left <  tau + dtau/2.
            && it->baryon_frac_l > 0.) {
            flag_left = 1;
        }

        int flag_right = 0;
        if (   it->tau_baryon_right >= tau - dtau/2.
            && it->tau_baryon_right <  tau + dtau/2.
            && it->baryon_frac_r > 0.) {
            flag_right = 1;
        }

        if (flag_left == 0 && flag_right == 0) continue;

        if (eta_s < it->eta_s_left - skip_dis_eta
                || eta_s > it->eta_s_right + skip_dis_eta) continue;

        double eta_frac_left = (
                (eta_s - it->eta_s_left)
                /std::max(Util::small_eps, it->eta_s_right - it->eta_s_left));
        double eta_frac_right = (
                (eta_s - it->eta_s_right)
                /std::max(Util::small_eps, it->eta_s_right - it->eta_s_left));
        eta_frac_left = std::max(0., std::min(1., eta_frac_left));
        eta_frac_right = std::max(0., std::min(1., eta_frac_right));

        const double x_perp_left = getStringTransverseCoord(
                                        it->x_pl, it->x_pr, eta_frac_left);
        const double x_perp_right = getStringTransverseCoord(
                                        it->x_pl, it->x_pr, eta_frac_right);
        const double x_dis_left  = x - x_perp_left;
        const double x_dis_right = x - x_perp_right;
        if (std::abs(x_dis_left) > skip_dis_x
                && std::abs(x_dis_right) > skip_dis_x) {
            continue;
        }

        const double y_perp_left = getStringTransverseCoord(
                                        it->y_pl, it->y_pr, eta_frac_left);
        const double y_perp_right = getStringTransverseCoord(
                                        it->y_pl, it->y_pr, eta_frac_right);
        const double y_dis_left  = y - y_perp_left;
        const double y_dis_right = y - y_perp_right;
        if (std::abs(y_dis_left) > skip_dis_x
                && std::abs(y_dis_right) > skip_dis_x) {
            continue;
        }

        double exp_eta_s_left = 0.0;
        if (flag_left == 1) {
            double eta_dis_left = std::abs(eta_s
                                           - it->eta_s_baryon_left);
            if (eta_dis_left < skip_dis_eta) {
                exp_eta_s_left = (exp(-eta_dis_left*eta_dis_left
                                      /(2.*sigma_eta*sigma_eta)));
            }
        }

        double exp_eta_s_right = 0.0;
        if (flag_right == 1) {
            double eta_dis_right = std::abs(eta_s
                                            - it->eta_s_baryon_right);
            if (eta_dis_right < skip_dis_eta) {
                exp_eta_s_right = (exp(-eta_dis_right*eta_dis_right
                                       /(2.*sigma_eta*sigma_eta)));
            }
        }

        double exp_xperp_l = exp(
            -(x_dis_left*x_dis_left + y_dis_left*y_dis_left)
            /(2.*sigma_x*sigma_x));
        double exp_xperp_r = exp(
            -(x_dis_right*x_dis_right + y_dis_right*y_dis_right)
            /(2.*sigma_x*sigma_x));

        double fsmear = exp_tau*(
                  exp_xperp_l*exp_eta_s_left*it->baryon_frac_l
                + exp_xperp_r*exp_eta_s_right*it->baryon_frac_r);
        if (fsmear > 0.) {
            double rapidity_local = (
                (  exp_eta_s_left*(it->baryon_frac_l)*(it->y_l_baryon)
                 + exp_eta_s_right*(it->baryon_frac_r)*(it->y_r_baryon))
                /(std::max(Util::small_eps,
                           (exp_eta_s_left*(it->baryon_frac_l)
                            + exp_eta_s_right*(it->baryon_frac_r)))
                 )
            );
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
            res += prefactor_etas*prefactor_prep*p_dot_u*fsmear;
        }
    }
    const double prefactor_tau  = 1./dtau;
    res *= prefactor_tau;
    return(res);
}


double HydroSourceStrings::getStringTransverseCoord(
            const double xl, const double xr, const double etaFrac) const {
    double xT = ((xl + xr)/2.
                 + stringTransverseShiftFrac_*(0.5 - etaFrac)*(xl - xr)/2.);
    return(xT);
}

