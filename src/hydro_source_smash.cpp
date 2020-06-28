// Copyright 2019 Chun Shen

#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <memory>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "hydro_source_smash.h"
#include "util.h"

using std::string;

HydroSourceSMASH::HydroSourceSMASH(const InitData &DATA_in) :
    DATA(DATA_in) {
    set_source_tau_min(100.0);
    set_source_tau_max(0.0);
    set_sigma_tau(0.1);
    set_sigma_x  (0.5);
    set_sigma_eta(0.2);
    parton_quench_factor = 1.;    // no diffusion current from the source
    int i_event = DATA.event_id_SMASH;
    read_in_SMASH_hadrons(i_event);
}


HydroSourceSMASH::~HydroSourceSMASH() {
    list_hadrons_.clear();
    list_hadrons_current_tau_.clear();
}


//! This function reads in the hadron information from the SMASH model
void HydroSourceSMASH::read_in_SMASH_hadrons(int i_event) {
    list_hadrons_.clear();
    string SMASH_filename = DATA.initName_SMASH;
    music_message << "hydro_source: "
                  << "read in SMASH parton list from " << SMASH_filename;
    music_message.flush("info");

    FILE *fin;
    fin = fopen(SMASH_filename.c_str(), "r");

    string text_string;

    if (fin == NULL) {
        music_message << "hydro_source::read_in_SMASH_hadrons: "
                      << "can not open the AMPT file: " << SMASH_filename;
        music_message.flush("error");
        exit(1);
    }

    int n_hadrons = 0;

    // reading the file header
    char line1_header[200];
    char line2_header[200];
    char line3_header[200];
    fgets(line1_header, 200, fin);
    fgets(line2_header, 200, fin);
    fgets(line3_header, 200, fin);

    // now we read in data
    for (int j_ev = 1; j_ev <= i_event; j_ev++) {
        // reading the event header
        char line_header_ev[200];
        fgets(line_header_ev, 200, fin);

        // read the event id and number of particles in the given event
        char entry1_dummy[10];
        char entry2_dummy[10];
        char entry3_dummy[10];
        int k_ev = 0;
        int n_had_init = 0;
        sscanf(line_header_ev, "%s %s %d %s %d", entry1_dummy, entry2_dummy, &k_ev, entry3_dummy, &n_had_init);

        char line_particle[500];
        bool end_of_event = false;
        while (!end_of_event && feof(fin) == 0) {
            fgets(line_particle, 500, fin);

            if (*line_particle == '#') {
                end_of_event = true;
                break;
            }

            double t;
            double x;
            double y;
            double z;
            double mass;
            double p0;
            double px;
            double py;
            double pz;
            int pdgid;
            int tag;
            int charge;

            sscanf(line_particle, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d",
                &t, &x, &y, &z, &mass, &p0, &px, &py, &pz, &pdgid, &tag, &charge);

            if (j_ev != i_event) {
                continue;
            }

            n_hadrons += 1;

            if (std::abs(t) <= std::abs(z)) {
                continue;
            }

            hadron new_hadron;
            new_hadron.pdgid = pdgid;
            new_hadron.tau = std::sqrt(t * t - z * z);
            new_hadron.x = x;
            new_hadron.y = y;
            new_hadron.eta_s = 0.5 * std::log((t + z) / (t - z));
            new_hadron.rapidity = 0.5 * std::log((p0 + pz) / (p0 - pz));
            new_hadron.E = p0;
            new_hadron.px = px;
            new_hadron.py = py;
            new_hadron.mass = mass;
            new_hadron.baryon_number = Util::get_baryon_number_from_pdg(pdgid);
            new_hadron.strangness = -Util::get_net_quark_from_pdg(pdgid, 3);
            new_hadron.electric_charge = charge;

            list_hadrons_.push_back(new_hadron);
            if (get_source_tau_max() < new_hadron.tau) {
                set_source_tau_max(new_hadron.tau);
            }
            if (get_source_tau_min() > new_hadron.tau) {
                set_source_tau_min(new_hadron.tau);
            }
        }
    }

    fclose(fin);

    music_message << "HydroSourceSMASH:: read in " << list_hadrons_.size() << "/"
                  << n_hadrons << " hadrons.";
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: tau_min = " << get_source_tau_min()
                  << " fm.";
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: tau_max = " << get_source_tau_max()
                  << " fm.";
    music_message.flush("info");
}


void HydroSourceSMASH::prepare_list_for_current_tau_frame(
                                                const double tau_local) {
    double dtau = DATA.delta_tau;
    list_hadrons_current_tau_.clear();
    int n_hadrons_all = list_hadrons_.size();
    for (int ipart = 0; ipart < n_hadrons_all; ipart++) {
        double tau_dis = list_hadrons_.at(ipart).tau - tau_local;
        if (tau_dis > 0. && tau_dis < dtau) {
            list_hadrons_current_tau_.push_back(list_hadrons_.at(ipart));
        }
    }
    music_message << "hydro_source: tau = " << tau_local
                  << " number of source: "
                  << list_hadrons_current_tau_.size();
    music_message.flush("info");
}


void HydroSourceSMASH::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s, 
    const FlowVec &u_mu, EnergyFlowVec &j_mu) const {
    j_mu = {0};
    if (list_hadrons_current_tau_.size() == 0) return;
    
    const double sigma_tau = get_sigma_tau();
    const double sigma_x = get_sigma_x();
    const double sigma_eta = get_sigma_eta();

    // flow velocity
    const double dtau = DATA.delta_tau;

    const double prefactor_prep = 1. / (M_PI * sigma_x * sigma_x);
    const double prefactor_tau = 1. / dtau;
    const double prefactor_etas = 1. / (sqrt(M_PI) * sigma_eta);
    const double n_sigma_skip = 5.;
    const double skip_dis_x = n_sigma_skip * sigma_x;
    const double skip_dis_eta = n_sigma_skip * sigma_eta;
    const double exp_tau = 1. / tau;

    // SMASH hadron sources
    double tau_dis_max = tau - get_source_tau_max();
    if (tau_dis_max < n_sigma_skip * sigma_tau) {
        int n_hadrons_now = list_hadrons_current_tau_.size();
        for (int ipart = 0; ipart < n_hadrons_now; ipart++) {
            double x_dis = x - list_hadrons_current_tau_.at(ipart).x;
            if (std::abs(x_dis) > skip_dis_x) continue;

            double y_dis = y - list_hadrons_current_tau_.at(ipart).y;
            if (std::abs(y_dis) > skip_dis_x) continue;

            double eta_s_dis = eta_s - list_hadrons_current_tau_.at(ipart).eta_s;
            if (std::abs(eta_s_dis) > skip_dis_eta) continue;

            double exp_xperp =
                exp(-(x_dis * x_dis + y_dis * y_dis) / (sigma_x * sigma_x));
            double exp_eta_s =
                exp(-eta_s_dis * eta_s_dis / (sigma_eta * sigma_eta));

            double f_smear = exp_tau * exp_xperp * exp_eta_s;
            double px_now = list_hadrons_current_tau_.at(ipart).px;
            double py_now = list_hadrons_current_tau_.at(ipart).py;
            double p_perp_sq = px_now * px_now + py_now * py_now;
            double mass_now = list_hadrons_current_tau_.at(ipart).mass;
            double m_perp = sqrt(mass_now * mass_now + p_perp_sq);
            j_mu[0] += m_perp * cosh(list_hadrons_current_tau_.at(ipart).rapidity - eta_s) * f_smear;
            j_mu[1] += px_now * f_smear;
            j_mu[2] += py_now * f_smear;
            j_mu[3] += m_perp * sinh(list_hadrons_current_tau_.at(ipart).rapidity - eta_s) * f_smear;
        }
        double norm = DATA.sFactor / Util::hbarc;     // 1/fm^4
        double prefactor = norm * prefactor_tau * prefactor_prep * prefactor_etas;
        j_mu[0] *= prefactor;
        j_mu[1] *= prefactor;
        j_mu[2] *= prefactor;
        j_mu[3] *= prefactor;
    }
}


double HydroSourceSMASH::get_hydro_rhob_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu) const {
    double res = 0.;
    if (list_hadrons_current_tau_.size() == 0) return(res);
    
    const double sigma_tau = get_sigma_tau();
    const double sigma_x = get_sigma_x();
    const double sigma_eta = get_sigma_eta();

    // flow velocity
    const double gamma_perp_flow = sqrt(1. + u_mu[1] * u_mu[1] + u_mu[2] * u_mu[2]);
    const double y_perp_flow = acosh(gamma_perp_flow);
    const double y_long_flow = asinh(u_mu[3] / gamma_perp_flow) + eta_s;
    const double sinh_y_perp_flow = sinh(y_perp_flow);
    const double dtau = DATA.delta_tau;

    const double exp_tau = 1.0 / tau;
    const double n_sigma_skip = 5.;
    const double prefactor_prep = 1. / (M_PI * sigma_x * sigma_x);
    const double prefactor_etas = 1. / (sqrt(M_PI) * sigma_eta);
    const double prefactor_tau = 1. / dtau;
    const double skip_dis_x = n_sigma_skip * sigma_x;
    const double skip_dis_eta = n_sigma_skip * sigma_eta;

    double tau_dis_max = tau - get_source_tau_max();
    if (tau_dis_max < n_sigma_skip * sigma_tau) {
        int n_hadrons_now = list_hadrons_current_tau_.size();
        for (int ipart = 0; ipart < n_hadrons_now; ipart++) {
            int bnumber = (int)list_hadrons_current_tau_.at(ipart).baryon_number;
            // mesons do not count towards the net baryon current.
            if (bnumber == 0) {
                continue;
            }

            // skip the evaluation if the strings is too far away in the
            // space-time grid
            double x_dis = x - list_hadrons_current_tau_.at(ipart).x;
            if (std::abs(x_dis) > skip_dis_x) continue;

            double y_dis = y - list_hadrons_current_tau_.at(ipart).y;
            if (std::abs(y_dis) > skip_dis_x) continue;

            double eta_s_dis = eta_s - list_hadrons_current_tau_.at(ipart).eta_s;
            if (std::abs(eta_s_dis) > skip_dis_eta) continue;

            double exp_xperp =
                exp(-(x_dis * x_dis + y_dis * y_dis) / (sigma_x * sigma_x));
            double exp_eta_s =
                exp(-eta_s_dis * eta_s_dis / (sigma_eta * sigma_eta));

            double f_smear = exp_tau * exp_xperp * exp_eta_s;
            double px_now = list_hadrons_current_tau_.at(ipart).px;
            double py_now = list_hadrons_current_tau_.at(ipart).py;
            double p_perp_sq = px_now * px_now + py_now * py_now;
            double mass_now = list_hadrons_current_tau_.at(ipart).mass;
            double u_perp = sqrt(p_perp_sq) / mass_now;
            double rapidity_perp = asinh(u_perp);

            double y_dump_long =
                (1. - parton_quench_factor) * list_hadrons_current_tau_.at(ipart).rapidity +
                parton_quench_factor * y_long_flow;
            double y_dump_perp =
                (1. - parton_quench_factor) * rapidity_perp +
                parton_quench_factor * y_perp_flow;
            double p_dot_u = u_mu[0] * (u_mu[0] -
                tanh(y_dump_perp) * sinh_y_perp_flow / cosh(y_dump_long - eta_s) -
                tanh(y_dump_long - eta_s) * u_mu[3]);
            res += p_dot_u * f_smear * (int)bnumber;
        }
        res *= prefactor_tau * prefactor_prep * prefactor_etas;
    }
    return res;
}
