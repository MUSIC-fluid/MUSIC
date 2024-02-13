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

HydroSourceSMASH::HydroSourceSMASH(InitData &DATA_in) :
    DATA(DATA_in) {
    set_source_tau_min(100.0);
    set_source_tau_max(0.0);
    set_sigma_tau(0.1);
    set_sigma_x  (0.5);
    set_sigma_eta(0.2);
    parton_quench_factor = 1.;    // no diffusion current from the source
    int i_event = DATA.event_id_SMASH_output;
    if (i_event < 1) {
        i_event = 1;
    }
    average_events_ = DATA.average_SMASH_events;
    int extended_output = DATA.extended_SMASH_output;
    int reject_spectators = DATA.reject_SMASH_spectators;
    read_in_SMASH_hadrons(i_event, extended_output, reject_spectators);
}


HydroSourceSMASH::~HydroSourceSMASH() {
    list_hadrons_.clear();
    list_hadrons_current_tau_.clear();
    list_spectators_.clear();
}


//! This function reads in the hadron information from the SMASH model
void HydroSourceSMASH::read_in_SMASH_hadrons(int i_event,
    int extended_output, int reject_spectators) {
    list_hadrons_.clear();
    list_spectators_.clear();
    string SMASH_filename = DATA.initName_SMASH;
    music_message << "hydro_source: "
                  << "read in SMASH parton list from " << SMASH_filename;
    music_message.flush("info");
    char commentChar ='#';

    // create a stream to file.
    std::ifstream SMASHinputFile(SMASH_filename);

    // Check if file opened properly
    if(!SMASHinputFile.is_open()){
        music_message << "hydro_source::read_in_SMASH_hadrons: "
                      << "can not open the SMASH IC file: " << SMASH_filename;
        music_message.flush("error");
        exit(1);
    }

    // Show useful info from the file header (three first lines)
    std::string line;
    std::string OscarVersion, SMASHVersion, dummy;

    //first line
    std::istringstream iss(line);
    std::getline(SMASHinputFile, line);
    (iss.str(line), iss.clear(), iss >> OscarVersion);
     
    // Ignore second line 
    std::getline(SMASHinputFile, line);

    // third line
    std::getline(SMASHinputFile, line);
    (iss.str(line), iss.clear(), iss >> dummy, iss >> SMASHVersion);

    music_message << "SMASH version: "
        << SMASHVersion
        << "Oscar version: "
        << OscarVersion;
    music_message.flush("info");

    // Define useful variables
    int n_hadrons = 0;
    int n_spectators = 0;

    baryon_total_ = 0;charge_total_ = 0;strangeness_total_ = 0;
    p0_total_ = 0.;px_total_ = 0.;py_total_ = 0.;pz_total_ = 0.;

    number_events_ = 0;

    // now we read in data
    for (int j_ev = 1; j_ev <= i_event; j_ev++) {
        // Ignore event header
        std::getline(SMASHinputFile, line);
        bool EndofEvent = false;
        bool event_of_interest = average_events_ || (j_ev == i_event);
        int n_hadrons_ev = 0; // Number of hadron in the current event
        int n_spectators_ev = 0;  // number of spectators in the current event
        while(!EndofEvent){
            std::getline(SMASHinputFile, line);
            // Check if end of the event
            if(line[0] == commentChar){EndofEvent = true;}
            // Avoid empty lines
            if(line.empty()) continue;
            // Only take on event at index i_event in SMASH Oscar
            // Or average over all events until i_event if 
            // average_events_ = true
            if(!event_of_interest) continue;
            n_hadrons_ev +=1;

            

            // Declare hadron properties
            double t, x, y, z, mass;
            double p0, px, py, pz;
            double form_time, xsecfac, time_last_coll;
            int pdgid, ID, tag, ncoll;
            int B, Q, S;
            int proc_id_origin, proc_type_origin;
            int pdg_mother1, pdg_mother2;

            if(extended_output){
                // t x y z mass p0 px py pz pdg ID charge 
                // ncoll form_time xsecfac proc_id_origin 
                // proc_type_origin time_last_coll 
                // pdg_mother1 pdg_mother2 baryon_number strangeness
                iss.str(line);
                iss.clear();
                iss >> t >> x >> y >> z >> mass;
                iss >> p0 >> px >> py >> pz; 
                iss >> pdgid >> ID >> Q >> ncoll;
                iss >> form_time >> xsecfac;
                iss >> proc_id_origin >> proc_type_origin; 
                iss >> time_last_coll;
                iss >> pdg_mother1 >> pdg_mother2;
                iss >> B >> S;
            }
            else{
                // t x y z mass p0 px py pz pdg ID charge 
                iss.str(line);
                iss.clear();
                iss >> t >> x >> y >> z >> mass;
                iss >> p0 >> px >> py >> pz; 
                iss >> pdgid >> ID >> Q;
                // Ncoll not consider here!
                ncoll = 1.0; // Consider all hadrons to be participants.
            }
            // Discard hadrons outside light-cone
            if (std::fabs(t) <= std::fabs(z)) continue;

            // Spectators are not considered.
            if(ncoll == 0.0){
                n_spectators_ev +=1;
                continue;
            }
            // cuts in rapidity.
            //double rapidity = 0.5 * std::log((p0 + pz) / (p0 - pz));
            //double eta_s =  0.5 * std::log((t + z) / (t - z));
            //if(eta_s > 2.0) continue;
            //if(eta_s < -2.0) continue;

            // Put info in a new hadron
            hadron new_hadron;
            new_hadron.pdgid = pdgid;
            new_hadron.ID = ID;
            new_hadron.tau = std::sqrt(t * t - z * z);
            new_hadron.x = x;
            new_hadron.y = y;
            new_hadron.eta_s = 0.5 * std::log((t + z) / (t - z));
            new_hadron.rapidity = 0.5 * std::log((p0 + pz) / (p0 - pz));
            new_hadron.E = p0;
            new_hadron.px = px;
            new_hadron.py = py;
            new_hadron.mass = mass;
            new_hadron.baryon_number = B;
            new_hadron.electric_charge = Q;
            new_hadron.strangeness = S;
            // Put hadron in hadron list.
            // Note, if average_events_ = 1
            // list_hadrons contains hadrons
            // for all event in IC file.
            list_hadrons_.push_back(new_hadron);

            // get/set source tau in the HydroSourceBase class
            if (get_source_tau_max() < new_hadron.tau) {
                    set_source_tau_max(new_hadron.tau);
                }
            if (get_source_tau_min() > new_hadron.tau) {
                    set_source_tau_min(new_hadron.tau);
                }

            // Check if hadronic input is within hydro grid
            bool out_of_range_x = std::fabs(new_hadron.x) > 0.5 * DATA.x_size;
            bool out_of_range_y = std::fabs(new_hadron.y) > 0.5 * DATA.y_size;
            bool out_of_range_eta = std::fabs(new_hadron.eta_s) > 0.5 * DATA.eta_size;

            if (out_of_range_x || out_of_range_y || out_of_range_eta) {
                music_message << "HydroSourceSMASH:: hadronic source is out of range.";
                music_message.flush("info");
                music_message << "HydroSourceSMASH::     pdgid = " << new_hadron.pdgid;
                music_message.flush("info");
                music_message << "HydroSourceSMASH::     (x, y) = (" << new_hadron.x << ", " << new_hadron.y << ") fm.";
                music_message.flush("info");
                music_message << "HydroSourceSMASH::     eta_s = " << new_hadron.eta_s;
                music_message.flush("info");
            }

            // Compute total conserved charge numbers 
            baryon_total_ += new_hadron.baryon_number;
            charge_total_ += new_hadron.electric_charge;
            strangeness_total_ += new_hadron.strangeness;

            // compute total momenta 
            p0_total_ += p0;px_total_ += px;
            py_total_ += py;pz_total_ += pz;
        }
        // Add number of participant hadrons in event in total number of hadrons.
        n_hadrons += n_hadrons_ev;
        // Add number of spectator hadrons in event in total number of hadrons.
        n_hadrons += n_spectators_ev;
        // Count event if it's an event of interest.
        if(event_of_interest){number_events_ += 1;}



    }// End of event loop.
    weight_event_ = 1. / (double)number_events_;
    music_message << "HydroSourceSMASH:: read in " << list_hadrons_.size() << "/"
                  << n_hadrons << " hadrons and " << n_spectators << " spectators from "
                  << number_events_ << " events.";
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: tau_min = " << get_source_tau_min()
                  << " fm.";
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: tau_max = " << get_source_tau_max()
                  << " fm.";
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: p0_total = " << p0_total_ << " GeV.";
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: px_total = " << px_total_ << " GeV.";
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: py_total = " << py_total_ << " GeV.";
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: pz_total = " << pz_total_ << " GeV.";
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: baryon_total = " << baryon_total_;
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: electric_charge_total = " << charge_total_;
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: strangeness_total = " << strangeness_total_;
    music_message.flush("info");
}


void HydroSourceSMASH::prepare_list_for_current_tau_frame(
                                                const double tau_local) {
    double dtau = DATA.delta_tau;
    list_hadrons_current_tau_.clear();
    int n_hadrons_all = list_hadrons_.size();
    for (int ipart = 0; ipart < n_hadrons_all; ipart++) {
        double tau_dis = list_hadrons_.at(ipart).tau - tau_local;
        if (tau_dis >= 0. && tau_dis < dtau) {
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
        j_mu[0] *= prefactor * weight_event_;
        j_mu[1] *= prefactor * weight_event_;
        j_mu[2] *= prefactor * weight_event_;
        j_mu[3] *= prefactor * weight_event_;
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
        res *= prefactor_tau * prefactor_prep * prefactor_etas * weight_event_;
    }
    return res;
}

double HydroSourceSMASH::get_hydro_rhoq_source(
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
            int qnumber = (int)list_hadrons_current_tau_.at(ipart).electric_charge;
            if (qnumber == 0) {
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
            res += p_dot_u * f_smear * (int)qnumber;
        }
        res *= prefactor_tau * prefactor_prep * prefactor_etas * weight_event_;
    }
    return res;
}

double HydroSourceSMASH::get_hydro_rhos_source(
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
            int snumber = (int)list_hadrons_current_tau_.at(ipart).strangeness;
            // mesons do not count towards the net baryon current.
            if (snumber == 0) {
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
            res += p_dot_u * f_smear * (int)snumber;
        }
        res *= prefactor_tau * prefactor_prep * prefactor_etas * weight_event_;
    }
    return res;
}
