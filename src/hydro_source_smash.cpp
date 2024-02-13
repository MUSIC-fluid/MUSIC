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

    //FILE *fin;
    //fin = fopen(SMASH_filename.c_str(), "r");

    // create a stream to file.
    std::ifstream SMASHinputFile(SMASH_filename);

    // Check if file opened properly
    if(!SMASHinputFile.is_open()){
        music_message << "hydro_source::read_in_SMASH_hadrons: "
                      << "can not open the SMASH IC file: " << SMASH_filename;
        music_message.flush("error");
        exit(1);
    }

    //if (fin == NULL) {
    //    music_message << "hydro_source::read_in_SMASH_hadrons: "
    //                  << "can not open the SMASH IC file: " << SMASH_filename;
    //    music_message.flush("error");
    //    exit(1);
    //}


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
    // reading the file header
    //char line1_header[200];
    //char line2_header[200];
    //char line3_header[200];
    //fgets(line1_header, 200, fin);
    //fgets(line2_header, 200, fin);
    //fgets(line3_header, 200, fin);

    // Define useful variables
    int n_hadrons = 0;
    int n_spectators = 0;

    baryon_total_ = 0;charge_total_ = 0;strangeness_total_ = 0;
    p0_total_ = 0.;px_total_ = 0.;py_total_ = 0.;pz_total_ = 0.;

    //char line_event_spec_ini[500];
    //char line_event_spec_fin[500];

    //if (reject_spectators != 0) {
    //    FILE *fout_spectators;
    //    fout_spectators = fopen("SMASH_spectators.oscar", "w");
    //    fprintf(fout_spectators, "#!OSCAR2013 SMASH_spectators t x y z mass p0 px py pz pdg ID charge\n");
    //    fprintf(fout_spectators, "# Units: fm fm fm fm GeV GeV GeV GeV GeV none none none\n");
    //    fprintf(fout_spectators, "%s", line3_header);
    //    fclose(fout_spectators);
    //}

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

            // SPECTATORS ARE NOT CONSIDERED.
            //n_spectators_ev +=1; Considering all hadrons to be participants
            
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
            }
            // Discard hadrons outside light-cone
            if (std::fabs(t) <= std::fabs(z)) continue;

            // Put info in a new hadron
            hadron new_hadron;
            new_hadron.pdgid = pdgid;
            new_hadron.ID = ID;
            new_hadron.tau = std::sqrt(t * t - z * z);
            new_hadron.x = x;
            new_hadron.y = y;
            new_hadron.eta_s = 0.5 * std::log((t + z) / (t - z));
            new_hadron.rapidity = 0.5 * std::log((p0 + pz) / (p0 - pz));
            new_hadron.pz = pz;
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

            // get/set source in the HydroSourceBase class
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
        // Add number of hadrons in event in total number of hadrons. 
        n_hadrons += n_hadrons_ev;
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

    //for (int j_ev = 1; j_ev <= i_event; j_ev++) {  // event loop
    //    // reading the event header
    //    char line_header_ev[200];
    //    fgets(line_header_ev, 200, fin);

    //    // read the event id and number of particles in the given event
    //    char entry1_dummy[10];
    //    char entry2_dummy[10];
    //    char entry3_dummy[10];
    //    int k_ev = 0;
    //    int n_had_init = 0;
    //    sscanf(line_header_ev, "%s %s %d %s %d",
    //        entry1_dummy, entry2_dummy, &k_ev, entry3_dummy, &n_had_init);

    //    int n_hadrons_ev = 0;     // number of hadrons in the current event
    //    int n_spectators_ev = 0;  // number of spectators in the current event

    //    char line_particle[500];
    //    bool end_of_event = false;
    //    bool event_of_interest = average_events_ != 0 || j_ev == i_event;
    //    while (!end_of_event && feof(fin) == 0) {  // particle loop
    //        fgets(line_particle, 500, fin);

    //        if (*line_particle == '#') {
    //            end_of_event = true;
    //            if (event_of_interest) {
    //                strcpy(line_event_spec_fin, line_particle);
    //            }

    //            break;
    //        }

    //        if (!event_of_interest) {
    //            continue;
    //        }

    //        double t;
    //        double x;
    //        double y;
    //        double z;
    //        double mass;
    //        double p0;
    //        double px;
    //        double py;
    //        double pz;
    //        int pdgid;
    //        int tag;
    //        int charge;

    //        int ncoll = 0;
    //        double form_time = 0.;
    //        double xsecfac = 1.;
    //        int proc_id_origin = 0;
    //        int proc_type_origin = 0;
    //        double time_last_coll = 0.;
    //        int pdg_mother1 = 0;
    //        int pdg_mother2 = 0;

    //        bool is_spectator_now = false;

    //        if (extended_output == 0) {
    //            sscanf(line_particle, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d",
    //                &t, &x, &y, &z, &mass, &p0, &px, &py, &pz, &pdgid, &tag, &charge);

    //            is_spectator_now =
    //                reject_spectators != 0 && std::abs(px) < 1e-4 && std::abs(py) < 1e-4;
    //        } else {
    //            sscanf(line_particle, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %lf %lf %d %d %lf %d %d",
    //                &t, &x, &y, &z, &mass, &p0, &px, &py, &pz, &pdgid, &tag, &charge,
    //                &ncoll, &form_time, &xsecfac, &proc_id_origin, &proc_type_origin, &time_last_coll,
    //                &pdg_mother1, &pdg_mother2);

    //            is_spectator_now = reject_spectators != 0 && proc_id_origin == 0;
    //        }

    //        if (is_spectator_now) {
    //            n_spectators_ev += 1;
    //        } else {
    //            n_hadrons_ev += 1;
    //        }

    //        if (std::abs(t) <= std::abs(z)) {
    //            continue;
    //        }

    //        hadron new_hadron;
    //        new_hadron.pdgid = pdgid;
    //        new_hadron.tau = std::sqrt(t * t - z * z);
    //        new_hadron.x = x;
    //        new_hadron.y = y;
    //        new_hadron.eta_s = 0.5 * std::log((t + z) / (t - z));
    //        new_hadron.rapidity = 0.5 * std::log((p0 + pz) / (p0 - pz));
    //        new_hadron.E = p0;
    //        new_hadron.px = px;
    //        new_hadron.py = py;
    //        new_hadron.mass = mass;
    //        new_hadron.baryon_number = Util::get_baryon_number_from_pdg(pdgid);
    //        new_hadron.strangness = -Util::get_net_quark_from_pdg(pdgid, 3);
    //        new_hadron.electric_charge = charge;

    //        if (is_spectator_now) {
    //            list_spectators_.push_back(new_hadron);
    //        } else {
    //            list_hadrons_.push_back(new_hadron);
    //            if (get_source_tau_max() < new_hadron.tau) {
    //                set_source_tau_max(new_hadron.tau);
    //            }
    //            if (get_source_tau_min() > new_hadron.tau) {
    //                set_source_tau_min(new_hadron.tau);
    //            }
    //        }

    //        bool out_of_range_x = fabs(new_hadron.x) > 0.5 * DATA.x_size;
    //        bool out_of_range_y = fabs(new_hadron.y) > 0.5 * DATA.y_size;
    //        bool out_of_range_eta = fabs(new_hadron.eta_s) > 0.5 * DATA.eta_size;
    //        if (out_of_range_x || out_of_range_y || out_of_range_eta) {
    //            music_message << "HydroSourceSMASH:: hadronic source is out of range.";
    //            music_message.flush("info");
    //            music_message << "HydroSourceSMASH::     pdgid = " << new_hadron.pdgid;
    //            music_message.flush("info");
    //            music_message << "HydroSourceSMASH::     (x, y) = (" << new_hadron.x << ", " << new_hadron.y << ") fm.";
    //            music_message.flush("info");
    //            music_message << "HydroSourceSMASH::     eta_s = " << new_hadron.eta_s;
    //            music_message.flush("info");
    //        }

    //        baryon_total_ += new_hadron.baryon_number;
    //        p0_total_ += p0;
    //        px_total_ += px;
    //        py_total_ += py;
    //        pz_total_ += pz;
    //    }  // end of the particle loop

    //    n_hadrons += n_hadrons_ev;
    //    n_spectators += n_spectators_ev;

    //    sprintf(line_event_spec_ini, "# event %d out %d", j_ev, n_spectators_ev);
    //    // print out the spectator list if necessary
    //    if (reject_spectators != 0 && event_of_interest) {
    //        FILE *fout_spectators;
    //        fout_spectators = fopen("SMASH_spectators.oscar", "a");
    //        fprintf(fout_spectators, "%s\n", line_event_spec_ini);
    //        int ipart_ini = n_spectators - n_spectators_ev;
    //        int ipart_fin = n_spectators - 1;
    //        for (int ipart = ipart_ini; ipart <= ipart_fin; ipart++) {
    //            double t = list_spectators_.at(ipart).tau * std::cosh(list_spectators_.at(ipart).eta_s);
    //            double x = list_spectators_.at(ipart).x;
    //            double y = list_spectators_.at(ipart).y;
    //            double z = list_spectators_.at(ipart).tau * std::sinh(list_spectators_.at(ipart).eta_s);
    //            double mass = list_spectators_.at(ipart).mass;
    //            double px = list_spectators_.at(ipart).px;
    //            double py = list_spectators_.at(ipart).py;
    //            double mT = std::sqrt(mass * mass + px * px + py * py);
    //            double p0 = mT * std::cosh(list_spectators_.at(ipart).rapidity);
    //            double pz = mT * std::sinh(list_spectators_.at(ipart).rapidity);
    //            int pdgid = list_spectators_.at(ipart).pdgid;
    //            int tag = ipart;
    //            int charge = list_spectators_.at(ipart).electric_charge;

    //            fprintf(fout_spectators, "%e  %e  %e  %e  %e  %e  %e  %e  %e  %d  %d  %d\n",
    //                t, x, y, z, mass, p0, px, py, pz, pdgid, tag, charge);
    //        }
    //        fprintf(fout_spectators, "%s", line_event_spec_fin);
    //        fclose(fout_spectators);
    //    }

    //    if (event_of_interest) {
    //        number_events_ += 1;
    //    }

    //    if (feof(fin) != 0) {
    //        break;
    //    }
    //}  // end of the event loop

    //weight_event_ = 1. / (double)number_events_;

    //fclose(fin);

    //music_message << "HydroSourceSMASH:: read in " << list_hadrons_.size() << "/"
    //              << n_hadrons << " hadrons and " << n_spectators << " spectators from "
    //              << number_events_ << " events.";
    //music_message.flush("info");
    //music_message << "HydroSourceSMASH:: tau_min = " << get_source_tau_min()
    //              << " fm.";
    //music_message.flush("info");
    //music_message << "HydroSourceSMASH:: tau_max = " << get_source_tau_max()
    //              << " fm.";
    //music_message.flush("info");
    //music_message << "HydroSourceSMASH:: p0_total = " << p0_total_ << " GeV.";
    //music_message.flush("info");
    //music_message << "HydroSourceSMASH:: px_total = " << px_total_ << " GeV.";
    //music_message.flush("info");
    //music_message << "HydroSourceSMASH:: py_total = " << py_total_ << " GeV.";
    //music_message.flush("info");
    //music_message << "HydroSourceSMASH:: pz_total = " << pz_total_ << " GeV.";
    //music_message.flush("info");
    //music_message << "HydroSourceSMASH:: baryon_total = " << baryon_total_;
    //music_message.flush("info");
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


double HydroSourceSMASH::covariant_smearing_kernel(const double x_diff, 
        const double y_diff, const double eta_diff, const double vx, 
        const double vy, const double veta, const double sigma, 
        const double tau) const {
    // Calculate the squared distance and scalar product
    // tau factors from Milne metric
    const double r_squared = x_diff * x_diff + y_diff * y_diff 
                        + eta_diff * eta_diff * tau * tau;
    const double ru_scalar = x_diff * vx + y_diff * vy + eta_diff * veta * tau * tau;

    // Compute and return the smoothing kernel
    return exp((-r_squared - ru_scalar * ru_scalar) / (sigma * sigma));
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
    const double prefactor_cov = 1. / pow(2.*M_PI*sigma_x*sigma_x, 3./2.);
    double val_smearing_kernel = 0.;

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

            // kinematic quantities
            const double mass = list_hadrons_current_tau_.at(ipart).mass;
            const double E = list_hadrons_current_tau_.at(ipart).E;
            const double px = list_hadrons_current_tau_.at(ipart).px;
            const double py = list_hadrons_current_tau_.at(ipart).py;
            const double pz = list_hadrons_current_tau_.at(ipart).pz;
            const double mT = sqrt(mass*mass + px*px + py*py);

            if (covariant_smearing_kernel_) {
                // computations for covariant smearing kernel
                const double vx = px / mass;
                const double vy = py / mass;
                const double vz = pz / mass;
                const double gamma = 1. / sqrt(1. - vx*vx - vy*vy - vz*vz);
                const double rapidity = 0.5 * std::log((E + pz) / (E - pz));
                const double peta = mT * sinh(rapidity - list_hadrons_current_tau_.at(ipart).eta_s);
                const double veta = peta / mass;

                val_smearing_kernel = exp_tau * gamma * covariant_smearing_kernel(x_dis,y_dis,eta_s_dis,vx,vy,veta,sigma_x,tau);
            } else {
                const double exp_xperp =
                    exp(-(x_dis * x_dis + y_dis * y_dis) / (sigma_x * sigma_x));
                const double exp_eta_s =
                    exp(-eta_s_dis * eta_s_dis / (sigma_eta * sigma_eta));

                val_smearing_kernel = exp_tau * exp_xperp * exp_eta_s;
            }

            j_mu[0] += mT * cosh(list_hadrons_current_tau_.at(ipart).rapidity - eta_s) * val_smearing_kernel;
            j_mu[1] += px * val_smearing_kernel;
            j_mu[2] += py * val_smearing_kernel;
            j_mu[3] += mT * sinh(list_hadrons_current_tau_.at(ipart).rapidity - eta_s) * val_smearing_kernel;
        }
        double norm = DATA.sFactor / Util::hbarc; // 1/fm^4
        double prefactor = norm;
        if (covariant_smearing_kernel_) {
            prefactor *= prefactor_tau * prefactor_cov;
        } else {
            prefactor *= prefactor_tau * prefactor_prep * prefactor_etas;
        }
        j_mu[0] *= prefactor * weight_event_;
        j_mu[1] *= prefactor * weight_event_;
        j_mu[2] *= prefactor * weight_event_;
        j_mu[3] *= prefactor * weight_event_;
    }
}

double HydroSourceSMASH::calculate_source(const double tau, const double x, 
                                        const double y, const double eta_s,
                                        const FlowVec &u_mu, 
                                        const int quantityType) const {
    double result = 0.0;

    if (list_hadrons_current_tau_.empty()) return result;

    const double sigma_tau = get_sigma_tau();
    const double sigma_x = get_sigma_x();
    const double sigma_eta = get_sigma_eta();

    // Extract flow velocity components
    const double gamma_perp_flow = sqrt(1.0 + u_mu[1] * u_mu[1] + u_mu[2] * u_mu[2]);
    const double y_perp_flow = acosh(gamma_perp_flow);
    const double y_long_flow = asinh(u_mu[3] / gamma_perp_flow) + eta_s;
    const double sinh_y_perp_flow = sinh(y_perp_flow);
    const double dtau = DATA.delta_tau;

    const double exp_tau = 1.0 / tau;
    const double prefactor_prep = 1.0 / (M_PI * sigma_x * sigma_x);
    const double prefactor_etas = 1.0 / (sqrt(M_PI) * sigma_eta);
    const double prefactor_tau = 1.0 / dtau;
    const double n_sigma_skip = 5.;
    const double skip_dis_x = n_sigma_skip * sigma_x;
    const double skip_dis_eta = n_sigma_skip * sigma_eta;

    double tau_dis_max = tau - get_source_tau_max();

    const double prefactor_cov = 1. / pow(2.*M_PI*sigma_x*sigma_x, 3./2.);
    double val_smearing_kernel = 0.;

    if (tau_dis_max < n_sigma_skip * sigma_tau) {
        int n_hadrons_now = list_hadrons_current_tau_.size();

        for (int ipart = 0; ipart < n_hadrons_now; ipart++) {
            // Extract hadron properties
            double quantity_now = 0.0;

            // Evaluate quantity based on type
            switch (quantityType) {
                case BARYON_NUMBER:
                    quantity_now = static_cast<double>(list_hadrons_current_tau_.at(ipart).baryon_number);
                    break;
                case ELECTRIC_CHARGE:
                    quantity_now = static_cast<double>(list_hadrons_current_tau_.at(ipart).electric_charge);
                    break;
                case STRANGENESS:
                    quantity_now = static_cast<double>(list_hadrons_current_tau_.at(ipart).strangeness);
                    break;
            }
            if (quantity_now == 0) {
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

            // kinematic quantities
            const double mass = list_hadrons_current_tau_.at(ipart).mass;
            const double E = list_hadrons_current_tau_.at(ipart).E;
            const double px = list_hadrons_current_tau_.at(ipart).px;
            const double py = list_hadrons_current_tau_.at(ipart).py;
            const double pz = list_hadrons_current_tau_.at(ipart).pz;
            const double mT = sqrt(mass*mass + px*px + py*py);

            if (covariant_smearing_kernel_) {
                // computations for covariant smearing kernel
                const double vx = px / mass;
                const double vy = py / mass;
                const double vz = pz / mass;
                const double gamma = 1. / sqrt(1. - vx*vx - vy*vy - vz*vz);
                const double rapidity = 0.5 * std::log((E + pz) / (E - pz));
                const double peta = mT * sinh(rapidity - list_hadrons_current_tau_.at(ipart).eta_s);
                const double veta = peta / mass;

                val_smearing_kernel = exp_tau * gamma * covariant_smearing_kernel(x_dis,y_dis,eta_s_dis,vx,vy,veta,sigma_x,tau);
            } else {
                const double exp_xperp =
                    exp(-(x_dis * x_dis + y_dis * y_dis) / (sigma_x * sigma_x));
                const double exp_eta_s =
                    exp(-eta_s_dis * eta_s_dis / (sigma_eta * sigma_eta));

                val_smearing_kernel = exp_tau * exp_xperp * exp_eta_s;
            }

            double p_perp_sq = px*px + py*py;
            double u_perp = sqrt(p_perp_sq) / mass;
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
            result += p_dot_u * val_smearing_kernel * (int)quantity_now;
        }
        if (covariant_smearing_kernel_) {
            result *= prefactor_tau * prefactor_cov * weight_event_;
        } else {
            result *= prefactor_tau * prefactor_prep * prefactor_etas * weight_event_;
        }
    }
    return result;
}

double HydroSourceSMASH::get_hydro_rhob_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu) const {
    QuantityType type = BARYON_NUMBER;
    return calculate_source(tau, x, y, eta_s, u_mu, type);
}

double HydroSourceSMASH::get_hydro_rhoq_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu) const {
    QuantityType type = ELECTRIC_CHARGE;
    return calculate_source(tau, x, y, eta_s, u_mu, type);
}

double HydroSourceSMASH::get_hydro_rhos_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu) const {
    QuantityType type = STRANGENESS;
    return calculate_source(tau, x, y, eta_s, u_mu, type);
}
