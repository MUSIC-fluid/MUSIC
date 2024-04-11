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
    set_sigma_eta(0.5);
    int i_event = DATA.event_id_SMASH_output;
    if (i_event < 1) {
        i_event = 1;
    }
    average_events_ = DATA.average_SMASH_events;
    int extended_output = DATA.extended_SMASH_output;
    int reject_spectators = DATA.reject_SMASH_spectators;
    read_in_SMASH_hadrons(i_event, extended_output, reject_spectators);
    set_covariant_smearing_kernel(DATA.use_cov_smearing);
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
                  << "read in SMASH hadron list from " << SMASH_filename;
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
            int pdgid, ID, ncoll;
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
            if(ncoll == 0){
                n_spectators_ev +=1;
                continue;
            }
            // cuts in rapidity.
            //double rapidity = 0.5 * std::log((p0 + pz) / (p0 - pz));
            //double eta_s =  0.5 * std::log((t + z) / (t - z));
            //if(eta_s > 1.0) continue;
            //if(eta_s < -1.0) continue;

            // Put info in a new hadron
            hadron new_hadron;
            new_hadron.pdgid = pdgid;
            new_hadron.ID = ID;
            new_hadron.t = t;
            new_hadron.tau = std::sqrt(t * t - z * z);
            new_hadron.x = x;
            new_hadron.y = y;
            new_hadron.z = z;
            new_hadron.eta_s = 0.5 * std::log((t + z) / (t - z));
            new_hadron.rapidity = 0.5 * std::log((p0 + pz) / (p0 - pz));
            new_hadron.E = p0;
            new_hadron.px = px;
            new_hadron.py = py;
            new_hadron.pz = pz;
            new_hadron.mass = mass;
            new_hadron.baryon_number = B;
            new_hadron.electric_charge = Q;
            new_hadron.strangeness = S;
            new_hadron.norm = 1.;
            // Put hadron in hadron list.
            // Note, if average_events_ = 1
            // list_hadrons contains hadrons
            // for all event in IC file.
            if (new_hadron.pdgid != 22) {
                list_hadrons_.push_back(new_hadron);
            }

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
                music_message << "HydroSourceSMASH::     pdgid = "
                              << new_hadron.pdgid;
                music_message.flush("info");
                music_message << "HydroSourceSMASH::     (x, y) = ("
                              << new_hadron.x << ", "
                              << new_hadron.y << ") fm.";
                music_message.flush("info");
                music_message << "HydroSourceSMASH::     eta_s = "
                              << new_hadron.eta_s;
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
    music_message << "HydroSourceSMASH:: read in "
                  << list_hadrons_.size() << "/"
                  << n_hadrons << " hadrons and "
                  << n_spectators << " spectators from "
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
    music_message << "HydroSourceSMASH:: electric_charge_total = "
                  << charge_total_;
    music_message.flush("info");
    music_message << "HydroSourceSMASH:: strangeness_total = "
                  << strangeness_total_;
    music_message.flush("info");
}


void HydroSourceSMASH::prepare_list_for_current_tau_frame(
                                                const double tau_local) {
    const double dtau = DATA.delta_tau;
    list_hadrons_current_tau_.clear();
    for (auto hadron_i : list_hadrons_) {
        if (hadron_i.tau >= (tau_local - 0.5*dtau)
            && hadron_i.tau < (tau_local + 0.5*dtau)) {
            list_hadrons_current_tau_.push_back(hadron_i);
        }
    }
    music_message << "hydro_source: tau = " << tau_local
                  << " number of sources: "
                  << list_hadrons_current_tau_.size();
    music_message.flush("info");

    // compute the norm for each source particle
    if (covariant_smearing_kernel_) {
        for (auto hadron_i: list_hadrons_current_tau_) {
            compute_covariant_norm(tau_local, &hadron_i);
            std::cout << "Norm = " << hadron_i.norm << std::endl;
        }
    } else {
        for (auto hadron_i: list_hadrons_current_tau_) {
            compute_norm(tau_local, &hadron_i);
            std::cout << "Norm = " << hadron_i.norm << std::endl;
        }
    }
}

void HydroSourceSMASH::compute_covariant_norm(const double tau,
                                                hadron* hadron_i) const {
    
    std::cout << "tau_MUSIC = " << tau << ", tau_hadron = " << hadron_i->tau << std::endl;
    double norm = 0.0;
    const double n_sigma_skip = 5.;
    const double sigma_x = get_sigma_x();
    const int nx = DATA.nx;
    const int ny = DATA.ny;
    const int neta = DATA.neta;
    const double dx = DATA.delta_x;
    const double dy = DATA.delta_y;
    const double deta = DATA.delta_eta;
    const double skip_dis_x = n_sigma_skip * sigma_x;
    const double skip_dis_eta = skip_dis_x / tau;
    const double cov_prefac = tau / (pow(M_PI,1.5)*sigma_x*sigma_x*sigma_x);

    #pragma omp parallel for collapse(3) schedule(guided) reduction(+: norm)
    for (int ieta = 0; ieta < neta; ieta++) {
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                const double eta_s = - DATA.eta_size/2. + ieta*deta;
                const double x = - DATA.x_size/2. + ix*dx;
                const double y = - DATA.y_size/2. + iy*dy;

                double x_dis = x - hadron_i->x;
                if (std::abs(x_dis) > skip_dis_x) continue;

                double y_dis = y - hadron_i->y;
                if (std::abs(y_dis) > skip_dis_x) continue;

                double eta_s_dis = eta_s - hadron_i->eta_s;
                if (std::abs(eta_s_dis) > skip_dis_eta) continue;

                double mass = hadron_i->mass;
                double px = hadron_i->px;
                double py = hadron_i->py;

                double ux = px / mass;
                double uy = py / mass;

                double mT = sqrt(mass * mass + px * px + py * py);
                double rapidity = hadron_i->rapidity;

                double peta = mT * sinh(rapidity - hadron_i->eta_s);
                double ueta = peta / mass;

                double gamma = mT * cosh(rapidity - hadron_i->eta_s) / mass;
                norm += (
                    gamma*covariant_smearing_kernel(x_dis, y_dis, eta_s_dis,
                                                    ux, uy, ueta,
                                                    sigma_x, tau)
                );
            }
        }
    }
    hadron_i->norm = (norm*cov_prefac
                      *(DATA.delta_x*DATA.delta_y*DATA.delta_eta)
                      );
}


void HydroSourceSMASH::compute_norm(const double tau, hadron* hadron_i) const {
    double norm = 0.0;
    const double n_sigma_skip = 5.;
    const double sigma_x = get_sigma_x();
    const double sigma_eta = get_sigma_eta();
    const int nx = DATA.nx;
    const int ny = DATA.ny;
    const int neta = DATA.neta;
    const double dx = DATA.delta_x;
    const double dy = DATA.delta_y;
    const double deta = DATA.delta_eta;
    const double skip_dis_x = n_sigma_skip * sigma_x;
    const double skip_dis_eta = n_sigma_skip * sigma_eta;
    const double prefactor_etas = 1. / (sqrt(M_PI) * sigma_eta);
    const double prefactor_prep = 1. / (M_PI * sigma_x * sigma_x);

    #pragma omp parallel for collapse(3) schedule(guided) reduction(+: norm)
    for (int ieta = 0; ieta < neta; ieta++) {
        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                const double eta_s = - DATA.eta_size/2. + ieta*deta;
                const double x = - DATA.x_size/2. + ix*dx;
                const double y = - DATA.y_size/2. + iy*dy;

                double x_dis = x - hadron_i->x;
                if (std::abs(x_dis) > skip_dis_x) continue;
                double y_dis = y - hadron_i->y;
                if (std::abs(y_dis) > skip_dis_x) continue;
                double eta_s_dis = eta_s - hadron_i->eta_s;
                if (std::abs(eta_s_dis) > skip_dis_eta) continue;

                const double exp_xperp = (
                    exp( - (x_dis * x_dis + y_dis * y_dis)/(sigma_x*sigma_x)));
                const double exp_eta_s = (
                    exp( - eta_s_dis*eta_s_dis/(sigma_eta*sigma_eta)));
                norm += exp_xperp*exp_eta_s;
            }
        }
    }
    hadron_i->norm = (norm*prefactor_etas*prefactor_prep
                      *(DATA.delta_x*DATA.delta_y*DATA.delta_eta));
}


double HydroSourceSMASH::covariant_smearing_kernel(
    const double x_diff, const double y_diff, const double eta_diff,
    const double ux, const double uy, const double ueta,
    const double sigma, const double tau) const {
    // Calculate the squared distance and scalar product r*u
    const double r_squared = (x_diff * x_diff + y_diff * y_diff
                              + eta_diff * eta_diff * tau * tau);
    const double ru_scalar = (x_diff * ux + y_diff * uy
                              + eta_diff * ueta * tau);
    // Compute and return the smoothing kernel
    return exp((-r_squared - ru_scalar * ru_scalar) / (sigma * sigma));
}


void HydroSourceSMASH::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s,
    const FlowVec &u_mu, EnergyFlowVec &j_mu) const {

    j_mu = {0};
    if (list_hadrons_current_tau_.size() == 0) return;

    const double sigma_x = get_sigma_x();
    const double sigma_eta = get_sigma_eta();
    const double dtau = DATA.delta_tau;
    const double prefactor_tau = 1. / dtau;
    const double n_sigma_skip = 5.;
    const double skip_dis_x = n_sigma_skip * sigma_x;

    double skip_dis_eta = n_sigma_skip * sigma_eta;
    if (covariant_smearing_kernel_) {
        skip_dis_eta = skip_dis_x / tau;
    }

    double val_smearing_kernel = 0.;
    const double cov_prefac = tau / (pow(M_PI,1.5)*sigma_x*sigma_x*sigma_x);
    const double prefactor_etas = 1. / (sqrt(M_PI) * sigma_eta);
    const double prefactor_prep = 1. / (M_PI * sigma_x * sigma_x);
    const double exp_tau = 1./tau;

    for (auto hadron_i: list_hadrons_current_tau_) {
        const double x_dis = x - hadron_i.x;
        if (std::abs(x_dis) > skip_dis_x) continue;
        const double y_dis = y - hadron_i.y;
        if (std::abs(y_dis) > skip_dis_x) continue;
        const double eta_s_dis = eta_s - hadron_i.eta_s;
        if (std::abs(eta_s_dis) > skip_dis_eta) continue;

        const double mass = hadron_i.mass;
        const double px = hadron_i.px;
        const double py = hadron_i.py;
        const double rapidity = hadron_i.rapidity;
        const double mT = sqrt(mass*mass + px*px + py*py);
        if (covariant_smearing_kernel_) {
            const double ux = px / mass;
            const double uy = py / mass;

            const double peta = mT * sinh(rapidity - hadron_i.eta_s);
            const double ueta = peta / mass;
            const double gamma = mT * cosh(rapidity - hadron_i.eta_s) / mass;

            val_smearing_kernel = (
                gamma*covariant_smearing_kernel(x_dis, y_dis, eta_s_dis,
                                                ux, uy, ueta,
                                                sigma_x, tau)/hadron_i.norm
            );
        } else {
            const double exp_xperp = (
                exp( - (x_dis*x_dis + y_dis*y_dis)/(sigma_x*sigma_x))
            );
            const double exp_eta_s = (
                exp( - eta_s_dis*eta_s_dis/(sigma_eta*sigma_eta))
            );

            val_smearing_kernel = exp_xperp*exp_eta_s/hadron_i.norm;
        }

        if (std::isinf(val_smearing_kernel)
            || std::isnan(val_smearing_kernel)) {
            std::cout << "px = " << hadron_i.px
                      << ", py = " << hadron_i.py
                      << ", pz = " << hadron_i.pz
                      << ", E = " << hadron_i.E
                      << std::endl;
            std::cout << "norm = " << hadron_i.norm << std::endl;
            exit(1);
        }

        j_mu[0] += mT*cosh(rapidity - eta_s)*val_smearing_kernel;
        j_mu[1] += px*val_smearing_kernel;
        j_mu[2] += py*val_smearing_kernel;
        j_mu[3] += mT*sinh(rapidity - eta_s)*val_smearing_kernel;
    }

    const double norm = DATA.sFactor / Util::hbarc; // 1/fm^4
    double prefactor = norm;
    if (covariant_smearing_kernel_) {
        prefactor *= exp_tau * prefactor_tau * cov_prefac;
    } else {
        prefactor *= exp_tau * prefactor_tau * prefactor_etas * prefactor_prep;
    }
    j_mu[0] *= prefactor * weight_event_;
    j_mu[1] *= prefactor * weight_event_;
    j_mu[2] *= prefactor * weight_event_;
    j_mu[3] *= prefactor * weight_event_;
}


double HydroSourceSMASH::calculate_source(
    const double tau, const double x, const double y, const double eta_s,
    const FlowVec &u_mu, const QuantityType quantityType) const {

    double result = 0.0;

    if (list_hadrons_current_tau_.empty()) return result;

    const double sigma_x = get_sigma_x();
    const double sigma_eta = get_sigma_eta();

    const double dtau = DATA.delta_tau;
    const double prefactor_tau = 1.0 / dtau;
    const double n_sigma_skip = 5.;
    const double skip_dis_x = n_sigma_skip * sigma_x;
    double skip_dis_eta = n_sigma_skip * sigma_eta;
    if (covariant_smearing_kernel_) {
        skip_dis_eta = skip_dis_x / tau;
    }

    const double cov_prefac = tau / (pow(M_PI,1.5)*sigma_x*sigma_x*sigma_x);
    const double prefactor_etas = 1./(sqrt(M_PI) * sigma_eta);
    const double prefactor_prep = 1./(M_PI * sigma_x * sigma_x);
    const double exp_tau = 1.0/tau;

    double val_smearing_kernel = 0.;
    for (auto hadron_i: list_hadrons_current_tau_) {
        int quantity_now = 0;
        if (quantityType == BARYON_NUMBER) {
            quantity_now = hadron_i.baryon_number;
        } else if (quantityType == ELECTRIC_CHARGE) {
            quantity_now = hadron_i.electric_charge;
        } else if (quantityType == STRANGENESS) {
            quantity_now = hadron_i.strangeness;
        }
        if (quantity_now == 0) continue;

        const double x_dis = x - hadron_i.x;
        if (std::abs(x_dis) > skip_dis_x) continue;
        const double y_dis = y - hadron_i.y;
        if (std::abs(y_dis) > skip_dis_x) continue;
        const double eta_s_dis = eta_s - hadron_i.eta_s;
        if (std::abs(eta_s_dis) > skip_dis_eta) continue;

        if (covariant_smearing_kernel_) {
            const double mass = hadron_i.mass;
            const double px = hadron_i.px;
            const double py = hadron_i.py;
            const double ux = px / mass;
            const double uy = py / mass;

            const double mT = sqrt(mass*mass + px*px + py*py);

            const double rapidity = hadron_i.rapidity;
            const double peta = mT*sinh(rapidity - hadron_i.eta_s);
            const double ueta = peta/mass;

            const double gamma = mT * cosh(rapidity - hadron_i.eta_s) / mass;
            val_smearing_kernel = (
                gamma*covariant_smearing_kernel(x_dis, y_dis, eta_s_dis,
                                                ux, uy, ueta,
                                                sigma_x, tau)/hadron_i.norm
            );
        } else {
            const double exp_xperp = (
                exp( - (x_dis*x_dis + y_dis * y_dis)/(sigma_x*sigma_x)));
            const double exp_eta_s = (
                exp( - eta_s_dis*eta_s_dis/(sigma_eta * sigma_eta)));

            val_smearing_kernel = exp_xperp*exp_eta_s/hadron_i.norm;
        }

        result += val_smearing_kernel * quantity_now;
    }

    if (covariant_smearing_kernel_) {
        result *= exp_tau * prefactor_tau * cov_prefac;
    } else {
        result *= exp_tau * prefactor_tau * prefactor_etas * prefactor_prep;
    }
    return (result*weight_event_);
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
