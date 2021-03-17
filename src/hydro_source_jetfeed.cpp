// Copyright 2019 Chun Shen

#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <memory>

#include "hydro_source_jetfeed.h"
#include "util.h"

using std::string;

HydroSourceJetFeed::HydroSourceJetFeed(const InitData &DATA_in) :
    DATA(DATA_in) {

    list_src_.clear();
    if (DATA.turn_on_source == 1) {
        import_from_file();
    }
}


HydroSourceJetFeed::~HydroSourceJetFeed() {
    list_src_.clear();
}


void HydroSourceJetFeed::import_from_file() {
    string JetFeed_filename = DATA.listName_source;
    music_message << "hydro_source: "
                  << "read in source list from " << JetFeed_filename;
    music_message.flush("info");

    FILE *fin;
    fin = fopen(JetFeed_filename.c_str(), "r");

    while (feof(fin) == 0) {
        SourceProfile new_src;
        new_src.import_from_file(fin);
        list_src_.push_back(new_src);
    }

    fclose(fin);
}


void HydroSourceJetFeed::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s, 
    const FlowVec &u_mu, EnergyFlowVec &j_mu) const {

    for (int imu = 0; imu < 4; imu++) {
        j_mu[imu] = 0.;
    }

    const int nsrc = list_src_.size();
    const double dtau = DATA.delta_tau;
    for (int isrc = 0; isrc < nsrc; isrc++) {
        double tau_src = list_src_.at(isrc).get_tau_src();

        if (tau_src < tau || tau_src >= tau + dtau) {
            continue;
        }

        std::vector<double> jvec_now;
        list_src_.at(isrc).get_source_jvec(eta_s, x, y, jvec_now);

        for (int imu = 0; imu < 4; imu++) {
            j_mu[imu] += tau_src * jvec_now.at(imu) / (tau * dtau);
        }
    }
}


double HydroSourceJetFeed::get_hydro_rhob_source(const double tau, const double x,
                                                 const double y, const double eta_s,
                                                 const FlowVec &u_mu) const {

    double res = 0.;
    if (DATA.turn_on_rhob == 0) {
        return 0.;
    }

    const int nsrc = list_src_.size();
    const double dtau = DATA.delta_tau;
    for (int isrc = 0; isrc < nsrc; isrc++) {
        double tau_src = list_src_.at(isrc).get_tau_src();

        if (tau_src < tau || tau_src >= tau + dtau) {
            continue;
        }

        std::vector<double> jvec_now;
        list_src_.at(isrc).get_source_jvec(eta_s, x, y, jvec_now);

        res += tau_src * jvec_now.at(4) / (tau * dtau);
    }

    return res;
}