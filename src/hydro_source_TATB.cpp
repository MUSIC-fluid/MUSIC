// Copyright 2019 Chun Shen

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <memory>

#include "hydro_source_TATB.h"
#include "util.h"

using std::string;


HydroSourceTATB::HydroSourceTATB(const InitData &DATA_in) :
    DATA_(DATA_in) {
    set_source_tau_min(DATA_.tau0);
    set_source_tau_max(DATA_.tau0);
    read_in_TATB();
}


HydroSourceTATB::~HydroSourceTATB() {
    profile_TA.clear();
    profile_TB.clear();
}


//! This function reads in the spatal information of the nuclear thickness
//! functions
void HydroSourceTATB::read_in_TATB() {
    music_message << "read in TA and TB from " << DATA_.initName_TA
                  << " and " << DATA_.initName_TB;
    music_message.flush("info");

    string text_string;

    std::ifstream TAfile(DATA_.initName_TA.c_str());
    if (!TAfile) {
        music_message << "hydro_source::read_in_TATB: "
                      << "can not open TA file: " << DATA_.initName_TA;
        music_message.flush("error");
        exit(1);
    }
    std::ifstream TBfile(DATA_.initName_TB.c_str());
    if (!TBfile) {
        music_message << "hydro_source::read_in_TATB: "
                      << "can not open TB file: " << DATA_.initName_TB;
        music_message.flush("error");
        exit(1);
    }

    const int nx = DATA_.nx;
    const int ny = DATA_.ny;
    double N_B = 0.0;
    for (int i = 0; i < nx; i++) {
        std::vector<double> TA_temp;
        std::vector<double> TB_temp;
        for (int j = 0; j < ny; j++) {
            double TA, TB;
            TAfile >> TA;
            TBfile >> TB;
            TA_temp.push_back(TA);
            TB_temp.push_back(TB);
            N_B += TA + TB;
        }
        profile_TA.push_back(TA_temp);
        profile_TB.push_back(TB_temp);
    }
    TAfile.close();
    TBfile.close();
    N_B *= DATA_.delta_x*DATA_.delta_y;
    double total_energy = DATA_.ecm/2.*N_B;
    music_message << "sqrt{s} = " << DATA_.ecm << " GeV, "
                  << "beam rapidity = " << DATA_.beam_rapidity << ", "
                  << "total energy = " << total_energy << " GeV, "
                  << "N_B = " << N_B;
    music_message.flush("info");

    music_message << "HydroSourceTATB: tau_min = " << get_source_tau_min()
                  << " fm/c.";
    music_message.flush("info");
    music_message << "HydroSourceTATB: tau_max = " << get_source_tau_max()
                  << " fm/c.";
    music_message.flush("info");
}


void HydroSourceTATB::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s,
    const FlowVec &u_mu, EnergyFlowVec &j_mu) const {
    j_mu = {0};
}


double HydroSourceTATB::get_hydro_rhob_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu) const {
    double res = 0.;
    return(res);
}


double HydroSourceTATB::eta_rhob_left_factor(const double eta) const {
    double eta_0       = -std::abs(DATA_.eta_rhob_0);
    double tau0        = DATA_.tau0;
    double delta_eta_1 = DATA_.eta_rhob_width_1;
    double delta_eta_2 = DATA_.eta_rhob_width_2;
    double norm        = 2./(sqrt(M_PI)*tau0*(delta_eta_1 + delta_eta_2));
    double exp_arg     = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_1;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_2;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}


double HydroSourceTATB::eta_rhob_right_factor(const double eta) const {
    double eta_0       = std::abs(DATA_.eta_rhob_0);
    double tau0        = DATA_.tau0;
    double delta_eta_1 = DATA_.eta_rhob_width_1;
    double delta_eta_2 = DATA_.eta_rhob_width_2;
    double norm        = 2./(sqrt(M_PI)*tau0*(delta_eta_1 + delta_eta_2));
    double exp_arg     = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_2;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_1;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}

