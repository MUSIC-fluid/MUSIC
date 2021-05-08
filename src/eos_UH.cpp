// Copyright 2020 @ Chun Shen

#include "eos_UH.h"
#include "util.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>

using std::stringstream;
using std::string;

EOS_UH::EOS_UH() {
    set_EOS_id(19);
    set_number_of_tables(0);
    set_eps_max(1e5);
}


EOS_UH::~EOS_UH() {
    int ntables = get_number_of_tables();
    for (int itable = 0; itable < ntables; itable++) {
        Util::mtx_free(mu_B_tb[itable],
                       nb_length[itable], e_length[itable]);
    }
    delete[] mu_B_tb;
}


void EOS_UH::initialize_eos() {
    // read the lattice EOS pressure, temperature, and 
    music_message.info("Using lattice EOS at finite muB from the UH Collaboration");
    music_message.info("reading EOS UH ...");

    stringstream slocalpath;
    slocalpath << "./EOS/UH/";

    string path = slocalpath.str();
    music_message << "from path " << path;
    music_message.flush("info");

    const int ntables = 6;
    set_number_of_tables(ntables);
    resize_table_info_arrays();

    string eos_file_string_array[6] = {"0", "1", "2", "3", "4", "5"};
    pressure_tb    = new double** [ntables];
    temperature_tb = new double** [ntables];
    mu_B_tb        = new double** [ntables];

    for (int itable = 0; itable < ntables; itable++) {
        std::ifstream eos_p(path + "UH_eos_p_"
                            + eos_file_string_array[itable] + ".dat");
        std::ifstream eos_T(path + "UH_eos_T_"
                            + eos_file_string_array[itable] + ".dat");
        std::ifstream eos_mub(path + "UH_eos_muB_"
                            + eos_file_string_array[itable] + ".dat");

        if (!eos_p) {
            music_message << "Can not found the EoS file! filename: "
                          << path + "UH_eos_p_" << eos_file_string_array[itable] + ".dat";
            music_message.flush("error");
            exit(1);
        }

        string dummy;
        std::getline(eos_p, dummy);
        stringstream ss;
        ss << dummy;
        // read the first line with general info:
        ss >> dummy >> dummy >> e_bounds[itable] >> dummy >> e_spacing[itable]
           >> dummy >> e_length[itable] >> dummy >> nb_bounds[itable]
           >> dummy >> nb_spacing[itable] >> dummy >> nb_length[itable];

        e_bounds[itable]  /= Util::hbarc;   // 1/fm^4
        e_spacing[itable] /= Util::hbarc;   // 1/fm^4

        // allocate memory for EOS arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
        mu_B_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                           e_length[itable]);

        // read pressure, temperature and chemical potential values
        for (int j = 0; j < e_length[itable]; j++) {
            for (int i = 0; i < nb_length[itable]; i++) {
                eos_p >> pressure_tb[itable][i][j];
                eos_T >> temperature_tb[itable][i][j];
                eos_mub >> mu_B_tb[itable][i][j];

                pressure_tb[itable][i][j]    /= Util::hbarc;    // 1/fm^4
                temperature_tb[itable][i][j] /= Util::hbarc;    // 1/fm
                mu_B_tb[itable][i][j]        /= Util::hbarc;    // 1/fm
            }
        }
    }

    double eps_max_in = e_bounds[5] + e_spacing[5]*e_length[5];
    set_eps_max(eps_max_in);

    music_message.info("Done reading EOS.");
}


double EOS_UH::p_e_func(double e, double rhob) const {
    return(get_dpOverde3(e, rhob));
}


double EOS_UH::p_rho_func(double e, double rhob) const {
    return(get_dpOverdrhob2(e, rhob));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_UH::get_temperature(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double T = interpolate2D(e, std::abs(rhob), table_idx,
                             temperature_tb);  // 1/fm
    return(T);
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_UH::get_pressure(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double f = interpolate2D(e, std::abs(rhob), table_idx, pressure_tb);
    return(f);
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_UH::get_muB(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double sign = rhob/(std::abs(rhob) + Util::small_eps);
    double mu = sign*interpolate2D(e, std::abs(rhob), table_idx,
                                   mu_B_tb);  // 1/fm
    return(mu);
}


double EOS_UH::get_s2e(double s, double rhob) const {
    double e = get_s2e_finite_rhob(s, rhob);
    return(e);
}
