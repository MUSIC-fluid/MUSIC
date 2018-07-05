// Copyright 2018 @ Chun Shen

#include "eos_neos3.h"
#include "util.h"

#include <sstream>
#include <fstream>

using std::stringstream;
using std::string;

EOS_neos3::EOS_neos3() {
    set_EOS_id(12);
    set_number_of_tables(0);
    set_eps_max(1e5);
}


EOS_neos3::~EOS_neos3() {
    int ntables = get_number_of_tables();
    for (int itable = 0; itable < ntables; itable++) {
        Util::mtx_free(mu_B_tb[itable],
                       nb_length[itable], e_length[itable]);
    }
}


void EOS_neos3::initialize_eos() {
    // read the lattice EOS pressure, temperature, and 
    music_message.info("Using lattice EOS from A. Monnai (up to mu_B^6)");
    music_message.info("reading EOS neos3 ...");
    
    stringstream slocalpath;
    slocalpath << "./EOS/neos_3/";

    string path = slocalpath.str();
    music_message << "from path " << path;
    music_message.flush("info");
    
    const int ntables = 7;
    set_number_of_tables(ntables);
    resize_table_info_arrays();

    string eos_file_string_array[7] = {"1", "2", "3", "4", "5", "6", "7"};
    pressure_tb    = new double** [ntables];
    temperature_tb = new double** [ntables];
    mu_B_tb        = new double** [ntables];

    for (int itable = 0; itable < ntables; itable++) {
        std::ifstream eos_p(path + "neos"
                            + eos_file_string_array[itable] + "_p.dat");
        std::ifstream eos_T(path + "neos"
                            + eos_file_string_array[itable] + "_t.dat");
        std::ifstream eos_muB(path + "neos"
                              + eos_file_string_array[itable] + "_mub.dat");
        // read the first two lines with general info:
        // first value of rhob, first value of epsilon
        // deltaRhob, deltaE, number of rhob points, number of epsilon points
        // the table size is
        // (number of rhob points + 1, number of epsilon points + 1)
        int N_e, N_rhob;
        eos_p >> nb_bounds[itable] >> e_bounds[itable];
        eos_p >> nb_spacing[itable] >> e_spacing[itable]
              >> N_rhob >> N_e;
        nb_length[itable] = N_rhob + 1;
        e_length[itable]  = N_e + 1;

        e_bounds[itable]  /= hbarc;   // 1/fm^4
        e_spacing[itable] /= hbarc;   // 1/fm^4

        // skip the header in T and mu_B files
        string dummy;
        std::getline(eos_T, dummy);
        std::getline(eos_T, dummy);
        std::getline(eos_muB, dummy);
        std::getline(eos_muB, dummy);

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
                eos_muB >> mu_B_tb[itable][i][j];
                
                pressure_tb[itable][i][j]    /= hbarc;    // 1/fm^4
                temperature_tb[itable][i][j] /= hbarc;    // 1/fm
                mu_B_tb[itable][i][j]        /= hbarc;    // 1/fm
            }
        }
    }
    
    //double eps_max_in = (e_bounds[6] + e_spacing[6]*e_length[6])/hbarc;
    double eps_max_in = e_bounds[6] + e_spacing[6]*e_length[6];
    set_eps_max(eps_max_in);

    music_message.info("Done reading EOS.");
}


double EOS_neos3::p_e_func(double e, double rhob) const {
    return(get_dpOverde3(e, rhob));
}


double EOS_neos3::p_rho_func(double e, double rhob) const {
    return(get_dpOverdrhob2(e, rhob));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_neos3::get_temperature(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double T = interpolate2D(e, std::abs(rhob), table_idx,
                             temperature_tb);  // 1/fm
    return(std::max(1e-15, T));
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_neos3::get_pressure(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double f = interpolate2D(e, std::abs(rhob), table_idx, pressure_tb);
    return(std::max(1e-15, f));
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_neos3::get_mu(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double sign = rhob/(std::abs(rhob) + 1e-15);
    double mu = sign*interpolate2D(e, std::abs(rhob), table_idx,
                                   mu_B_tb);  // 1/fm
    return(mu);
}


double EOS_neos3::get_s2e(double s, double rhob) const {
    double e = get_s2e_finite_rhob(s, rhob);
    return(e);
}
