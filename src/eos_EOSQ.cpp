// Copyright 2018 @ Chun Shen

#include "eos_EOSQ.h"
#include "util.h"

#include <sstream>
#include <fstream>
#include <cmath>

using std::stringstream;
using std::string;

EOS_eosQ::EOS_eosQ() {
    set_EOS_id(1);
    set_number_of_tables(0);
    set_eps_max(1e5);
    set_flag_muB(true);
    set_flag_muS(false);
    set_flag_muC(false);
}


EOS_eosQ::~EOS_eosQ() {
    int ntables = get_number_of_tables();
    for (int itable = 0; itable < ntables; itable++) {
        Util::mtx_free(mu_B_tb[itable],
                       nb_length[itable], e_length[itable]);
    }
}


void EOS_eosQ::initialize_eos() {
    // read the lattice EOS pressure, temperature, and 
    music_message.info("Using EOS-Q from AZHYDRO");
    music_message.info("reading EOS eosQ ...");
    
    auto envPath = get_hydro_env_path();
    music_message << "from path " << envPath.c_str() << "/EOS";
    music_message.flush("info");
    
    const int ntables = 2;
    set_number_of_tables(ntables);
    resize_table_info_arrays();

    string eos_file_string_array[2] = {"1", "2"};
    pressure_tb    = new double** [ntables];
    temperature_tb = new double** [ntables];
    mu_B_tb        = new double** [ntables];

    for (int itable = 0; itable < ntables; itable++) {
        std::ifstream eos_p(envPath + "/EOS/EOS-Q/aa"
                            + eos_file_string_array[itable] + "_p.dat");
        std::ifstream eos_T(envPath + "/EOS/EOS-Q/aa"
                            + eos_file_string_array[itable] + "_t.dat");
        std::ifstream eos_mub(envPath + "/EOS/EOS-Q/aa"
                              + eos_file_string_array[itable] + "_mb.dat");
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

        e_bounds[itable]  /= Util::hbarc;   // 1/fm^4
        e_spacing[itable] /= Util::hbarc;   // 1/fm^4

        // skip the header in T and mu_B files
        string dummy;
        std::getline(eos_T, dummy);
        std::getline(eos_T, dummy);
        std::getline(eos_mub, dummy);
        std::getline(eos_mub, dummy);

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
    
    //double eps_max_in = (e_bounds[1] + e_spacing[1]*e_length[1])/hbarc;
    double eps_max_in = e_bounds[1] + e_spacing[1]*e_length[1];
    set_eps_max(eps_max_in);

    music_message.info("Done reading EOS.");
}


double EOS_eosQ::p_e_func(double e, double rhob) const {
    return(get_dpOverde3(e, rhob));
}


double EOS_eosQ::p_rho_func(double e, double rhob) const {
    return(get_dpOverdrhob2(e, rhob));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_eosQ::get_temperature(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double T = interpolate2D(e, std::abs(rhob), table_idx,
                             temperature_tb);  // 1/fm
    return(std::max(Util::small_eps, T));
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_eosQ::get_pressure(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double f = interpolate2D(e, std::abs(rhob), table_idx, pressure_tb);
    return(std::max(Util::small_eps, f));
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_eosQ::get_mu(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double sign = rhob/(std::abs(rhob) + Util::small_eps);
    double mu = sign*interpolate2D(e, std::abs(rhob), table_idx,
                                   mu_B_tb);  // 1/fm
    return(mu);
}


double EOS_eosQ::get_s2e(double s, double rhob) const {
    double e = get_s2e_finite_rhob(s, rhob);
    return(e);
}
