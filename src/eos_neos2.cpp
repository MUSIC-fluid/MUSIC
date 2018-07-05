// Copyright 2018 @ Chun Shen

#include "eos_neos2.h"
#include "util.h"

#include <sstream>
#include <fstream>

using std::stringstream;
using std::string;

EOS_neos2::EOS_neos2() {
    set_EOS_id(10);
    set_number_of_tables(0);
    set_eps_max(1e5);
}


void EOS_neos2::initialize_eos() {
    // read the lattice EOS pressure, temperature, and 
    music_message.info("Using lattice EOS at finite muB from A. Monnai");
    music_message.info("reading EOS neos2 ...");
    
    stringstream slocalpath;
    slocalpath << "./EOS/neos_2/";

    string path = slocalpath.str();
    music_message << "from path " << path;
    music_message.flush("info");
    
    const int ntables = 7;
    set_number_of_tables(ntables);
    resize_table_info_arrays();

    string eos_file_string_array[7] = {"0a", "0b", "0c", "1a", "2", "3", "4"};
    pressure_tb    = new double** [ntables];
    temperature_tb = new double** [ntables];
    mu_B_tb        = new double** [ntables];

    for (int itable = 0; itable < ntables; itable++) {
        std::ifstream eos_p(path + "neos" + eos_file_string_array[itable]
                            + "_p.dat");
        std::ifstream eos_T(path + "neos" + eos_file_string_array[itable]
                            + "_t.dat");
        std::ifstream eos_mub(path + "neos" + eos_file_string_array[itable]
                              + "_mb.dat");
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
            }
        }
    }
    
    double eps_max_in = (e_bounds[6] + e_spacing[6]*e_length[6])/hbarc;
    set_eps_max(eps_max_in);

    music_message.info("Done reading EOS.");
}


double EOS_neos2::get_cs2(double e, double rhob) const {
    double f = calculate_velocity_of_sound_sq(e, rhob);
    return(f);
}
    

double EOS_neos2::p_e_func(double e, double rhob) const {
    return(get_dpOverde3(e, rhob));
}


double EOS_neos2::p_rho_func(double e, double rhob) const {
    return(get_dpOverdrhob2(e, rhob));
}


double EOS_neos2::get_dpOverdrhob2(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double deltaRhob = nb_spacing[table_idx];
    //double rhob_max = nb_bounds[table_idx] + nb_length[table_idx]*deltaRhob;
    
    double rhobLeft  = rhob - deltaRhob*0.5;
    double rhobRight = rhob + deltaRhob*0.5;

    double pL = get_pressure(e, rhobLeft);      // 1/fm^4
    double pR = get_pressure(e, rhobRight);     // 1/fm^4
      
    double dpdrho = (pR - pL)/(rhobRight - rhobLeft);  // 1/fm
    return (dpdrho);   // in 1/fm
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_neos2::get_temperature(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double T = interpolate2D(e, std::abs(rhob), table_idx,
                             temperature_tb)/hbarc;  // 1/fm
    return(std::max(1e-15, T));
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_neos2::get_pressure(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double f = interpolate2D(e, std::abs(rhob), table_idx, pressure_tb)/hbarc;
    return(std::max(1e-15, f));
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_neos2::get_mu(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double sign = rhob/(std::abs(rhob) + 1e-15);
    double mu = sign*interpolate2D(e, std::abs(rhob), table_idx,
                                   mu_B_tb)/hbarc;  // 1/fm
    return(mu);
}


double EOS_neos2::get_s2e(double s, double rhob) const {
    double e = get_s2e_finite_rhob(s, rhob);
    return(e);
}
