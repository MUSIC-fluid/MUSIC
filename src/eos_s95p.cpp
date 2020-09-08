// Copyright 2018 @ Chun Shen

#include "eos_s95p.h"
#include "util.h"

#include <sstream>
#include <fstream>

using std::stringstream;
using std::string;

EOS_s95p::EOS_s95p(const int eos_id_in) : eos_id(eos_id_in) {
    set_EOS_id(eos_id);
    set_number_of_tables(0);
    set_eps_max(1e5);
    set_flag_muB(false);
    set_flag_muS(false);
    set_flag_muC(false);
}


void EOS_s95p::initialize_eos() {
    // read the lattice EOS pressure, temperature, and 
    music_message.info("reading EOS s95p ...");
    
    auto envPath = get_hydro_env_path();
    stringstream spath;
    spath << envPath;
    if (eos_id == 2) {
        music_message.info("Using lattice EOS from Huovinen/Petreczky");
        spath << "/EOS/s95p-v1/";
    } else if (eos_id == 3) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) "
                      << "chem. f.o. at 150 MeV";
        music_message.flush("info");
        spath << "/EOS/s95p-PCE-v1/";
    } else if (eos_id == 4) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) "
                      << "chem. f.o. at 155 MeV";
        spath << "/EOS/s95p-PCE155/";
    } else if (eos_id == 5) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) "
                      << "chem. f.o. at 160 MeV";
        spath << "/EOS/s95p-PCE160/";
    } else if (eos_id == 6) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) chem. f.o. "
                       << "at 165 MeV";
        spath << "/EOS/s95p-PCE165-v0/";
    } else if (eos_id == 7) {
        music_message.info(
            "Using lattice EOS from Huovinen/Petreczky s95p-v1.2 (for UrQMD)");
        spath << "/EOS/s95p-v1.2/";
    }
    
    music_message << "from path " << spath.str();
    music_message.flush("info");

    if (eos_id == 2) {
        spath << "s95p-v1_";
    } else if (eos_id == 7) {
        spath << "s95p-v1.2_";
    }
    
    const int ntables = 7;
    set_number_of_tables(ntables);
    resize_table_info_arrays();
    
    string eos_file_string_array[7] = {"1", "2", "3", "4", "5", "6", "7"};
    pressure_tb    = new double** [ntables];
    temperature_tb = new double** [ntables];
    for (int itable = 0; itable < ntables; itable++) {
        std::ifstream eos_d(spath.str() + "dens"
                            + eos_file_string_array[itable] + ".dat");
        std::ifstream eos_T(spath.str() + "par"
                            + eos_file_string_array[itable] + ".dat");
        // read the first two lines with general info:
        // lowest value of epsilon
        // deltaEpsilon, number of epsilon steps (i.e. # of lines)
        eos_d >> e_bounds[itable];
        eos_d >> e_spacing[itable] >> e_length[itable];
        
        e_bounds[itable]  /= Util::hbarc;   // 1/fm^4
        e_spacing[itable] /= Util::hbarc;   // 1/fm^4
        
        // skip the header in T file
        string dummy;
        std::getline(eos_T, dummy);
        std::getline(eos_T, dummy);
        
        // no rho_b dependence at the moment
        nb_length[itable] = 1;
        
        // allocate memory for pressure arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);

        // read pressure, temperature and chemical potential values
        // files have it backwards, so I start with maximum j and count down
        int i = 0;
        double d_dummy;
        for (int j = e_length[itable] - 1; j >= 0; j--) {
            eos_d >> d_dummy;
            eos_d >> pressure_tb[itable][i][j];
            eos_d >> d_dummy >> dummy >> dummy;;
            eos_T >> temperature_tb[itable][i][j] >> dummy >> dummy;
                
            pressure_tb[itable][i][j]    /= Util::hbarc;    // 1/fm^4
            temperature_tb[itable][i][j] /= Util::hbarc;    // 1/fm
        }
    }

    //double eps_max_in = (e_bounds[6] + e_spacing[6]*e_length[6])/hbarc;
    double eps_max_in = e_bounds[6] + e_spacing[6]*e_length[6];
    set_eps_max(eps_max_in);

    music_message.info("Done reading EOS.");
}


double EOS_s95p::p_e_func(double e, double rhob) const {
    return(get_dpOverde3(e, rhob));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_s95p::get_temperature(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double T = interpolate1D(e, table_idx, temperature_tb);  // 1/fm
    return(std::max(Util::small_eps, T));
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_s95p::get_pressure(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double f = interpolate1D(e, table_idx, pressure_tb);  // 1/fm^4
    return(std::max(Util::small_eps, f));
}


double EOS_s95p::get_s2e(double s, double rhob) const {
    double e = get_s2e_finite_rhob(s, 0.0);
    return(e);
}

double EOS_s95p::get_T2e(double T, double rhob) const {
    double e = get_T2e_finite_rhob(T, 0.0);
    return(e);
}
