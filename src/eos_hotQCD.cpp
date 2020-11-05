// Copyright 2018 @ Chun Shen

#include "eos_hotQCD.h"
#include "util.h"

#include <sstream>
#include <fstream>
#include <cmath>

using std::stringstream;
using std::string;

EOS_hotQCD::EOS_hotQCD(const int eos_id_in) : eos_id(eos_id_in) {
    set_EOS_id(eos_id);
    set_number_of_tables(0);
    set_eps_max(1e5);
    set_flag_muB(false);
    set_flag_muS(false);
    set_flag_muC(false);
}


void EOS_hotQCD::initialize_eos() {
    // read the lattice EOS pressure, temperature, and 
    music_message.info("reading EOS hotQCD ...");

    auto envPath = get_hydro_env_path();
    stringstream slocalpath;
    slocalpath << envPath << "/EOS/hotQCD";

    string path = slocalpath.str();
    music_message << "from path " << path;
    music_message.flush("info");

    set_number_of_tables(1);
    resize_table_info_arrays();

    int ntables = get_number_of_tables();

    pressure_tb    = new double** [ntables];
    temperature_tb = new double** [ntables];
    for (int itable = 0; itable < ntables; itable++) {
        std::ifstream eos_file;
        if (eos_id == 9) {
            eos_file.open(path + "/hrg_hotqcd_eos_binary.dat",
                          std::ios::binary);
        } else if (eos_id == 91) {
            eos_file.open(path + "/hrg_hotqcd_eos_SMASH_binary.dat",
                          std::ios::binary);
        }

        if (!eos_file) {
            music_message.error("Can not find the EoS file.");
            exit(1);
        }

        e_length[itable]  = 100000;
        nb_length[itable] = 1;
        // allocate memory for pressure arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
        double temp;
        for (int ii = 0; ii < e_length[itable]; ii++) {
            eos_file.read((char*)&temp, sizeof(double));  // e
            temp /= Util::hbarc;      // 1/fm^4
            if (ii == 0) e_bounds[itable] = temp;
            if (ii == 1) e_spacing[itable] = temp - e_bounds[itable];
            if (ii == e_length[itable] - 1) set_eps_max(temp);

            eos_file.read((char*)&temp, sizeof(double));  // P
            pressure_tb[itable][0][ii] = temp/Util::hbarc;      // 1/fm^4

            eos_file.read((char*)&temp, sizeof(double));  // s

            eos_file.read((char*)&temp, sizeof(double));  // T
            temp /= Util::hbarc;   // 1/fm
            // store T^5
            // scaling with T^5 seems to be the best to get a smooth curve
            temperature_tb[itable][0][ii] = pow(temp, 5);
        }
    }
    music_message.info("Done reading EOS.");
}


double EOS_hotQCD::p_e_func(double e, double rhob) const {
    return(get_dpOverde3(e, rhob));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_hotQCD::get_temperature(double e, double rhob) const {
    double T5 = interpolate1D(e, 0, temperature_tb);  // e/T^5
    double T = pow(T5, 0.2);  // 1/fm
    return(std::max(Util::small_eps, T));
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_hotQCD::get_pressure(double e, double rhob) const {
    double f = interpolate1D(e, 0, pressure_tb);  // 1/fm^4
    return(std::max(Util::small_eps, f));
}


double EOS_hotQCD::get_s2e(double s, double rhob) const {
    double e = get_s2e_finite_rhob(s, 0.0);
    return(e);
}

double EOS_hotQCD::get_T2e(double T, double rhob) const {
    double e = get_T2e_finite_rhob(T, 0.0);
    return(e);
}
