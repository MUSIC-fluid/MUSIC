// Copyright 2018 @ Chun Shen

#include "eos_neos.h"
#include "util.h"

#include <sstream>
#include <fstream>
#include <cmath>

using std::stringstream;
using std::string;

EOS_neos::EOS_neos(const int eos_id_in) : eos_id(eos_id_in) {
    set_EOS_id(eos_id);
    set_number_of_tables(0);
    set_eps_max(1e5);
    set_flag_muB(true);
    set_flag_muS(false);
    set_flag_muC(false);
}


EOS_neos::~EOS_neos() {
    int ntables = get_number_of_tables();
    for (int itable = 0; itable < ntables; itable++) {
        Util::mtx_free(mu_B_tb[itable],
                       nb_length[itable], e_length[itable]);
        if (get_flag_muS()) {
            Util::mtx_free(mu_S_tb[itable],
                           nb_length[itable], e_length[itable]);
        }
        if (get_flag_muC()) {
            Util::mtx_free(mu_C_tb[itable],
                           nb_length[itable], e_length[itable]);
        }
    }
    if (ntables > 0) {
        delete [] mu_B_tb;
        if (get_flag_muS()) {
            delete [] mu_S_tb;
        }
        if (get_flag_muC()) {
            delete [] mu_C_tb;
        }
    }
}


void EOS_neos::initialize_eos() {
    // read the lattice EOS pressure, temperature, and 
    music_message.info("Using lattice EOS at finite muB from A. Monnai");

    auto envPath = get_hydro_env_path();
    stringstream spath;
    spath << envPath;

    bool flag_muS = false;
    bool flag_muC = false;
    string eos_file_string_array[7];
    if (eos_id == 10) {
        music_message.info("reading EOS neos ...");
        spath << "/EOS/neos_2/";
        string string_tmp[] = {"0a", "0b", "0c", "1a", "2", "3", "4"};
        std::copy(std::begin(string_tmp), std::end(string_tmp),
                  std::begin(eos_file_string_array));
    } else if (eos_id == 11) {
        music_message.info("reading EOS neos3 ...");
        spath << "/EOS/neos_3/";
        string string_tmp[] = {"1", "2", "3", "4", "5", "6", "7"};
        std::copy(std::begin(string_tmp), std::end(string_tmp),
                  std::begin(eos_file_string_array));
    } else if (eos_id == 12) {
        music_message.info("reading EOS neos_b ...");
        spath << "/EOS/neos_b/";
        string string_tmp[] = {"1", "2", "3", "4", "5", "6", "7"};
        std::copy(std::begin(string_tmp), std::end(string_tmp),
                  std::begin(eos_file_string_array));
    } else if (eos_id == 13) {
        music_message.info("reading EOS neos_bs ...");
        spath << "/EOS/neos_bs/";
        string string_tmp[] = {"1s", "2s", "3s", "4s", "5s", "6s", "7s"};
        std::copy(std::begin(string_tmp), std::end(string_tmp),
                  std::begin(eos_file_string_array));
        flag_muS = true;
        set_flag_muS(flag_muS);
    } else if (eos_id == 14) {
        music_message.info("reading EOS neos_bqs ...");
        spath << "/EOS/neos_bqs/";
        string string_tmp[] = {"1qs", "2qs", "3qs", "4qs", "5qs", "6qs", "7qs"};
        std::copy(std::begin(string_tmp), std::end(string_tmp),
                  std::begin(eos_file_string_array));
        flag_muS = true;
        flag_muC = true;
        set_flag_muS(flag_muS);
        set_flag_muC(flag_muC);
    }

    string path = spath.str();
    music_message << "from path " << path;
    music_message.flush("info");

    const int ntables = 7;
    set_number_of_tables(ntables);
    resize_table_info_arrays();

    pressure_tb    = new double** [ntables];
    temperature_tb = new double** [ntables];
    mu_B_tb        = new double** [ntables];
    if (flag_muS) {
        mu_S_tb = new double** [ntables];
    }
    if (flag_muC) {
        mu_C_tb = new double** [ntables];
    }

    for (int itable = 0; itable < ntables; itable++) {
        std::ifstream eos_p(path + "neos" + eos_file_string_array[itable]
                            + "_p.dat");
        if (!eos_p.is_open()) {
            music_message << "Can not open EOS files: "
                          << (path + "neos" + eos_file_string_array[itable]
                              + "_p.dat");
            music_message.flush("error");
            exit(1);
        }
        std::ifstream eos_T(path + "neos" + eos_file_string_array[itable]
                            + "_t.dat");
        std::ifstream eos_mub(path + "neos" + eos_file_string_array[itable]
                              + "_mub.dat");
        std::ifstream eos_muS;
        std::ifstream eos_muC;
        if (flag_muS) {
            eos_muS.open(path + "neos" + eos_file_string_array[itable]
                         + "_mus.dat");
            if (!eos_muS.is_open()) {
                music_message << "Can not open EOS files: "
                              << (path + "neos" + eos_file_string_array[itable]
                                  + "_mus.dat");
                music_message.flush("error");
                exit(1);
            }
        }
        if (flag_muC) {
            eos_muC.open(path + "neos" + eos_file_string_array[itable]
                         + "_muq.dat");
            if (!eos_muC.is_open()) {
                music_message << "Can not open EOS files: "
                              << (path + "neos" + eos_file_string_array[itable]
                                  + "_muc.dat");
                music_message.flush("error");
                exit(1);
            }
        }
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
        if (flag_muS) {
            std::getline(eos_muS, dummy);
            std::getline(eos_muS, dummy);
        }
        if (flag_muC) {
            std::getline(eos_muC, dummy);
            std::getline(eos_muC, dummy);
        }

        // allocate memory for EOS arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
        mu_B_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                           e_length[itable]);
        if (flag_muS) {
            mu_S_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        }
        if (flag_muC) {
            mu_C_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        }

        // read pressure, temperature and chemical potential values
        for (int j = 0; j < e_length[itable]; j++) {
            for (int i = 0; i < nb_length[itable]; i++) {
                eos_p >> pressure_tb[itable][i][j];
                eos_T >> temperature_tb[itable][i][j];
                eos_mub >> mu_B_tb[itable][i][j];

                if (flag_muS) {
                    eos_muS >> mu_S_tb[itable][i][j];
                    mu_S_tb[itable][i][j] /= Util::hbarc;    // 1/fm
                }
                if (flag_muC) {
                    eos_muC >> mu_C_tb[itable][i][j];
                    mu_C_tb[itable][i][j] /= Util::hbarc;    // 1/fm
                }

                pressure_tb[itable][i][j]    /= Util::hbarc;    // 1/fm^4
                temperature_tb[itable][i][j] /= Util::hbarc;    // 1/fm
                temperature_tb[itable][i][j] = pow(temperature_tb[itable][i][j], 5.);    // 1/fm^5
                mu_B_tb[itable][i][j]        /= Util::hbarc;    // 1/fm
            }
        }
    }

    //double eps_max_in = (e_bounds[6] + e_spacing[6]*e_length[6])/hbarc;
    double eps_max_in = e_bounds[6] + e_spacing[6]*e_length[6];
    set_eps_max(eps_max_in);

    music_message.info("Done reading EOS.");
}


double EOS_neos::p_e_func(double e, double rhob) const {
    return(get_dpOverde3(e, rhob));
}


double EOS_neos::p_rho_func(double e, double rhob) const {
    return(get_dpOverdrhob2(e, rhob));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_neos::get_temperature(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double T5 = interpolate2D(e, std::abs(rhob), table_idx,
                              temperature_tb);  // 1/fm^5
    T5 = std::max(Util::small_eps, T5);
    double T = pow(T5, 0.2);  // 1/fm
    return(T);
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_neos::get_pressure(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double f = interpolate2D(e, std::abs(rhob), table_idx, pressure_tb);
    f = std::max(Util::small_eps, f);
    return(f);
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_neos::get_muB(double e, double rhob) const {
    int table_idx = get_table_idx(e);
    double sign = rhob/(std::abs(rhob) + Util::small_eps);
    double mu = sign*interpolate2D(e, std::abs(rhob), table_idx,
                                   mu_B_tb);  // 1/fm
    return(mu);
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_neos::get_muS(double e, double rhob) const {
    if (!get_flag_muS()) return(0.0);
    int table_idx = get_table_idx(e);
    double sign = rhob/(std::abs(rhob) + Util::small_eps);
    double mu = sign*interpolate2D(e, std::abs(rhob), table_idx,
                                   mu_S_tb);  // 1/fm
    return(mu);
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_neos::get_muC(double e, double rhob) const {
    if (!get_flag_muC()) return(0.0);
    int table_idx = get_table_idx(e);
    double sign = rhob/(std::abs(rhob) + Util::small_eps);
    double mu = sign*interpolate2D(e, std::abs(rhob), table_idx,
                                   mu_C_tb);  // 1/fm
    return(mu);
}


double EOS_neos::get_s2e(double s, double rhob) const {
    double e = get_s2e_finite_rhob(s, rhob);
    return(e);
}
