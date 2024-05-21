// Copyright 2024 @ Chun Shen

#include "eos_1DGenerator.h"
#include "util.h"

#include <sstream>
#include <fstream>
#include <cmath>

using std::stringstream;
using std::string;

EOS_1DGenerator::EOS_1DGenerator(const int eos_id_in) : eos_id(eos_id_in) {
    set_EOS_id(eos_id);
    set_number_of_tables(0);
    set_eps_max(1e5);
    set_flag_muB(false);
    set_flag_muS(false);
    set_flag_muQ(false);
}


void EOS_1DGenerator::initialize_eos() {
    // read the lattice EOS pressure, temperature, and 
    music_message.info("reading EOS 1D Generator ...");

    auto envPath = get_hydro_env_path();
    stringstream slocalpath;
    slocalpath << envPath << "/EOS/EoS_1DGen.bin";

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
        eos_file.open(path, std::ios::binary);

        if (!eos_file) {
            music_message.error("Can not find the EoS file.");
            exit(1);
        }

        e_length[itable]  = 200;
        nb_length[itable] = 1;
        // allocate memory for pressure arrays
        pressure_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                               e_length[itable]);
        temperature_tb[itable] = Util::mtx_malloc(nb_length[itable],
                                                  e_length[itable]);
        double temp;
        for (int ii = 0; ii < e_length[itable]; ii++) {
            eos_file.read((char*)&temp, sizeof(double));  // e^{1/4} (GeV)
            temp /= Util::hbarc;      // 1/fm
            if (ii == 0) e_bounds[itable] = temp;
            if (ii == 1) e_spacing[itable] = temp - e_bounds[itable];
            if (ii == e_length[itable] - 1) set_eps_max(temp);

            eos_file.read((char*)&temp, sizeof(double));  // P (GeV^4)
            pressure_tb[itable][0][ii] = temp/pow(Util::hbarc, 4);  // 1/fm^4

            eos_file.read((char*)&temp, sizeof(double));  // T (GeV)
            temp /= Util::hbarc;   // 1/fm
            temperature_tb[itable][0][ii] = temp;
        }
    }
    music_message.info("Done reading EOS.");
}


double EOS_1DGenerator::p_e_func(double e, double rhob) const {
    return(get_dpOverde3(e, rhob));
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_1DGenerator::get_temperature(double e, double rhob) const {
    double etilde = pow(e, 0.25);
    double T = interpolate1D(etilde, 0, temperature_tb);  // 1/fm
    return(std::max(Util::small_eps, T));
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_1DGenerator::get_pressure(double e, double rhob) const {
    double etilde = pow(e, 0.25);
    double f = interpolate1D(etilde, 0, pressure_tb);  // 1/fm^4
    return(std::max(Util::small_eps, f));
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
void EOS_1DGenerator::get_pressure_with_gradients(double e, double rhob,
        double &p, double &dpde, double &dpdrhob, double &cs2) const {
    double etilde = pow(e, 0.25);
    interpolate1D_with_weighted_gradient(etilde, 0, pressure_tb, p, dpde);
    p = std::max(Util::small_eps, p);           // [1/fm^4]
    dpdrhob = 0.;
    //cs2 = std::max(0.01, std::min(1./3., dpde));
    cs2 = std::max(0.01, dpde);
}


double EOS_1DGenerator::get_s2e(double s, double rhob) const {
    double e = get_s2e_finite_rhob(s, 0.0);
    return(e);
}

double EOS_1DGenerator::get_T2e(double T_in_GeV, double rhob) const {
    double e = get_T2e_finite_rhob(T_in_GeV, 0.0);
    return(e);
}


void EOS_1DGenerator::interpolate1D_with_weighted_gradient(
    double e, int table_idx, double ***table, double &p, double &dpde) const {
// This is a generic linear interpolation routine for EOS at zero mu_B
// it assumes the class has already read in
//        P(e), T(e), s(e)
// as one-dimensional arrays on an equally spacing lattice grid
// units: e is in 1/fm^4
// This function uses the interpolation points to compute the local
// derivatives dP/de.
    const double e0 = e_bounds[table_idx];
    const double delta_e = e_spacing[table_idx];
    const int N_e = e_length[table_idx];

    // compute the indices
    int idx_e = static_cast<int>((e - e0)/delta_e);

    // treatment for overflow, use the last two points to do extrapolation
    idx_e = std::min(N_e - 2, std::max(0, idx_e));

    if (e < e0) {
        // check underflow
        dpde = table[table_idx][0][0]/e0;
        p = dpde*e;
        dpde = 5./4.*table[table_idx][0][0]/pow(e0, 5)*e;
    } else {
        double e1 = e0 + idx_e*delta_e;
        double e2 = e1 + delta_e;
        double temp1 = table[table_idx][0][idx_e];
        double temp2 = table[table_idx][0][idx_e + 1];
        p = temp1 + (temp2 - temp1)/delta_e*(e - e1);
        dpde = 5./4.*(temp2 - temp1)/(pow(e2, 5) - pow(e1, 5))*e;
    }
}
