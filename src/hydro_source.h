// Copyright 2017 Chun Shen
#ifndef SRC_HYDRO_SOURCE_H_
#define SRC_HYDRO_SOURCE_H_

#include <vector>
#include "./data.h"

//! This data structure contains a QCD string object
struct QCD_string {
    double norm;              // normalization for the string energy
    double delta_E;           // the energy difference between
                              // before and after the collisions [GeV]
    double tau_form;
    double x_perp, y_perp;    // transverse position of the string
    double eta_s_left, eta_s_right;
    double y_l, y_r;          // rapidity of the two ends of the string
};


//! This data structure stores parton information
struct parton {
    double tau, x, y, eta_s;
    double rapidity;
    double baryon_number;     //!< nucleon = 1, quark = 1/3
    // this data structure can be extended to include charge and strangeness
    // quantum number
};

//! This is a class that feeds source currents to hydrodynamics
class hydro_source {
 private:
    InitData *DATA_ptr;
    double volume;
    double sigma_tau, sigma_x, sigma_eta;
    std::vector<QCD_string> QCD_strings_list;
    std::vector<parton> parton_list;

 public:
    hydro_source(InitData *DATA_in);
    ~hydro_source();

    //! this function returns the energy source term J^\mu at a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source(double tau, double x, double y, double eta_s,
                                 double *u_mu, double *j_mu);
    //! this function returns the net baryon density source term rho
    //! at a given point (tau, x, y, eta_s)
    double get_hydro_rhob_source(double tau, double x, double y, double eta_s);

    void get_hydro_energy_source_before_tau(
        double tau, double x, double y, double eta_s, double *j_mu);
    double get_hydro_rhob_source_before_tau(double tau, double x, double y,
                                            double eta_s);

    //! This function reads in the spatial positions of the QCD strings
    //! and partons
    void read_in_QCD_strings_and_partons();
};

#endif  // SRC_HYDRO_SOURCE_H_
