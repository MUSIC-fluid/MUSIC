// Copyright 2017 Chun Shen
#ifndef SRC_HYDRO_SOURCE_H_
#define SRC_HYDRO_SOURCE_H_

#include <array>
#include <vector>
#include <memory>

#include "data.h"
#include "pretty_ostream.h"
#include "data_struct.h"


//! This data structure contains a QCD string object
struct QCD_string {
    double norm;              // normalization for the string energy
    double delta_E;           // the energy difference between
                              // before and after the collisions [GeV]
    double tau_form;
    double tau_start, eta_s_start;
    double tau_0, eta_s_0;
    double x_perp, y_perp;    // transverse position of the string
    double tau_end_left, tau_end_right;
    double eta_s_left, eta_s_right;
    double y_l, y_r;          // rapidity of the two ends of the string
    double frac_l, frac_r;
    double y_l_i, y_r_i;
};


//! This data structure stores parton information
struct parton {
    double tau, x, y, eta_s;
    double rapidity;
    double rapidity_perp;
    double E, px, py;
    double mass;
    double baryon_number;
    double strangness;
    double electric_charge;
};


//! This is a class that feeds source currents to hydrodynamics
class hydro_source {
 private:
    const InitData &DATA;
    pretty_ostream music_message;
    double volume;
    double string_quench_factor;
    double parton_quench_factor;
    double sigma_tau, sigma_x, sigma_eta;
    int string_dump_mode;
    double source_tau_max;
    double source_tau_min;
    std::vector<std::shared_ptr<QCD_string>> QCD_strings_list;
    std::vector<std::shared_ptr<QCD_string>> QCD_strings_list_current_tau;
    std::vector<std::shared_ptr<QCD_string>> QCD_strings_baryon_list_current_tau;
    std::vector<std::shared_ptr<parton>> parton_list;
    std::vector<std::shared_ptr<parton>> parton_list_current_tau;

 public:
    hydro_source() = default;
    hydro_source(const InitData &DATA_in);
    ~hydro_source();

    //! this function returns the energy source term J^\mu at a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source(double tau, double x, double y, double eta_s,
                                 FlowVec &u_mu, EnergyFlowVec &j_mu);

    //! this function returns the net baryon density source term rho
    //! at a given point (tau, x, y, eta_s)
    double get_hydro_rhob_source(double tau, double x, double y, double eta_s,
                                 FlowVec &u_mu);

    //! this function returns the energy source term J^\mu up to a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source_before_tau(
        double tau, double x, double y, double eta_s, double *j_mu);
    
    //! this function returns the net baryon density source term rho
    //! up to a given point (tau, x, y, eta_s)
    double get_hydro_rhob_source_before_tau(double tau, double x, double y,
                                            double eta_s);

    //! This function reads in the spatal information of the strings
    //! and partons which are produced from the MC-Glauber-LEXUS model
    void read_in_QCD_strings_and_partons();

    //! This function reads in the partons information from the AMPT model
    void read_in_AMPT_partons();

    //! Get the minimum and maximum tau for the source term
    double get_source_tau_min() const {return(source_tau_min);}
    double get_source_tau_max() const {return(source_tau_max);}

    void prepare_list_for_current_tau_frame(double tau_local);
    void compute_norm_for_strings();
};

#endif  // SRC_HYDRO_SOURCE_H_
