// Copyright 2019 Chun Shen

#ifndef SRC_HYDRO_SOURCE_STRINGS_H_
#define SRC_HYDRO_SOURCE_STRINGS_H_

#include <vector>
#include <memory>
#include "hydro_source_base.h"

//! This data structure contains a QCD string object
struct QCD_string {
    double norm;              // normalization for the string energy
    double E_remnant_norm_L, E_remnant_norm_R;
    double m_over_sigma;      // m/sigma [fm] sigma is the string tension

    double mass;
    double tau_form;
    double sigma_x, sigma_eta;
    double tau_start, eta_s_start;
    double tau_0, eta_s_0;
    double x_perp, y_perp;    // transverse position of the string
    double x_pl, y_pl, x_pr, y_pr;
    double tau_end_left, tau_end_right;
    double eta_s_left, eta_s_right;
    double y_l, y_r;          // rapidity of the two ends of the string
    double remnant_l, remnant_r;
    double y_l_i, y_r_i;
    double tau_baryon_left, tau_baryon_right;
    double eta_s_baryon_left, eta_s_baryon_right;
    double y_l_baryon, y_r_baryon;
    double baryon_frac_l, baryon_frac_r;
};


class HydroSourceStrings : public HydroSourceBase {
 private:
    const InitData &DATA;
    int string_dump_mode;
    double string_quench_factor;
    double parton_quench_factor;
    std::vector<std::shared_ptr<QCD_string>> QCD_strings_list;
    std::vector<std::shared_ptr<QCD_string>> QCD_strings_list_current_tau;
    std::vector<std::shared_ptr<QCD_string>> QCD_strings_remnant_list_current_tau;
    std::vector<std::shared_ptr<QCD_string>> QCD_strings_baryon_list_current_tau;

 public:
    HydroSourceStrings() = default;
    HydroSourceStrings(const InitData &DATA_in);
    ~HydroSourceStrings();

    //! This function reads in the spatal information of the strings
    //! and partons which are produced from the MC-Glauber-LEXUS model
    void read_in_QCD_strings_and_partons();

    //! this function returns the energy source term J^\mu at a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu, EnergyFlowVec &j_mu) const ;

    //! this function returns the net baryon density source term rho
    //! at a given point (tau, x, y, eta_s)
    double get_hydro_rhob_source(const double tau, const double x,
                                 const double y, const double eta_s,
                                 const FlowVec &u_mu) const ;

    void prepare_list_for_current_tau_frame(const double tau_local);
    void compute_norm_for_strings(const double total_energy);
};

#endif  // SRC_HYDRO_SOURCE_STRINGS_H_
