// Copyright 2020 Chun Shen

#ifndef SRC_HYDRO_SOURCE_TATB_H_
#define SRC_HYDRO_SOURCE_TATB_H_

#include <vector>
#include <memory>
#include "hydro_source_base.h"


class HydroSourceTATB : public HydroSourceBase {
 private:
    const InitData &DATA_;
    double yL_frac_;
    double tau_source;
    double TA_, TB_;
    std::vector<std::vector<double>> profile_TA;
    std::vector<std::vector<double>> profile_TB;

 public:
    HydroSourceTATB() = default;
    HydroSourceTATB(const InitData &DATA_in);
    ~HydroSourceTATB();

    //! This function reads in the spatal information of the nuclear thickness
    //! functions
    void read_in_TATB();

    double eta_rhob_left_factor(const double eta) const;
    double eta_rhob_right_factor(const double eta) const;
    double energy_eta_profile_normalisation(
        const double y_CM, const double eta_0, const double sigma_eta) const;
    double eta_profile_plateau(
        const double eta, const double eta_0, const double sigma_eta) const;

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
    double eta_profile_plateau_frag(
        const double eta, const double eta_0, const double sigma_eta) const;
    double energy_eta_profile_normalisation_numerical(
        const double y_CM, const double eta_0, const double sigma_eta) const;
};

#endif  // SRC_HYDRO_SOURCE_TATB_H_
