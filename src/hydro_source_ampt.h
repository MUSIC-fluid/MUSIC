// Copyright 2019 Chun Shen

#ifndef SRC_HYDRO_SOURCE_AMPT_H_
#define SRC_HYDRO_SOURCE_AMPT_H_

#include <vector>
#include <memory>
#include "hydro_source_base.h"

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


class HydroSourceAMPT : public HydroSourceBase {
 private:
    const InitData &DATA;
    double parton_quench_factor;
    std::vector<std::shared_ptr<parton>> parton_list;
    std::vector<std::shared_ptr<parton>> parton_list_current_tau;

 public:
    HydroSourceAMPT() = default;
    HydroSourceAMPT(const InitData &DATA_in);
    ~HydroSourceAMPT();
    
    //! This function reads in the partons information from the AMPT model
    void read_in_AMPT_partons();
    
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
};

#endif  // SRC_HYDRO_SOURCE_AMPT_H_

