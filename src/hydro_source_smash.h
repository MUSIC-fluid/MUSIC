// Copyright 2019 Chun Shen

#ifndef SRC_HYDRO_SOURCE_SMASH_H_
#define SRC_HYDRO_SOURCE_SMASH_H_

#include <vector>
#include <memory>
#include "hydro_source_base.h"

//! This data structure stores hadron information
struct hadron {
    int pdgid;
    double tau, x, y, eta_s;
    double rapidity;
    double rapidity_perp;
    double E, px, py;
    double mass;
    int baryon_number;
    int strangness;
    int electric_charge;
};

class HydroSourceSMASH : public HydroSourceBase {
 private:
    const InitData &DATA;
    double parton_quench_factor;
    std::vector<hadron> list_hadrons_;
    std::vector<hadron> list_hadrons_current_tau_;
    std::vector<hadron> list_spectators_;

    int average_events_;
    int number_events_;
    double weight_event_;

    int baryon_total_;
    double p0_total_;
    double px_total_;
    double py_total_;
    double pz_total_;

 public:
    HydroSourceSMASH() = default;
    HydroSourceSMASH(const InitData &DATA_in);
    ~HydroSourceSMASH();
    
    //! This function reads in the hadron information from the SMASH model
    void read_in_SMASH_hadrons(int i_event,
        int extended_output, int reject_spectators);
    
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

#endif  // SRC_HYDRO_SOURCE_SMASH_H_

