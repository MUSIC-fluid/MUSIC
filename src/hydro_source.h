// Copyright 2019 Chun Shen

#ifndef SRC_HYDRO_SOURCE_H_
#define SRC_HYDRO_SOURCE_H_

#include "data.h"
#include "hydro_source_base.h"
#include <memory>

class HydroSource {
 private:
    const InitData &DATA;
    std::unique_ptr<HydroSourceBase> hydro_source_ptr;

 public:
    HydroSource() = default;
    HydroSource(const InitData &DATA_in);

    ~HydroSource() {};
    
    //! this function returns the energy source term J^\mu at a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu, EnergyFlowVec &j_mu) const {
        hydro_source_ptr->get_hydro_energy_source(
                                    tau, x, y, eta_s, u_mu, j_mu);
    }

    //! this function returns the net baryon density source term rho
    //! at a given point (tau, x, y, eta_s)
    double get_hydro_rhob_source(const double tau, const double x,
                                 const double y, const double eta_s,
                                 const FlowVec &u_mu) const {
        return(hydro_source_ptr->get_hydro_rhob_source(tau, x, y, eta_s, u_mu));
    }
    
    //! this function returns the energy source term J^\mu up to a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source_before_tau(
        const double tau, const double x, const double y, const double eta_s,
        EnergyFlowVec &j_mu) const {
        hydro_source_ptr->get_hydro_energy_source_before_tau(
                                            tau, x, y, eta_s, j_mu);
    }

    //! this function returns the net baryon density source term rho
    //! up to a given point (tau, x, y, eta_s)
    double get_hydro_rhob_source_before_tau(const double tau, const double x,
                                            const double y,
                                            const double eta_s) const {
        return(hydro_source_ptr->get_hydro_rhob_source_before_tau(
                                                    tau, x, y, eta_s));
    }


    void prepare_list_for_current_tau_frame(const double tau_local) {
        hydro_source_ptr->prepare_list_for_current_tau_frame(tau_local);
    }

    //! Get the minimum and maximum tau for the source term
    double get_source_tau_min() const {
        return(hydro_source_ptr->get_source_tau_min());
    }

    double get_source_tau_max() const {
        return(hydro_source_ptr->get_source_tau_max());
    }
};

#endif  // SRC_HYDRO_SOURCE_H_
