// Copyright 2019 Chun Shen
#ifndef SRC_HYDRO_SOURCE_BASE_H_
#define SRC_HYDRO_SOURCE_BASE_H_

#include "data.h"
#include "pretty_ostream.h"
#include "data_struct.h"

class HydroSourceBase {
 private:
    double source_tau_max;
    double source_tau_min;
    double sigma_tau, sigma_x, sigma_eta;

 public:
    pretty_ostream music_message;

    HydroSourceBase() = default;
    virtual ~HydroSourceBase() {}

    virtual int get_number_of_sources() const {return(0);}

    void set_sigma_tau(double sigma_tau_in) {sigma_tau = sigma_tau_in;}
    void set_sigma_x  (double sigma_x_in  ) {sigma_x   = sigma_x_in  ;}
    void set_sigma_eta(double sigma_eta_in) {sigma_eta = sigma_eta_in;}
    double get_sigma_tau() const {return(sigma_tau);}
    double get_sigma_x  () const {return(sigma_x  );}
    double get_sigma_eta() const {return(sigma_eta);}

    //! Set the minimum and maximum tau for the source term
    void set_source_tau_min(double tau_in) {source_tau_min = tau_in;}
    void set_source_tau_max(double tau_in) {source_tau_max = tau_in;}

    //! Get the minimum and maximum tau for the source term
    double get_source_tau_min() const {return(source_tau_min);}
    double get_source_tau_max() const {return(source_tau_max);}

    //! this function returns the energy source term J^\mu at a given point
    //! (tau, x, y, eta_s)
    virtual void get_hydro_energy_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu, EnergyFlowVec &j_mu) const {
        j_mu = {0.0};
    }

    //! this function returns the net baryon density source term rho
    //! at a given point (tau, x, y, eta_s)
    virtual double get_hydro_rhob_source(const double tau, const double x,
                                 const double y, const double eta_s,
                                 const FlowVec &u_mu) const {
        return(0.0);
    }

    //! this function returns the energy source term J^\mu up to a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source_before_tau(
        const double tau, const double x, const double y, const double eta_s,
        EnergyFlowVec &j_mu) const;

    //! this function returns the net baryon density source term rho
    //! up to a given point (tau, x, y, eta_s)
    double get_hydro_rhob_source_before_tau(const double tau, const double x,
                                            const double y,
                                            const double eta_s) const;

    virtual void prepare_list_for_current_tau_frame(const double tau_local) {}
};

#endif  // SRC_HYDRO_SOURCE_BASE_H_
