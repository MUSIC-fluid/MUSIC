// Copyright 2019 Chun Shen

#ifndef SRC_HYDRO_SOURCE_SMASH_H_
#define SRC_HYDRO_SOURCE_SMASH_H_

#include <vector>
#include <memory>
#include "hydro_source_base.h"

//! This data structure stores hadron information
struct hadron {
    int pdgid;
    int ID;
    double tau, x, y, eta_s, t, z;
    double rapidity;
    double E, px, py, pz;
    double mass;
    int baryon_number;
    int electric_charge;
    int strangeness;
};

class HydroSourceSMASH : public HydroSourceBase {
 private:
    InitData &DATA;
    double parton_quench_factor;
    std::vector<hadron> list_hadrons_;
    std::vector<hadron> list_hadrons_current_tau_;
    std::vector<hadron> list_spectators_;
    int covariant_smearing_kernel_;

    int average_events_;
    int number_events_;
    double weight_event_;

    int baryon_total_;
    int charge_total_;
    int strangeness_total_;
    double p0_total_;
    double px_total_;
    double py_total_;
    double pz_total_;
    double current_covariant_norm_;

    enum QuantityType { BARYON_NUMBER, ELECTRIC_CHARGE, STRANGENESS, ENERGY_SOURCE };

 public:
    HydroSourceSMASH() = default;
    HydroSourceSMASH(InitData &DATA_in);
    ~HydroSourceSMASH();
    
    //! This function reads in the hadron information from the SMASH model
    void read_in_SMASH_hadrons(int i_event,
        int extended_output, int reject_spectators);
    
    //! set the covariant smearing kernel to true or false
    void set_covariant_smearing_kernel(int cov_kernel) { 
        covariant_smearing_kernel_ = cov_kernel;
    }

    //! Compute the norm for the covariant smearing kernel
    void compute_covariant_norm(double tau);

    //! defines a covariant formulation of a smearing kernel
    double covariant_smearing_kernel(const double x_diff, const double y_diff, 
        const double z_diff, const double ux, const double uy, 
        const double ueta, const double sigma) const;

    //! this function returns the energy source term J^\mu at a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu, EnergyFlowVec &j_mu) const ;

    //! this function computes the source for a given QuantityType at point
    //! (tau, x, y, eta_s)
    double calculate_source(const double tau, const double x, 
                        const double y, const double eta_s,
                        const FlowVec &u_mu, const int quantityType) const ;

    //! this function returns the net baryon density source term rho
    //! at a given point (tau, x, y, eta_s)
    double get_hydro_rhob_source(const double tau, const double x,
                                 const double y, const double eta_s,
                                 const FlowVec &u_mu) const ;
    double get_hydro_rhoq_source(const double tau, const double x,
                                 const double y, const double eta_s,
                                 const FlowVec &u_mu) const ;
    double get_hydro_rhos_source(const double tau, const double x,
                                 const double y, const double eta_s,
                                 const FlowVec &u_mu) const ;

    void prepare_list_for_current_tau_frame(const double tau_local);
};

#endif  // SRC_HYDRO_SOURCE_SMASH_H_

