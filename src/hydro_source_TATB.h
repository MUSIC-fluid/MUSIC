// Copyright 2020 Chun Shen

#ifndef SRC_HYDRO_SOURCE_TATB_H_
#define SRC_HYDRO_SOURCE_TATB_H_

#include <memory>
#include <vector>

#include "hydro_source_base.h"

struct participant {
    double x, y;
    int dir;
    double e;
};

class HydroSourceTATB : public HydroSourceBase {
  private:
    InitData &DATA_;
    double yL_frac_;
    double tau_source;
    double TA_, TB_, TATB_;
    std::vector<std::vector<double>> profile_TA;
    std::vector<std::vector<double>> profile_TB;
    double ybeam_, tanhYbeam_, cosh2Ybeam_;
    // grid on which profile_TA/TB are stored (the initial hydro grid). It is
    // kept fixed even if beastMode later changes the hydro grid.
    double gridXmin_, gridYmin_, gridDX_, gridDY_;

    bool profileIndexInRange(const int ix, const int iy) const {
        return (
            ix >= 0 && iy >= 0 && ix < static_cast<int>(profile_TA.size())
            && iy < static_cast<int>(profile_TA[ix].size()));
    }
    int gridNX_, gridNY_;
    double preEqFlowFactor_;
    double eta0_, eta_m_, sigma_eta_;
    double beta_;

  public:
    //! Integrals of exp(k eta) cosh(eta) and exp(k eta) sinh(eta) over the
    //! plateau [-eta_m, eta_m], finite also at k = +-1 (|ln a| = 2 eta_m)
    static void plateauExpIntegrals(
        const double k, const double eta_m, double &I_cosh, double &I_sinh);

    HydroSourceTATB() = delete;
    HydroSourceTATB(InitData &DATA_in);
    ~HydroSourceTATB();

    //! This function reads in the spatal information of the nuclear thickness
    //! functions
    void read_in_TATB();

    void read_in_participants_and_compute_TATB();
    void computeTATB(const double x_0, const double y_0, const int dir);

    double eta_rhob_left_factor(const double eta) const;
    double eta_rhob_right_factor(const double eta) const;
    double energy_eta_profile_normalisation(
        const double y_CM, const double eta_0, const double sigma_eta) const;
    double eta_profile_plateau(
        const double eta, const double eta_0, const double sigma_eta) const;
    double compute_yL(
        const double TA, const double TB, const double ycm, const double eta_0,
        const double sigma_eta, const double eta_m) const;
    double eta_envelope_f(const double eta, const double eta_m) const;
    double energy_eta_profile_normalisation_tilted(
        const double TA, const double TB, const double eta_0,
        const double eta_m, const double sigma_eta, const double ycm,
        const double M, const double yL) const;
    void compute_energy_density(
        const double TA, const double TB, const double eta_s, double &epsilon,
        double &y_L) const;
    //! this function returns the energy source term J^\mu at a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu, EnergyFlowVec &j_mu) const;

    //! this function returns the net baryon density source term rho
    //! at a given point (tau, x, y, eta_s)
    double get_hydro_rhob_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu) const;
    double eta_profile_plateau_frag(
        const double eta, const double eta_0, const double sigma_eta) const;
    double energy_eta_profile_normalisation_numerical(
        const double y_CM, const double eta_0, const double sigma_eta) const;
};

#endif  // SRC_HYDRO_SOURCE_TATB_H_
