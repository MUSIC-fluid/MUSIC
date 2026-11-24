// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_IDEALGAS_H_
#define SRC_EOS_IDEALGAS_H_

#include <vector>

#include "eos_base.h"
#include "util.h"

//! Free massless quark-gluon gas with Nc colors and Nf flavors,
//!   P(T, mu_B) = a T^4 + Nc Nf (mu_q^2 T^2 / 6 + mu_q^4 / (12 pi^2)),
//! with mu_q = mu_B / 3 and a = pi^2 / 90 (2 (Nc^2 - 1) + 7/2 Nc Nf).
//! It is conformal, P = e / 3, and thermodynamically consistent at any n_B,
//! so s = (e + P - mu_B n_B) / T = dP/dT > 0. At e and n_B the state is
//! found in x = mu_B / T, which is capped at muB_over_T_max_; beyond it
//! (cold, baryon-rich cells) the state of the cap is used, as for the
//! tables (see EOS_base::get_rhob_in_table).
class EOS_idealgas : public EOS_base {
  private:
    double Nc;
    double Nf;
    static constexpr double muB_over_T_max_ = 10.;

    double coeff_a() const;
    double coeff_b() const { return (Nc * Nf); }
    // P = T^4 D(x), n_B = T^3 N(x), s = T^3 S(x), with dD/dx = N
    double D_of_x(double x) const;
    double N_of_x(double x) const;
    double S_of_x(double x) const;
    //! mu_B / T >= 0 at (e, |n_B|), capped at muB_over_T_max_
    double get_muB_over_T(double e, double rhob) const;
    double T_of_e_x(double e, double x) const;
    //! x = mu_B / T at rho = (n_B / e^{3/4}) / (its value at the cap) < 1,
    //! by Newton from the tabulated guess
    double x_of_rho(double rho, double x_guess) const;
    // x on a uniform grid in rho, the starting point of x_of_rho
    std::vector<double> x_table_;
    double r_cap_ = 0.;  // n_B / e^{3/4} at the cap
    void build_x_table();

  public:
    EOS_idealgas();
    ~EOS_idealgas() {}

    void set_Nc_and_Nf(double Nc_in, double Nf_in) {
        Nc = Nc_in;
        Nf = Nf_in;
        build_x_table();
    }

    void initialize_eos();
    double get_cs2(double e, double rhob) const { return (1. / 3.); }
    double p_rho_func(double e, double rhob) const { return (0.0); }
    double get_rhob_in_table(double e, double rhob) const override;
    // one solve for x per call
    double get_entropy(double e, double rhob) const override;
    void getThermalVariables(
        double e, double rhob, std::vector<double> &thermalVec) const override;
    double p_e_func(double e, double rhob) const { return (1. / 3.); }
    double get_temperature(double e, double rhob) const;
    double get_muB(double e, double rhob) const;
    double get_muS(double e, double rhob) const { return (0.0); }
    double get_pressure(double e, double rhob) const { return (1. / 3. * e); }
    double get_s2e(double s, double rhob) const;
    double get_T2e(double T_in_GeV, double rhob) const;
    virtual void get_pressure_with_gradients(
        double epsilon, double rhob, double &p, double &dpde, double &dpdrhob,
        double &cs2) const {
        dpde = 1 / 3.;
        p = dpde * epsilon;
        dpdrhob = 0.;
        cs2 = dpde;
    }

    void check_eos() const { check_eos_no_muB(); }
};

#endif  // SRC_EOS_IDEALGAS_H_
