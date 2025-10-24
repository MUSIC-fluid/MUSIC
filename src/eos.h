// Copyright 2018 @ Chun Shen
#ifndef SRC_EOS_H_
#define SRC_EOS_H_

#include <memory>
#include <vector>
#include "eos_base.h"

//! Wrapper class for the equation of state
class EOS {
  private:
    const int eos_id;
    std::unique_ptr<EOS_base> eos_ptr;

  public:
    EOS() = delete;
    EOS(const int eos_id_in);
    ~EOS() {};

    // Core thermodynamics
    double get_pressure(double e, double rhob) const {
        return eos_ptr->get_pressure(e, rhob);
    }
    void get_pressure_with_gradients(
        double e, double rhob, double &p, double &dpde, double &dpdrhob, double &cs2) const {
        eos_ptr->get_pressure_with_gradients(e, rhob, p, dpde, dpdrhob, cs2);
    }
    double get_temperature(double e, double rhob) const {
        return eos_ptr->get_temperature(e, rhob);
    }
    double get_entropy(double e, double rhob) const {
        return eos_ptr->get_entropy(e, rhob);
    }
    double get_cs2(double e, double rhob) const {
        return eos_ptr->get_cs2(e, rhob);
    }
    double get_dpde(double e, double rhob) const {
        return eos_ptr->p_e_func(e, rhob);
    }
    double get_dpdrhob(double e, double rhob) const {
        return eos_ptr->p_rho_func(e, rhob);
    }
    double get_muB(double e, double rhob) const {
        return eos_ptr->get_muB(e, rhob);
    }
    double get_muS(double e, double rhob) const {
        return eos_ptr->get_muS(e, rhob);
    }
    double get_muQ(double e, double rhob) const {
        return eos_ptr->get_muQ(e, rhob);
    }
    double get_s2e(double s, double rhob) const {
        return eos_ptr->get_s2e(s, rhob);
    }
    double get_T2e(double T_in_GeV, double rhob) const {
        return eos_ptr->get_T2e(T_in_GeV, rhob);
    }
    void getThermalVariables(double e, double rhob, std::vector<double> &thermalVec) const {
        eos_ptr->getThermalVariables(e, rhob, thermalVec);
    }

    double get_eps_max() const { return eos_ptr->get_eps_max(); }
    void check_eos()     const { return eos_ptr->check_eos(); }

    // ---- Chemistry helpers ----
    // IQCD(e) := I(e)|_{muB=0} = ε - 3P(e, rhob=0)
    inline double get_IQCD(double e) const {
        return e - 3.0 * get_pressure(e, 0.0);
    }
    // I0(e): baseline (chem-eq) reference; 0 for now until you wire a model
    inline double get_I0(double /*e*/) const {
        return 0.0;
    }
    // Full interaction measure at arbitrary rhob (general utility)
    inline double get_I(double e, double rhob) const {
        return e - 3.0 * get_pressure(e, rhob);
    }
    inline double get_dIde(double e, double rhob) const {
        // dI/de = 1 − 3 * (dP/de)
        return 1.0 - 3.0 * get_dpde(e, rhob);
    }
};

#endif  // SRC_EOS_H_
