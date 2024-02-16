// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_H_
#define SRC_EOS_H_

#include "eos_base.h"
#include <memory>
#include <vector>

//! This is a wrapper class for the equation of state
class EOS {
 private:
    const int eos_id;

    std::unique_ptr<EOS_base> eos_ptr;

 public:
    EOS() = default;
    EOS(const int eos_id_in);

    ~EOS() {};

    // functions to call the function pointers
    double get_pressure(double e, double rhob,
                        double rhoq=0.0, double rhos=0.0) const {
        return(eos_ptr->get_pressure(e, rhob));
    }

    void get_pressure_with_gradients(double epsilon, double rhob, double rhoq,
            double rhos, double &p, double &dpde, double &dpdrhob,
            double &dpdrhoq, double &dpdrhos) const {
        eos_ptr->get_pressure_with_gradients(
                epsilon, rhob, rhoq, rhos,
                p, dpde, dpdrhob, dpdrhoq, dpdrhos);
    }

    void get_pressure_with_gradients_and_cs2(
            double epsilon, double rhob, double rhoq,
            double rhos, double &p, double &dpde, double &dpdrhob,
            double &dpdrhoq, double &dpdrhos, double &cs2) const {
        eos_ptr->get_pressure_with_gradients_and_cs2(
                epsilon, rhob, rhoq, rhos,
                p, dpde, dpdrhob, dpdrhoq, dpdrhos, cs2);
    }

    double get_temperature(double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(eos_ptr->get_temperature(e, rhob, rhoq, rhos));}
    double get_entropy    (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(eos_ptr->get_entropy(e, rhob, rhoq, rhos));}
    double get_cs2        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(eos_ptr->get_cs2(e, rhob, rhoq, rhos));}
    double get_dpde       (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(eos_ptr->p_e_func(e, rhob, rhoq, rhos));}
    double get_dpdrhob    (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(eos_ptr->p_rho_func(e, rhob));}
    double get_muB        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(eos_ptr->get_muB(e, rhob, rhoq, rhos));}
    double get_muS        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(eos_ptr->get_muS(e, rhob, rhoq, rhos));}
    double get_muQ        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(eos_ptr->get_muQ(e, rhob, rhoq, rhos));}
    double get_s2e        (double s, double rhob, double rhoq=0.0, double rhos=0.0) const {return(eos_ptr->get_s2e(s, rhob, rhoq, rhos));}
    double get_T2e        (double T_in_GeV, double rhob, double rhoq=0.0, double rhos=0.0) const {return(eos_ptr->get_T2e(T_in_GeV, rhob, rhoq, rhos));}

    void getThermalVariables(double e, double rhob, double rhoq, double rhos,
                             std::vector<double> &thermalVec) const {
        eos_ptr->getThermalVariables(e, rhob, rhoq, rhos, thermalVec);
    }

    double get_eps_max() const {return(eos_ptr->get_eps_max());}
    void   check_eos()   const {return(eos_ptr->check_eos());}
};

#endif  // SRC_EOS_H_
