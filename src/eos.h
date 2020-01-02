// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_H_
#define SRC_EOS_H_

#include "eos_base.h"
#include <memory>

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
    double get_pressure   (double e, double rhob) const {return(eos_ptr->get_pressure(e, rhob));}
    double get_temperature(double e, double rhob) const {return(eos_ptr->get_temperature(e, rhob));}
    double get_entropy    (double e, double rhob) const {return(eos_ptr->get_entropy(e, rhob));}
    double get_cs2        (double e, double rhob) const {return(eos_ptr->get_cs2(e, rhob));}
    double get_dpde       (double e, double rhob) const {return(eos_ptr->p_e_func(e, rhob));}
    double get_dpdrhob    (double e, double rhob) const {return(eos_ptr->p_rho_func(e, rhob));}
    double get_muB        (double e, double rhob) const {return(eos_ptr->get_muB(e, rhob));}
    double get_muS        (double e, double rhob) const {return(eos_ptr->get_muS(e, rhob));}
    double get_muC        (double e, double rhob) const {return(eos_ptr->get_muC(e, rhob));}
    double get_s2e        (double s, double rhob) const {return(eos_ptr->get_s2e(s, rhob));}
    double get_T2e        (double T, double rhob) const {return(eos_ptr->get_T2e(T, rhob));}

    double get_eps_max() const {return(eos_ptr->get_eps_max());}
    void   check_eos()   const {return(eos_ptr->check_eos());}
};

#endif  // SRC_EOS_H_
