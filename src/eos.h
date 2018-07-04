// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_H_
#define SRC_EOS_H_

#include "eos_idealgas.h"
#include "eos_hotQCD.h"

#include <iostream>

//! This is a wrapper class for the equation of state
class EOS {
 private:
    const int eos_id; EOS_idealgas eos_ideal;
    EOS_hotQCD eos_HQCD;

    // function pointers
    double (EOS::*pressure_ptr)   (double e, double rhob) const;
    double (EOS::*temperature_ptr)(double e, double rhob) const;
    double (EOS::*entropy_ptr)    (double e, double rhob) const;
    double (EOS::*cs2_ptr)        (double e, double rhob) const;
    double (EOS::*dpde_ptr)       (double e, double rhob) const;
    double (EOS::*dpdrhob_ptr)    (double e, double rhob) const;
    double (EOS::*muB_ptr)        (double e, double rhob) const;
    double (EOS::*muS_ptr)        (double e, double rhob) const;
    double (EOS::*s2e_ptr)        (double e, double rhob) const;
    double (EOS::*get_eps_max_ptr)() const;
    void   (EOS::*check_eos_ptr)  () const;

 public:
    EOS() = default;
    EOS(const int eos_id_in);

    ~EOS() {};

    // functions from the ideal gas EOS
    double get_pressure_idealgas   (double e, double rhob) const {return(eos_ideal.get_pressure(e, rhob));}
    double get_temperature_idealgas(double e, double rhob) const {return(eos_ideal.get_temperature(e, rhob));}
    double get_entropy_idealgas    (double e, double rhob) const {return(eos_ideal.get_entropy(e, rhob));}
    double get_cs2_idealgas        (double e, double rhob) const {return(eos_ideal.get_cs2(e, rhob));}
    double get_dpde_idealgas       (double e, double rhob) const {return(eos_ideal.p_e_func(e, rhob));}
    double get_dpdrhob_idealgas    (double e, double rhob) const {return(eos_ideal.p_rho_func(e, rhob));}
    double get_muB_idealgas        (double e, double rhob) const {return(eos_ideal.get_mu(e, rhob));}
    double get_muS_idealgas        (double e, double rhob) const {return(eos_ideal.get_muS(e, rhob));}
    double get_s2e_idealgas        (double s, double rhob) const {return(eos_ideal.get_s2e(s, rhob));}
    double get_eps_max_idealgas    () const {return(eos_ideal.get_eps_max());}
    void   check_eos_idealgas      () const {return(eos_ideal.check_eos());}
    
    // functions from the hot QCD EOS
    double get_pressure_hotQCD   (double e, double rhob) const {return(eos_HQCD.get_pressure(e, rhob));}
    double get_temperature_hotQCD(double e, double rhob) const {return(eos_HQCD.get_temperature(e, rhob));}
    double get_entropy_hotQCD    (double e, double rhob) const {return(eos_HQCD.get_entropy(e, rhob));}
    double get_cs2_hotQCD        (double e, double rhob) const {return(eos_HQCD.get_cs2(e, rhob));}
    double get_dpde_hotQCD       (double e, double rhob) const {return(eos_HQCD.p_e_func(e, rhob));}
    double get_dpdrhob_hotQCD    (double e, double rhob) const {return(eos_HQCD.p_rho_func(e, rhob));}
    double get_muB_hotQCD        (double e, double rhob) const {return(eos_HQCD.get_mu(e, rhob));}
    double get_muS_hotQCD        (double e, double rhob) const {return(eos_HQCD.get_muS(e, rhob));}
    double get_s2e_hotQCD        (double s, double rhob) const {return(eos_HQCD.get_s2e(s, rhob));}
    double get_eps_max_hotQCD    () const {return(eos_HQCD.get_eps_max());}
    void   check_eos_hotQCD      () const {return(eos_HQCD.check_eos());}


    // functions to call the function pointers
    double get_pressure   (double e, double rhob) const {return((this->*(pressure_ptr))(e, rhob));}
    double get_temperature(double e, double rhob) const {return((this->*(temperature_ptr))(e, rhob));}
    double get_entropy    (double e, double rhob) const {return((this->*(entropy_ptr))(e, rhob));}
    double get_cs2        (double e, double rhob) const {return((this->*(cs2_ptr))(e, rhob));}
    double get_dpde       (double e, double rhob) const {return((this->*(dpde_ptr))(e, rhob));}
    double get_dpdrhob    (double e, double rhob) const {return((this->*(dpdrhob_ptr))(e, rhob));}
    double get_muB        (double e, double rhob) const {return((this->*(muB_ptr))(e, rhob));}
    double get_muS        (double e, double rhob) const {return((this->*(muS_ptr))(e, rhob));}
    double get_s2e        (double s, double rhob) const {return((this->*(s2e_ptr))(s, rhob));}

    double get_eps_max() const {return((this->*(get_eps_max_ptr))());}
    void   check_eos()   const {return((this->*(check_eos_ptr))());}
};

#endif  // SRC_EOS_H_
