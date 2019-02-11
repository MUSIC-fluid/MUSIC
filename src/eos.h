// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_H_
#define SRC_EOS_H_

#include "eos_idealgas.h"
#include "eos_EOSQ.h"
#include "eos_s95p.h"
#include "eos_WB.h"
#include "eos_hotQCD.h"
#include "eos_best.h"
#include "eos_neos.h"

//! This is a wrapper class for the equation of state
class EOS {
 private:
    const int eos_id;

    EOS_idealgas ideal;
    EOS_eosQ eosQ;
    EOS_s95p s95p;
    EOS_WB WB;
    EOS_hotQCD hotQCD;
    EOS_BEST best;
    EOS_neos neos;

    // function pointers
    double (EOS::*pressure_ptr)   (double e, double rhob) const;
    double (EOS::*temperature_ptr)(double e, double rhob) const;
    double (EOS::*entropy_ptr)    (double e, double rhob) const;
    double (EOS::*cs2_ptr)        (double e, double rhob) const;
    double (EOS::*dpde_ptr)       (double e, double rhob) const;
    double (EOS::*dpdrhob_ptr)    (double e, double rhob) const;
    double (EOS::*muB_ptr)        (double e, double rhob) const;
    double (EOS::*muS_ptr)        (double e, double rhob) const;
    double (EOS::*muC_ptr)        (double e, double rhob) const;
    double (EOS::*s2e_ptr)        (double s, double rhob) const;
    double (EOS::*T2e_ptr)        (double T, double rhob) const;
    double (EOS::*get_eps_max_ptr)() const;
    void   (EOS::*check_eos_ptr)  () const;

 public:
    EOS() = default;
    EOS(const int eos_id_in);

    ~EOS() {};

    // functions from the ideal gas EOS
    double get_pressure_idealgas   (double e, double rhob) const {return(ideal.get_pressure(e, rhob));}
    double get_temperature_idealgas(double e, double rhob) const {return(ideal.get_temperature(e, rhob));}
    double get_entropy_idealgas    (double e, double rhob) const {return(ideal.get_entropy(e, rhob));}
    double get_cs2_idealgas        (double e, double rhob) const {return(ideal.get_cs2(e, rhob));}
    double get_dpde_idealgas       (double e, double rhob) const {return(ideal.p_e_func(e, rhob));}
    double get_dpdrhob_idealgas    (double e, double rhob) const {return(ideal.p_rho_func(e, rhob));}
    double get_muB_idealgas        (double e, double rhob) const {return(ideal.get_muB(e, rhob));}
    double get_muS_idealgas        (double e, double rhob) const {return(ideal.get_muS(e, rhob));}
    double get_muC_idealgas        (double e, double rhob) const {return(ideal.get_muC(e, rhob));}
    double get_s2e_idealgas        (double s, double rhob) const {return(ideal.get_s2e(s, rhob));}
    double get_T2e_idealgas        (double T, double rhob) const {return(ideal.get_T2e(T, rhob));}
    double get_eps_max_idealgas    () const {return(ideal.get_eps_max());}
    void   check_eos_idealgas      () const {return(ideal.check_eos());}
    
    // functions from the EOSQ
    double get_pressure_eosQ   (double e, double rhob) const {return(eosQ.get_pressure(e, rhob));}
    double get_temperature_eosQ(double e, double rhob) const {return(eosQ.get_temperature(e, rhob));}
    double get_entropy_eosQ    (double e, double rhob) const {return(eosQ.get_entropy(e, rhob));}
    double get_cs2_eosQ        (double e, double rhob) const {return(eosQ.get_cs2(e, rhob));}
    double get_dpde_eosQ       (double e, double rhob) const {return(eosQ.p_e_func(e, rhob));}
    double get_dpdrhob_eosQ    (double e, double rhob) const {return(eosQ.p_rho_func(e, rhob));}
    double get_muB_eosQ        (double e, double rhob) const {return(eosQ.get_muB(e, rhob));}
    double get_muS_eosQ        (double e, double rhob) const {return(eosQ.get_muS(e, rhob));}
    double get_muC_eosQ        (double e, double rhob) const {return(eosQ.get_muC(e, rhob));}
    double get_s2e_eosQ        (double s, double rhob) const {return(eosQ.get_s2e(s, rhob));}
    double get_T2e_eosQ        (double T, double rhob) const {return(eosQ.get_T2e(T, rhob));}
    double get_eps_max_eosQ    () const {return(eosQ.get_eps_max());}
    void   check_eos_eosQ      () const {return(eosQ.check_eos());}
    
    // functions from the hot QCD EOS
    double get_pressure_hotQCD   (double e, double rhob) const {return(hotQCD.get_pressure(e, rhob));}
    double get_temperature_hotQCD(double e, double rhob) const {return(hotQCD.get_temperature(e, rhob));}
    double get_entropy_hotQCD    (double e, double rhob) const {return(hotQCD.get_entropy(e, rhob));}
    double get_cs2_hotQCD        (double e, double rhob) const {return(hotQCD.get_cs2(e, rhob));}
    double get_dpde_hotQCD       (double e, double rhob) const {return(hotQCD.p_e_func(e, rhob));}
    double get_dpdrhob_hotQCD    (double e, double rhob) const {return(hotQCD.p_rho_func(e, rhob));}
    double get_muB_hotQCD        (double e, double rhob) const {return(hotQCD.get_muB(e, rhob));}
    double get_muS_hotQCD        (double e, double rhob) const {return(hotQCD.get_muS(e, rhob));}
    double get_muC_hotQCD        (double e, double rhob) const {return(hotQCD.get_muC(e, rhob));}
    double get_s2e_hotQCD        (double s, double rhob) const {return(hotQCD.get_s2e(s, rhob));}
    double get_T2e_hotQCD        (double T, double rhob) const {return(hotQCD.get_T2e(T, rhob));}
    double get_eps_max_hotQCD    () const {return(hotQCD.get_eps_max());}
    void   check_eos_hotQCD      () const {return(hotQCD.check_eos());}
    
    // functions from the s95p EOS
    double get_pressure_s95p   (double e, double rhob) const {return(s95p.get_pressure(e, rhob));}
    double get_temperature_s95p(double e, double rhob) const {return(s95p.get_temperature(e, rhob));}
    double get_entropy_s95p    (double e, double rhob) const {return(s95p.get_entropy(e, rhob));}
    double get_cs2_s95p        (double e, double rhob) const {return(s95p.get_cs2(e, rhob));}
    double get_dpde_s95p       (double e, double rhob) const {return(s95p.p_e_func(e, rhob));}
    double get_dpdrhob_s95p    (double e, double rhob) const {return(s95p.p_rho_func(e, rhob));}
    double get_muB_s95p        (double e, double rhob) const {return(s95p.get_muB(e, rhob));}
    double get_muS_s95p        (double e, double rhob) const {return(s95p.get_muS(e, rhob));}
    double get_muC_s95p        (double e, double rhob) const {return(s95p.get_muC(e, rhob));}
    double get_s2e_s95p        (double s, double rhob) const {return(s95p.get_s2e(s, rhob));}
    double get_T2e_s95p        (double T, double rhob) const {return(s95p.get_T2e(T, rhob));}
    double get_eps_max_s95p    () const {return(s95p.get_eps_max());}
    void   check_eos_s95p      () const {return(s95p.check_eos());}
    
    // functions from the WB EOS
    double get_pressure_WB   (double e, double rhob) const {return(WB.get_pressure(e, rhob));}
    double get_temperature_WB(double e, double rhob) const {return(WB.get_temperature(e, rhob));}
    double get_entropy_WB    (double e, double rhob) const {return(WB.get_entropy(e, rhob));}
    double get_cs2_WB        (double e, double rhob) const {return(WB.get_cs2(e, rhob));}
    double get_dpde_WB       (double e, double rhob) const {return(WB.p_e_func(e, rhob));}
    double get_dpdrhob_WB    (double e, double rhob) const {return(WB.p_rho_func(e, rhob));}
    double get_muB_WB        (double e, double rhob) const {return(WB.get_muB(e, rhob));}
    double get_muS_WB        (double e, double rhob) const {return(WB.get_muS(e, rhob));}
    double get_muC_WB        (double e, double rhob) const {return(WB.get_muC(e, rhob));}
    double get_s2e_WB        (double s, double rhob) const {return(WB.get_s2e(s, rhob));}
    double get_T2e_WB        (double T, double rhob) const {return(WB.get_T2e(T, rhob));}
    double get_eps_max_WB    () const {return(WB.get_eps_max());}
    void   check_eos_WB      () const {return(WB.check_eos());}

    // functions from the neos EOS
    double get_pressure_neos   (double e, double rhob) const {return(neos.get_pressure(e, rhob));}
    double get_temperature_neos(double e, double rhob) const {return(neos.get_temperature(e, rhob));}
    double get_entropy_neos    (double e, double rhob) const {return(neos.get_entropy(e, rhob));}
    double get_cs2_neos        (double e, double rhob) const {return(neos.get_cs2(e, rhob));}
    double get_dpde_neos       (double e, double rhob) const {return(neos.p_e_func(e, rhob));}
    double get_dpdrhob_neos    (double e, double rhob) const {return(neos.p_rho_func(e, rhob));}
    double get_muB_neos        (double e, double rhob) const {return(neos.get_muB(e, rhob));}
    double get_muS_neos        (double e, double rhob) const {return(neos.get_muS(e, rhob));}
    double get_muC_neos        (double e, double rhob) const {return(neos.get_muC(e, rhob));}
    double get_s2e_neos        (double s, double rhob) const {return(neos.get_s2e(s, rhob));}
    double get_T2e_neos        (double T, double rhob) const {return(neos.get_T2e(T, rhob));}
    double get_eps_max_neos    () const {return(neos.get_eps_max());}
    void   check_eos_neos      () const {return(neos.check_eos());}
    
    // functions from the BEST EOS
    double get_pressure_best   (double e, double rhob) const {return(best.get_pressure(e, rhob));}
    double get_temperature_best(double e, double rhob) const {return(best.get_temperature(e, rhob));}
    double get_entropy_best    (double e, double rhob) const {return(best.get_entropy(e, rhob));}
    double get_cs2_best        (double e, double rhob) const {return(best.get_cs2(e, rhob));}
    double get_dpde_best       (double e, double rhob) const {return(best.p_e_func(e, rhob));}
    double get_dpdrhob_best    (double e, double rhob) const {return(best.p_rho_func(e, rhob));}
    double get_muB_best        (double e, double rhob) const {return(best.get_mu(e, rhob));}
    double get_muS_best        (double e, double rhob) const {return(best.get_muS(e, rhob));}
    double get_muC_best        (double e, double rhob) const {return(best.get_muC(e, rhob));}
    double get_s2e_best        (double s, double rhob) const {return(best.get_s2e(s, rhob));}
    double get_T2e_best        (double T, double rhob) const {return(best.get_T2e(T, rhob));}
    double get_eps_max_best    () const {return(best.get_eps_max());}
    void   check_eos_best      () const {return(best.check_eos());}

    // functions to call the function pointers
    double get_pressure   (double e, double rhob) const {return((this->*(pressure_ptr))(e, rhob));}
    double get_temperature(double e, double rhob) const {return((this->*(temperature_ptr))(e, rhob));}
    double get_entropy    (double e, double rhob) const {return((this->*(entropy_ptr))(e, rhob));}
    double get_cs2        (double e, double rhob) const {return((this->*(cs2_ptr))(e, rhob));}
    double get_dpde       (double e, double rhob) const {return((this->*(dpde_ptr))(e, rhob));}
    double get_dpdrhob    (double e, double rhob) const {return((this->*(dpdrhob_ptr))(e, rhob));}
    double get_muB        (double e, double rhob) const {return((this->*(muB_ptr))(e, rhob));}
    double get_muS        (double e, double rhob) const {return((this->*(muS_ptr))(e, rhob));}
    double get_muC        (double e, double rhob) const {return((this->*(muC_ptr))(e, rhob));}
    double get_s2e        (double s, double rhob) const {return((this->*(s2e_ptr))(s, rhob));}
    double get_T2e        (double T, double rhob) const {return((this->*(T2e_ptr))(T, rhob));}

    double get_eps_max() const {return((this->*(get_eps_max_ptr))());}
    void   check_eos()   const {return((this->*(check_eos_ptr))());}
};

#endif  // SRC_EOS_H_
