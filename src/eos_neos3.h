// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_NEOS3_H_
#define SRC_EOS_NEOS3_H_

#include "eos_base.h"

class EOS_neos3 : public EOS_base {
 private:
   
 public:
    EOS_neos3();
    ~EOS_neos3() {}
    
    void initialize_eos();
    double get_cs2        (double e, double rhob) const;
    double p_rho_func     (double e, double rhob) const;
    double p_e_func       (double e, double rhob) const;
    double get_temperature(double e, double rhob) const;
    double get_mu         (double e, double rhob) const;
    double get_muS        (double e, double rhob) const {return(0.0);}
    double get_pressure   (double e, double rhob) const;
    double get_s2e        (double s, double rhob) const;
    double get_dpOverdrhob2(double e, double rhob) const;

    void check_eos() const {check_eos_with_finite_muB();}
};

#endif  // SRC_EOS_NEOS2_H_
