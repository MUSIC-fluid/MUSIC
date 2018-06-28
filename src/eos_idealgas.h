// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_IDEALGAS_H_
#define SRC_EOS_IDEALGAS_H_

#include "eos_base.h"

class EOS_idealgas : public EOS_base {
 private:
     double Nc;
     double Nf;
   
 public:
    EOS_idealgas();
    ~EOS_idealgas() {}
    
    double get_cs2        (double e, double rhob) const {return(1./3.);}
    double p_rho_func     (double e, double rhob) const {return(0.0);}
    double p_e_func       (double e, double rhob) const {return(1./3.);}
    double get_temperature(double e, double rhob) const;
    double get_mu         (double e, double rhob) const {return(0.0);}
    double get_muS        (double e, double rhob) const {return(0.0);}
    double get_pressure   (double e, double rhob) const {return(1./3.*e);}
    double get_s2e        (double s, double rhob) const;

    void check_eos() const {check_eos_no_muB();}
};

#endif  // SRC_EOS_IDEALGAS_H_
