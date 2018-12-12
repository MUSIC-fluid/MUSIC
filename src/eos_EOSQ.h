// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_EOSQ_H_
#define SRC_EOS_EOSQ_H_

#include "eos_base.h"

class EOS_eosQ : public EOS_base {
 private:
   
 public:
    EOS_eosQ();
    ~EOS_eosQ();
    
    void initialize_eos();
    double p_rho_func     (double e, double rhob) const;
    double p_e_func       (double e, double rhob) const;
    double get_temperature(double e, double rhob) const;
    double get_mu         (double e, double rhob) const;
    double get_pressure   (double e, double rhob) const;
    double get_s2e        (double s, double rhob) const;

    void check_eos() const {check_eos_with_finite_muB();}
};

#endif  // SRC_EOS_EOSQ_H_
