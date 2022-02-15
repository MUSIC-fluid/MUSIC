// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_s95p_H_
#define SRC_EOS_s95p_H_

#include "eos_base.h"

class EOS_s95p : public EOS_base {
 private:
    const int eos_id;
   
 public:
    EOS_s95p(const int eos_id_in);
    
    void initialize_eos();
    double p_e_func       (double e, double rhob) const;
    double get_temperature(double e, double rhob) const;
    double get_pressure   (double e, double rhob) const;
    double get_s2e        (double s, double rhob) const;
    double get_T2e        (double T_in_GeV, double rhob) const;

    void check_eos() const {check_eos_no_muB();}
};

#endif  // SRC_EOS_s95p_H_
