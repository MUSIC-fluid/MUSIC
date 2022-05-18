// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_neos_H_
#define SRC_EOS_neos_H_

#include "eos_base.h"

class EOS_neos : public EOS_base {
 private:
    const int eos_id;

 public:
    EOS_neos(const int eos_id_in);
    ~EOS_neos();

    void initialize_eos();
    double p_rho_func     (double e, double rhob) const;
    double p_e_func       (double e, double rhob) const;
    double get_temperature(double e, double rhob) const;
    double get_muB        (double e, double rhob) const;
    double get_muS        (double e, double rhob) const;
    double get_muC        (double e, double rhob) const;
    double get_pressure   (double e, double rhob) const;
    double get_s2e        (double s, double rhob) const;

    void check_eos() const {
        check_eos_with_finite_muB();
        outputMutable();
    }
};

#endif  // SRC_EOS_neos_H_
