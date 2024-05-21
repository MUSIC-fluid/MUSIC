// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_1DGenerator_H_
#define SRC_EOS_1DGenerator_H_

#include "eos_base.h"

class EOS_1DGenerator : public EOS_base {
 private:
    const int eos_id;

 public:
    EOS_1DGenerator(const int eos_id_in);

    void initialize_eos();
    double p_e_func       (double e, double rhob) const;
    double get_temperature(double e, double rhob) const;
    double get_pressure   (double e, double rhob) const;
    double get_s2e        (double s, double rhob) const;
    double get_T2e        (double T, double rhob) const;
    void get_pressure_with_gradients(double e, double rhob,
            double &p, double &dpde, double &dpdrhob, double &cs2) const;

    void check_eos() const {check_eos_no_muB();}
};

#endif  // SRC_EOS_hotQCD_H_
