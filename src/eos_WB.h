// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_WB_H_
#define SRC_EOS_WB_H_

#include "eos_base.h"

class EOS_WB : public EOS_base {
 private:

 public:
    EOS_WB();
    ~EOS_WB() {}

    void initialize_eos();
    double get_cs2        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double p_e_func       (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_temperature(double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_pressure   (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_s2e        (double s, double rhob, double rhoq=0.0, double rhos=0.0) const;

    void check_eos() const {check_eos_no_muB();}
};

#endif  // SRC_EOS_WB_H_
