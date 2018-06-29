// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_WRAPPER_H_
#define SRC_EOS_WRAPPER_H_

#include "eos_idealgas.h"
#include "eos_hotQCD.h"

class EOS_wrapper {
 private:
    double (EOS_wrapper::*pressure_ptr)(double e, double rhob);

    EOS_idealgas eos_ig;
    EOS_hotQCD eos_HQ;

 public:
    EOS_wrapper() {
        pressure_ptr = &EOS_wrapper::get_pressure_idealgas;
    }

    ~EOS_wrapper() {};

    double get_pressure_idealgas(double e, double rhob) {
        return(eos_ig.get_pressure(e, rhob));
    }

    double get_pressure(double e, double rhob) {
        return((this->*(pressure_ptr))(e, rhob));
    }


};

#endif  // SRC_EOS_WRAPPER_H_
