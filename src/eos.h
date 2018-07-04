// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_H_
#define SRC_EOS_H_

#include "eos_idealgas.h"
#include "eos_hotQCD.h"


//! This is a wrapper class for the equation of state
class EOS {
 private:
    double (EOS::*pressure_ptr)(double e, double rhob);

    EOS_idealgas eos_ideal;
    EOS_hotQCD eos_HQCD;

 public:
    EOS() {
        pressure_ptr = &EOS::get_pressure_idealgas;
    }

    ~EOS() {};

    double get_pressure_idealgas(double e, double rhob) {
        return(eos_ideal.get_pressure(e, rhob));
    }

    double get_pressure(double e, double rhob) {
        return((this->*(pressure_ptr))(e, rhob));
    }


};

#endif  // SRC_EOS_H_
