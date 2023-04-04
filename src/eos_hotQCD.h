// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_hotQCD_H_
#define SRC_EOS_hotQCD_H_

#include "eos_base.h"

class EOS_hotQCD : public EOS_base {
 private:
    const int eos_id;

 public:
    EOS_hotQCD(const int eos_id_in);

    void initialize_eos();
    double p_e_func       (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_temperature(double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_pressure   (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_s2e        (double s, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_T2e        (double T, double rhob, double rhoq=0.0, double rhos=0.0) const;

    void get_pressure_with_gradients(double epsilon, double rhob, double rhoq, 
		    double rhos, double &p, double &dpde, double &dpdrhob, 
		    double &dpdrhoq, double &dpdrhos, double &cs2) const;
 
    void check_eos() const {check_eos_no_muB();}
};

#endif  // SRC_EOS_hotQCD_H_
