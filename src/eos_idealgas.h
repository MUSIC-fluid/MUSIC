// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_IDEALGAS_H_
#define SRC_EOS_IDEALGAS_H_

#include "eos_base.h"
#include "util.h"

class EOS_idealgas : public EOS_base {
 private:
     double Nc;
     double Nf;

 public:
    EOS_idealgas();
    ~EOS_idealgas() {}

    void set_Nc_and_Nf(double Nc_in, double Nf_in) {Nc = Nc_in; Nf = Nf_in;}

    void   initialize_eos();
    double get_cs2        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(1./3.);}
    double p_rho_func     (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(0.0);}
    double p_e_func       (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(1./3.);}
    double get_temperature(double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muB        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muQ        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muS        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_pressure   (double e, double rhob, double rhoq=0.0, double rhos=0.0) const {return(1./3.*e);}
    double get_s2e        (double s, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_T2e        (double T_in_GeV, double rhob, double rhoq=0.0, double rhos=0.0) const;

    void get_pressure_with_gradients(
            double epsilon, double rhob, double rhoq,
            double rhos, double &p, double &dpde, double &dpdrhob,
            double &dpdrhoq, double &dpdrhos) const {
        dpde = 1/3.;
        p = dpde*epsilon;
        dpdrhob = 0.;
        dpdrhoq = 0.;
        dpdrhos = 0.;
    }

    void check_eos() const {check_eos_no_muB();}
};

#endif  // SRC_EOS_IDEALGAS_H_
