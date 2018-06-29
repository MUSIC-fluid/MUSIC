// Copyright 2018 @ Chun Shen
#include "eos_idealgas.h"

#include <cmath>

EOS_idealgas::EOS_idealgas() {
    set_EOS_id(0);
    set_number_of_tables(0);
    set_eps_max(1e5);
    Nc = 3.;
    Nf = 2.5;
}

double EOS_idealgas::get_temperature(double eps, double rhob) const {
    return pow(90.0/M_PI/M_PI*(eps/3.0)/(2*(Nc*Nc-1)+7./2*Nc*Nf), .25);
}

double EOS_idealgas::get_s2e(double s, double rhob) const {
    return(3./4.*s*pow(3.*s/4./(M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0), 1./3.));  // in 1/fm^4
}


