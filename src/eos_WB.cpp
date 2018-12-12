// Copyright 2018 @ Chun Shen

#include "eos_WB.h"
#include "util.h"

#include <sstream>
#include <fstream>

using std::stringstream;
using std::string;

EOS_WB::EOS_WB() {
    set_EOS_id(8);
    set_number_of_tables(0);
    set_eps_max(1e5);
    set_flag_muB(false);
    set_flag_muS(false);
    set_flag_muC(false);
}


void EOS_WB::initialize_eos() {
    // read the lattice EOS pressure, temperature, and 
    music_message.info("Using lattice EOS parameterization from WB");
}


double EOS_WB::get_cs2(double e, double rhob) const {
    double f = calculate_velocity_of_sound_sq(e, rhob);
    return(f);
}
    

double EOS_WB::p_e_func(double e_local, double rhob) const {
    double cs2_local;
    double e1 = e_local;
	double e2 = e1*e1;
	double e3 = e2*e1;
	double e4 = e3*e1;
	double e5 = e4*e1;
	double e6 = e5*e1;
	double e7 = e6*e1;
	double e8 = e7*e1;
	double e9 = e8*e1;
	double e10 = e9*e1;
	double e11 = e10*e1;
	double e12 = e11*e1;
	double e13 = e12*e1;
	cs2_local = ((5.191934309650155e-32 + 4.123605749683891e-23*e1
                 + 3.1955868410879504e-16*e2 + 1.4170364808063119e-10*e3
                 + 6.087136671592452e-6*e4 + 0.02969737949090831*e5
                 + 15.382615282179595*e6 + 460.6487249985994*e7
                 + 1612.4245252438795*e8 + 275.0492627924299*e9
                 + 58.60283714484669*e10 + 6.504847576502024*e11
                 + 0.03009027913262399*e12 + 8.189430244031285e-6*e13)
		        /(1.4637868900982493e-30 + 6.716598285341542e-22*e1
                  + 3.5477700458515908e-15*e2 + 1.1225580509306008e-9*e3
                  + 0.00003551782901018317*e4 + 0.13653226327408863*e5
                  + 60.85769171450653*e6 + 1800.5461219450308*e7
                  + 15190.225535036281*e8 + 590.2572000057821*e9
                  + 293.99144775704605*e10 + 21.461303090563028*e11
                  + 0.09301685073435291*e12 + 0.000024810902623582917*e13));
    return(cs2_local);
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_WB::get_temperature(double e_local, double rhob) const {
    double temperature;
    double e1 = e_local;
	double e2 = e1*e1;
	double e3 = e2*e1;
	double e4 = e3*e1;
	double e5 = e4*e1;
	double e6 = e5*e1;
	double e7 = e6*e1;
	double e8 = e7*e1;
	double e9 = e8*e1;
	double e10 = e9*e1;
	double e11 = e10*e1;
	temperature = ((1.510073201405604e-29 + 8.014062800678687e-18*e1
                    + 2.4954778310451065e-10*e2 + 0.000063810382643387*e3
                    + 0.4873490574161924*e4 + 207.48582344326206*e5
                    + 6686.07424325115*e6 + 14109.766109389702*e7
                    + 1471.6180520527757*e8 + 14.055788949565482*e9
                    + 0.015421252394182246*e10 + 1.5780479034557783e-6*e11)
                   /(7.558667139355393e-28 + 1.3686372302041508e-16*e1
                     + 2.998130743142826e-9*e2 + 0.0005036835870305458*e3
                     + 2.316902328874072*e4 + 578.0778724946719*e5
                     + 11179.193315394154*e6 + 17965.67607192861*e7
                     + 1051.0730543534657*e8 + 5.916312075925817*e9
                     + 0.003778342768228011*e10 + 1.8472801679382593e-7*e11));
    return(temperature);
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_WB::get_pressure(double e_local, double rhob) const {
    double p;
    double e1 = e_local;
    double e2 = e1*e_local;
    double e3 = e2*e_local;
    double e4 = e3*e_local;
    double e5 = e4*e_local;
    double e6 = e5*e_local;
    double e7 = e6*e_local;
    double e8 = e7*e_local;
    double e9 = e8*e_local;
    double e10 = e9*e_local;
    double e11 = e10*e_local;
    double e12 = e11*e_local;
	
	p = ((  1.9531729608963267e-11*e12 + 3.1188455176941583e-7*e11
          + 0.0009417586777847889*e10 + 0.7158279081255019*e9
          + 141.5073484468774*e8 + 6340.448389300905*e7
          + 41913.439282708554*e6 + 334334.4309240126*e5
          + 1.6357487344679043e6*e4 + 3.1729694865420084e6*e3
          + 1.077580993288114e6*e2 + 9737.845799644809*e1
          - 0.25181736420168666)
         /(  3.2581066229887368e-18*e12 + 5.928138360995685e-11*e11
           + 9.601103399348206e-7*e10 + 0.002962497695527404*e9
           + 2.3405487982094204*e8 + 499.04919730607065*e7
           + 26452.34905933697*e6 + 278581.2989342773*e5
           + 1.7851642641834426e6*e4 + 1.3512402226067686e7*e3
           + 2.0931169138134286e7*e2 + 4.0574329080826794e6*e1
           + 45829.44617893836));
    p = std::max(1e-16, p);
    return(p);
}


double EOS_WB::get_s2e(double s, double rhob) const {
    double e = get_s2e_finite_rhob(s, 0.0);
    return(e);
}
