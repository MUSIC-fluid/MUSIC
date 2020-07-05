// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale

#include <cmath>
#include "util.h"
#include "transport_coeffs.h"

using Util::hbarc;

TransportCoeffs::TransportCoeffs(const EOS &eosIn, const InitData &Data_in)
    : DATA(Data_in), eos(eosIn) {
    shear_relax_time_factor_ = DATA.shear_relax_time_factor;
    bulk_relax_time_factor_  = DATA.bulk_relax_time_factor;
}


double TransportCoeffs::get_eta_over_s(double T) const {
    double eta_over_s;
    if (DATA.T_dependent_shear_to_s == 1) {
        eta_over_s = get_temperature_dependent_eta_over_s_default(T);
    } else if (DATA.T_dependent_bulk_to_s == 2) {
        eta_over_s = get_temperature_dependent_eta_over_s_duke(T);
    } else if (DATA.T_dependent_bulk_to_s == 3) {
        eta_over_s = get_temperature_dependent_eta_over_s_sims(T);
    } else {
        eta_over_s = DATA.shear_to_s;
    }
    return eta_over_s;
}


double TransportCoeffs::get_temperature_dependent_eta_over_s_default(double T) const {
    double Ttr = 0.18/hbarc;  // phase transition temperature
    double Tfrac = T/Ttr;
    double eta_over_s;
    if (T < Ttr) {
        eta_over_s = (DATA.shear_to_s + 0.0594*(1. - Tfrac)
                      + 0.544*(1. - Tfrac*Tfrac));
    } else {
        eta_over_s = (DATA.shear_to_s + 0.288*(Tfrac - 1.)
                      + 0.0818*(Tfrac*Tfrac - 1.));
    }
    return(eta_over_s);
}


double TransportCoeffs::get_temperature_dependent_eta_over_s_duke(
                                                double T_in_fm) const {
    double T_in_GeV = T_in_fm*hbarc;
    double Ttr_in_GeV = 0.154;
    double Tfrac = T_in_GeV/Ttr_in_GeV;

    double eta_over_s = (DATA.eta_over_s_min
                         + (DATA.eta_over_s_slope)*(T_in_GeV - Ttr_in_GeV)
                           *pow(Tfrac,DATA.eta_over_s_curv));
    return eta_over_s;
}


double TransportCoeffs::get_temperature_dependent_eta_over_s_sims(
                                                double T_in_fm) const {
    double T_in_GeV = T_in_fm*hbarc;
    double T_kink_in_GeV = DATA.eta_over_s_T_kink_in_GeV;
    double low_T_slope = DATA.eta_over_s_low_T_slope_in_GeV;
    double high_T_slope = DATA.eta_over_s_high_T_slope_in_GeV;
    double eta_over_s_at_kink = DATA.eta_over_s_at_kink;

    const double eta_over_s_min = 1e-6;
    double eta_over_s;
    if (T_in_GeV < T_kink_in_GeV) {
        eta_over_s = (eta_over_s_at_kink
                      + low_T_slope*(T_in_GeV - T_kink_in_GeV));
    } else {
        eta_over_s = (eta_over_s_at_kink
                      + high_T_slope*(T_in_GeV - T_kink_in_GeV));
    }
    eta_over_s = std::max(eta_over_s,eta_over_s_min);
    return(eta_over_s);
}


double TransportCoeffs::get_zeta_over_s(double T) const {
    double zeta_over_s = 0.;;
    if (DATA.T_dependent_bulk_to_s == 2) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_duke(T);
    } else if (DATA.T_dependent_bulk_to_s == 3) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_sims(T);
    } else if (DATA.T_dependent_bulk_to_s == 1) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_default(T);
    }
    zeta_over_s = std::max(0., zeta_over_s);
    return zeta_over_s;
}


//! Cauchy distribution
double TransportCoeffs::get_temperature_dependent_zeta_over_s_duke(
                                                double T_in_fm) const {
    const double A=DATA.bulk_viscosity_normalisation;
    const double G=DATA.bulk_viscosity_width_in_GeV;
    const double Tpeak_in_GeV=DATA.bulk_viscosity_peak_in_GeV;
    const double T_in_GeV=T_in_fm*hbarc;
    const double diff_ratio=(T_in_GeV-Tpeak_in_GeV)/G;

    //const double T_delta=(T_in_GeV*T_in_GeV)/(Tpeak_in_GeV*Tpeak_in_GeV)-1;
    //double  bulk_over_sden=A*(G*G)/(T_delta*T_delta+G*G);
    return A/(1+diff_ratio*diff_ratio);
}


//! Skewed Cauchy distribution
double TransportCoeffs::get_temperature_dependent_zeta_over_s_sims(
                                                    double T_in_fm) const {
    const double T_in_GeV=T_in_fm*hbarc;

    const double max=DATA.zeta_over_s_max;
    const double width=DATA.zeta_over_s_width_in_GeV;
    const double T_peak_in_GeV=DATA.zeta_over_s_T_peak_in_GeV;
    const double lambda=DATA.zeta_over_s_lambda_asymm;
    const double diff=T_in_GeV-T_peak_in_GeV;
    const double sign=(diff > 0) - (diff < 0);
    const double diff_ratio=(diff)/(width*(lambda*sign+1));

    return max/(1+diff_ratio*diff_ratio);
}


//! T dependent bulk viscosity from Gabriel
double TransportCoeffs::get_temperature_dependent_zeta_over_s_default(
                                                double temperature) const {
    /////////////////////////////////////////////
    //           Parametrization 1             //
    /////////////////////////////////////////////
    double Ttr=0.18/0.1973;
    double dummy=temperature/Ttr;
    double A1=-13.77, A2=27.55, A3=13.45;
    double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;

    double bulk = A1*dummy*dummy + A2*dummy - A3;
    if (temperature < 0.995*Ttr) {
        bulk = (lambda3*exp((dummy-1)/sigma3)
                + lambda4*exp((dummy-1)/sigma4) + 0.03);
    }
    if (temperature > 1.05*Ttr) {
        bulk = (lambda1*exp(-(dummy-1)/sigma1)
                + lambda2*exp(-(dummy-1)/sigma2) + 0.001);
    }

    /////////////////////////////////////////////
    //           Parametrization 2             //
    /////////////////////////////////////////////
    //double Ttr=0.18/0.1973;
    //double dummy=temperature/Ttr;
    //double A1=-79.53, A2=159.067, A3=79.04;
    //double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    //double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;

    //bulk = A1*dummy*dummy + A2*dummy - A3;

    //if (temperature < 0.997*Ttr) {
    //    bulk = (lambda3*exp((dummy-1)/sigma3)
    //            + lambda4*exp((dummy-1)/sigma4) + 0.03);
    //}
    //if (temperature > 1.04*Ttr) {
    //    bulk = (lambda1*exp(-(dummy-1)/sigma1)
    //            + lambda2*exp(-(dummy-1)/sigma2) + 0.001);
    //}

    ////////////////////////////////////////////
    //           Parametrization 3            //
    ////////////////////////////////////////////
    //double Ttr=0.18/0.1973;
    //double dummy=temperature/Ttr;
    //double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    //double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;

    //if (temperature<0.99945*Ttr) {
    //    bulk = (lambda3*exp((dummy-1)/sigma3)
    //            + lambda4*exp((dummy-1)/sigma4) + 0.03);
    //}
    //if (temperature>0.99945*Ttr) {
    //    bulk = 0.901*exp(14.5*(1.0-dummy)) + 0.061/dummy/dummy;
    //}

    return(bulk);
}


