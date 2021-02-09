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


double TransportCoeffs::get_eta_over_s(const double T, const double muB) const {
    // inputs T [1/fm], muB [1/fm]
    // outputs \eta/s
    double eta_over_s = DATA.shear_to_s;
    if (DATA.T_dependent_shear_to_s == 1) {
        eta_over_s = get_temperature_dependent_eta_over_s_default(T);
    } else if (DATA.T_dependent_shear_to_s == 2) {
        eta_over_s = get_temperature_dependent_eta_over_s_duke(T);
    } else if (DATA.T_dependent_shear_to_s == 3) {
        eta_over_s = get_temperature_dependent_eta_over_s_sims(T);
    } else if (DATA.T_dependent_shear_to_s == 11) {
        eta_over_s *= get_temperature_dependence_shear_profile(T);
    } else {
        eta_over_s = DATA.shear_to_s;
    }
    if (DATA.muB_dependent_shear_to_s == 10) {
        eta_over_s *= get_muB_dependence_shear_profile(muB);
    }
    return eta_over_s;
}


double TransportCoeffs::get_temperature_dependence_shear_profile(
                                            const double T_in_fm) const {
    const double T_in_GeV = T_in_fm*hbarc;
    const double Tc = 0.165;
    double f_T = 1.0;
    if (T_in_GeV < Tc) {
        const double Tslope = 1.2;
        const double Tlow = 0.1;
        f_T += Tslope*(Tc - T_in_GeV)/(Tc - Tlow);
    } else {
        const double Tslope2 = 0.0;
        const double Thigh = 0.4;
        f_T += Tslope2*(T_in_GeV - Tc)/(Thigh - Tc);
    }
    return(f_T);
}


double TransportCoeffs::get_muB_dependence_shear_profile(
                                                const double muB_in_fm) const {
    const double muB_in_GeV = muB_in_fm*hbarc;
    const double alpha = 0.8;
    const double muB_slope = 0.9;
    const double muB_scale = 0.6;
    double f_muB = 1. + muB_slope*pow(muB_in_GeV/muB_scale, alpha);
    return(f_muB);
}


double TransportCoeffs::get_temperature_dependent_eta_over_s_default(
                                                    const double T) const {
    const double Ttr = 0.18/hbarc;  // phase transition temperature
    const double Tfrac = T/Ttr;
    double eta_over_s;
    if (Tfrac < 1.) {
        eta_over_s = (DATA.shear_to_s + 0.0594*(1. - Tfrac)
                      + 0.544*(1. - Tfrac*Tfrac));
    } else {
        eta_over_s = (DATA.shear_to_s + 0.288*(Tfrac - 1.)
                      + 0.0818*(Tfrac*Tfrac - 1.));
    }
    return(eta_over_s);
}


double TransportCoeffs::get_temperature_dependent_eta_over_s_duke(
                                            const double T_in_fm) const {
    double T_in_GeV = T_in_fm*hbarc;
    double Ttr_in_GeV = 0.154;
    double Tfrac = T_in_GeV/Ttr_in_GeV;

    double eta_over_s = (DATA.shear_2_min
                         + (DATA.shear_2_slope)*(T_in_GeV - Ttr_in_GeV)
                           *pow(Tfrac,DATA.shear_2_curv));
    return eta_over_s;
}


double TransportCoeffs::get_temperature_dependent_eta_over_s_sims(
                                            const double T_in_fm) const {
    double T_in_GeV = T_in_fm*hbarc;
    double T_kink_in_GeV = DATA.shear_3_T_kink_in_GeV;
    double low_T_slope = DATA.shear_3_low_T_slope_in_GeV;
    double high_T_slope = DATA.shear_3_high_T_slope_in_GeV;
    double eta_over_s_at_kink = DATA.shear_3_at_kink;

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


double TransportCoeffs::get_zeta_over_s(const double T) const {
    // input T [1/fm]
    double zeta_over_s = 0.;
    if (DATA.T_dependent_bulk_to_s == 2) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_duke(T);
    } else if (DATA.T_dependent_bulk_to_s == 3) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_sims(T);
    } else if (DATA.T_dependent_bulk_to_s == 1) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_default(T);
    } else if (DATA.T_dependent_bulk_to_s == 7) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_bigbroadP(T);
    } else if (DATA.T_dependent_bulk_to_s == 8) {
        // latest param. for IPGlasma + MUSIC + UrQMD
        const double peak_norm = 0.13;
        zeta_over_s = get_temperature_dependent_zeta_over_s_AsymGaussian(
                                                                T, peak_norm);
    } else if (DATA.T_dependent_bulk_to_s == 9) {
        // latest param. for IPGlasma + KoMPoST + MUSIC + UrQMD
        const double peak_norm = 0.175;
        zeta_over_s = get_temperature_dependent_zeta_over_s_AsymGaussian(
                                                                T, peak_norm);
    }
    zeta_over_s = std::max(0., zeta_over_s);
    return zeta_over_s;
}


//! Cauchy distribution
double TransportCoeffs::get_temperature_dependent_zeta_over_s_duke(
                                            const double T_in_fm) const {
    const double A = DATA.bulk_2_normalisation;
    const double G = DATA.bulk_2_width_in_GeV;
    const double Tpeak_in_GeV = DATA.bulk_2_peak_in_GeV;
    const double T_in_GeV = T_in_fm*hbarc;
    const double diff_ratio = (T_in_GeV-Tpeak_in_GeV)/G;

    //const double T_delta=(T_in_GeV*T_in_GeV)/(Tpeak_in_GeV*Tpeak_in_GeV)-1;
    //double  bulk_over_sden=A*(G*G)/(T_delta*T_delta+G*G);
    return A/(1+diff_ratio*diff_ratio);
}


//! Skewed Cauchy distribution
double TransportCoeffs::get_temperature_dependent_zeta_over_s_sims(
                                                const double T_in_fm) const {
    const double T_in_GeV = T_in_fm*hbarc;

    const double max = DATA.bulk_3_max;
    const double width = DATA.bulk_3_width_in_GeV;
    const double T_peak_in_GeV = DATA.bulk_3_T_peak_in_GeV;
    const double lambda = DATA.bulk_3_lambda_asymm;
    const double diff = T_in_GeV-T_peak_in_GeV;
    const double sign = (diff > 0) - (diff < 0);
    const double diff_ratio = (diff)/(width*(lambda*sign+1));

    return max/(1+diff_ratio*diff_ratio);
}


//! T dependent bulk viscosity from Gabriel
double TransportCoeffs::get_temperature_dependent_zeta_over_s_default(
                                            const double T_in_fm) const {
    /////////////////////////////////////////////
    //           Parametrization 1             //
    /////////////////////////////////////////////
    // T dependent bulk viscosity from Gabriel
    // used in arXiv: 1502.01675 and 1704.04216
    const double T_in_GeV = T_in_fm*hbarc;
    double Ttr = 0.18;
    double dummy = T_in_GeV/Ttr;
    double A1=-13.77, A2=27.55, A3=13.45;
    double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;

    double bulk = A1*dummy*dummy + A2*dummy - A3;
    if (T_in_GeV < 0.995*Ttr) {
        bulk = (lambda3*exp((dummy-1)/sigma3)
                + lambda4*exp((dummy-1)/sigma4) + 0.03);
    }
    if (T_in_GeV > 1.05*Ttr) {
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


double TransportCoeffs::get_temperature_dependent_zeta_over_s_bigbroadP(
                                                const double T_in_fm) const {
    // used in arXiv: 1901.04378 and 1908.06212
    const double T_in_GeV = T_in_fm*hbarc;
    const double B_norm = 0.24;
    const double B_width = 1.5;
    const double Tpeak = 0.165;
    const double Ttilde = (T_in_GeV/Tpeak - 1.)/B_width;
    double bulk = B_norm/(Ttilde*Ttilde + 1.);
    if (T_in_GeV < Tpeak) {
        const double Tdiff = (T_in_GeV - Tpeak)/0.01;
        bulk = B_norm*exp(-Tdiff*Tdiff);
    }
    return(bulk);
}


double TransportCoeffs::get_temperature_dependent_zeta_over_s_AsymGaussian(
                            const double T_in_fm, const double norm) const {
    const double T_in_GeV = T_in_fm*hbarc;
    const double B_norm = norm;
    const double B_width1 = 0.01;
    const double B_width2 = 0.12;
    const double Tpeak = 0.160;
    double Tdiff = T_in_GeV - Tpeak;
    if (Tdiff > 0.) {
        Tdiff = Tdiff/B_width2;
    } else {
        Tdiff = Tdiff/B_width1;
    }
    double bulk = B_norm*exp(-Tdiff*Tdiff);
    return(bulk);
}
