// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale

#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

#include "util.h"
#include "transport_coeffs.h"

using Util::hbarc;

TransportCoeffs::TransportCoeffs(const InitData &Data_in)
    : DATA(Data_in),
      shear_T_(DATA.T_dependent_shear_to_s),
      shear_muB_(DATA.muB_dependent_shear_to_s),
      bulk_T_(DATA.T_dependent_bulk_to_s) {
    shear_relax_time_factor_ = DATA.shear_relax_time_factor;
    bulk_relax_time_factor_  = DATA.bulk_relax_time_factor;

    if (shear_T_ == 23) {
        read_in_shear_from_file();
    }
    if (bulk_T_ == 23) {
        read_in_bulk_from_file();
    }
}


void TransportCoeffs::read_in_shear_from_file() {
    std::string path = "./EOS/shear_1DGen.bin";
    std::ifstream shearFile;
    shearFile.open(path, std::ios::binary);
    if (!shearFile) {
        std::cout << "Can not find the shear viscosity file: "
                  << path << std::endl;
        exit(1);
    }
    TArr_.resize(100);
    shearArr_.resize(100);

    for (int i = 0; i < 100; i++) {
        shearFile.read((char*)&TArr_[i], sizeof(double));
        shearFile.read((char*)&shearArr_[i], sizeof(double));
    }

    shearFile.close();
}


void TransportCoeffs::read_in_bulk_from_file() {
    std::string path = "./EOS/bulk_1DGen.bin";
    std::ifstream bulkFile;
    bulkFile.open(path, std::ios::binary);
    if (!bulkFile) {
        std::cout << "Can not find the bulk viscosity file: "
                  << path << std::endl;
        exit(1);
    }
    TArr_.resize(100);
    bulkArr_.resize(100);

    for (int i = 0; i < 100; i++) {
        bulkFile.read((char*)&TArr_[i], sizeof(double));
        bulkFile.read((char*)&bulkArr_[i], sizeof(double));
    }

    bulkFile.close();
}


double TransportCoeffs::get_eta_over_s(const double T, const double muB) const {
    // inputs T [1/fm], muB [1/fm]
    // outputs \eta/s
    double eta_over_s = DATA.shear_to_s;

    if (shear_T_ == 0) {
        eta_over_s = DATA.shear_to_s;
    } else if (shear_T_ == 23) {
        double T_in_GeV = T*hbarc;
        int Tidx = static_cast<int>(
            (T_in_GeV - TArr_[0])/(TArr_[1] - TArr_[0]));
        Tidx = std::max(0, std::min(static_cast<int>(TArr_.size() - 2), Tidx));
        double Tfrac = (T_in_GeV - TArr_[Tidx])/(TArr_[1] - TArr_[0]);
        Tfrac = std::max(0., std::min(1., Tfrac));
        eta_over_s = (1. - Tfrac)*shearArr_[Tidx] + Tfrac*shearArr_[Tidx + 1];
    } else if (shear_T_ == 3) {
        eta_over_s = get_temperature_dependent_eta_over_s_sims(T);
    } else if (shear_T_ == 1) {
        eta_over_s = get_temperature_dependent_eta_over_s_default(T);
    } else if (shear_T_ == 2) {
        eta_over_s = get_temperature_dependent_eta_over_s_duke(T);
    } else if (shear_T_ == 11) {
        eta_over_s *= get_temperature_dependence_shear_profile(T);
    } else {
        eta_over_s = DATA.shear_to_s;
    }

    if (shear_muB_ == 7) {
        eta_over_s *= get_muB_dependence_shear_piecewise(muB);
    } else if (shear_muB_ == 10) {
        eta_over_s *= get_muB_dependence_shear_profile(muB);
    }
    eta_over_s = std::max(0., eta_over_s);
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


double TransportCoeffs::get_muB_dependence_shear_piecewise(
                                                const double muB_in_fm) const {
    const double muB_in_GeV = std::abs(muB_in_fm)*hbarc;
    const double f0 = 1.0;
    const double f1 = DATA.shear_muBf0p2;
    const double f2 = DATA.shear_muBf0p4;
    double f_muB = f0;
    if (muB_in_GeV < 0.2) {
        f_muB = f0 + (f1 - f0)/0.2*muB_in_GeV;
    } else if (muB_in_GeV < 0.4) {
        f_muB = f1 + (f2 - f1)/0.2*(muB_in_GeV - 0.2);
    } else {
        f_muB = f2;
    }
    return(f_muB);
}


double TransportCoeffs::get_muB_dependence_shear_profile(
                                                const double muB_in_fm) const {
    const double muB_in_GeV = muB_in_fm*hbarc;
    const double alpha = DATA.shear_muBDep_alpha;
    const double muB_slope = DATA.shear_muBDep_slope;
    const double muB_scale = DATA.shear_muBDep_scale;
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

    const double eta_over_s_min = 1e-3;
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


double TransportCoeffs::get_zeta_over_s(const double T,
                                        const double mu_B) const {
    // input T [1/fm], mu_B [1/fm]
    double muB_in_GeV = mu_B*hbarc;
    double zeta_over_s = 0.;
    if (bulk_T_ == 23) {
        double T_in_GeV = T*hbarc;
        int Tidx = static_cast<int>(
            (T_in_GeV - TArr_[0])/(TArr_[1] - TArr_[0]));
        Tidx = std::max(0, std::min(static_cast<int>(TArr_.size() - 2), Tidx));
        double Tfrac = (T_in_GeV - TArr_[Tidx])/(TArr_[1] - TArr_[0]);
        Tfrac = std::max(0., std::min(1., Tfrac));
        zeta_over_s = (1. - Tfrac)*bulkArr_[Tidx] + Tfrac*bulkArr_[Tidx + 1];
    } else if (bulk_T_ == 10) {
        // param. for 3D-Glauber + MUSIC + UrQMD
        double peak_norm = DATA.bulk_10_max;
        if (muB_in_GeV < 0.2) {
            peak_norm = (DATA.bulk_10_max
                         + (DATA.bulk_10_max_muB0p2 - DATA.bulk_10_max)
                           /0.2*muB_in_GeV);
        } else if (muB_in_GeV < 0.4) {
            peak_norm = (DATA.bulk_10_max_muB0p2
                         + (DATA.bulk_10_max_muB0p4 - DATA.bulk_10_max_muB0p2)
                           /0.2*(muB_in_GeV - 0.2));
        } else {
            peak_norm = DATA.bulk_10_max_muB0p4;
        }
        const double Tpeak = (DATA.bulk_10_Tpeak
                              + DATA.bulk_10_Tpeak_muBcurv
                                *muB_in_GeV*muB_in_GeV);       // GeV
        const double B_width1 = DATA.bulk_10_width_low;        // GeV
        const double B_width2 = DATA.bulk_10_width_high;       // GeV
        zeta_over_s = get_temperature_dependent_zeta_over_s_AsymGaussian(
                            T, peak_norm, B_width1, B_width2, Tpeak);
    } else if (bulk_T_ == 8) {
        // latest param. for IPGlasma + MUSIC + UrQMD
        // Phys.Rev.C 102 (2020) 4, 044905, e-Print: 2005.14682 [nucl-th]
        const double peak_norm = 0.13;
        const double B_width1 = 0.01;
        const double B_width2 = 0.12;
        const double Tpeak = 0.160;
        zeta_over_s = get_temperature_dependent_zeta_over_s_AsymGaussian(
                            T, peak_norm, B_width1, B_width2, Tpeak);
    } else if (bulk_T_ == 1) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_default(T);
    } else if (bulk_T_ == 2) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_duke(T);
    } else if (bulk_T_ == 3) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_sims(T);
    } else if (bulk_T_ == 7) {
        zeta_over_s = get_temperature_dependent_zeta_over_s_bigbroadP(T);
    } else if (bulk_T_ == 9) {
        // latest param. for IPGlasma + KoMPoST + MUSIC + UrQMD
        // Phys.Rev.C 105 (2022) 1, 014909, e-Print: 2106.11216 [nucl-th]
        const double peak_norm = 0.175;
        const double B_width1 = 0.01;
        const double B_width2 = 0.12;
        const double Tpeak = 0.160;
        zeta_over_s = get_temperature_dependent_zeta_over_s_AsymGaussian(
                            T, peak_norm, B_width1, B_width2, Tpeak);
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
        const double T_in_fm, const double B_norm, const double B_width1,
        const double B_width2, const double Tpeak) const {
    const double T_in_GeV = T_in_fm*hbarc;
    double Tdiff = T_in_GeV - Tpeak;        // GeV
    if (Tdiff > 0.) {
        Tdiff = Tdiff/B_width2;
    } else {
        Tdiff = Tdiff/B_width1;
    }
    double bulk = B_norm*exp(-Tdiff*Tdiff);
    return(bulk);
}
