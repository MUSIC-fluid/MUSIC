// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_TRANSPORT_H_
#define SRC_TRANSPORT_H_

#include "data.h"
#include "eos.h"

class TransportCoeffs {
 private:
    const InitData &DATA;
    const EOS &eos;
    double shear_relax_time_factor_;
    double bulk_relax_time_factor_;

    const int shear_T_;
    const int shear_muB_;
    const int bulk_T_;

 public:
    TransportCoeffs(const EOS &eosIn, const InitData &DATA_in);

    double get_eta_over_s(const double T, const double muB) const;
    double get_zeta_over_s(const double T, const double muB) const;

    double get_temperature_dependent_eta_over_s_default(const double T) const;
    double get_temperature_dependent_zeta_over_s_default(const double T) const;

    double get_temperature_dependent_eta_over_s_duke(const double T) const;
    double get_temperature_dependent_zeta_over_s_duke(const double T) const;

    double get_temperature_dependent_eta_over_s_sims(const double T) const;
    double get_temperature_dependent_zeta_over_s_sims(const double T) const;

    double get_muB_dependence_shear_profile(const double muB) const;
    double get_muB_dependence_shear_piecewise(const double muB) const;
    double get_temperature_dependence_shear_profile(const double T) const;

    double get_temperature_dependent_zeta_over_s_bigbroadP(
                                                    const double T) const;
    double get_temperature_dependent_zeta_over_s_AsymGaussian(
        const double T_in_fm, const double B_norm, const double B_width1,
        const double B_width2, const double Tpeak) const;

    // Second-order coefficients
    double get_shear_relax_time_factor() const {
        return(shear_relax_time_factor_);
    }

    double get_bulk_relax_time_factor() const {
        return(bulk_relax_time_factor_);
    }

    // shear
    double get_tau_pipi_coeff() const { return(10./7.); }
    double get_delta_pipi_coeff() const { return(4./3.); }
    double get_phi7_coeff() const { return(9./70.); }
    double get_lambda_pibulkPi_coeff() const { return(6./5.); }

    // bulk
    double get_lambda_bulkPipi_coeff() const { return(8./5.); }
    double get_delta_bulkPibulkPi_coeff() const { return(2./3.); }
    double get_tau_bulkPibulkPi_coeff() const { return(0.); }  // unknown


    // net baryon diffusion

    // from conformal kinetic theory
    double get_delta_qq_coeff() const { return(1.0); }

    // from 14-momentum massless
    double get_lambda_qq_coeff() const { return(3./5.); }

    double get_l_qpi_coeff() const { return(0.); }
    double get_lambda_qpi_coeff() const { return(0.); }

};

#endif  // SRC_TRANSPORT_H_
