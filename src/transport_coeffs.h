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

 public:
    TransportCoeffs(const EOS &eosIn, const InitData &DATA_in);

    double get_eta_over_s(double T) const;
    double get_zeta_over_s(double T) const;

    double get_temperature_dependent_eta_over_s_default(double T) const;
    double get_temperature_dependent_zeta_over_s_default(double T) const;

    double get_temperature_dependent_eta_over_s_duke(double T) const;
    double get_temperature_dependent_zeta_over_s_duke(double T) const;

    double get_temperature_dependent_eta_over_s_sims(double T) const;
    double get_temperature_dependent_zeta_over_s_sims(double T) const;

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
    double get_lambda_piPi_coeff() const { return(6./5.); }

    // bulk
    double get_lambda_Pipi_coeff() const { return(8./5.); }
    double get_delta_PiPi_coeff() const { return(2./3.); }
    double get_tau_PiPi_coeff() const { return(0.); }  // unknown


    // net baryon diffusion

    // from conformal kinetic theory
    double get_delta_qq_coeff() const { return(1.0); }

    // from 14-momentum massless
    double get_lambda_qq_coeff() const { return(3./5.); }

    double get_l_qpi_coeff() const { return(0.); }
    double get_lambda_qpi_coeff() const { return(0.); }

};

#endif  // SRC_TRANSPORT_H_
