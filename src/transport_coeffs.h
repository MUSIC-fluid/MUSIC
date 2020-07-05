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

    double get_shear_relax_time_factor() const {
        return(shear_relax_time_factor_);
    }

    double get_bulk_relax_time_factor() const {
        return(bulk_relax_time_factor_);
    }
};

#endif  // SRC_TRANSPORT_H_
