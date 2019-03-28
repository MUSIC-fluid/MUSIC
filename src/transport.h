// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_TRANSPORT_H_
#define SRC_TRANSPORT_H_

#include "util.h"
#include "data.h"
#include "eos.h"
#include "pretty_ostream.h"

class Transport {
 private:
    const InitData &DATA;
    const EOS &eos;

 public:
    Transport(const EOS &eosIn, const InitData &DATA_in);

    double get_eta_over_s(double T);
    double get_zeta_over_s(double T);

    double get_temperature_dependent_eta_over_s_default(double T);
    double get_temperature_dependent_zeta_over_s_default(double T);

    double get_temperature_dependent_eta_over_s_duke(double T);
    double get_temperature_dependent_zeta_over_s_duke(double T);

    double get_temperature_dependent_eta_over_s_sims(double T);
    double get_temperature_dependent_zeta_over_s_sims(double T);
};

#endif  // SRC_TRANSPORT_H_
