// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_1DGenerator_H_
#define SRC_EOS_1DGenerator_H_

#include "eos_base.h"

class EOS_1DGenerator : public EOS_base {
  private:
    const int eos_id;

  public:
    EOS_1DGenerator(const int eos_id_in);

    void initialize_eos();
    double p_e_func(
        double e, double rhob, double rhoq = 0., double rhos = 0.) const;
    double get_temperature(
        double e, double rhob, double rhoq = 0., double rhos = 0.) const;
    double get_pressure(
        double e, double rhob, double rhoq = 0., double rhos = 0.) const;
    double get_s2e(
        double s, double rhob, double rhoq = 0., double rhos = 0.) const;
    double get_T2e(
        double T, double rhob, double rhoq = 0., double rhos = 0.) const;
    void get_pressure_with_gradients(
        double e, double rhob, double rhoq, double rhos, double &p,
        double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos) const;
    // c_s^2 = dP/de, clamped to [0.01, 0.5] in get_pressure_with_gradients
    void get_pressure_with_gradients_and_cs2(
        double e, double rhob, double rhoq, double rhos, double &p,
        double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos,
        double &cs2) const;
    void interpolate1D_with_weighted_gradient(
        double e, int table_idx, const EosTable &table, double &p,
        double &dpde) const;

    void check_eos() const { check_eos_no_muB(); }
};

#endif  // SRC_EOS_hotQCD_H_
