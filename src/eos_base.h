// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_BASE_H_
#define SRC_EOS_BASE_H_

#include <string>
#include <vector>

#include "data_struct.h"
#include "eos_table.h"
#include "pretty_ostream.h"

class EOS_base {
  private:
    // initialized: a default-constructed EOS_base must destruct safely
    int whichEOS = 0;
    int number_of_tables = 0;
    double eps_max = 0.;
    bool flag_muB = false;
    bool flag_muS = false;
    bool flag_muQ = false;

  public:
    pretty_ostream music_message;
    std::vector<double> nb_bounds;
    std::vector<double> e_bounds;

    std::vector<double> nb_spacing;
    std::vector<double> e_spacing;

    std::vector<int> nb_length;
    std::vector<int> e_length;

    EosTable pressure_tb;
    EosTable temperature_tb;
    EosTable mu_B_tb;
    EosTable mu_S_tb;
    EosTable mu_Q_tb;

    EOS_base() = default;
    virtual ~EOS_base();

    std::string get_hydro_env_path() const;

    void set_number_of_tables(int ntables) { number_of_tables = ntables; }
    int get_number_of_tables() const { return (number_of_tables); }
    void resize_table_info_arrays();
    bool check_table_seams() const;

    void set_EOS_id(int eos_id) { whichEOS = eos_id; }
    int get_EOS_id() const { return (whichEOS); }

    void set_flag_muB(bool flag_muB_in) { flag_muB = flag_muB_in; }
    bool get_flag_muB() const { return (flag_muB); }
    void set_flag_muS(bool flag_muS_in) { flag_muS = flag_muS_in; }
    bool get_flag_muS() const { return (flag_muS); }
    void set_flag_muQ(bool flag_muQ_in) { flag_muQ = flag_muQ_in; }
    bool get_flag_muQ() const { return (flag_muQ); }

    // returns maximum local energy density of the EoS table
    // in the unit of [1/fm^4]
    void set_eps_max(double eps_max_in) { eps_max = eps_max_in; }
    double get_eps_max() const { return (eps_max); }

    double interpolate1D(double e, int table_idx, const EosTable &table) const;
    void interpolate1D_with_gradients(
        double e, int table_idx, const EosTable &table, double &p,
        double &dpde) const;
    double interpolate2D(
        const double e, const double rhob, const int table_idx,
        const EosTable &table) const;
    void interpolate2D_with_gradients(
        const double e, const double rhob, const int table_idx,
        const EosTable &table, double &p, double &dpde, double &dpdrhob) const;

    int get_table_idx(double e) const;

    // n_B [1/fm^3] at which each e column (nb_edge[table][ie]) of a
    // (e, n_B) table stops: beyond it mu_B is capped at the table maximum
    // or the table runs out. Empty for EoS without n_B tables.
    std::vector<std::vector<double>> nb_edge;
    void find_baryon_table_edges();
    double get_rhob_edge(double e) const;
    //! n_B at which the tabulated state is evaluated, with the sign of rhob
    virtual double get_rhob_in_table(double e, double rhob) const;

    //! true if the EoS takes n_Q and n_S as independent inputs (4D EoS);
    //! otherwise n_Q and n_S follow from (e, n_B), e.g. n_Q = 0.4 n_B
    virtual bool has_QS_densities() const { return (false); }
    //! n_B, n_Q and n_S of the state the EoS evaluates at (e, n_B, n_Q, n_S)
    void get_densities_in_table(
        double e, double rhob, double rhoq, double rhos, double &nB, double &nQ,
        double &nS) const;

    virtual double get_entropy(
        double epsilon, double rhob, double rhoq = 0.0,
        double rhos = 0.0) const;
    virtual void getThermalVariables(
        double epsilon, double rhob, double rhoq, double rhos,
        ThermalVec &thermalVec) const;

    double calculate_velocity_of_sound_sq(
        double e, double rhob, double rhoq = 0.0, double rhos = 0.0) const;
    double get_dpOverde3(
        double e, double rhob, double rhoq = 0.0, double rhos = 0.0) const;
    double get_dpOverdrhob2(
        double e, double rhob, double rhoq = 0.0, double rhos = 0.0) const;
    double get_s2e_finite_rhob(
        double s, double rhob, double rhoq = 0.0, double rhos = 0.0) const;
    double get_T2e_finite_rhob(
        const double T, const double rhob, double rhoq = 0.0,
        double rhos = 0.0) const;
    void map_TmuB2erhoB(
        const double T, const double muB, double &e, double &rhob) const;

    virtual void initialize_eos() {}
    virtual void initialize_eos(int eos_id_in) {}
    virtual double get_cs2(
        double e, double rhob, double rhoq = 0., double rhos = 0.) const;
    virtual double p_rho_func(
        double e, double rhob, double rhoq = 0., double rhos = 0.) const {
        return (0.0);
    }
    virtual double p_e_func(
        double e, double rhob, double rhoq = 0., double rhos = 0.) const {
        return (0.0);
    }
    virtual double get_temperature(
        double epsilon, double rhob, double rhoq = 0.0,
        double rhos = 0.0) const {
        return (0.0);
    }
    virtual double get_muB(
        double epsilon, double rhob, double rhoq = 0.0,
        double rhos = 0.0) const {
        return (0.0);
    }
    virtual double get_muS(
        double epsilon, double rhob, double rhoq = 0.0,
        double rhos = 0.0) const {
        return (0.0);
    }
    virtual double get_muQ(
        double epsilon, double rhob, double rhoq = 0.0,
        double rhos = 0.0) const {
        return (0.0);
    }
    virtual double get_rhoS(
        double epsilon, double rhob, double rhoq = 0.0,
        double rhos = 0.0) const {
        return (0.0);
    }
    virtual double get_rhoQ(
        double epsilon, double rhob, double rhoq = 0.0,
        double rhos = 0.0) const {
        return (0.4 * rhob);
    }
    virtual double get_pressure(
        double epsilon, double rhob, double rhoq = 0.0,
        double rhos = 0.0) const {
        return (0.0);
    }
    virtual double get_s2e(
        double s, double rhob, double rhoq = 0.0, double rhos = 0.0) const {
        return (0.0);
    }
    virtual double get_T2e(
        double T_in_GeV, double rhob, double rhoq = 0.0,
        double rhos = 0.0) const {
        return (0.0);
    }
    virtual void check_eos() const {}

    virtual void get_pressure_with_gradients(
        double epsilon, double rhob, double rhoq, double rhos, double &p,
        double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos) const {
        p = get_pressure(epsilon, rhob);
        dpde = get_dpOverde3(epsilon, rhob);
        dpdrhob = get_dpOverdrhob2(epsilon, rhob);
        dpdrhoq = 0;
        dpdrhos = 0;
    }

    virtual void get_pressure_with_gradients_and_cs2(
        double epsilon, double rhob, double rhoq, double rhos, double &p,
        double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos,
        double &cs2) const {
        get_pressure_with_gradients(
            epsilon, rhob, rhoq, rhos, p, dpde, dpdrhob, dpdrhoq, dpdrhos);
        cs2 = dpde
              + (dpdrhob * rhob + dpdrhos * rhos + dpdrhoq * rhoq)
                    / (epsilon + p);
        cs2 = std::max(0.01, std::min(1. / 3., cs2));
    }

    void check_eos_with_finite_muB() const;
    void check_4D_eos() const;
    void check_eos_no_muB() const;
    void outputMutable() const;
};

#endif  // SRC_EOS_BASE_H_
