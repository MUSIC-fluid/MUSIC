// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_EOS_H_
#define SRC_EOS_H_

#include <iostream>
#include <vector>
#include <string>

#include "util.h"
#include "data.h"
#include "pretty_ostream.h"

class EOS {
 private:
    const InitData& parameters_ptr;
    pretty_ostream music_message;

    int number_of_tables;
    std::vector<double> nb_bounds;
    std::vector<double> e_bounds;

    std::vector<double> nb_spacing;
    std::vector<double> e_spacing;

    std::vector<int> nb_length;
    std::vector<int> e_length;

    double ***pressure_tb;
    double ***temperature_tb;
    double ***mu_B_tb;
    double ***mu_S_tb;

    int whichEOS;
    double eps_max;

 public:
    EOS() = default;
    EOS(const InitData &para_in);  // constructor
    ~EOS();  // destructor
    void initialize_eos();
    std::string get_hydro_env_path() const;
    void resize_table_info_arrays();
    void init_eos();                 // for whichEOS=1

    void init_eos_s95p(int selector);
    void init_eos10();   // for EOS at finite mu_B from A. M.
    void init_eos11();   // foe EoS at finite mu_B from Pasi
    void init_eos12();   // for EOS at finite mu_B from A. M.

    // returns maximum local energy density of the EoS table
    // in the unit of [1/fm^4]
    double get_eps_max() const { return(eps_max); }

    // for EOS at finite mu_B
    int get_table_idx(double e) const;
    double interpolate1D(double e, int table_idx, double ***table) const;
    double interpolate2D(double e, double rhob, int table_idx, double ***table) const;

    double get_cs2(double e, double rhob) const;
    double calculate_velocity_of_sound_sq(double e, double rhob) const;
    double get_rhob_from_mub   (double e, double mub) const;
    double get_dpOverde_WB     (double e) const;
    double get_dpOverde3       (double e, double rhob) const;
    double get_dpOverdrhob2    (double e, double rhob) const;
    double p_rho_func          (double e, double rhob) const;
    double p_e_func            (double e, double rhob) const;
    double T_from_eps_ideal_gas(double eps) const;
    double get_entropy         (double epsilon, double rhob) const;
    double get_temperature     (double epsilon, double rhob) const;
    double get_temperature_WB  (double e_local) const;
    double get_mu              (double epsilon, double rhob) const;
    double get_muS             (double epsilon, double rhob) const;
    double get_pressure        (double epsilon, double rhob) const;
    double get_pressure_WB     (double e_local) const;
    double s2e_ideal_gas(double s) const;
    double get_s2e(double s, double rhob) const;
    double get_s2e_finite_rhob(double s, double rhob) const;
    void check_eos() const;
    void check_eos_with_finite_muB() const;
    void check_eos_no_muB() const;
    void output_eos_matrix(int ne, int nrhob, double** matrix_ptr,
                           std::string filename) const;
};

#endif  // SRC_EOS_H_
