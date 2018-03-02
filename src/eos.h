// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_EOS_H_
#define SRC_EOS_H_

#include <iostream>

#include "util.h"
#include "data.h"
#include "pretty_ostream.h"

class EOS {
 private:
    const InitData& parameters_ptr;
    pretty_ostream music_message;

    double BNP1, EPP1;            // start value for \mu_B and epsilon
    double BNP2, EPP2;            // start value for \mu_B and epsilon
    double BNP3, EPP3;            // start value for \mu_B and epsilon
    double BNP4, EPP4;            // start value for \mu_B and epsilon
    double BNP5, EPP5;            // start value for \mu_B and epsilon
    double BNP6, EPP6;            // start value for \mu_B and epsilon
    double BNP7, EPP7;            // start value for \mu_B and epsilon

    double deltaBNP1, deltaEPP1;  // step size for \mu_B and epsilon
    double deltaBNP2, deltaEPP2;  // step size for \mu_B and epsilon
    double deltaBNP3, deltaEPP3;  // step size for \mu_B and epsilon
    double deltaBNP4, deltaEPP4;  // step size for \mu_B and epsilon
    double deltaBNP5, deltaEPP5;  // step size for \mu_B and epsilon
    double deltaBNP6, deltaEPP6;  // step size for \mu_B and epsilon
    double deltaBNP7, deltaEPP7;  // step size for \mu_B and epsilon

    int NBNP1, NEPP1;             // number of entries for \mu_B and epsilon
    int NBNP2, NEPP2;             // number of entries for \mu_B and epsilon
    int NBNP3, NEPP3;             // number of entries for \mu_B and epsilon
    int NBNP4, NEPP4;             // number of entries for \mu_B and epsilon
    int NBNP5, NEPP5;             // number of entries for \mu_B and epsilon
    int NBNP6, NEPP6;             // number of entries for \mu_B and epsilon
    int NBNP7, NEPP7;             // number of entries for \mu_B and epsilon

    double **pressure1;
    double **pressure2;
    double **pressure3;
    double **pressure4;
    double **pressure5;
    double **pressure6;
    double **pressure7;
    double **temperature1;
    double **temperature2;
    double **temperature3;
    double **temperature4;
    double **temperature5;
    double **temperature6;
    double **temperature7;
    double **QGPfraction1;
    double **QGPfraction2;
    double **QGPfraction3;
    double **QGPfraction4;
    double **QGPfraction5;
    double **QGPfraction6;
    double **QGPfraction7;
    double **entropyDensity1;
    double **entropyDensity2;
    double **entropyDensity3;
    double **entropyDensity4;
    double **entropyDensity5;
    double **entropyDensity6;
    double **entropyDensity7;
    double **mu1;
    double **mu2;
    double **mu3;
    double **mu4;
    double **mu5;
    double **mu6;
    double **mu7;
    double **mus1;
    double **mus2;
    double **mus3;
    double **mus4;
    double **mus5;
    double **mus6;
    double **mus7;
    double **cs2_1;
    double **cs2_2;
    double **cs2_3;
    double **cs2_4;
    double **cs2_5;
    double **cs2_6;
    double **cs2_7;

    double * eps_list_rho0, * s_list_rho0;
    int s_list_rho0_length;
    int whichEOS;
    double eps_max;


 public:
    EOS() = default;
    EOS(const InitData &para_in);  // constructor
    ~EOS();  // destructor
    void initialize_eos();
    void init_eos0();                // for whichEOS=0
    void init_eos();                 // for whichEOS=1
    void init_eos2();                // for whichEOS=2
    void init_eos3(int selector);    // for whichEOS=3 (PCE 150 MeV),
                                     // whichEOS=4 (PCE 155 MeV),
                                     // whichEOS=5 (PCE 160 MeV),
                                     // whichEOS=6 (PCE 165 MeV)
    void init_eos7();                // for whichEOS=7 s95p-v1.2 (for UrQMD)
    void init_eos10(int selector);   // for EOS at finite mu_B from A. M.
    void init_eos11(int selector);   // foe EoS at finite mu_B from Pasi
    void init_eos12(int selector);   // for EOS at finite mu_B from A. M.

    // returns maximum local energy density of the EoS table
    // in the unit of [1/fm^4]
    double get_eps_max() const { return(eps_max); }

    void checkForReadError(FILE *file, const char* name) const;
    double interpolate_pressure(double e, double rhob) const;  // for whichEOS == 1
    double interpolate(double e, double rhob, int selector) const;

    // for whichEOS == 2
    double interpolate2(double e, double rhob, int selector) const;

    // for EOS at finite mu_B
    double interpolate2D(double e, double rhob, int selector) const;

    double get_cs2(double e, double rhob) const;
    double calculate_velocity_of_sound_sq(double e, double rhob) const;
    void fill_cs2_matrix(double e0, double de, int ne,
                         double rhob0, double drhob, int nrhob,
                         double** cs2_ptr);
    void build_velocity_of_sound_sq_matrix();
    double get_rhob_from_mub   (double e, double mub) const;
    double get_dpOverde        (double e, double rhob) const;
    double get_dpOverde2       (double e, double rhob) const;
    double get_dpOverde_WB     (double e) const;
    double get_dpOverde3       (double e, double rhob) const;
    double get_dpOverdrhob     (double e, double rhob) const;
    double get_dpOverdrhob2    (double e, double rhob) const;
#pragma omp declare simd
    double p_rho_func          (double e, double rhob) const;
#pragma omp declare simd
    double p_e_func            (double e, double rhob) const;
    double T_from_eps_ideal_gas(double eps) const;
    double get_entropy         (double epsilon, double rhob) const;
    double get_temperature     (double epsilon, double rhob) const;
    double get_temperature_WB  (double e_local) const;
    double get_mu              (double epsilon, double rhob) const;
    double get_muS             (double epsilon, double rhob) const;
#pragma omp declare simd
    double get_pressure        (double epsilon, double rhob) const;
    double get_pressure_WB     (double e_local) const;
    double ssolve              (double e, double rhob, double s) const;
    double Tsolve              (double e, double rhob, double T);
    double findRoot(double (EOS::*function)(double, double, double),
                    double rhob, double s, double e1, double e2, double eacc);
    double s2e_ideal_gas(double s) const;
    double get_s2e(double s, double rhob) const;
    double get_s2e_finite_rhob(double s, double rhob) const;
    void check_eos() const;
    void check_eos_with_finite_muB() const;
    void check_eos_no_muB() const;
    void output_eos_matrix(int ne, int nrhob, double** matrix_ptr,
                           string filename) const;
};

#endif  // SRC_EOS_H_
