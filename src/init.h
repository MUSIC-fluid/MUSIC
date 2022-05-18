// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_INIT_H_
#define SRC_INIT_H_

#include <stdio.h>
#include <vector>
#include <cmath>
#include <memory>

#include "data.h"
#include "cell.h"
#include "grid.h"
#include "eos.h"
#include "hydro_source_base.h"
#include "pretty_ostream.h"

class Init {
 private:
    InitData &DATA;
    const EOS &eos;
    std::weak_ptr<HydroSourceBase> hydro_source_terms_ptr;
    pretty_ostream music_message;

    // support for JETSCAPE
    std::vector<double> jetscape_initial_energy_density;
    std::vector<double> jetscape_initial_pressure;
    std::vector<double> jetscape_initial_u_tau;
    std::vector<double> jetscape_initial_u_x;
    std::vector<double> jetscape_initial_u_y;
    std::vector<double> jetscape_initial_u_eta;
    std::vector<double> jetscape_initial_pi_00;
    std::vector<double> jetscape_initial_pi_01;
    std::vector<double> jetscape_initial_pi_02;
    std::vector<double> jetscape_initial_pi_03;
    std::vector<double> jetscape_initial_pi_11;
    std::vector<double> jetscape_initial_pi_12;
    std::vector<double> jetscape_initial_pi_13;
    std::vector<double> jetscape_initial_pi_22;
    std::vector<double> jetscape_initial_pi_23;
    std::vector<double> jetscape_initial_pi_33;
    std::vector<double> jetscape_initial_bulk_pi;

 public:
    Init(const EOS &eos, InitData &DATA_in,
         std::shared_ptr<HydroSourceBase> hydro_source_ptr_in);

    void InitArena(SCGrid &arena_prev, SCGrid &arena_current,
                   SCGrid &arena_future);
    void InitTJb  (SCGrid &arena_prev, SCGrid &arena_current);
    void print_num_of_threads();

    void initial_Gubser_XY               (int ieta, SCGrid &arena_prev, SCGrid &arena_current);
    void initial_1p1D_eta                (SCGrid &arena_prev, SCGrid &arena_current);
    void initial_IPGlasma_XY             (int ieta, SCGrid &arena_prev, SCGrid &arena_current);
    void initial_IPGlasma_XY_with_pi     (int ieta, SCGrid &arena_prev, SCGrid &arena_current);
    void initial_with_zero_XY            (int ieta, SCGrid &arena_prev, SCGrid &arena_current);
    void initial_AMPT_XY                 (int ieta, SCGrid &arena_prev, SCGrid &arena_current);
    void initial_MCGlb_with_rhob         (SCGrid &arena_prev, SCGrid &arena_current);
    void initial_UMN_with_rhob           (SCGrid &arena_prev, SCGrid & arena_current);
    void initial_with_jetscape           (int ieta, SCGrid &arena_prev, SCGrid &arena_current);

    void get_jetscape_preequilibrium_vectors(
        std::vector<double> e_in, std::vector<double> P_in,
        std::vector<double> u_tau_in, std::vector<double> u_x_in,
        std::vector<double> u_y_in,   std::vector<double> u_eta_in,
        std::vector<double> pi_00_in, std::vector<double> pi_01_in,
        std::vector<double> pi_02_in, std::vector<double> pi_03_in,
        std::vector<double> pi_11_in, std::vector<double> pi_12_in,
        std::vector<double> pi_13_in, std::vector<double> pi_22_in,
        std::vector<double> pi_23_in, std::vector<double> pi_33_in,
        std::vector<double> Bulk_pi_in);
    void clean_up_jetscape_arrays();

    double eta_profile_plateau(const double eta, const double eta_0,
                               const double sigma_eta) const;
    double energy_eta_profile_normalisation(
        const double y_CM, const double eta_0, const double sigma_eta) const;
    double Pz_eta_profile_normalisation(
        const double eta_0, const double sigma_eta) const;
    double eta_rhob_profile_normalisation  (const double eta) const;
    double eta_profile_left_factor         (const double eta) const;
    double eta_profile_right_factor        (const double eta) const;
    double eta_rhob_left_factor            (const double eta) const;
    double eta_rhob_right_factor           (const double eta) const;
    void   output_initial_density_profiles (SCGrid &arena);
};

#endif  // SRC_INIT_H_
