// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_DISSIPATIVE_H_
#define SRC_DISSIPATIVE_H_

#include <array>
#include "util.h"
#include "cell.h"
#include "data.h"
#include "fields.h"
#include "transport_coeffs.h"
#include "minmod.h"
#include "pretty_ostream.h"

class Diss {
 private:
    const InitData &DATA;
    const EOS &eos;
    const Minmod minmod;
    TransportCoeffs transport_coeffs_;
    int map_2d_idx_to_1d(int a, int b) {
        static const int index_map[5][4] = {{0,   1,  2,  3},
                                            {1,   4,  5,  6},
                                            {2,   5,  7,  8},
                                            {3,   6,  8,  9},
                                            {10, 11, 12, 13}};
        return index_map[a][b];
    }

    pretty_ostream music_message;

 public:
    Diss(const EOS &eosIn, const InitData &DATA_in);
    void MakeWSource(const double tau,
                     const int ix, const int iy, const int ieta, TJbVec &dwmn,
                     Fields &arenaFieldsCurr, Fields &arenaFieldsPrev,
                     const int fieldIdx);

    void Make_uWSource(const double tau, const Cell_small &grid_pt,
                       const Cell_small &grid_pt_prev, const int rk_flag,
                       const double theta_local, const DumuVec &a_local,
                       const VelocityShearVec &sigma_1d,
                       const VorticityVec &omega_1d,
                       std::array<double, 5> &sourceTerms);

    void Make_uWRHS(const double tau, Fields &arena,
                    const int fieldIdx,
                    const int ix, const int iy, const int ieta,
                    std::array<double, 9> &w_rhs,
                    const double theta_local, const DumuVec &a_local);

    double Make_uPiSource(const double tau, const Cell_small &grid_pt,
                          const Cell_small &grid_pt_prev, const int rk_flag,
                          const double theta_local,
                          const VelocityShearVec &sigma_1d);

    void Make_uqSource(const double tau, const Cell_small &grid_pt,
                       const Cell_small &grid_pt_prev,
                       const int rk_flag, const double theta_local,
                       const DumuVec &a_local,
                       const VelocityShearVec &sigma_1d,
                       const VorticityVec &omega_1d,
                       const DmuMuBoverTVec &baryon_diffusion_vec,
                       std::array<double, 3> &sourceTerms);

    void output_kappa_T_and_muB_dependence();
    void output_kappa_along_const_sovernB();
    void output_eta_over_s_T_and_muB_dependence();
    void output_eta_over_s_along_const_sovernB();
    void output_zeta_over_s_T_and_muB_dependence();
    void output_zeta_over_s_along_const_sovernB();
};

#endif  // SRC_DISSIPATIVE_H_
