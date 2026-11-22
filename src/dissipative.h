// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_DISSIPATIVE_H_
#define SRC_DISSIPATIVE_H_

#include <array>
#include <vector>

#include "cell.h"
#include "data.h"
#include "eos.h"
#include "fields.h"
#include "minmod.h"
#include "pretty_ostream.h"
#include "transport_coeffs.h"
#include "util.h"

class Diss {
  private:
    const InitData &DATA;
    const EOS &eos;
    const Minmod minmod;
    TransportCoeffs transport_coeffs_;
    int map_2d_idx_to_1d(int a, int b) {
        static const int index_map[5][4] = {
            {0, 1, 2, 3},
            {1, 4, 5, 6},
            {2, 5, 7, 8},
            {3, 6, 8, 9},
            {10, 11, 12, 13}};
        return index_map[a][b];
    }

    const int echo_level_ = DATA.echo_level;
    pretty_ostream music_message;

    //! the face fluxes of one direction (computeWFluxDivergence)
    std::vector<double> wFaceFluxes_;
    //! the face averages of tau W^{i alpha} and of the bulk part of one
    //! direction (computeWSourceDivergence)
    std::vector<double> wSourceFaces_;
    //! cosh(deta/2)/deta and sinh(deta/2)/deta of MakeWSource
    double wSourceCoshDeta_ = 0.;
    double wSourceSinhDeta_ = 0.5;

    //! the number of W components that Make_uWRHS evolves
    int nWComponents() const;

  public:
    Diss(const EOS &eosIn, const InitData &DATA_in);
    //! d_m (tau W^{m alpha}) + geometric terms of one cell, with the
    //! neighbouring cells read for every cell (two reads of each face).
    //! computeWSourceDivergence + MakeWSourceFromDivergence give the same
    //! result bit for bit
    void MakeWSource(
        const double tau, const int ix, const int iy, const int ieta,
        TJbVec &dwmn, Fields &arenaFieldsCurr, Fields &arenaFieldsPrev,
        const int fieldIdx);

    //! values per cell of computeWSourceDivergence: d_m (tau W^{m alpha})
    //! for alpha < 5 without the eta faces of alpha = 0, 3, then
    //! tau W^{eta alpha} (+ bulk) at eta +- deta/2 for alpha = 0, 3
    static constexpr int kWSourceComponents = 9;

    //! the part of MakeWSource that reads the neighbouring cells, for every
    //! cell of arenaCurr, with each face average computed once. Component
    //! k of cell fieldIdx is at wSource[k * arenaCurr.nStorage() + fieldIdx]
    void computeWSourceDivergence(
        const double tau, const Fields &arenaCurr, const Fields &arenaPrev,
        std::vector<double> &wSource);

    //! MakeWSource of cell fieldIdx from computeWSourceDivergence
    void MakeWSourceFromDivergence(
        const std::vector<double> &wSource, const int nStorage,
        const int fieldIdx, TJbVec &dwmn) const;

    void Make_uWSource(
        const double tau, const Cell_small &grid_pt, const double theta_local,
        const DumuVec &a_local, const VelocityShearVec &sigma_1d,
        const VorticityVec &omega_1d, const ThermalVec &thermalVec,
        std::array<double, 5> &sourceTerms);

    //! components of the W flux divergence, see computeWFluxDivergence
    static constexpr int kWComponents = 10;

    //! the KT flux part of Make_uWRHS for every cell of arena:
    //! -delta_tau (H_{i+1/2} - H_{i-1/2}) / delta_i of u^i W, summed over
    //! the directions, with each face flux computed once. Component idx of
    //! cell fieldIdx is at wFluxDiv[idx * arena.nStorage() + fieldIdx].
    void computeWFluxDivergence(
        const double tau, const Fields &arena, std::vector<double> &wFluxDiv);

    //! w_rhs: the spatial part of the W update; wFluxDiv from
    //! computeWFluxDivergence on the same arena and tau
    void Make_uWRHS(
        const double tau, Fields &arena, const int fieldIdx,
        const std::vector<double> &wFluxDiv, std::array<double, 10> &w_rhs,
        const double theta_local, const DumuVec &a_local);

    double Make_uPiSource(
        const double tau, const Cell_small &grid_pt, const double theta_local,
        const VelocityShearVec &sigma_1d, const ThermalVec &thermalVec);

    double Make_uPiChemSource(
        const double tau, const Cell_small &grid_pt, const double theta_local,
        const VelocityShearVec &sigma_1d, const ThermalVec &thermalVec);

    void Make_uqSource(
        const double tau, const Cell_small &grid_pt, const double theta_local,
        const DumuVec &a_local, const VelocityShearVec &sigma_1d,
        const VorticityVec &omega_1d,
        const DmuMuBoverTVec &baryon_diffusion_vec,
        const ThermalVec &thermalVec, std::array<double, 3> &sourceTerms);

    double D_muB_over_T_ideal(
        const double epsilon, const double rhob, const double pressure,
        const double theta) const;

    void computeInverseReynoldsNumbers(
        const double enthalpy, const Cell_small &grid_pt, double &R_shear,
        double &R_bulk) const;

    void output_kappa_T_and_muB_dependence();
    void output_kappa_along_const_sovernB();
    void output_eta_over_s_T_and_muB_dependence();
    void output_eta_over_s_along_const_sovernB();
    void output_zeta_over_s_T_and_muB_dependence();
    void output_zeta_over_s_along_const_sovernB();
};

#endif  // SRC_DISSIPATIVE_H_
