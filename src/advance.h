// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_ADVANCE_H_
#define SRC_ADVANCE_H_

#include <memory>

#include "cell.h"
#include "data.h"
#include "dissipative.h"
#include "fields.h"
#include "hydro_source_base.h"
#include "minmod.h"
#include "pretty_ostream.h"
#include "reconst.h"
#include "transport_coeffs.h"
#include "u_derivative.h"

class Advance {
  private:
    const InitData &DATA;
    const EOS &eos;
    std::shared_ptr<HydroSourceBase> hydro_source_terms_ptr;

    Diss diss_helper;
    Minmod minmod;
    Reconst reconst_helper;
    pretty_ostream music_message;
    TransportCoeffs transport_coeffs_;

    bool flag_add_hydro_source;

    //! P(e, n_B, n_Q, n_S) of arenaFieldsCurr, filled at the start of each
    //! AdvanceIt so that the stencils do not call the EoS for every neighbour
    std::vector<double> pressureCurr_;
    //! the slice whose pressure is in pressureCurr_ (nullptr: unknown)
    const Fields *pressureCurrOwner_ = nullptr;
    //! in the second RK stage: the pressure of arenaFieldsPrev, which was
    //! the current slice of the first stage
    std::vector<double> pressurePrev_;
    //! mu_B / T of arenaFieldsCurr and arenaFieldsPrev (baryon diffusion),
    //! filled at the start of each AdvanceIt so that the derivatives of
    //! mu_B / T do not call the EoS for every neighbour. In the second RK
    //! stage muBOverTPrev_ is the first stage's muBOverTCurr_
    std::vector<double> muBOverTCurr_;
    std::vector<double> muBOverTPrev_;
    //! the slice whose mu_B / T is in muBOverTCurr_ (nullptr: unknown)
    const Fields *muBOverTCurrOwner_ = nullptr;
    //! the KT flux part of the W update for every cell of the current
    //! stage (Diss::computeWFluxDivergence)
    std::vector<double> wFluxDiv_;
    //! d_m (tau W^{m alpha}) of every cell of the current stage
    //! (Diss::computeWSourceDivergence)
    std::vector<double> wSource_;

    //! KT fluxes through the faces of arenaFieldsCurr, filled at the start
    //! of each AdvanceIt (computeFaceFluxes)
    FaceFluxes faceFluxes_;

    //! vacuumSkipEnergy: per storage index of arenaFieldsCurr, ghost cells
    //! included, 1 if e < vacuumSkipEnergy (vacuum_), and for the cells of
    //! the grid 1 if the cell and the cells up to two away along each axis
    //! are vacuum (frozen_: its update is skipped)
    std::vector<unsigned char> vacuum_;
    std::vector<unsigned char> frozen_;
    //! the cells whose mu_B / T is read: all but the frozen cells whose
    //! neighbours are frozen too (U_derivative reads one cell away)
    std::vector<unsigned char> needMuBOverT_;
    //! cell updates and faces of the run, and how many were skipped
    long long nCellUpdates_ = 0, nCellsFrozen_ = 0;
    long long nFaces_ = 0, nFacesSkipped_ = 0;

    //! fills vacuum_ and frozen_ for arena; returns the frozen cells
    long long markVacuum(const Fields &arena);
    //! the viscous update of a frozen cell (see AdvanceIt)
    void frozenViscousUpdate(
        const double tau, const int rk_flag, const Fields &arenaFieldsPrev,
        const Fields &arenaFieldsCurr, Fields &arenaFieldsNext,
        const int fieldIdx) const;

    //! KT flux through one face (see computeFaceFluxes)
    void computeFaceFlux(
        const double tau, const Fields &arena,
        const std::vector<double> &pressure, const int direction, const int fx,
        const int fy, const int feta, TJbVec &flux);

    // Local fixed point that puts the end-of-step viscous part V_{n+1} into
    // the time derivative of the second RK stage (see
    // stage2ViscousFixedPoint): at most stage2MaxIterations_ updates,
    // stopped when the shift changes by less than stage2Tolerance_ times
    // |tau T^{tau tau}| (|tau J^tau| for the baryon current).
    static constexpr int stage2MaxIterations_ = 16;
    static constexpr double stage2Tolerance_ = 1e-8;

    void stage2ViscousFixedPoint(
        const double tau, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
        const Cell_small &cellCurrFull, Fields &arenaFieldsNext,
        const TJbVec &qiFinal, const double theta_local, const DumuVec &a_local,
        const VelocityShearVec &sigma_local, const VorticityVec &omega_local,
        const DmuMuBoverTVec &baryon_diffusion_vector,
        const std::array<double, 10> &w_rhs, const ThermalVec &thermalVec,
        const int ieta, const int ix, const int iy, const int fieldIdx);

    //! the viscous part of T^{tau alpha} (alpha < 4) and J^tau (alpha = 4)
    //! that enters the time derivative in Diss::MakeWSource
    double viscousTimeComponent(
        const Fields &arena, const int fieldIdx, const int alpha) const;

    //! reconstruct e, n and u from tau T^{tau mu} and store them
    //! pressureGuess: the pressure of guess (see Reconst::ReconstIt_shell)
    void reconstructAndStore(
        const double tau_next, const TJbVec &qi, const ReconstCell &guess,
        Fields &arenaFieldsNext, const int fieldIdx,
        const double pressureGuess = -1.);

  public:
    Advance(
        const EOS &eosIn, const InitData &DATA_in,
        std::shared_ptr<HydroSourceBase> hydro_source_ptr_in);
    //! reports the skipped vacuum work (vacuumSkipEnergy > 0)
    ~Advance();

    //! source terms enter the update (the device path has none yet)
    bool hasHydroSource() const { return (flag_add_hydro_source); }
    const Reconst &getReconst() const { return (reconst_helper); }

    //! currGhostsValid: the ghost cells of arenaFieldsCurr are already up
    //! to date (Evolve keeps them so), else AdvanceIt updates them first.
    //! AdvanceIt always updates the ghost cells of arenaFieldsNext.
    //! currPressureValid: preparePressure(arenaFieldsCurr) was called and
    //! arenaFieldsCurr has not changed since
    void AdvanceIt(
        const double tau_init, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
        Fields &arenaFieldsNext, const int rk_flag,
        const bool currGhostsValid = false,
        const bool currPressureValid = false);

    //! the pressure of every cell of arena (ghost cells included), kept for
    //! the next AdvanceIt on arena (currPressureValid)
    const std::vector<double> &preparePressure(const Fields &arena) {
        computePressure(arena, pressureCurr_);
        pressureCurrOwner_ = &arena;
        return (pressureCurr_);
    }

    void FirstRKStepT(
        const double tau, const double x_local, const double y_local,
        const double eta_s_local, const int ix, const int iy, const int ieta,
        const int rk_flag, const int fieldIdx, Fields &arenaFieldsCurr,
        Fields &arenaFieldsNext, Fields &arenaFieldsPrev,
        const std::vector<double> &pressureCurr, TJbVec &qiFinal,
        const std::vector<double> *pressurePrev = nullptr);

    //! fills pressure[idx] = P(e, n_B, n_Q, n_S) of every cell of arena
    void computePressure(
        const Fields &arena, std::vector<double> &pressure) const;

    //! fills alpha[idx] = mu_B / max(T, small_eps) of every cell of arena,
    //! ghost cells included, as U_derivative evaluates it
    //! need (optional): only the cells with need[idx] != 0
    void computeMuBOverT(
        const Fields &arena, std::vector<double> &alpha,
        const std::vector<unsigned char> *need = nullptr) const;

    void FirstRKStepW(
        const double tau_it, Fields &arenaFieldsPrev, const Cell_small &grid_c,
        Fields &arenaFieldsNext, const int rk_flag, const double theta_local,
        const DumuVec &a_local, const VelocityShearVec &sigma_local,
        const VorticityVec &omega_local,
        const DmuMuBoverTVec &baryon_diffusion_vector,
        const std::array<double, 10> &w_rhs, const ThermalVec &thermalVec,
        const int ieta, const int ix, const int iy, const int fieldIdx);

    void QuestRevertResummedTransCoeff(Cell_small &grid_pt);
    void QuestRevert(
        Cell_small &grid_pt, const int ieta, const int ix, const int iy);
    void QuestRevert_qmu(
        Cell_small &grid_pt, const int ieta, const int ix, const int iy);

    //! fills fluxes with the KT fluxes of arena at tau through every face;
    //! pressure: computePressure of arena
    //! vacuum (optional, vacuum_): faces whose four stencil cells are
    //! vacuum get zero flux
    void computeFaceFluxes(
        const double tau, const Fields &arena,
        const std::vector<double> &pressure, FaceFluxes &fluxes,
        const std::vector<unsigned char> *vacuum = nullptr);

    //! adds the flux divergence (and the Milne geometric terms) of the
    //! cell (ix, iy, ieta) times delta_tau to qi; fluxes: computeFaceFluxes
    void MakeDeltaQI(
        const FaceFluxes &fluxes, const int ix, const int iy, const int ieta,
        TJbVec &qi);

    double MaxSpeed(
        const double tau, const int direc, const ReconstCell &grid_p,
        double &pressure);
    //! the KT flux through one face with the half-way states limited in
    //! the primitive variables (face_reconstruction = 1); cLL, cL, cR, cRR
    //! are the four cells along direction around the face
    void computeFaceFluxPrimitive(
        const double tau, const int direction, const ReconstCell &cLL,
        const ReconstCell &cL, const ReconstCell &cR, const ReconstCell &cRR,
        TJbVec &flux);

    //! MaxSpeed with the EoS at grid_p already known
    double MaxSpeed(
        const double tau, const int direc, const ReconstCell &grid_p,
        const ReconstEos &eosState);

    double get_TJb(
        const ReconstCell &grid_p, const int mu, const int nu,
        const double pressure);

    void solveEigenvaluesWmunu(Cell_small &grid_pt);
    void nCausalityConstraints(Cell_small &grid_pt);
};

#endif  // SRC_ADVANCE_H_
