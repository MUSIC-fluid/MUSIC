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

    // Local fixed point that puts the end-of-step viscous part V_{n+1} into
    // the time derivative of the second RK stage (see
    // stage2ViscousFixedPoint): at most stage2MaxIterations_ updates,
    // stopped when the shift changes by less than stage2Tolerance_ times
    // |tau T^{tau tau}| (|tau J^tau| for the baryon current).
    static constexpr int stage2MaxIterations_ = 16;
    static constexpr double stage2Tolerance_ = 1e-8;

    void stage2ViscousFixedPoint(
        const double tau, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
        Fields &arenaFieldsNext, const TJbVec &qiFinal,
        const double theta_local, const DumuVec &a_local,
        const VelocityShearVec &sigma_local, const VorticityVec &omega_local,
        const DmuMuBoverTVec &baryon_diffusion_vector,
        const std::array<double, 10> &w_rhs, const int ieta, const int ix,
        const int iy, const int fieldIdx);

    //! the viscous part of T^{tau alpha} (alpha < 4) and J^tau (alpha = 4)
    //! that enters the time derivative in Diss::MakeWSource
    double viscousTimeComponent(
        const Fields &arena, const int fieldIdx, const int alpha) const;

    //! reconstruct e, n and u from tau T^{tau mu} and store them
    void reconstructAndStore(
        const double tau_next, const TJbVec &qi, const ReconstCell &guess,
        Fields &arenaFieldsNext, const int fieldIdx);

  public:
    Advance(
        const EOS &eosIn, const InitData &DATA_in,
        std::shared_ptr<HydroSourceBase> hydro_source_ptr_in);

    void AdvanceIt(
        const double tau_init, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
        Fields &arenaFieldsNext, const int rk_flag);

    void FirstRKStepT(
        const double tau, const double x_local, const double y_local,
        const double eta_s_local, const int ix, const int iy, const int ieta,
        const int rk_flag, const int fieldIdx, Fields &arenaFieldsCurr,
        Fields &arenaFieldsNext, Fields &arenaFieldsPrev, TJbVec &qiFinal);

    void FirstRKStepW(
        const double tau_it, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
        Fields &arenaFieldsNext, const int rk_flag, const double theta_local,
        const DumuVec &a_local, const VelocityShearVec &sigma_local,
        const VorticityVec &omega_local,
        const DmuMuBoverTVec &baryon_diffusion_vector,
        const std::array<double, 10> &w_rhs, const int ieta, const int ix,
        const int iy, const int fieldIdx);

    void QuestRevertResummedTransCoeff(Cell_small &grid_pt);
    void QuestRevert(
        Cell_small &grid_pt, const int ieta, const int ix, const int iy);
    void QuestRevert_qmu(
        Cell_small &grid_pt, const int ieta, const int ix, const int iy);

    void MakeDeltaQI(
        const double tau, Fields &arenaFieldsCurr, const int ix, const int iy,
        const int ieta, TJbVec &qi, const int rk_flag);

    double MaxSpeed(
        const double tau, const int direc, const ReconstCell &grid_p,
        double &pressure);

    double get_TJb(
        const ReconstCell &grid_p, const int mu, const int nu,
        const double pressure);

    void solveEigenvaluesWmunu(Cell_small &grid_pt);
    void nCausalityConstraints(Cell_small &grid_pt);
};

#endif  // SRC_ADVANCE_H_
