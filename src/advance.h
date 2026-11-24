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

    // Deferred correction of the viscous time-derivative term (see
    // FirstRKStepT). Valid only if the grid and delta_tau are unchanged
    // since the last full RK2 step.
    bool viscousCorrectionValid_ = false;
    double viscousCorrectionTau_ = 0.;  // tau_{n+1} + delta_tau
    // the correction is applied where |V_{n+1} - V*| <= threshold * |V_{n+1}
    // - V_n|, i.e. where V is resolved smoothly in time. The threshold is 4
    // where partial_tau u comes from the ideal equations of motion and 0.25
    // where it comes from the backward time difference, blended with the
    // same weight (U_derivative::idealEoMWeight). NOTE: use 0.25 whenever
    // the backward difference is used (Du_from_ideal_EoM = 0, or dense
    // cells): it makes V noisy in dilute, fast cells and the correction
    // then feeds back on the flow (lumpy events blow up with 4).
    static constexpr double viscousCorrectionSmoothnessBackward_ = 0.25;
    static constexpr double viscousCorrectionSmoothnessIdealEoM_ = 4.;
    // the J^tau correction (q^tau) has its own guard: it reaches the flow
    // only through rho_B in the EoS, and with 0.25 it was skipped in the
    // cold, baryon-rich cells between the lumps (N_B drifted by up to
    // 1.6e-3; 2e-4 with 4)
    static constexpr double viscousCorrectionSmoothnessBaryon_ = 4.;

  public:
    //! call when the grid or delta_tau changes between time steps
    void invalidateViscousCorrection() { viscousCorrectionValid_ = false; }

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
        Fields &arenaFieldsNext, Fields &arenaFieldsPrev);

    void FirstRKStepW(
        const double tau_it, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
        Fields &arenaFieldsNext, const int rk_flag, const double theta_local,
        const DumuVec &a_local, const VelocityShearVec &sigma_local,
        const VorticityVec &omega_local,
        const DmuMuBoverTVec &baryon_diffusion_vector, const int ieta,
        const int ix, const int iy, const int fieldIdx);

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
