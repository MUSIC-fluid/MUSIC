// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_ADVANCE_H_
#define SRC_ADVANCE_H_

#include <memory>
#include "data.h"
#include "cell.h"
#include "fields.h"
#include "dissipative.h"
#include "minmod.h"
#include "u_derivative.h"
#include "reconst.h"
#include "hydro_source_base.h"
#include "pretty_ostream.h"

class Advance {
 private:
    const InitData &DATA;
    const EOS &eos;
    std::shared_ptr<HydroSourceBase> hydro_source_terms_ptr;

    Diss diss_helper;
    Minmod minmod;
    Reconst reconst_helper;
    pretty_ostream music_message;

    bool flag_add_hydro_source;

 public:
    Advance(const EOS &eosIn, const InitData &DATA_in,
            std::shared_ptr<HydroSourceBase> hydro_source_ptr_in);

    void AdvanceIt(const double tau_init, Fields &arenaFieldsPrev,
                   Fields &arenaFieldsCurr, Fields &arenaFieldsNext,
                   const int rk_flag);

    void FirstRKStepT(const double tau, const double x_local,
                      const double y_local, const double eta_s_local,
                      const int ix, const int iy,
                      const int ieta, const int rk_flag,
                      const int fieldIdx, Fields &arenaFieldsCurr,
                      Fields &arenaFieldsNext, Fields &arenaFieldsPrev);

    void FirstRKStepW(const double tau_it, Fields &arenaFieldsPrev,
                      Fields &arenaFieldsCurr, Fields &arenaFieldsNext,
                      const int rk_flag, const double theta_local,
                      const DumuVec &a_local,
                      const VelocityShearVec &sigma_local,
                      const VorticityVec &omega_local,
                      const DmuMuBoverTVec &baryon_diffusion_vector,
                      const int ieta, const int ix, const int iy,
                      const int fieldIdx);

    void QuestRevert(const double tau, Cell_small &grid_pt,
                     const int ieta, const int ix, const int iy);
    void QuestRevert_qmu(const double tau, Cell_small &grid_pt,
                         const int ieta, const int ix, const int iy);

    void MakeDeltaQI(const double tau, Fields &arenaFieldsCurr,
                     const int ix, const int iy, const int ieta, TJbVec &qi,
                     const int rk_flag);

    double MaxSpeed(const double tau, const int direc,
                    const ReconstCell &grid_p, const double pressure);

    double get_TJb(const ReconstCell &grid_p, const int mu, const int nu,
                   const double pressure);
};

#endif  // SRC_ADVANCE_H_
