// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_EVOLVE_H_
#define SRC_EVOLVE_H_

#include <memory>
#include <vector>
#include "util.h"
#include "data.h"
#include "cell.h"
#include "grid_info.h"
#include "eos.h"
#include "advance.h"
#include "fields.h"
#include "hydro_source_base.h"
#include "pretty_ostream.h"
#include "HydroinfoMUSIC.h"
#include "surfaceCell.h"

// this is a control class for the hydrodynamic evolution
class Evolve {
 private:
    const EOS &eos;        // declare EOS object
    InitData &DATA;
    std::shared_ptr<HydroSourceBase> hydro_source_terms_ptr;

    Cell_info grid_info;
    Advance advance;
    pretty_ostream music_message;

    // simulation information
    int rk_order;

    // information about freeze-out surface
    // (only used when freezeout_method == 4)
    std::vector<double> epsFO_list;

    std::vector<SurfaceCell> surfaceCellVec_;
    std::vector<double> FO_nBvsEta_;
    std::vector<double> FO_nQvsEta_;
    std::vector<double> FO_nSvsEta_;

 public:
    Evolve(const EOS &eos, InitData &DATA_in,
           std::shared_ptr<HydroSourceBase> hydro_source_ptr_in);

    ~Evolve() {clearSurfaceCellVector();}

    void clearSurfaceCellVector() {surfaceCellVec_.clear();}
    int get_number_of_surface_cells() const {return(surfaceCellVec_.size());}
    void get_surface_cell_with_index(const int i, SurfaceCell &cell_i) {
        cell_i = surfaceCellVec_[i];
    }

    int EvolveIt(Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
                 Fields &arenaFieldsNext, HydroinfoMUSIC &hydro_info_ptr);

    int EvolveOneTimeStep(const int itau, Fields &arenaFieldsPrev,
                          Fields &arenaFieldsCurr, Fields &arenaFieldsNext,
                          Fields &freezeoutFieldPrev,
                          Fields &freezeoutFieldCurr,
                          HydroinfoMUSIC &hydro_info_ptr);

    void AdvanceRK(double tau, Fields* &fpPrev, Fields* &fpCurr,
                   Fields* &fpNext);

    int FreezeOut_equal_tau_Surface(double tau, Fields &arena_current);
    void FreezeOut_equal_tau_Surface_XY(double tau,
                                        int ieta, Fields &arena_current,
                                        int thread_id, double epsFO);
    int FindFreezeOutSurface_Cornelius(double tau,
        Fields &arena_prev, Fields &arena_current,
        Fields &arena_freezeout_prev, Fields &arena_freezeout);

    int FindFreezeOutSurface_Cornelius_XY(double tau, int ieta,
                                          Fields &arena_prev,
                                          Fields &arena_current,
                                          Fields &arena_freezeout_prev,
                                          Fields &arena_freezeout,
                                          int thread_id, double epsFO);
    int FindFreezeOutSurface_boostinvariant_Cornelius(
                double tau, Fields &arena_current, Fields &arena_freezeout);

    void store_previous_step_for_freezeout(Fields &arenaCurr,
                                           Fields &arenaFreeze);
    void regulate_qmu(const FlowVec u, const double q[],
                      double q_regulated[]) const;
    void regulate_Wmunu(const FlowVec u, const double Wmunu[4][4],
                        double Wmunu_regulated[4][4]) const;

    void initialize_freezeout_surface_info();

    Cell_small four_dimension_linear_interpolation(
        double* lattice_spacing, double fraction[2][4], Cell_small**** cube);
    Cell_small three_dimension_linear_interpolation(
        double* lattice_spacing, double fraction[2][3], Cell_small*** cube);
    Cell_aux four_dimension_linear_interpolation(
        double* lattice_spacing, double fraction[2][4], Cell_aux**** cube);
    Cell_aux three_dimension_linear_interpolation(
        double* lattice_spacing, double fraction[2][3], Cell_aux*** cube);
};

#endif  // SRC_EVOLVE_H_

