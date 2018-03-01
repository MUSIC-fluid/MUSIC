// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_INIT_H_
#define SRC_INIT_H_

#include <stdio.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "eos.h"
#include "hydro_source.h"
#include "pretty_ostream.h"

class Init {
 private:
    InitData &DATA;
    const EOS &eos;
    hydro_source *hydro_source_ptr;
    pretty_ostream music_message;

 public:
    Init(const EOS &eos, InitData &DATA_in, hydro_source *hydro_source_in);

    void   InitArena                       (Grid &arena);
    int    InitTJb                         (Grid &arena);
    void   initial_Gubser_XY               (int ieta, Grid &arena);
    void   initial_1p1D_eta                (Grid &arena);
    void   initial_IPGlasma_XY             (int ieta, Grid &arena);
    void   initial_IPGlasma_XY_with_pi     (int ieta, Grid &arena);
    void   initial_MCGlb_with_rhob_XY      (int ieta, Grid &arena);
    void   initial_MCGlbLEXUS_with_rhob_XY (int ieta, Grid &arena);
    void   initial_AMPT_XY                 (int ieta, Grid &arena);
    void   initial_UMN_with_rhob           (Grid &arena);
    
    double eta_profile_normalisation       (double eta);
    double eta_rhob_profile_normalisation  (double eta);
    double eta_profile_left_factor         (double eta);
    double eta_profile_right_factor        (double eta);
    double eta_rhob_left_factor            (double eta);
    double eta_rhob_right_factor           (double eta);
    void   output_initial_density_profiles (Grid &arena);
};

#endif  // SRC_INIT_H_
