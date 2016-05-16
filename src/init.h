// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_INIT_H_
#define SRC_INIT_H_

#include <stdio.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "./data.h"
#include "./grid.h"


class Init {
 private:
    EOS *eos;
    Util * util;

 public:
    Init(EOS *eos);  // constructor
    ~Init();  // destructor

    void InitArena(InitData *DATA, Grid ****arena, Grid ****Lneighbor,
                   Grid ****Rneighbor, int size, int rank);
    void LinkNeighbors(InitData *DATA, Grid ****arena, int size, int rank);
    int InitTJb(InitData *DATA, Grid ****arena, int size, int rank);

    double eta_profile_normalisation(InitData *DATA, double eta);
    double eta_rhob_profile_normalisation(InitData *DATA, double eta);
    double eta_profile_left_factor(InitData *Data, double eta);
    double eta_profile_right_factor(InitData *Data, double eta);
    double eta_rhob_left_factor(InitData *Data, double eta);
    double eta_rhob_right_factor(InitData *Data, double eta);
};

#endif  // SRC_INIT_H_
