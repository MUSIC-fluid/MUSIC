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
#include "./glauber.h"


class Init {
 private:
    Random *random;
    EOS *eos;
    Util * util;
    Glauber *glauber;

    // list of x and y coordinates of nucleons in nucleus A
    vector<ReturnValue> nucleusA;
    // list of x and y coordinates of nucleons in nucleus B
    vector<ReturnValue> nucleusB;

 public:
    Init(EOS *eos, Glauber* glauber);  // constructor
    ~Init();  // destructor

    void sampleTA();
    void InitArena(InitData *DATA, Grid ****arena, Grid ****Lneighbor,
                   Grid ****Rneighbor, int size, int rank);
    void LinkNeighbors(InitData *DATA, Grid ****arena, int size, int rank);
    int InitTJb(InitData *DATA, Grid ****arena, Grid ****Lneighbor,
                Grid ****Rneighbor, int size, int rank);

    // The following two functions return T_A and T_B.
    // Normalization: \int r T_A(r) dr dphi = A
    double TATarget(InitData *DATA, double r);
    double TAProjectile(InitData *DATA, double r);

    double eta_profile_normalisation(InitData *DATA, double eta);
    double eta_rhob_profile_normalisation(InitData *DATA, double eta);
    double eta_profile_left_factor(InitData *Data, double eta);
    double eta_profile_right_factor(InitData *Data, double eta);
    double eta_rhob_left_factor(InitData *Data, double eta);
    double eta_rhob_right_factor(InitData *Data, double eta);
};

#endif  // SRC_INIT_H_
