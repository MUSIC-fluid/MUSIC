// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_MINMOD_H_
#define SRC_MINMOD_H_

#include "data.h"

class Minmod {
 private:
    double theta_flux;

 public:
    Minmod(InitData* DATA);
    double minmod_dx(double up1, double u, double um1);
};

#endif  // SRC_MINMOD_H_
