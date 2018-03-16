// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef _SRC_CELL_H_
#define _SRC_CELL_H_

#include "data_struct.h"
#include <array>

class Cell_small {
 public:
    double epsilon = 0;
    double rhob    = 0;
    FlowVec u;

    ViscousVec Wmunu;
    double pi_b    = 0.;
};

#endif  // SRC_GRID_H_
