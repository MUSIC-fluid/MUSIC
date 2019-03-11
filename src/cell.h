// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef _SRC_CELL_H_
#define _SRC_CELL_H_

#include "data_struct.h"
#include <vector>

class Cell_small {
 public:
    double epsilon = 0;
    double rhob    = 0;
    FlowVec u;

    ViscousVec Wmunu;
    double pi_b    = 0.;

    std::vector<double> phi_Q;  //!< for critical slow modes
    
    //! array to store the renormalization to the EoS and transport
    //! coefficients from out-of-equilibrium phi_Q fields
    std::vector<double> critical_renormalization;
};

#endif  // SRC_CELL_H_
