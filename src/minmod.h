// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_MINMOD_H_
#define SRC_MINMOD_H_

#include "data.h"

class Minmod {
 private:
    double theta_flux;

 public:
    Minmod(InitData* DATA);

    double minmod_dx(double up1, double u, double um1) {
        const double diffup   = (up1 - u)*theta_flux;
        const double diffdown = (u - um1)*theta_flux;
        const double diffmid  = (up1 - um1)*0.5;

        if ( (diffup > 0.0) && (diffdown > 0.0) && (diffmid > 0.0) ) {
            return std::min(std::min(diffdown, diffmid), diffup);
        } else if ( (diffup < 0.0) && (diffdown < 0.0) && (diffmid < 0.0) ) {
            return std::max(std::max(diffdown, diffmid), diffup);
        } else {
          return 0.0;
        }
    }/* minmod_dx */
};

#endif  // SRC_MINMOD_H_
