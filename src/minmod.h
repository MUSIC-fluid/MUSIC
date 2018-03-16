// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_MINMOD_H_
#define SRC_MINMOD_H_

#include "data.h"
#include "iostream"
class Minmod {
 private:
    const double theta_flux;

 public:
    Minmod(const InitData &DATA);
    Minmod(double theta_in);
    double get_theta() const {return(theta_flux);}

    /* double minmod_dx(const double up1, const double u, const double um1) { */
    /*     const double diffup   = (up1 - u)*theta_flux; */
    /*     const double diffdown = (u - um1)*theta_flux; */
    /*     const double diffmid  = (up1 - um1)*0.5; */

    /*     if ( (diffup > 0.0) && (diffdown > 0.0) && (diffmid > 0.0) ) { */
    /*         return std::min(std::min(diffdown, diffmid), diffup); */
    /*     } else if ( (diffup < 0.0) && (diffdown < 0.0) && (diffmid < 0.0) ) { */
    /*         return std::max(std::max(diffdown, diffmid), diffup); */
    /*     } else { */
    /*       return 0.0; */
    /*     } */
    /* }/\* minmod_dx *\/ */

    double minmod_dx(const double up1, const double u, const double um1) const {
        const double diffup   = (up1 - u)*theta_flux;
        const double diffdown = (u - um1)*theta_flux;
        const double diffmid  = (up1 - um1)*0.5;
	if(diffup==0.)
	  return 0.;
	else
	  return diffup*std::max(0., std::min(1.,std::min(diffdown/diffup, diffmid/diffup)));
    }/* minmod_dx */

};

#endif  // SRC_MINMOD_H_
