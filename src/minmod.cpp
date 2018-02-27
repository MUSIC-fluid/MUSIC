#include "util.h"
#include "minmod.h"

Minmod::Minmod(InitData *DATA) {
    theta_flux = DATA->minmod_theta;
}

double Minmod::minmod_dx(double up1, double u, double um1) {
    double diffup   = (up1 - u)*theta_flux;
    double diffdown = (u - um1)*theta_flux;
    double diffmid  = (up1 - um1)*0.5;

    if ( (diffup > 0.0) && (diffdown > 0.0) && (diffmid > 0.0) ) {
        return std::min(std::min(diffdown, diffmid), diffup);
    } else if ( (diffup < 0.0) && (diffdown < 0.0) && (diffmid < 0.0) ) {
        return std::max(std::max(diffdown, diffmid), diffup);
    } else {
      return 0.0;
    }
}/* minmod_dx */
