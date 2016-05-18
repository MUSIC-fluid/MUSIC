#include "util.h"
#include "minmod.h"

Minmod::Minmod(InitData *DATA) {
    theta_flux = DATA->minmod_theta;
}

// destructor
Minmod::~Minmod() {
}

double Minmod::minmod_dx(double up1, double u, double um1) {
    double diffup = (up1 - u)*theta_flux;
    double diffdown = (u - um1)*theta_flux;
    double diffmid = (up1 - um1)*0.5;

    double tempf;
    if ( (diffup > 0.0) && (diffdown > 0.0) && (diffmid > 0.0) ) {
        tempf = mini(diffdown, diffmid);
        return mini(diffup, tempf);
    } else if ( (diffup < 0.0) && (diffdown < 0.0) && (diffmid < 0.0) ) {
        tempf = maxi(diffdown, diffmid);
        return maxi(diffup, tempf);
    } else {
      return 0.0;
    }
}/* minmod_dx */


double Minmod::minmod_theta_dx(double up1, double u, double um1, double theta) {
    double diffup = (up1 - u)*theta;
    double diffdown = (u - um1)*theta;
    double diffmid = (up1 - um1)*0.5;

    double tempf;
    if ( (diffup > 0.0) && (diffdown > 0.0) && (diffmid > 0.0) ) {
        tempf = mini(diffdown, diffmid);
        return mini(diffup, tempf);
    } else if ( (diffup < 0.0) && (diffdown < 0.0) && (diffmid < 0.0) ) {
        tempf = maxi(diffdown, diffmid);
        return maxi(diffup, tempf);
    } else {
        return 0.0;
    }
}/* minmod_dx */
