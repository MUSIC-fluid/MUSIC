#include "util.h"
#include "minmod.h"

Minmod::Minmod()
{
}

// destructor
Minmod::~Minmod()
{
}

double Minmod::minmod_dx(double up1, double u, double um1, InitData *DATA)
{
 double theta, diffup, diffdown, diffmid;
 double tempf;
 
 theta = DATA->minmod_theta;
 
 diffup = up1 - u;
 diffup *= theta;

 diffdown = u - um1;
 diffdown *= theta;

 diffmid = up1 - um1;
 diffmid *= 0.5;

 if( (diffup > 0.0) && (diffdown > 0.0) && (diffmid > 0.0) )
  {
   tempf = mini(diffdown, diffmid);
   return mini(diffup, tempf);
  }
 else if( (diffup < 0.0) && (diffdown < 0.0) && (diffmid < 0.0) )
  {
   tempf = maxi(diffdown, diffmid);
   return maxi(diffup, tempf);
  }
 else return 0.0;

}/* minmod_dx */


double Minmod::minmod_theta_dx(double up1, double u, double um1, double theta)
{
 double diffup, diffdown, diffmid;
 double tempf;
 
 diffup = up1 - u;
 diffup *= theta;

 diffdown = u - um1;
 diffdown *= theta;

 diffmid = up1 - um1;
 diffmid *= 0.5;

 if( (diffup > 0.0) && (diffdown > 0.0) && (diffmid > 0.0) )
  {
   tempf = mini(diffdown, diffmid);
   return mini(diffup, tempf);
  }
 else if( (diffup < 0.0) && (diffdown < 0.0) && (diffmid < 0.0) )
  {
   tempf = maxi(diffdown, diffmid);
   return maxi(diffup, tempf);
  }
 else return 0.0;

}/* minmod_dx */

