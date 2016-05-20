#ifndef U_DERIVATIVE_H_
#define U_DERIVATIVE_H_

#include "util.h"
#include "data.h"
#include "grid.h"
#include <string.h>
#include <iostream>

class U_derivative {
 private:
     Minmod *minmod;
     // Sangyong Nov 18 2014: added EOS *eos;
     EOS *eos;
  
 public:
     // Sangyong Nov 18 2014: added EOS *eos in the argument
     U_derivative(EOS *eosIn, InitData* DATA_in);  // constructor
     ~U_derivative();
     int MakedU(double tau, InitData *DATA, Grid ***arena, int rk_flag);
     int MakeDSpatial(double tau, InitData *DATA, Grid *grid_pt, int rk_flag);
     int MakeDTau(double tau, InitData *DATA, Grid *grid_pt, int rk_flag);
};
#endif
