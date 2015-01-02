#ifndef U_DERIVATIVE_H
#define U_DERIVATIVE_H

#include "util.h"
#include "data.h"
#include "grid.h"
#include <string.h>
#include <iostream>

class U_derivative{

 private:
  Minmod *minmod;
  // Sangyong Nov 18 2014: added EOS *eos;
  EOS *eos;
  
 public:
  // Sangyong Nov 18 2014: added EOS *eos in the argument
  U_derivative(EOS *eosIn, InitData* DATA_in);//constructor
  ~U_derivative();
  int UpdateDSpatial(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		     Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank);
  int MakedU(double tau, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int rk_flag, int size, int rank);
  int MakeDSpatial(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		   Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank);
  int MakeDTau(double tau, InitData *DATA, Grid *grid_pt, int rk_flag);
  int UpdateDTau(double tau, InitData *DATA, Grid *grid_pt, int rk_flag);
  int UpdatedU(double tau, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int rk_flag, int size, int rank);
};
#endif
