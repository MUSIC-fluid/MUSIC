#ifndef U_DERIVATIVE_H_
#define U_DERIVATIVE_H_

#include "util.h"
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "data_struct.h"
#include <string.h>
#include <iostream>

class U_derivative {
 private:
     const InitData &DATA;
     const EOS &eos; // Sangyong Nov 18 2014: added EOS *eos;
     Minmod minmod;
  
 public:
    // Sangyong Nov 18 2014: added EOS *eos in the argument
    U_derivative(const EOS &eosIn, const InitData &DATA_in);  // constructor
    int MakedU(double tau, const InitData &DATA, Grid &arena, int rk_flag);

    //! this function returns the expansion rate on the grid
    double calculate_expansion_rate(double tau, Grid &arena, int ieta, int ix, int iy, int rk_flag);

    //! this function returns Du^\mu
    void calculate_Du_supmu(double tau, Grid &arena, int ieta, int ix, int iy, int rk_flag, DumuVec &a);

    //! This funciton returns the velocity shear tensor sigma^\mu\nu
    void calculate_velocity_shear_tensor(double tau, Grid &arena, int ieta, int ix, int iy, int rk_flag, DumuVec &a_local, VelocityShearVec &sigma);
    int MakeDSpatial(double tau, const InitData &DATA, Grid &arena, int ix, int iy, int ieta, int rk_flag);
    int MakeDTau(double tau, const InitData &DATA, Cell *grid_pt, int rk_flag);
};
#endif
