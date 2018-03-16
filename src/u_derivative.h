#ifndef SRC_U_DERIVATIVE_H_
#define SRC_U_DERIVATIVE_H_

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
     const EOS &eos;
     Minmod minmod;
     dUsupMat dUsup;

 public:
    U_derivative(const InitData &DATA_in, const EOS &eosIn);
    void MakedU(double tau, SCGrid &arena_prev, SCGrid &arena_current,
                int ix, int iy, int ieta);

    //! this function returns the expansion rate on the grid
    double calculate_expansion_rate(double tau, SCGrid &arena,
                                    int ieta, int ix, int iy);

    //! this function returns Du^\mu
    void calculate_Du_supmu(double tau, SCGrid &arena, int ieta, int ix, int iy,
                            DumuVec &a);

    //! this function returns the vector D^\mu(\mu_B/T)
    void get_DmuMuBoverTVec(DmuMuBoverTVec &vec);

    //! This funciton returns the velocity shear tensor sigma^\mu\nu
    void calculate_velocity_shear_tensor(
        double tau, SCGrid &arena, int ieta, int ix, int iy,
        DumuVec &a_local, VelocityShearVec &sigma);
    int MakeDSpatial(double tau, SCGrid &arena, int ix, int iy, int ieta);
    int MakeDTau(double tau, Cell_small *grid_pt_prev, Cell_small *grid_pt);
};

#endif
