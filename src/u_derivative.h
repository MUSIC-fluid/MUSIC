#ifndef SRC_U_DERIVATIVE_H_
#define SRC_U_DERIVATIVE_H_

#include "util.h"
#include "data.h"
#include "cell.h"
#include "grid.h"
#include "minmod.h"
#include "data_struct.h"
#include <string.h>
#include <iostream>

class U_derivative {
 private:
     const InitData &DATA;
     const EOS &eos;
     const double T_tol = 1e-5;
     Minmod minmod;
     dUsupMat dUsup;
     Mat4x4 dUoverTsup;
     Mat4x4 dUTsup;

 public:
    U_derivative(const InitData &DATA_in, const EOS &eosIn);
    void MakedU(const double tau, SCGrid &arena_prev, SCGrid &arena_current,
                const int ix, const int iy, const int ieta);

    //! this function returns the expansion rate on the grid
    double calculate_expansion_rate(double tau, SCGrid &arena,
                                    int ieta, int ix, int iy);

    //! this function returns Du^\mu
    void calculate_Du_supmu(const double tau, SCGrid &arena, const int ieta,
                            const int ix, const int iy, DumuVec &a);

    //! this function returns the vector D^\mu(\mu_B/T)
    void get_DmuMuBoverTVec(DmuMuBoverTVec &vec);

    //! this function computes the kinetic vorticity
    void calculate_kinetic_vorticity_with_spatial_projector(
        const double tau, SCGrid &arena,
        const int ieta, const int ix, const int iy, const DumuVec &a_local,
        VorticityVec &omega);

    //! this function computes the kinetic vorticity without spatial projection
    void calculate_kinetic_vorticity_no_spatial_projection(
        const double tau, SCGrid &arena,
        const int ieta, const int ix, const int iy, VorticityVec &omega);

    //! this function computes the thermal vorticity
    void calculate_thermal_vorticity(const double tau, SCGrid &arena,
        const int ieta, const int ix, const int iy, VorticityVec &omega);

    //! this function computes the T-vorticity
    void calculate_T_vorticity(const double tau, SCGrid &arena,
        const int ieta, const int ix, const int iy, VorticityVec &omega);

    //! This funciton returns the velocity shear tensor sigma^\mu\nu
    void calculate_velocity_shear_tensor(
        const double tau, SCGrid &arena, const int ieta, const int ix,
        const int iy, const DumuVec &a_local, VelocityShearVec &sigma);

    int MakeDSpatial(const double tau, SCGrid &arena, const int ix,
                     const int iy, const int ieta);
    int MakeDTau(const double tau, const Cell_small *grid_pt_prev,
                 const Cell_small *grid_pt);

    //! This is a shell function to compute all 4 kinds of vorticity tensors
    void compute_vorticity_shell(
        const double tau, SCGrid &arena_prev, SCGrid &arena_curr,
        const int ieta, const int ix, const int iy, const double eta,
        VorticityVec &omega_local_k, VorticityVec &omega_local_knoSP,
        VorticityVec &omega_local_th, VorticityVec &omega_local_T);

    //! This function transforms the vorticity tensor from tau-eta to tz
    VorticityVec transform_vorticity_to_tz(const VorticityVec omega_Mline,
                                           const double eta);
};

#endif
