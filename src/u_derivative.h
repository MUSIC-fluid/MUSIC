#ifndef SRC_U_DERIVATIVE_H_
#define SRC_U_DERIVATIVE_H_

#include <string.h>

#include <iostream>

#include "cell.h"
#include "data.h"
#include "data_struct.h"
#include "eos.h"
#include "fields.h"
#include "minmod.h"
#include "util.h"

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
    void MakedU(
        const double tau, Fields &arenaFieldsPrev, Fields &arenaFieldsCurr,
        const int fieldIdx, const int ix, const int iy, const int ieta);

    //! this function returns the expansion rate on the grid
    double calculate_expansion_rate(
        const double tau, Fields &arena, const int fieldIdx);

    //! this function returns Du^\mu
    void calculate_Du_supmu(
        const double tau, Fields &arena, const int fieldIdx, DumuVec &a);

    //! this function returns the vector D^\mu(\mu_B/T)
    void get_DmuMuBoverTVec(DmuMuBoverTVec &vec);

    //! this function computes the kinetic vorticity
    void calculate_kinetic_vorticity_with_spatial_projector(
        const double tau, Fields &arena, const int fieldIdx,
        const DumuVec &a_local, VorticityVec &omega);

    //! this function computes the kinetic vorticity without spatial projection
    void calculate_kinetic_vorticity_no_spatial_projection(
        const double tau, Fields &arena, const int fieldIdx,
        VorticityVec &omega);

    //! this function computes the thermal vorticity
    void calculate_thermal_vorticity(
        const double tau, Fields &arena, const int fieldIdx,
        VorticityVec &omega);

    //! this function computes the T-vorticity
    void calculate_T_vorticity(
        const double tau, Fields &arena, const int fieldIdx,
        VorticityVec &omega);

    //! This funciton returns the velocity shear tensor sigma^\mu\nu
    void calculate_velocity_shear_tensor(
        const double tau, Fields &arena, const int fieldIdx,
        const double theta_u_local, const DumuVec &a_local,
        VelocityShearVec &sigma);

    int MakeDSpatial(
        const double tau, Fields &arena, const int fieldIdx, const int ix,
        const int iy, const int ieta);

    int MakeDTau(
        const double tau, const Fields &arenaFieldsPrev,
        const Fields &arenaFieldsCurr, const int fieldIdx);

    //! this function computes the thermal shear tensor
    void calculate_thermal_shear_tensor(
        const double tau, Fields &arena, const int fieldIdx,
        VelocityShearVec &sigma_th);

    //! This is a shell function to compute all 4 kinds of vorticity tensors
    void compute_vorticity_shell(
        const double tau, Fields &arena_prev, Fields &arena_curr,
        const int ieta, const int ix, const int iy, const int fieldIdx,
        const double eta, VorticityVec &omega_local_k,
        VorticityVec &omega_local_knoSP, VorticityVec &omega_local_th,
        VorticityVec &omega_local_T, VelocityShearVec &sigma_local,
        DmuMuBoverTVec &DbetaMu);

    //! This function transforms the vorticity tensor from tau-eta to tz
    VorticityVec transform_vorticity_to_tz(
        const VorticityVec omega_Mline, const double eta);
    DmuMuBoverTVec transform_vector_to_tz(
        const DmuMuBoverTVec vec_Mline, const double eta);
    VelocityShearVec transform_SigmaMuNu_to_tz(
        const VelocityShearVec sigma_Mline, const double eta);
};

#endif
