// Copyright 2019 Chun Shen
#ifndef SRC_HYDRO_SOURCE_BULLET_H_
#define SRC_HYDRO_SOURCE_BULLET_H_

#include "hydro_source_base.h"
#include "grid.h"
#include "util.h"

#include <array>
#include <cmath>
#include <sstream>
#include <iostream>

class HydroSourceBullet : public HydroSourceBase{
 private:
    std::array<double, 3> r_;   ///< Position of the bullet
    std::array<double, 4> pmu_; ///< 4-momentum of the thermalized jet
    double tau_;                ///< Instant in which it p^\mu is deposited
    const InitData &DATA_;      ///< MUSIC initialization parameters

    std::array<double, 4> pmu_frac_;



    //Auxiliary functions
    GridT<int> TagPointsForFillEllipse(double x0, double y0, double z0);

    //Convert space coordinate to grid coordinates
    int get_ix(double x) const;
    int get_iy(double y) const;
    int get_ieta(double eta) const;

    //Spatial Rapidity <-> z coordinate conversion
    double get_spatial_rapidity(double tau, double z) const;
    double get_z(double tau, double etaS) const;

    //Convert grid coordinates to space coordinates
    double get_x(int ix) const;
    double get_y(int iy) const;
    double get_eta(int ieta) const;

    //Stores points flaged for pmu insertion
    GridT<int> point_flagged;

 public:

    //Constructors
    HydroSourceBullet(const InitData &DATA_in);
    
    //! this function returns the energy source term J^\mu at a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu, EnergyFlowVec &j_mu) ;

};

#endif  // SRC_HYDRO_SOURCE_BASE_H_
