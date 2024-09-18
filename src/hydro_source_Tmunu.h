#ifndef SRC_HYDRO_SOURCE_TMUNU_H_
#define SRC_HYDRO_SOURCE_TMUNU_H_

#include "data.h"
#include "hydro_source_base.h"
#include "eos.h"
#include "cell.h"
#include <vector>
#include <array>
#include <memory>
#include <iostream>
#include <cstdlib> // for exit function

class HydroSourceTmunu : public HydroSourceBase {
 private:
    InitData &DATA;

    std::vector<std::vector<std::vector<double>>> energy_density_;
    std::vector<std::vector<std::vector<std::array<double, 4>>>> velocity_;
    std::vector<std::vector<std::vector<std::array<double, 10>>>> shear_tensor_;

    void readIPGevent();

 public:
    HydroSourceTmunu() = delete; 
    HydroSourceTmunu(InitData &DATA_in);
    //~HydroSourceTmunu() = default; 
    void get_hydro_energy_source(const double tau, const double x, const double y, const double eta_s,
                                 const std::array<double, 4>& u_mu, std::array<double, 4>& j_mu) const;
};
#endif
