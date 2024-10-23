// hydro_source_Tmunu.h

#ifndef HYDRO_SOURCE_TMUNU_H
#define HYDRO_SOURCE_TMUNU_H

#include "hydro_source_base.h"  // Include the base class header
#include "data_struct.h"
#include <vector>
#include <array>

class HydroSourceTmunu : public HydroSourceBase {  // Inherit from HydroSourceBase
public:
    // Constructor
    HydroSourceTmunu(InitData &DATA_in);

    // Override the virtual function from HydroSourceBase
    void get_hydro_energy_source(
        const double tau, const double x, const double y, const double eta_s,
        const std::array<double, 4>& u_mu, std::array<double, 4>& j_mu) const;

private:
    // Add your private members and functions here
    InitData &DATA;

    double tau_source_;
    // Data containers
    std::vector<std::vector<double>> energy_density_;
    std::vector<std::vector<std::array<double, 4>>> velocity_;
    std::vector<std::vector<std::array<double, 10>>> shear_tensor_;

    // Function to read the IP-Glasma event file and initialize data
    void readIPGevent();
};

#endif  // HYDRO_SOURCE_TMUNU_H
