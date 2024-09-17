#include "hydro_source_Tmunu.h"
#include "util.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <array>

using Util::hbarc;

// Constructor
HydroSourceTmunu::HydroSourceTmunu(InitData &DATA_in, const EOS &eos_in)
    : DATA(DATA_in), eos(eos_in) {
    // Read the IP-Glasma event file
    readIPGevent();
}

// Function to read the IP-Glasma event file and initialize data
void HydroSourceTmunu::readIPGevent() {
    // Open file
    std::string IPG_filename = DATA.initName;  
    std::ifstream IPGinputFile(IPG_filename);

    if (!IPGinputFile.is_open()) {
        std::cerr << "Error: Cannot open the IP-Glasma file: " << IPG_filename << std::endl;
        std::exit(1);
    }

    // Step1: Initialize containers for energy density, velocity, and shear tensor
    energy_density_ = std::vector<std::vector<std::vector<double>>>(
        DATA.nx, std::vector<std::vector<double>>(DATA.ny, std::vector<double>(DATA.neta, 0.0)));
    
    velocity_ = std::vector<std::vector<std::vector<std::array<double, 4>>>>(
        DATA.nx, std::vector<std::vector<std::array<double, 4>>>(
            DATA.ny, std::vector<std::array<double, 4>>(DATA.neta)));
    
    shear_tensor_ = std::vector<std::vector<std::vector<std::array<double, 10>>>>(
        DATA.nx, std::vector<std::vector<std::array<double, 10>>>(
            DATA.ny, std::vector<std::array<double, 10>>(DATA.neta)));

    std::string line;
    // Step2: Reading IP-Glasma data
    int ix, iy, ieta;
    double pitautau, pitaux, pitauy, pitaueta;
    double pixx*DATA.sFactor, pixy*DATA.sFactor, pixeta*DATA.sFactor, piyy*DATA.sFactor, piyeta*DATA.sFactor, pietaeta*DATA.sFactor;
    double e*DATA.sFactor/hbarc, utau, ux, uy, ueta;
    double tau = DATA.tau0;

    while (std::getline(IPGinputFile, line)) {
        std::istringstream iss(line);

        iss >> ix >> iy >> ieta >> e >> utau >> ux >> uy >> ueta
            >> pitautau >> pitaux >> pitauy >> pitaueta
            >> pixx*DATA.sFactor >> pixy >> pixeta >> piyy >> piyeta >> pietaeta;


       ueta = ueta * tau;

        
        pitaueta *= tau;
        pietaeta *= tau * tau*Data.sFactor;
        pixeta *= tau*Data.sFactor;
        piyeta *= tau*Data.sFactor;

        energy_density_[ix][iy][ieta] = e;
        velocity_[ix][iy][ieta] = {utau, ux, uy, ueta};
        shear_tensor_[ix][iy][ieta] = {pitautau, pitaux, pitauy, pitaueta*tau*Data.sFactor,
                                       pixx*DATA.sFactor, pixy*DATA.sFactor, pixeta*DATA.sFactor, piyy*DATA.sFactor, piyeta*DATA.sFactor, pietaeta}; //multiply Data.sFactor
    }

    IPGinputFile.close();
}

// Function to calculate energy-momentum source term
void HydroSourceTmunu::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s,
    const std::array<double, 4>& u_mu, std::array<double, 4>& j_mu) const {

    j_mu = {0};


    // Grid indices from spatial coordinates
    int ix = static_cast<int>((x + DATA.x_size / 2.0) / DATA.delta_x);
    int iy = static_cast<int>((y + DATA.y_size / 2.0) / DATA.delta_y);
    int ieta = static_cast<int>((eta_s + DATA.eta_size / 2.0) / DATA.delta_eta);

    // Bounds for grid indices
    if (ix < 0 || ix >= DATA.nx || iy < 0 || iy >= DATA.ny || ieta < 0 || ieta >= DATA.neta) {
        std::cerr << "Out of bounds grid access in HydroSourceTmunu::get_hydro_energy_source" << std::endl;
        std::exit(1);
    }

    // Read the energy density, velocity, and viscous tensor from the IP-Glasma data
    double epsilon = energy_density_[ix][iy][ieta];  // Energy density e
    const auto &Wmunu = shear_tensor_[ix][iy][ieta]; // Shear viscous tensor pi^{mu nu}
    const auto &u = velocity_[ix][iy][ieta];         // Fluid velocity u^mu

    // Compute j^0
    j_mu[0] = epsilon * (u[0] * u[0]) + (epsilon / 3) * (1 + u[0] * u[0]) + Wmunu[0];  // \pi^{\tau \tau}

    // Compute spatial components
    j_mu[1] = epsilon * u[0] * u[1] + (epsilon / 3) * u[0] * u[1] + Wmunu[1];  // \pi^{\tau x}
    j_mu[2] = epsilon * u[0] * u[2] + (epsilon / 3) * u[0] * u[2] + Wmunu[2];  // \pi^{\tau y}
    j_mu[3] = epsilon * u[0] * u[3] + (epsilon / 3) * u[0] * u[3] + Wmunu[3];  // \pi^{\tau \eta}


    

    // Unit conversion factor
    const double prefactors = 1.0 / (DATA.delta_tau * Util::hbarc);
    j_mu[0] *= prefactors;
    j_mu[1] *= prefactors;
    j_mu[2] *= prefactors;
    j_mu[3] *= prefactors;
}
