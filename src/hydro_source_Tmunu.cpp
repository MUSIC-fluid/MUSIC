#include "hydro_source_Tmunu.h"
#include "util.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <array>
#include <omp.h>  // Include for OpenMP

using Util::hbarc;

// Constructor
HydroSourceTmunu::HydroSourceTmunu(InitData &DATA_in)
    : DATA(DATA_in) {
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
    double tau0 = DATA.tau0;
    std::ifstream profile(DATA.initName.c_str());  // Fixed the ifstream

    std::string dummy;
    // read the information line
    std::getline(profile, dummy);

    const int nx = DATA.nx;  // Changed to correct variable from DATA
    const int ny = DATA.ny;  // Changed to correct variable from DATA
    std::vector<double> temp_profile_ed(nx * ny, 0.0);
    std::vector<double> temp_profile_utau(nx * ny, 0.0);
    std::vector<double> temp_profile_ux(nx * ny, 0.0);
    std::vector<double> temp_profile_uy(nx * ny, 0.0);
    std::vector<double> temp_profile_ueta(nx * ny, 0.0);
    std::vector<double> temp_profile_pitautau(nx * ny, 0.0);
    std::vector<double> temp_profile_pitaux(nx * ny, 0.0);
    std::vector<double> temp_profile_pitauy(nx * ny, 0.0);
    std::vector<double> temp_profile_pitaueta(nx * ny, 0.0);
    std::vector<double> temp_profile_pixx(nx * ny, 0.0);
    std::vector<double> temp_profile_pixy(nx * ny, 0.0);
    std::vector<double> temp_profile_pixeta(nx * ny, 0.0);
    std::vector<double> temp_profile_piyy(nx * ny, 0.0);
    std::vector<double> temp_profile_piyeta(nx * ny, 0.0);
    std::vector<double> temp_profile_pietaeta(nx * ny, 0.0);

    // Step2: Reading IP-Glasma data
    int ix, iy, ieta;

    double e, utau, ux, uy, ueta;
    double density, dummy1, dummy2, dummy3;
    double pixx, pixy, pixeta, piyy, piyeta, pietaeta, pitautau, pitaux, pitauy, pitaueta;

    for (ix = 0; ix < nx; ix++) {
        for (iy = 0; iy < ny; iy++) {
            int idx = iy + ix * ny;
            std::getline(profile, dummy);
            std::stringstream ss(dummy);  // Fixed the issue with ss and iss
            ss >> ix >> iy >> ieta >> e >> utau >> ux >> uy >> ueta
               >> pitautau >> pitaux >> pitauy >> pitaueta
               >> pixx >> pixy >> pixeta >> piyy >> piyeta >> pietaeta;

            ueta = ueta * tau0;
            temp_profile_ed[idx] = density * DATA.sFactor / hbarc;  // 1/fm^4
            temp_profile_ux[idx] = ux;
            temp_profile_uy[idx] = uy;
            temp_profile_ueta[idx] = ueta;
            temp_profile_utau[idx] = sqrt(1. + ux * ux + uy * uy + ueta * ueta);

            // Viscous components already in 1/fm^4 from the IP-Glasma output
            double visFactor = DATA.sFactor * DATA.preEqVisFactor;
            temp_profile_pixx[idx] = pixx * visFactor;
            temp_profile_pixy[idx] = pixy * visFactor;
            temp_profile_pixeta[idx] = pixeta * tau0 * visFactor;
            temp_profile_piyy[idx] = piyy * visFactor;
            temp_profile_piyeta[idx] = piyeta * tau0 * visFactor;

            utau = temp_profile_utau[idx];
            temp_profile_pietaeta[idx] = (
                (2. * (ux * uy * temp_profile_pixy[idx]
                       + ux * ueta * temp_profile_pixeta[idx]
                       + uy * ueta * temp_profile_piyeta[idx])
                 - (utau * utau - ux * ux) * temp_profile_pixx[idx]
                 - (utau * utau - uy * uy) * temp_profile_piyy[idx])
                / (utau * utau - ueta * ueta));
            temp_profile_pitaux[idx] = (1. / utau
                                        * (temp_profile_pixx[idx] * ux
                                           + temp_profile_pixy[idx] * uy
                                           + temp_profile_pixeta[idx] * ueta));
            temp_profile_pitauy[idx] = (1. / utau
                                        * (temp_profile_pixy[idx] * ux
                                           + temp_profile_piyy[idx] * uy
                                           + temp_profile_piyeta[idx] * ueta));
            temp_profile_pitaueta[idx] = (1. / utau
                                          * (temp_profile_pixeta[idx] * ux
                                             + temp_profile_piyeta[idx] * uy
                                             + temp_profile_pietaeta[idx] * ueta));
            temp_profile_pitautau[idx] = (1. / utau
                                          * (temp_profile_pitaux[idx] * ux
                                             + temp_profile_pitauy[idx] * uy
                                             + temp_profile_pitaueta[idx] * ueta));

            if (ix == 0 && iy == 0) {
                DATA.x_size = -dummy2 * 2;
                DATA.y_size = -dummy3 * 2;
                if (omp_get_thread_num() == 0) {
                    music_message << "eta_size=" << DATA.eta_size
                                  << ", x_size=" << DATA.x_size
                                  << ", y_size=" << DATA.y_size;
                    music_message.flush("info");
                }
            }
        }

        // Store the data in the containers
        double tau = DATA.tau0;  // Fixed tau declaration
        shear_tensor_[ix][iy][ieta] = {pitautau, pitaux, pitauy, pitaueta * tau * DATA.sFactor,
                                       pixx * DATA.sFactor, pixy * DATA.sFactor, pixeta * DATA.sFactor, piyy * DATA.sFactor, piyeta * DATA.sFactor, pietaeta}; // multiply by DATA.sFactor
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
