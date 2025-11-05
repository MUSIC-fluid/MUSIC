#include "hydro_source_Tmunu.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "util.h"

// Constructor
HydroSourceTmunu::HydroSourceTmunu(InitData &DATA_in) : DATA(DATA_in) {
    // Read the IP-Glasma event file
    readIPGevent();
    tau_source_ = DATA.tau0;
    DATA.tau0 = DATA.tau0 - DATA.delta_tau;
    set_source_tau_min(tau_source_);
    set_source_tau_max(tau_source_);
}

// Function to read the IP-Glasma event file and initialize data
void HydroSourceTmunu::readIPGevent() {
    // Open file
    std::string IPG_filename = DATA.initName;
    std::ifstream profile(IPG_filename);

    if (!profile.is_open()) {
        std::cerr << "Error: Cannot open the IP-Glasma file: " << IPG_filename
                  << std::endl;
        std::exit(1);
    }

    // Read the information line
    std::string dummy;
    std::getline(profile, dummy);
    std::istringstream ss1(dummy);

    int nx, ny, neta;
    double deta, dx, dy, dummy2;
    // read the first line with general info
    ss1 >> dummy >> dummy >> dummy2 >> dummy >> neta >> dummy >> nx >> dummy
        >> ny >> dummy >> deta >> dummy >> dx >> dummy >> dy;

    DATA.nx = nx;
    DATA.ny = ny;
    DATA.delta_x = dx;
    DATA.delta_y = dy;

    if (nx <= 0 || ny <= 0) {
        std::cerr << "Error: Invalid grid size nx = " << nx << ", ny = " << ny
                  << std::endl;
        std::exit(1);
    }

    // Initialize containers
    energy_density_.resize(nx, std::vector<double>(ny, 0.0));
    velocity_.resize(nx, std::vector<std::array<double, 4>>(ny));
    shear_tensor_.resize(nx, std::vector<std::array<double, 10>>(ny));

    double tau0 = DATA.tau0;

    // Parallelize the loops if desired
    // #pragma omp parallel for collapse(2)
    for (int ix = 0; ix < nx; ix++) {
        for (int iy = 0; iy < ny; iy++) {
            double dummy1, dummy2, dummy3;
            double e, utau, ux, uy, ueta;
            double pitautau, pitaux, pitauy, pitaueta;
            double pixx, pixy, pixeta, piyy, piyeta, pietaeta;

            if (profile.eof()) {
                std::cerr << "Error: Unexpected end of file while reading data."
                          << std::endl;
                std::exit(1);
            }

            std::getline(profile, dummy);
            std::istringstream ss(dummy);
            // Read data from the file
            if (!(ss >> dummy1 >> dummy2 >> dummy3 >> e >> utau >> ux >> uy
                  >> ueta >> pitautau >> pitaux >> pitauy >> pitaueta >> pixx
                  >> pixy >> pixeta >> piyy >> piyeta >> pietaeta)) {
                std::cerr << "Error: Failed to parse data line: " << dummy
                          << std::endl;
                std::exit(1);
            }

            // Set x_size and y_size based on the data
            if (ix == 0 && iy == 0) {
                DATA.x_size = -dummy2 * 2;
                DATA.y_size = -dummy3 * 2;
                DATA.eta_size = DATA.delta_eta * DATA.neta;
                music_message << "eta_size=" << DATA.eta_size
                              << ", x_size=" << DATA.x_size
                              << ", y_size=" << DATA.y_size;
                music_message.flush("info");
            }

            // Adjust ueta and shear tensor components
            ueta *= tau0;
            utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
            double visFactor = DATA.sFactor * DATA.preEqVisFactor;
            if (DATA.FlagResumTransportCoeff) {
                visFactor = DATA.sFactor;
            }
            pixx = pixx * visFactor;
            pixy = pixy * visFactor;
            pixeta = pixeta * visFactor;
            piyy = piyy * visFactor;
            piyeta = piyeta * visFactor;
            pietaeta =
                ((2.
                      * (ux * uy * pixy + ux * ueta * pixeta
                         + uy * ueta * piyeta)
                  - (utau * utau - ux * ux) * pixx
                  - (utau * utau - uy * uy) * piyy)
                 / (utau * utau - ueta * ueta));
            pitaux = (pixx * ux + pixy * uy + pixeta * ueta) / utau;
            pitauy = (pixy * ux + piyy * uy + piyeta * ueta) / utau;
            pitaueta = (pixeta * ux + piyeta * uy + pietaeta * ueta) / utau;
            pitautau = (pitaux * ux + pitauy * uy + pitaueta * ueta) / utau;

            // Store the values in the containers
            energy_density_[ix][iy] = e * DATA.sFactor / Util::hbarc;  // 1/fm^4
            velocity_[ix][iy] = {utau, ux, uy, ueta};
            shear_tensor_[ix][iy] = {pitautau, pitaux,  pitauy, pitaueta,
                                     pixx,     pixy,    pixeta, piyy,
                                     piyeta,   pietaeta};  // 1/fm^4
        }
    }
    profile.close();
    std::cout << "neta = " << DATA.neta << ", eta_size = " << DATA.eta_size
              << ", delta_eta = " << DATA.delta_eta << std::endl;
    std::cout << "nx = " << DATA.nx << ", x_size = " << DATA.x_size
              << ", delta_x = " << DATA.delta_x << std::endl;
    std::cout << "ny = " << DATA.ny << ", x_size = " << DATA.y_size
              << ", delta_y = " << DATA.delta_y << std::endl;
}

// Function to calculate energy-momentum source terms
void HydroSourceTmunu::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s,
    const std::array<double, 4> &u_mu, std::array<double, 4> &j_mu) const {
    const double dtau = DATA.delta_tau;
    j_mu = {0};

    if (std::abs((tau - tau_source_)) > 0.5 * dtau) return;

    // Compute indices from positions
    const int ix =
        static_cast<int>((x + DATA.x_size / 2.) / DATA.delta_x + 0.1);
    const int iy =
        static_cast<int>((y + DATA.y_size / 2.) / DATA.delta_y + 0.1);

    // Bounds checking
    if (ix < 0 || ix >= DATA.nx || iy < 0 || iy >= DATA.ny) {
        std::cerr
            << "Error: Index out of bounds in get_hydro_energy_source: ix = "
            << ix << ", iy = " << iy << std::endl;
        std::cerr << "x = " << x << ", y = " << y << std::endl;
        std::cerr << "x_size = " << DATA.x_size << ", y_size = " << DATA.y_size
                  << std::endl;
        std::cerr << "delta_x = " << DATA.delta_x
                  << ", delta_y = " << DATA.delta_y << std::endl;
        return;
    }

    // Access data
    double epsilon = energy_density_[ix][iy];
    const auto &Wmunu = shear_tensor_[ix][iy];
    const auto &u = velocity_[ix][iy];

    // Compute j^0 and j^i  (energy-momentum source term)
    j_mu[0] = 4. / 3. * epsilon * u[0] * u[0] - (epsilon / 3.) + Wmunu[0];
    j_mu[1] = 4. / 3. * epsilon * u[0] * u[1] + Wmunu[1];
    j_mu[2] = 4. / 3. * epsilon * u[0] * u[2] + Wmunu[2];
    j_mu[3] = 4. / 3. * epsilon * u[0] * u[3] + Wmunu[3];

    const double prefactors = 1.0 / dtau;
    for (int i = 0; i < 4; ++i) {
        j_mu[i] *= prefactors;  // 1/fm^5
    }
}
