// Integration test driver: the hydro history kept in memory
// (store_hydro_info_in_memory = 1) must reproduce the evolution file
// (output_evolution_data = 1, text) at every stored grid point.
//
// usage: integration_hydro_memory.e music_input
// The input must set store_hydro_info_in_memory 1, output_evolution_data 1
// and outputBinaryEvolution 0. Prints a summary and returns 0 on success.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "data_struct.h"
#include "music.h"

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " music_input" << std::endl;
        return 2;
    }
    MUSIC music(argv[1]);
    music.initialize_hydro();
    music.run_hydro();

    const InitData &DATA = music.getParameters();
    if (DATA.store_hydro_info_in_memory != 1 || DATA.outputEvolutionData != 1
        || DATA.outputBinaryEvolution != 0) {
        std::cerr << "the input must set store_hydro_info_in_memory 1, "
                  << "output_evolution_data 1, outputBinaryEvolution 0"
                  << std::endl;
        return 2;
    }

    std::ifstream evo("evolution_xyeta.dat");
    if (!evo.is_open()) {
        std::cerr << "evolution_xyeta.dat not found" << std::endl;
        return 2;
    }

    const int nx = DATA.nx, ny = DATA.ny, neta = DATA.neta;
    const int sx = DATA.output_evolution_every_N_x;
    const int sy = DATA.output_evolution_every_N_y;
    const int seta = DATA.output_evolution_every_N_eta;
    const int nslice = music.get_hydro_tau_max() > 0
                           ? static_cast<int>(
                                 (music.get_hydro_tau_max()
                                  - music.get_hydro_tau0())
                                     / music.get_hydro_dtau()
                                 + 0.5)
                                 + 1
                           : 0;

    long npoints = 0;
    double maxdT = 0., maxdv = 0., Tscale = 0.;
    // the last stored slice is skipped: the interpolation needs tau + dtau
    for (int itau = 0; itau < nslice - 1; itau++) {
        const double tau =
            music.get_hydro_tau0() + itau * music.get_hydro_dtau();
        for (int ieta = 0; ieta < neta; ieta += seta) {
            const double eta =
                DATA.boost_invariant
                    ? 0.
                    : ieta * DATA.delta_eta - DATA.eta_size / 2.;
            for (int iy = 0; iy < ny; iy += sy) {
                const double y = iy * DATA.delta_y - DATA.y_size / 2.;
                for (int ix = 0; ix < nx; ix += sx) {
                    const double x = ix * DATA.delta_x - DATA.x_size / 2.;
                    double T, muB, vx, vy, vz;
                    if (!(evo >> T >> muB >> vx >> vy >> vz)) {
                        std::cerr << "evolution_xyeta.dat ended early at "
                                  << "itau = " << itau << std::endl;
                        return 1;
                    }
                    fluidCell cell;
                    music.get_hydro_info(
                        x, y, tau * std::sinh(eta), tau * std::cosh(eta),
                        &cell);
                    // std::max drops NaNs, so a non-finite value on either
                    // side would otherwise leave the errors at zero
                    const double vals[8] = {
                        T,       vx,      vy,     vz, cell.temperature,
                        cell.vx, cell.vy, cell.vz};
                    for (const double v : vals) {
                        if (!std::isfinite(v)) {
                            std::cerr << "non-finite value at itau = " << itau
                                      << ", ix = " << ix << ", iy = " << iy
                                      << ", ieta = " << ieta << std::endl;
                            return 1;
                        }
                    }
                    Tscale = std::max(Tscale, T);
                    maxdT = std::max(maxdT, std::abs(cell.temperature - T));
                    maxdv = std::max(
                        maxdv, std::max(
                                   std::abs(cell.vx - vx),
                                   std::max(
                                       std::abs(cell.vy - vy),
                                       std::abs(cell.vz - vz))));
                    npoints++;
                }
            }
        }
    }

    const bool ok = npoints > 0 && maxdT <= 1e-5 * Tscale && maxdv <= 1e-5;
    std::printf(
        "compared %ld points in %d slices: max |dT| = %.3e GeV "
        "(T_max = %.4f GeV), max |dv| = %.3e -> %s\n",
        npoints, nslice - 1, maxdT, Tscale, maxdv, ok ? "OK" : "FAILED");
    return ok ? 0 : 1;
}
