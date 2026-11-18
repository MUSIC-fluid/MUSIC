// Unit tests for the bug fixes on chun_dev (Nov 2026).
// Each test case names the fix it protects; most of them fail on the code
// before the fix.

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "HydroinfoMUSIC.h"
#include "advance.h"
#include "cell.h"
#include "data.h"
#include "dissipative.h"
#include "doctest.h"
#include "eos.h"
#include "eos_base.h"
#include "evolve.h"
#include "fields.h"
#include "grid_info.h"
#include "hydro_source_TATB.h"
#include "hydro_source_Tmunu.h"
#include "hydro_source_ampt.h"
#include "music.h"
#include "pretty_ostream.h"
#include "read_in_parameters.h"
#include "reconst.h"
#include "util.h"

namespace {

//! a fresh, empty directory for the files a test writes
std::string makeTempDir() {
    const char *base = std::getenv("TMPDIR");
    std::string tmpl = std::string(base ? base : "/tmp") + "/music_ut_XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    REQUIRE(mkdtemp(buf.data()) != nullptr);
    return std::string(buf.data());
}

std::string writeFile(
    const std::string &dir, const std::string &name,
    const std::string &content) {
    const std::string path = dir + "/" + name;
    std::ofstream f(path.c_str());
    f << content;
    return path;
}

bool fileExists(const std::string &path) {
    struct stat st;
    return (stat(path.c_str(), &st) == 0);
}

//! MUSIC's default parameters, read from a minimal input file
InitData defaultParameters(const std::string &extraLines = "") {
    const std::string dir = makeTempDir();
    const std::string input = writeFile(
        dir, "music_input", "echo_level 0\n" + extraLines + "EndOfData\n");
    return ReadInParameters::read_in_parameters(input);
}

//! run f() inside directory dir, then go back
template <typename F>
void inDirectory(const std::string &dir, F f) {
    char cwd[4096];
    REQUIRE(getcwd(cwd, sizeof(cwd)) != nullptr);
    REQUIRE(chdir(dir.c_str()) == 0);
    f();
    REQUIRE(chdir(cwd) == 0);
}

}  // namespace

// ---------------------------------------------------------------------------
// item 2: in-memory hydro history (HydroinfoMUSIC)
// ---------------------------------------------------------------------------

namespace {
//! 3x3x3 grid, dx = dy = deta = 1, two stored time slices tau = 1, 1.1
InitData hydroInfoParameters() {
    InitData DATA = defaultParameters();
    DATA.boost_invariant = false;
    DATA.tau0 = 1.0;
    DATA.delta_tau = 0.1;
    DATA.output_evolution_every_N_timesteps = 1;
    DATA.nx = 3;
    DATA.ny = 3;
    DATA.neta = 3;
    DATA.delta_x = 1.;
    DATA.delta_y = 1.;
    DATA.delta_eta = 1.;
    DATA.x_size = 2.;
    DATA.y_size = 2.;
    DATA.eta_size = 2.;
    DATA.output_evolution_every_N_x = 1;
    DATA.output_evolution_every_N_y = 1;
    DATA.output_evolution_every_N_eta = 1;
    return DATA;
}
}  // namespace

TEST_CASE("HydroinfoMUSIC: stored u^mu give the right Cartesian velocity") {
    // item 2: grid_info stored v instead of u^x, u^y, tau u^eta
    InitData DATA = hydroInfoParameters();
    EOS eos(0);
    Cell_info grid_info(DATA, eos);
    HydroinfoMUSIC hydro_info;
    hydro_info.set_grid_infomatioin(DATA);

    const double ux = 0.5, uy = -0.2, ueta = 0.3;  // ueta = tau u^eta
    Fields arena(3, 3, 3);
    for (unsigned int i = 0; i < arena.e_.size(); i++) {
        arena.e_[i] = 1.0 / Util::hbarc;
        arena.u_[1][i] = ux;
        arena.u_[2][i] = uy;
        arena.u_[3][i] = ueta;
        arena.u_[0][i] = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
    }
    grid_info.OutputEvolutionDataXYEta_memory(arena, 1.0, hydro_info);
    grid_info.OutputEvolutionDataXYEta_memory(arena, 1.1, hydro_info);

    const double tau = 1.05, eta = 0.4;
    fluidCell cell;
    hydro_info.getHydroValues(
        0.2, -0.3, tau * sinh(eta), tau * cosh(eta), &cell);

    const double utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
    const double ut = utau * cosh(eta) + ueta * sinh(eta);
    const double uz = utau * sinh(eta) + ueta * cosh(eta);
    CHECK(cell.vx == doctest::Approx(ux / ut).epsilon(1e-5));
    CHECK(cell.vy == doctest::Approx(uy / ut).epsilon(1e-5));
    CHECK(cell.vz == doctest::Approx(uz / ut).epsilon(1e-5));
}

TEST_CASE("HydroinfoMUSIC: eta interpolation has no half-cell offset") {
    // item 2: etafrac was shifted by +0.5 cell
    InitData DATA = hydroInfoParameters();
    HydroinfoMUSIC hydro_info;
    hydro_info.set_grid_infomatioin(DATA);
    // ed = 1 + eta on the eta grid -1, 0, 1 (x outer, y, eta inner)
    for (double tau : {1.0, 1.1}) {
        for (int ix = 0; ix < 3; ix++) {
            for (int iy = 0; iy < 3; iy++) {
                for (int ieta = 0; ieta < 3; ieta++) {
                    const double eta = ieta - 1.;
                    hydro_info.dump_ideal_info_to_memory(
                        tau, eta, (1. + eta) / Util::hbarc, 0., 0., 0., 0., 0.,
                        0.);
                }
            }
        }
    }
    for (double eta : {-0.7, 0.3, 0.6}) {
        fluidCell cell;
        const double tau = 1.0;
        hydro_info.getHydroValues(
            0., 0., tau * sinh(eta), tau * cosh(eta), &cell);
        CHECK(cell.ed == doctest::Approx(1. + eta).epsilon(1e-5));
    }
}

TEST_CASE("HydroinfoMUSIC: tau after the last slice is out of range") {
    // item 2: itau > itaumax allowed itau == itaumax (read past the end)
    InitData DATA = hydroInfoParameters();
    HydroinfoMUSIC hydro_info;
    hydro_info.set_grid_infomatioin(DATA);
    for (double tau : {1.0, 1.1}) {
        for (int i = 0; i < 27; i++) {
            hydro_info.dump_ideal_info_to_memory(
                tau, (i % 3) - 1., 2. / Util::hbarc, 0., 0., 0., 0., 0., 0.);
        }
    }
    fluidCell cell;
    hydro_info.getHydroValues(0., 0., 0., 1.15, &cell);  // last slice
    CHECK(cell.ed == doctest::Approx(2.).epsilon(1e-5));
    hydro_info.getHydroValues(0., 0., 0., 1.25, &cell);  // past the end
    CHECK(cell.ed == 0.);
}

// ---------------------------------------------------------------------------
// item 3: cleanup of old *.dat outputs before a run
// ---------------------------------------------------------------------------

TEST_CASE("removeOldOutputDatFiles keeps the input files") {
    const std::string dir = makeTempDir();
    mkdir((dir + "/initial").c_str(), 0755);
    for (const char *name :
         {"eps_freeze_list_s95p_v1.dat", "myinit.dat", "surface_eps_0.18.dat",
          "evolution_xyeta.dat", "notes.txt", "initial/initial_TA.dat",
          "initial/other.dat"}) {
        writeFile(dir, name, "x\n");
    }
    InitData DATA;
    DATA.initName = "./myinit.dat";
    DATA.initName_TA = "initial/initial_TA.dat";
    DATA.freeze_list_filename = "eps_freeze_list_s95p_v1.dat";
    inDirectory(
        dir, [&] { MUSICRunControl::removeOldOutputDatFiles(DATA); });

    CHECK(fileExists(dir + "/eps_freeze_list_s95p_v1.dat"));
    CHECK(fileExists(dir + "/myinit.dat"));
    CHECK(fileExists(dir + "/notes.txt"));
    CHECK(fileExists(dir + "/initial/initial_TA.dat"));
    CHECK(fileExists(dir + "/initial/other.dat"));
    CHECK_FALSE(fileExists(dir + "/surface_eps_0.18.dat"));
    CHECK_FALSE(fileExists(dir + "/evolution_xyeta.dat"));
}

// ---------------------------------------------------------------------------
// item 4: EoS interpolation across a table seam
// ---------------------------------------------------------------------------

TEST_CASE("EOS_base::interpolate2D is exact for a linear EoS at a seam") {
    // two tables with different nb grids; P(e, nb) = 2 e + 3 nb is linear,
    // so bilinear interpolation must reproduce it exactly
    auto P = [](double e, double nb) { return (2. * e + 3. * nb); };
    EOS_base eos;
    eos.set_number_of_tables(2);
    eos.resize_table_info_arrays();
    const double e0[2] = {0., 4.}, de[2] = {1., 2.};
    const double dnb[2] = {1., 0.5};
    const int Ne[2] = {4, 4}, Nnb[2] = {4, 8};
    eos.pressure_tb = new double **[2];
    eos.temperature_tb = new double **[2];
    for (int t = 0; t < 2; t++) {
        eos.e_bounds[t] = e0[t];
        eos.e_spacing[t] = de[t];
        eos.e_length[t] = Ne[t];
        eos.nb_bounds[t] = 0.;
        eos.nb_spacing[t] = dnb[t];
        eos.nb_length[t] = Nnb[t];
        eos.pressure_tb[t] = Util::mtx_malloc(Nnb[t], Ne[t]);
        eos.temperature_tb[t] = Util::mtx_malloc(Nnb[t], Ne[t]);
        for (int i = 0; i < Nnb[t]; i++) {
            for (int j = 0; j < Ne[t]; j++) {
                eos.pressure_tb[t][i][j] = P(e0[t] + j * de[t], i * dnb[t]);
            }
        }
    }
    // e between the last point of table 0 (e = 3) and table 1 (e = 4)
    for (double nb : {0.2, 1.3, 2.7}) {
        CHECK(
            eos.interpolate2D(3.5, nb, 0, eos.pressure_tb)
            == doctest::Approx(P(3.5, nb)));
    }
}

// ---------------------------------------------------------------------------
// item 5: AMPT net-baryon source
// ---------------------------------------------------------------------------

TEST_CASE("HydroSourceAMPT: rhob source is weighted by baryon number") {
    // one u quark and one anti-u quark at the same point
    const std::string dir = makeTempDir();
    const std::string partons =
        "1 0 2\n"
        "2 0 0 0 0.3 0 0 0 0.55\n"
        "-2 0 0 0 0.3 0 0 0 0.55\n";
    const std::string quarkOnly =
        "1 0 1\n"
        "2 0 0 0 0.3 0 0 0 0.55\n";
    InitData DATA = defaultParameters();
    DATA.delta_tau = 0.1;
    const double tau = 0.5;
    const FlowVec u = {1., 0., 0., 0.};

    DATA.initName_AMPT = writeFile(dir, "pair.dat", partons);
    HydroSourceAMPT pair(DATA);
    pair.prepare_list_for_current_tau_frame(tau);
    CHECK(
        pair.get_hydro_rhob_source(tau, 0., 0., 0., u)
        == doctest::Approx(0.));

    DATA.initName_AMPT = writeFile(dir, "quark.dat", quarkOnly);
    HydroSourceAMPT quark(DATA);
    quark.prepare_list_for_current_tau_frame(tau);
    // sigma_x = 0.5 fm, sigma_eta = 0.2 (constructor defaults)
    const double expected =
        (1. / 3.) / tau
        / (DATA.delta_tau * M_PI * 0.25 * sqrt(M_PI) * 0.2);
    CHECK(
        quark.get_hydro_rhob_source(tau, 0., 0., 0., u)
        == doctest::Approx(expected));
}

// ---------------------------------------------------------------------------
// item 6: thread-safe logger
// ---------------------------------------------------------------------------

TEST_CASE("pretty_ostream: concurrent messages stay whole") {
    std::ostringstream captured;
    std::streambuf *old = std::cout.rdbuf(captured.rdbuf());
    pretty_ostream logger;
    const int nThreads = 8, nMessages = 200;
    std::vector<std::thread> threads;
    for (int t = 0; t < nThreads; t++) {
        threads.emplace_back([&logger, t] {
            for (int i = 0; i < 200; i++) {
                logger << "thread " << t << " msg " << i << " end";
                logger.flush(i % 2 ? "warning" : "info");
            }
        });
    }
    for (auto &th : threads) th.join();
    std::cout.rdbuf(old);

    std::istringstream lines(captured.str());
    std::string line;
    std::set<std::string> seen;
    int nLines = 0;
    while (std::getline(lines, line)) {
        nLines++;
        const auto pos = line.find("thread ");
        REQUIRE(pos != std::string::npos);
        const auto end = line.find(" end", pos);
        REQUIRE(end != std::string::npos);
        seen.insert(line.substr(pos, end - pos));
    }
    CHECK(nLines == nThreads * nMessages);
    CHECK(seen.size() == static_cast<size_t>(nThreads * nMessages));
}

// ---------------------------------------------------------------------------
// items 7, 8: resetDtau and the re-run decision
// ---------------------------------------------------------------------------

TEST_CASE("resetDtau is always set from reset_dtau_use_CFL_condition") {
    // item 7: resetDtau was only set when the initial CFL check fired;
    // Delta_Tau = 0.01 passes that check
    CHECK(defaultParameters("Delta_Tau 0.01\n").resetDtau == true);
    CHECK(
        defaultParameters("Delta_Tau 0.01\nreset_dtau_use_CFL_condition 0\n")
            .resetDtau
        == false);
}

TEST_CASE("reRunEnlargesGrid only for setups that enlarge the grid") {
    // item 8
    InitData DATA;
    DATA.dynamicGridDetermination = false;
    for (int profile : {11, 111}) {
        DATA.Initial_profile = profile;
        CHECK(MUSICRunControl::reRunEnlargesGrid(DATA));
    }
    for (int profile : {13, 131}) {
        DATA.Initial_profile = profile;
        DATA.dynamicGridDetermination = false;
        CHECK_FALSE(MUSICRunControl::reRunEnlargesGrid(DATA));
        DATA.dynamicGridDetermination = true;
        CHECK(MUSICRunControl::reRunEnlargesGrid(DATA));
    }
    for (int profile : {8, 9, 42, 112}) {
        DATA.Initial_profile = profile;
        CHECK_FALSE(MUSICRunControl::reRunEnlargesGrid(DATA));
    }
}

// ---------------------------------------------------------------------------
// item 9: IP-Glasma Tmunu source, tau0 scaling of pi^{x eta}, pi^{y eta}
// ---------------------------------------------------------------------------

TEST_CASE("HydroSourceTmunu: pi^{x eta} is scaled by tau0 like u^eta") {
    const std::string dir = makeTempDir();
    // header: tau etamax nx ny deta dx dy, then nx*ny lines of
    // x y z e utau ux uy ueta pi^{tau tau} ... pi^{eta eta}
    const double e = 2.0, ux = 0.3, uy = 0., uetaIn = 0.2, pixetaIn = 0.05;
    std::ostringstream file;
    file << "# tau_in_fm 0.4 etamax= 1 xmax= 2 ymax= 2 deta= 0 dx= 1 dy= 1\n";
    for (double x : {-0.5, 0.5}) {
        for (double y : {-0.5, 0.5}) {
            file << x << " " << y << " 0 " << e << " 1 " << ux << " " << uy
                 << " " << uetaIn << " 0 0 0 0 0 0 " << pixetaIn
                 << " 0 0 0\n";
        }
    }
    InitData DATA = defaultParameters();
    DATA.initName = writeFile(dir, "ipglasma.dat", file.str());
    DATA.tau0 = 0.4;
    DATA.delta_tau = 0.02;
    DATA.sFactor = 1.;
    DATA.preEqVisFactor = 1.;
    DATA.FlagResumTransportCoeff = false;
    HydroSourceTmunu source(DATA);

    EnergyFlowVec j_mu;
    source.get_hydro_energy_source(0.4, -0.5, -0.5, 0., {1, 0, 0, 0}, j_mu);

    // expected pi^{tau eta} with tau*pi^{x eta}, as in init.cpp
    const double tau0 = 0.4;
    const double ueta = uetaIn * tau0;
    const double utau = sqrt(1. + ux * ux + uy * uy + ueta * ueta);
    const double pixeta = pixetaIn * tau0;
    const double pietaeta =
        2. * ux * ueta * pixeta / (utau * utau - ueta * ueta);
    const double pitaueta = (pixeta * ux + pietaeta * ueta) / utau;
    const double eps = e / Util::hbarc;
    const double expected =
        (4. / 3. * eps * utau * ueta + pitaueta) / DATA.delta_tau;
    CHECK(j_mu[3] == doctest::Approx(expected).epsilon(1e-10));
}

// ---------------------------------------------------------------------------
// item 11: Pi_chem enters T^{mu nu} without bulk viscosity
// ---------------------------------------------------------------------------

TEST_CASE("Diss::MakeWSource couples Pi_chem with turn_on_bulk = 0") {
    InitData DATA = defaultParameters();
    DATA.boost_invariant = true;
    DATA.delta_tau = 0.02;
    EOS eos(0);

    const double Pi = 0.1, ux = 0.4, tau = 1.0;
    Fields curr(1, 1, 1), prev(1, 1, 1);
    for (Fields *f : {&curr, &prev}) {
        f->e_[0] = 1.;
        f->u_[0][0] = sqrt(1. + ux * ux);
        f->u_[1][0] = ux;
        f->piBulkChem_[0] = Pi;
    }

    // uniform, static Pi: d_tau term gives Pi (u^0 u^alpha + g^{0 alpha}),
    // the boost-invariant eta flux gives Pi to alpha = 0
    DATA.turn_on_bulk = 0;
    DATA.turn_on_bulk_chem = 1;
    Diss withChem(eos, DATA);
    TJbVec dwmn = {0.};
    withChem.MakeWSource(tau, 0, 0, 0, dwmn, curr, prev, 0);
    CHECK(dwmn[0] == doctest::Approx(Pi * (ux * ux + 1.)));
    CHECK(dwmn[1] == doctest::Approx(Pi * ux * sqrt(1. + ux * ux)));

    DATA.turn_on_bulk_chem = 0;
    Diss without(eos, DATA);
    dwmn = {0.};
    without.MakeWSource(tau, 0, 0, 0, dwmn, curr, prev, 0);
    CHECK(dwmn[0] == doctest::Approx(0.));
}

// ---------------------------------------------------------------------------
// item 13: bulk source at T = 0
// ---------------------------------------------------------------------------

TEST_CASE("Diss::Make_uPiSource is finite at T = 0") {
    InitData DATA = defaultParameters();
    EOS eos(0);
    Diss diss(eos, DATA);
    Cell_small cell;
    cell.epsilon = 1e-3;
    std::vector<double> thermalVec(13, 0.);
    thermalVec[0] = 1e-3;        // e
    thermalVec[2] = 1e-3 / 3.;   // P
    thermalVec[5] = 0.2;         // cs^2
    thermalVec[6] = 0.;          // T
    VelocityShearVec sigma = {0.};
    const double res = diss.Make_uPiSource(1.0, cell, 0.1, sigma, thermalVec);
    CHECK(std::isfinite(res));
}

// ---------------------------------------------------------------------------
// item 14: tilted TATB source, removable singularity
// ---------------------------------------------------------------------------

namespace {
double plateauQuadrature(double k, double eta_m, bool useCosh) {
    const int n = 20000;  // Simpson
    const double h = 2. * eta_m / n;
    double sum = 0.;
    for (int i = 0; i <= n; i++) {
        const double eta = -eta_m + i * h;
        const double f = exp(k * eta) * (useCosh ? cosh(eta) : sinh(eta));
        sum += f * ((i == 0 || i == n) ? 1. : (i % 2 ? 4. : 2.));
    }
    return (sum * h / 3.);
}
}  // namespace

TEST_CASE("HydroSourceTATB::plateauExpIntegrals matches quadrature") {
    for (double eta_m : {0.5, 2.0}) {
        // k = +-1 is |ln(TA/TB)| = 2 eta_m, where the old closed form blew up
        for (double k : {-3., -1., -1. + 1e-9, -0.3, 0., 0.7, 1., 2.5}) {
            double Ic, Is;
            HydroSourceTATB::plateauExpIntegrals(k, eta_m, Ic, Is);
            CHECK(std::isfinite(Ic));
            CHECK(std::isfinite(Is));
            CHECK(
                Ic
                == doctest::Approx(plateauQuadrature(k, eta_m, true))
                       .epsilon(1e-8));
            CHECK(
                Is
                == doctest::Approx(plateauQuadrature(k, eta_m, false))
                       .epsilon(1e-8));
        }
    }
    double Ic, Is;
    HydroSourceTATB::plateauExpIntegrals(1., 0., Ic, Is);
    CHECK(Ic == 0.);
    CHECK(Is == 0.);
}

// ---------------------------------------------------------------------------
// minor: Rbulk = |Pi|/(e + P) in the viscous regulator
// ---------------------------------------------------------------------------

TEST_CASE("Advance::QuestRevertResummedTransCoeff uses e + P for Rbulk") {
    InitData DATA = defaultParameters();
    EOS eos(0);  // ideal gas, P = e/3
    Advance advance(eos, DATA, nullptr);

    // |Pi|/(e+P) = 0.375 -> 1.5 * R < 1, no regulation
    // (with |Pi|/(e+Pi) = 1 it was cut by 1/1.5)
    Cell_small cell;
    cell.epsilon = 1.;
    cell.pi_b = -0.5;
    advance.QuestRevertResummedTransCoeff(cell);
    CHECK(cell.pi_b == doctest::Approx(-0.5));

    // |Pi|/(e+P) = 0.75 -> scaled by (1 - 1e-6)/1.125
    cell.pi_b = -1.;
    advance.QuestRevertResummedTransCoeff(cell);
    CHECK(cell.pi_b == doctest::Approx(-(1. - 1e-6) / 1.125));
}

// ---------------------------------------------------------------------------
// minor: freeze-out list reader
// ---------------------------------------------------------------------------

TEST_CASE("Evolve reads the last freeze-out entry without a newline") {
    const std::string dir = makeTempDir();
    const std::string list =
        writeFile(dir, "fo.dat", "# e\n0.1 0 0 0 0 0 0\n0.2 0 0 0 0 0 0\n"
                                 "0.3 0 0 0 0 0 0");
    InitData DATA = defaultParameters();
    DATA.freezeOutMethod = 4;
    DATA.useEpsFO = 1;
    DATA.freeze_eps_flag = 1;
    DATA.freeze_list_filename = list;
    DATA.N_freeze_out = 1;  // the list length wins
    EOS eos(0);
    Evolve evolve(eos, DATA, nullptr);
    const auto &eps = evolve.getFreezeOutEnergyDensities();
    REQUIRE(eps.size() == 3);
    CHECK(eps[2] == doctest::Approx(0.3));
}

// ---------------------------------------------------------------------------
// reconstruction reverts: conservation of the KT update
// ---------------------------------------------------------------------------

TEST_CASE("Reconst caps a superluminal T^{0mu} and keeps T00 and J0") {
    // |M| > T00 has no solution; reverting to the old cell did not conserve
    // energy and net baryon number
    EOS eos(0);
    Reconst reconst(eos, 0, 0);
    ReconstCell guess;
    guess.e = 0.1;
    guess.rhob = 0.01;
    guess.u = {1., 0., 0., 0.};

    const double tau = 1.0;
    const TJbVec tauq = {1.0, 0.8, 0.0, 0.9, 0.05};  // |M| = 1.2 > T00
    const ReconstCell sol = reconst.ReconstIt_shell(tau, tauq, guess);

    const double P = eos.get_pressure(sol.e, sol.rhob);
    const double T00 = (sol.e + P) * sol.u[0] * sol.u[0] - P;
    CHECK(T00 == doctest::Approx(1.0).epsilon(1e-10));
    CHECK(sol.rhob * sol.u[0] == doctest::Approx(0.05).epsilon(1e-10));
    // the momentum keeps its direction
    CHECK(sol.u[3] / sol.u[1] == doctest::Approx(0.9 / 0.8).epsilon(1e-10));
    CHECK(sol.e > 0.);
}

TEST_CASE("Advance::MakeDeltaQI conserves P_x and N_B at a fast/rest interface") {
    // the half-way state on the neighbour's side was reconstructed with the
    // own cell as the guess; when its solution was rejected the two
    // neighbours used different interface fluxes
    InitData DATA = defaultParameters("Include_Rhob_Yes_1_No_0 1\n");
    DATA.boost_invariant = true;
    DATA.nx = 12;
    DATA.ny = 1;
    DATA.neta = 1;
    DATA.delta_x = 0.1;
    DATA.delta_y = 0.1;
    DATA.delta_eta = 0.1;
    DATA.delta_tau = 0.01;
    DATA.x_size = DATA.delta_x * (DATA.nx - 1);
    DATA.y_size = 0.;
    DATA.eta_size = 0.;
    EOS eos(0);
    Advance advance(eos, DATA, nullptr);

    Fields arena(DATA.nx, 1, 1);
    for (int ix = 0; ix < DATA.nx; ix++) {
        const int idx = arena.getFieldIdx(ix, 0, 0);
        const bool moving = (ix >= 5 && ix <= 6);
        const double ux = moving ? 300. : 0.;
        arena.e_[idx] = 1.;
        arena.rhob_[idx] = 0.1;
        arena.u_[0][idx] = sqrt(1. + ux * ux);
        arena.u_[1][idx] = ux;
        arena.u_[2][idx] = 0.;
        arena.u_[3][idx] = 0.;
    }

    // cells at rest at both edges: no flux through the boundary, so the
    // updates must sum to zero
    const double tau = 1.0;
    // (tau T^{tau x} and J^tau are pure fluxes in x; the energy row also
    // carries the Milne geometric source)
    double sumPx = 0., sumNB = 0., scalePx = 0., scaleNB = 0.;
    for (int ix = 0; ix < DATA.nx; ix++) {
        // MakeDeltaQI adds the flux update to the cell's own tau T^{tau mu}
        const ReconstCell cell = arena.getCellIdeal(ix, 0, 0);
        const double P = eos.get_pressure(cell.e, cell.rhob);
        TJbVec q0, qi;
        for (int alpha = 0; alpha < 5; alpha++) {
            q0[alpha] = tau * advance.get_TJb(cell, alpha, 0, P);
            qi[alpha] = q0[alpha];
        }
        advance.MakeDeltaQI(tau, arena, ix, 0, 0, qi, 0);
        sumPx += qi[1] - q0[1];
        sumNB += qi[4] - q0[4];
        scalePx = std::max(scalePx, std::abs(qi[1] - q0[1]));
        scaleNB = std::max(scaleNB, std::abs(qi[4] - q0[4]));
    }
    REQUIRE(scalePx > 0.);
    REQUIRE(scaleNB > 0.);
    CHECK(std::abs(sumPx) < 1e-12 * scalePx);
    CHECK(std::abs(sumNB) < 1e-12 * scaleNB);
}
