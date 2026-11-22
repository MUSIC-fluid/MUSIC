// Unit tests for the bug fixes on chun_dev (Nov 2026).
// Each test case names the fix it protects; most of them fail on the code
// before the fix.

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <array>
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
#include "eos_idealgas.h"
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
#include "u_derivative.h"
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
    for (int i = 0; i < arena.nStorage(); i++) {
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
    inDirectory(dir, [&] { MUSICRunControl::removeOldOutputDatFiles(DATA); });

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
    eos.pressure_tb.resize(2);
    eos.temperature_tb.resize(2);
    for (int t = 0; t < 2; t++) {
        eos.e_bounds[t] = e0[t];
        eos.e_spacing[t] = de[t];
        eos.e_length[t] = Ne[t];
        eos.nb_bounds[t] = 0.;
        eos.nb_spacing[t] = dnb[t];
        eos.nb_length[t] = Nnb[t];
        eos.pressure_tb.allocate(t, Nnb[t], Ne[t]);
        eos.temperature_tb.allocate(t, Nnb[t], Ne[t]);
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

namespace {
// a one-table (e, n_B) EoS whose columns stop at n_edge = (3 + j) 0.1,
// where mu_B reaches its cap. Beyond it the rows repeat the edge state (as
// neos_2 does) or, with continued = true, only mu_B stays capped while P
// keeps changing (as neos_bqs does)
class CappedTableEOS : public EOS_base {
  public:
    static constexpr int Ne = 4, Nnb = 10;
    static constexpr double de = 1., dnb = 0.1;
    static double P(double e, double nb) { return (e / 3. + nb); }
    static double mu(double nb) { return (10. * nb); }
    explicit CappedTableEOS(const bool continued = false) {
        set_flag_muB(true);
        set_flag_muS(false);
        set_flag_muQ(false);
        set_number_of_tables(1);
        resize_table_info_arrays();
        e_bounds[0] = 1.;
        e_spacing[0] = de;
        e_length[0] = Ne;
        nb_bounds[0] = 0.;
        nb_spacing[0] = dnb;
        nb_length[0] = Nnb;
        pressure_tb.resize(1);
        temperature_tb.resize(1);
        mu_B_tb.resize(1);
        pressure_tb.allocate(0, Nnb, Ne);
        temperature_tb.allocate(0, Nnb, Ne);
        mu_B_tb.allocate(0, Nnb, Ne);
        for (int j = 0; j < Ne; j++) {
            for (int i = 0; i < Nnb; i++) {
                const double nb = std::min(i, 3 + j) * dnb;
                pressure_tb[0][i][j] = P(1. + j * de, continued ? i * dnb : nb);
                temperature_tb[0][i][j] = 1.;
                mu_B_tb[0][i][j] = mu(nb);
            }
        }
        find_baryon_table_edges();
    }
    double get_pressure(
        double e, double rhob, double rhoq = 0.,
        double rhos = 0.) const override {
        return (interpolate2D(e, std::abs(rhob), 0, pressure_tb));
    }
    double get_temperature(
        double e, double rhob, double rhoq = 0.,
        double rhos = 0.) const override {
        return (interpolate2D(e, std::abs(rhob), 0, temperature_tb));
    }
    double get_muB(double e, double rhob, double rhoq = 0., double rhos = 0.)
        const override {
        // signed, as in EOS_neos
        return (
            std::copysign(interpolate2D(e, std::abs(rhob), 0, mu_B_tb), rhob));
    }
};
}  // namespace

TEST_CASE("EOS_base: s is evaluated at the edge of an (e, n_B) table") {
    CappedTableEOS eos;
    using C = CappedTableEOS;
    for (int j = 0; j < C::Ne; j++) {
        CHECK(eos.nb_edge[0][j] == doctest::Approx((3 + j) * C::dnb));
    }
    CHECK(eos.get_rhob_edge(1.5) == doctest::Approx(0.35));
    for (int j = 0; j < C::Ne; j++) {
        const double e = 1. + j * C::de;
        const double nEdge = (3 + j) * C::dnb;
        // inside the table: s = (e + P - mu n) / T
        const double nIn = 0.5 * nEdge;
        CHECK(
            eos.get_entropy(e, nIn)
            == doctest::Approx(e + C::P(e, nIn) - C::mu(nIn) * nIn));
        // beyond the edge T, P and mu are those of (e, n_edge), and so is s;
        // with the actual n it would be negative (at n = 0.9 for j = 0)
        const double sEdge = e + C::P(e, nEdge) - C::mu(nEdge) * nEdge;
        for (double nb : {nEdge + 0.05, 0.9, 1.5, -0.9}) {
            CHECK(eos.get_entropy(e, nb) == doctest::Approx(sEdge));
            ThermalVec thermalVec;
            eos.getThermalVariables(e, nb, 0., 0., thermalVec);
            CHECK(thermalVec[12] == doctest::Approx(sEdge));
        }
        CHECK(sEdge > 0.);
    }
    CHECK(1. + C::P(1., 0.3) - C::mu(0.3) * 0.9 < 0.);
    // mu_B capped while P keeps changing (neos_bqs): the edge is where the
    // cap starts, and only the mu n term uses it
    CappedTableEOS bqs(true);
    for (int j = 0; j < C::Ne; j++) {
        const double e = 1. + j * C::de;
        const double nEdge = (3 + j) * C::dnb;
        CHECK(bqs.nb_edge[0][j] == doctest::Approx(nEdge));
        const double nb = 0.85;
        CHECK(
            bqs.get_entropy(e, nb)
            == doctest::Approx(e + C::P(e, nb) - C::mu(nEdge) * nEdge));
    }
    // EoS without n_B tables are not touched
    EOS_base plain;
    CHECK(plain.get_rhob_in_table(1., 123.) == 123.);
}

TEST_CASE("EOS_idealgas is a thermodynamically consistent free QGP") {
    EOS_idealgas eos;
    const double Nc = 3., Nf = 2.5;
    const double a = M_PI * M_PI / 90. * (2. * (Nc * Nc - 1.) + 3.5 * Nc * Nf);
    // P(T, mu_B) of the free massless quark-gluon gas
    auto P = [&](const double T, const double muB) {
        const double mq = muB / 3.;
        return (
            a * pow(T, 4)
            + Nc * Nf
                  * (mq * mq * T * T / 6. + pow(mq, 4) / (12. * M_PI * M_PI)));
    };
    int nChecked = 0;
    for (double e : {1e-3, 0.05, 1., 20., 500.}) {         // [1/fm^4]
        for (double n : {0., 1e-3, 0.03, 0.3, 3., 30.}) {  // [1/fm^3]
            const double T = eos.get_temperature(e, n);
            const double mu = eos.get_muB(e, n);
            const double s = eos.get_entropy(e, n);
            const double nIn = eos.get_rhob_in_table(e, n);
            CHECK(s > 0.);
            CHECK(mu / T <= 10. + 1e-12);
            CHECK(nIn <= n);
            // the state at (e, n_B) solves e = 3P and n_B = dP/dmu_B
            const double h = 1e-6 * (std::abs(mu) + T);
            const double nOfState = (P(T, mu + h) - P(T, mu - h)) / (2. * h);
            const double sOfState = (P(T + h, mu) - P(T - h, mu)) / (2. * h);
            CHECK(3. * P(T, mu) == doctest::Approx(e).epsilon(1e-10));
            CHECK(nOfState == doctest::Approx(nIn).epsilon(1e-6));
            CHECK(s == doctest::Approx(sOfState).epsilon(1e-6));
            CHECK(eos.get_pressure(e, n) == doctest::Approx(e / 3.));
            // inverse functions
            CHECK(eos.get_s2e(s, nIn) == doctest::Approx(e).epsilon(1e-9));
            CHECK(
                eos.get_T2e(T * Util::hbarc, nIn)
                == doctest::Approx(e).epsilon(1e-9));
            ThermalVec thermalVec;
            eos.getThermalVariables(e, n, 0., 0., thermalVec);
            CHECK(thermalVec[12] == doctest::Approx(s));
            if (nIn == n) nChecked++;
        }
        // n_B = 0 is the old Stefan-Boltzmann state; mu_B is odd in n_B
        CHECK(
            eos.get_temperature(e, 0.)
            == doctest::Approx(pow(e / (3. * a), 0.25)));
        CHECK(eos.get_muB(e, -0.3) == doctest::Approx(-eos.get_muB(e, 0.3)));
    }
    CHECK(nChecked > 15);  // most points are below the mu_B / T cap
    // the old mu_B = 5 n_B / T^2 with n_B-independent T gave s < 0 here
    const double e = 1., n = 0.3, T0 = pow(e / (3. * a), 0.25);
    CHECK((e + e / 3. - 5. * n / (T0 * T0) * n) / T0 < 0.);
    CHECK(eos.get_entropy(e, n) > 0.);
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
        pair.get_hydro_rhob_source(tau, 0., 0., 0., u) == doctest::Approx(0.));

    DATA.initName_AMPT = writeFile(dir, "quark.dat", quarkOnly);
    HydroSourceAMPT quark(DATA);
    quark.prepare_list_for_current_tau_frame(tau);
    // sigma_x = 0.5 fm, sigma_eta = 0.2 (constructor defaults)
    const double expected =
        (1. / 3.) / tau / (DATA.delta_tau * M_PI * 0.25 * sqrt(M_PI) * 0.2);
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
    for (int profile : {11, 111, 112, 113}) {
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
    for (int profile : {8, 9, 42}) {
        DATA.Initial_profile = profile;
        CHECK_FALSE(MUSICRunControl::reRunEnlargesGrid(DATA));
    }
}

// ---------------------------------------------------------------------------
// item 9: IP-Glasma Tmunu source, tau0 scaling of pi^{x eta}, pi^{y eta}
// ---------------------------------------------------------------------------

TEST_CASE("HydroSourceTATB: zero outside the stored TA, TB profile") {
    // a re-run pads the hydro grid; the cells next to the stored profile
    // must not repeat its edge column (int cast of -0.9 gives 0)
    const std::string dir = makeTempDir();
    std::ostringstream profile;
    for (int i = 0; i < 9; i++) profile << "1.0\n";
    InitData DATA = defaultParameters();
    DATA.Initial_profile = 112;
    DATA.initName_TA = writeFile(dir, "TA.dat", profile.str());
    DATA.initName_TB = writeFile(dir, "TB.dat", profile.str());
    DATA.nx = DATA.ny = 3;
    DATA.x_size = DATA.y_size = 2.;
    DATA.delta_x = DATA.delta_y = 1.;
    DATA.delta_tau = 0.02;
    DATA.ecm = 19.6;
    DATA.beam_rapidity = acosh(DATA.ecm / (2. * 0.938));
    DATA.tau0 = 1.;
    HydroSourceTATB source(DATA);
    const double tau = source.get_source_tau_min();

    EnergyFlowVec j_mu;
    source.get_hydro_energy_source(tau, -1., 0., 0., {1, 0, 0, 0}, j_mu);
    CHECK(j_mu[0] > 0.);  // edge of the profile
    for (double x : {-2., -1.5, 2., 3.}) {
        source.get_hydro_energy_source(tau, x, 0., 0., {1, 0, 0, 0}, j_mu);
        CHECK(j_mu[0] == 0.);
        CHECK(source.get_hydro_rhob_source(tau, x, 0., 0., {1, 0, 0, 0}) == 0.);
    }
    source.get_hydro_energy_source(tau, 0., -2., 0., {1, 0, 0, 0}, j_mu);
    CHECK(j_mu[0] == 0.);
}

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
                 << " " << uetaIn << " 0 0 0 0 0 0 " << pixetaIn << " 0 0 0\n";
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

TEST_CASE(
    "Diss::computeWSourceDivergence and computeWDivergences match the "
    "per-cell MakeWSource and computeWFluxDivergence bit for bit") {
    EOS eos(0);
    for (int bulk = 0; bulk < 2; bulk++) {
        InitData DATA = defaultParameters();
        DATA.boost_invariant = false;
        DATA.delta_tau = 0.02;
        DATA.delta_x = 0.3;
        DATA.delta_y = 0.25;
        DATA.delta_eta = 0.2;
        DATA.turn_on_bulk = bulk;
        DATA.turn_on_bulk_chem = bulk;
        DATA.turn_on_diff = bulk;
        // more than one block of x cells in the fused pass
        const int nx = 37, ny = 4, neta = 3;
        Fields curr(nx, ny, neta), prev(nx, ny, neta);
        unsigned int seed = 12345u;
        auto rnd = [&seed]() {
            seed = seed * 1103515245u + 12345u;
            return ((seed >> 8) & 0xffff) / 65536.0 - 0.5;
        };
        for (Fields *f : {&curr, &prev}) {
            for (int ieta = 0; ieta < neta; ieta++)
                for (int iy = 0; iy < ny; iy++)
                    for (int ix = 0; ix < nx; ix++) {
                        const int I = f->getFieldIdx(ix, iy, ieta);
                        f->e_[I] = 1. + rnd();
                        const double ux = rnd(), uy = rnd(), ue = rnd();
                        f->u_[1][I] = ux;
                        f->u_[2][I] = uy;
                        f->u_[3][I] = ue;
                        f->u_[0][I] = sqrt(1. + ux * ux + uy * uy + ue * ue);
                        for (int k = 0; k < 14; k++) f->Wmunu_[k][I] = rnd();
                        f->piBulk_[I] = rnd();
                        f->piBulkChem_[I] = rnd();
                    }
            f->updateGhosts();
        }
        const double tau = 1.3;
        Diss diss(eos, DATA);
        std::vector<double> wSource;
        diss.computeWSourceDivergence(tau, curr, prev, wSource);
        std::vector<double> wFluxDiv, fusedFlux, fusedSource;
        diss.computeWFluxDivergence(tau, curr, wFluxDiv);
        diss.computeWDivergences(tau, curr, prev, true, fusedFlux, fusedSource);
        // the W components of a run with bulk and diffusion on or off
        const int nComp = (bulk == 1 ? 10 : 5);
        int fusedMismatches = 0;
        for (int ieta = 0; ieta < neta; ieta++)
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++) {
                    const size_t I = curr.getFieldIdx(ix, iy, ieta);
                    const size_t n = curr.nStorage();
                    for (int k = 0; k < Diss::kWSourceComponents; k++) {
                        if (wSource[k * n + I] != fusedSource[k * n + I]) {
                            fusedMismatches++;
                        }
                    }
                    for (int k = 0; k < nComp; k++) {
                        if (wFluxDiv[k * n + I] != fusedFlux[k * n + I]) {
                            fusedMismatches++;
                        }
                    }
                }
        CHECK(fusedMismatches == 0);
        int mismatches = 0;
        for (int ieta = 0; ieta < neta; ieta++)
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++) {
                    const int I = curr.getFieldIdx(ix, iy, ieta);
                    TJbVec ref = {0.}, got = {0.};
                    diss.MakeWSource(tau, ix, iy, ieta, ref, curr, prev, I);
                    diss.MakeWSourceFromDivergence(
                        wSource, curr.nStorage(), I, got);
                    for (int a = 0; a < 5; a++) {
                        if (ref[a] != got[a]) mismatches++;
                    }
                }
        CHECK(mismatches == 0);
    }
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
    ThermalVec thermalVec = {0.};
    thermalVec[0] = 1e-3;       // e
    thermalVec[2] = 1e-3 / 3.;  // P
    thermalVec[5] = 0.2;        // cs^2
    thermalVec[6] = 0.;         // T
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
    const std::string list = writeFile(
        dir, "fo.dat",
        "# e\n0.1 0 0 0 0 0 0\n0.2 0 0 0 0 0 0\n"
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
    CHECK(sol.e > 0.);
    // T^{tau eta} is kept (it enters the Cartesian E and P_z); only the
    // transverse momentum is reduced, to |M| = 0.999 T00
    const double h = (sol.e + P) * sol.u[0];
    CHECK(h * sol.u[3] == doctest::Approx(0.9).epsilon(1e-8));
    const double M = h * std::sqrt(sol.u[1] * sol.u[1] + sol.u[3] * sol.u[3]);
    CHECK(M == doctest::Approx(0.999).epsilon(1e-8));
    CHECK(sol.u[2] == 0.);
}

TEST_CASE("Reconst caps a superluminal T^{tau eta} by scaling all of M") {
    // T^{tau eta} alone above T00: the transverse momentum cannot absorb it
    EOS eos(0);
    Reconst reconst(eos, 0, 0);
    ReconstCell guess;
    guess.e = 0.1;
    guess.rhob = 0.01;
    guess.u = {1., 0., 0., 0.};
    const TJbVec tauq = {1.0, 0.3, 0.0, 1.1, 0.05};
    const ReconstCell sol = reconst.ReconstIt_shell(1.0, tauq, guess);
    const double P = eos.get_pressure(sol.e, sol.rhob);
    const double h = (sol.e + P) * sol.u[0];
    CHECK(h * sol.u[0] - P == doctest::Approx(1.0).epsilon(1e-10));
    CHECK(sol.u[3] / sol.u[1] == doctest::Approx(1.1 / 0.3).epsilon(1e-8));
    CHECK(
        h * std::sqrt(sol.u[1] * sol.u[1] + sol.u[3] * sol.u[3])
        == doctest::Approx(0.999).epsilon(1e-8));
}

TEST_CASE("U_derivative: Du from the ideal equations of motion") {
    // a fluid at rest in a linear pressure gradient starts to accelerate,
    // Du^x = -partial_x P / (e + P); the backward time difference of u
    // (u unchanged since the previous step) gives 0. The ideal form is used
    // only in dilute cells: e ~ 0.3 GeV/fm^3 here, so with the default
    // e_dilute = 0.01 GeV/fm^3 the cell keeps the backward difference.
    struct Config {
        std::string input;
        bool ideal;
    };
    const Config configs[] = {
        {"Du_from_ideal_EoM 0\n", false},
        {"Du_from_ideal_EoM 1\nDu_from_ideal_EoM_e_dilute 10.\n", true},
        {"Du_from_ideal_EoM 1\n", false},
    };
    for (const auto &config : configs) {
        InitData DATA = defaultParameters(config.input);
        DATA.boost_invariant = true;
        DATA.nx = 11;
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
        Fields arena(DATA.nx, 1, 1);
        for (int ix = 0; ix < DATA.nx; ix++) {
            const int idx = arena.getFieldIdx(ix, 0, 0);
            arena.e_[idx] = 1. + 0.1 * ix;
            arena.rhob_[idx] = 0.;
            arena.u_[0][idx] = 1.;
            for (int i = 1; i < 4; i++) arena.u_[i][idx] = 0.;
        }
        arena.updateGhosts();
        U_derivative uder(DATA, eos);
        const double tau = 1.0;
        const int ix = 5;
        const int idx = arena.getFieldIdx(ix, 0, 0);
        uder.MakedU(tau, arena, arena, idx, ix, 0, 0);
        DumuVec a;
        uder.calculate_Du_supmu(tau, arena, idx, a);
        const double e = arena.e_[idx];
        const double P = eos.get_pressure(e, 0.);
        const double dPdx =
            (eos.get_pressure(e + 0.1, 0.) - eos.get_pressure(e - 0.1, 0.))
            / (2. * DATA.delta_x);
        const double expected = (config.ideal ? -dPdx / (e + P) : 0.);
        CHECK(a[1] == doctest::Approx(expected).epsilon(1e-6));
        CHECK(a[2] == doctest::Approx(0.));
        CHECK(a[3] == doctest::Approx(0.));
        CHECK(
            uder.calculate_expansion_rate(tau, arena, idx)
            == doctest::Approx(1. / tau).epsilon(1e-10));
    }
}

TEST_CASE(
    "Advance::MakeDeltaQI conserves P_x and N_B at a fast/rest interface") {
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
    arena.updateGhosts();
    std::vector<double> pressure;
    advance.computePressure(arena, pressure);
    FaceFluxes fluxes;
    advance.computeFaceFluxes(tau, arena, pressure, fluxes);
    for (int ix = 0; ix < DATA.nx; ix++) {
        // MakeDeltaQI adds the flux update to the cell's own tau T^{tau mu}
        const ReconstCell cell = arena.getCellIdeal(ix, 0, 0);
        const double P = eos.get_pressure(cell.e, cell.rhob);
        TJbVec q0, qi;
        for (int alpha = 0; alpha < 5; alpha++) {
            q0[alpha] = tau * advance.get_TJb(cell, alpha, 0, P);
            qi[alpha] = q0[alpha];
        }
        advance.MakeDeltaQI(fluxes, ix, 0, 0, qi);
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

// ---------------------------------------------------------------------------
// PR #38 review: EoS table seams, conservative coarsening, freeze-out
// ---------------------------------------------------------------------------

namespace {
//! two tables with nb spacings 1 and 0.5 and P(e, nb) = e + nb^2; the seam
//! is between e = 3 (last column of table 0) and e = 4 (table 1)
struct SeamEOS {
    EOS_base eos;
    static double P(double e, double nb) { return (e + nb * nb); }
    SeamEOS() {
        eos.set_number_of_tables(2);
        eos.resize_table_info_arrays();
        const double e0[2] = {0., 4.}, de[2] = {1., 2.};
        const double dnb[2] = {1., 0.5};
        const int Ne[2] = {4, 4}, Nnb[2] = {4, 8};
        eos.pressure_tb.resize(2);
        eos.temperature_tb.resize(2);
        for (int t = 0; t < 2; t++) {
            eos.e_bounds[t] = e0[t];
            eos.e_spacing[t] = de[t];
            eos.e_length[t] = Ne[t];
            eos.nb_bounds[t] = 0.;
            eos.nb_spacing[t] = dnb[t];
            eos.nb_length[t] = Nnb[t];
            eos.pressure_tb.allocate(t, Nnb[t], Ne[t]);
            eos.temperature_tb.allocate(t, Nnb[t], Ne[t]);
            for (int i = 0; i < Nnb[t]; i++) {
                for (int j = 0; j < Ne[t]; j++) {
                    eos.pressure_tb[t][i][j] = P(e0[t] + j * de[t], i * dnb[t]);
                }
            }
        }
    }
};
}  // namespace

TEST_CASE("EOS_base::interpolate2D is continuous in e at a table seam") {
    // the next table was interpolated at the current table's nb nodes and
    // then again in nb: at nb = 1.3 it gave 4 + 1.90 just below e = 4 but
    // 4 + 1.75 at it
    SeamEOS s;
    for (double nb : {0.2, 1.3, 2.7}) {
        const double below =
            s.eos.interpolate2D(4. - 1e-10, nb, 0, s.eos.pressure_tb);
        const double at = s.eos.interpolate2D(4., nb, 1, s.eos.pressure_tb);
        CHECK(below == doctest::Approx(at).epsilon(1e-8));
    }
    CHECK(
        s.eos.interpolate2D(4., 1.3, 1, s.eos.pressure_tb)
        == doctest::Approx(4. + 1.75));
}

TEST_CASE("EOS_base::get_rhob_edge interpolates across a table seam") {
    SeamEOS s;
    s.eos.nb_edge = {{0.3, 0.4, 0.5, 0.6}, {1.0, 1.2, 1.4, 1.6}};
    // halfway between e = 3 (edge 0.6) and e = 4 (edge 1.0)
    CHECK(s.eos.get_rhob_edge(3.5) == doctest::Approx(0.8));
    CHECK(
        s.eos.get_rhob_edge(4. - 1e-10)
        == doctest::Approx(s.eos.get_rhob_edge(4.)).epsilon(1e-8));
    // inside a table and past the last table it is unchanged
    CHECK(s.eos.get_rhob_edge(1.5) == doctest::Approx(0.45));
    CHECK(s.eos.get_rhob_edge(20.) == doctest::Approx(1.6));
}

TEST_CASE("EOS_base::check_table_seams accepts only seams the blend handles") {
    // table 0 has columns e = 0, 1, 2, 3; the blend across the seam runs
    // from e = 3 to e = 3 + 1
    SeamEOS s;
    CHECK(s.eos.check_table_seams());  // next table at 4: one spacing
    s.eos.e_bounds[1] = 3.;            // shares the last column (NEoS)
    CHECK(s.eos.check_table_seams());
    s.eos.e_bounds[1] = 2.;  // overlaps by one column (s95p tables 5, 6)
    CHECK(s.eos.check_table_seams());
    s.eos.e_bounds[1] = 4.5;  // 1.5 spacings: the blend overshoots
    CHECK_FALSE(s.eos.check_table_seams());
    s.eos.e_bounds[1] = 3.5;  // 0.5 spacings: the blend stops short
    CHECK_FALSE(s.eos.check_table_seams());
    s.eos.e_bounds[1] = 0.;  // not increasing
    CHECK_FALSE(s.eos.check_table_seams());
    // rounding in an ASCII e spacing, amplified by the column count
    s.eos.e_bounds[1] = 4. + 5e-4;
    CHECK(s.eos.check_table_seams());
    s.eos.e_bounds[1] = 4. + 5e-3;
    CHECK_FALSE(s.eos.check_table_seams());
    // tolerance grows with the column count: 1e-6 * 2000 = 2e-3 spacings
    s.eos.e_length[0] = 2001;
    s.eos.e_bounds[1] = 2000. + 1.5e-3;  // last column at 2000
    CHECK(s.eos.check_table_seams());
    s.eos.e_length[0] = 4;
    // 1D tables (s95p) use interpolate1D, which never blends across seams,
    // but they still need increasing bounds and a positive spacing
    const int nbLen[2] = {s.eos.nb_length[0], s.eos.nb_length[1]};
    s.eos.nb_length[0] = 1;
    s.eos.nb_length[1] = 1;
    s.eos.e_bounds[1] = 4.5;
    CHECK(s.eos.check_table_seams());
    s.eos.e_bounds[1] = 0.;
    CHECK_FALSE(s.eos.check_table_seams());
    s.eos.e_bounds[1] = 4.;
    s.eos.e_spacing[0] = 0.;
    CHECK_FALSE(s.eos.check_table_seams());
    s.eos.e_spacing[0] = 1.;
    // the destructor frees nb_length rows
    s.eos.nb_length[0] = nbLen[0];
    s.eos.nb_length[1] = nbLen[1];
}

TEST_CASE("Evolve::coarseGrainAndEnlargeGrid conserves T^{tau mu}, J^tau") {
    // two cells e = 1, u^x = +-1 with the same rest-frame shear; averaging
    // e and u gave a rest cell with e = 1 (T^{tau tau} 14/3 -> 2) and
    // pi^{tau tau} != 0
    InitData DATA = defaultParameters("Include_Rhob_Yes_1_No_0 1\n");
    DATA.boost_invariant = true;
    DATA.nx = 2;
    DATA.ny = 1;
    DATA.neta = 1;
    DATA.delta_x = 1.;
    DATA.delta_y = 1.;
    DATA.x_size = 1.;
    DATA.y_size = 0.;
    EOS eos(0);
    Evolve evolve(eos, DATA, nullptr);
    Fields arena(2, 1, 1);
    const double a = 0.02;
    const double sq2 = std::sqrt(2.);
    for (int ix = 0; ix < 2; ix++) {
        const int idx = arena.getFieldIdx(ix, 0, 0);
        const double s = (ix == 0 ? 1. : -1.);
        arena.e_[idx] = 1.;
        arena.rhob_[idx] = 0.1;
        arena.u_[0][idx] = sq2;
        arena.u_[1][idx] = s;
        arena.u_[2][idx] = 0.;
        arena.u_[3][idx] = 0.;
        for (int k = 0; k < 14; k++) arena.Wmunu_[k][idx] = 0.;
        // rest-frame pi = diag(0, 2a, -a, -a) boosted along x
        arena.Wmunu_[0][idx] = 2. * a;
        arena.Wmunu_[1][idx] = s * 2. * sq2 * a;
        arena.Wmunu_[4][idx] = 4. * a;
        arena.Wmunu_[7][idx] = -a;
        arena.Wmunu_[9][idx] = -a;
        arena.piBulk_[idx] = 0.;
        arena.piBulkChem_[idx] = 0.;
    }
    auto content = [&](const Fields &f, double dx, double out[5]) {
        for (int m = 0; m < 5; m++) out[m] = 0.;
        for (int ix = 0; ix < f.nX(); ix++) {
            const int idx = f.getFieldIdx(ix, 0, 0);
            const double e = f.e_[idx];
            const double P = eos.get_pressure(e, f.rhob_[idx]);
            for (int m = 0; m < 4; m++) {
                const double g = (m == 0 ? -1. : 0.);
                out[m] += dx
                          * ((e + P) * f.u_[0][idx] * f.u_[m][idx] + P * g
                             + f.Wmunu_[m][idx]);
            }
            out[4] += dx * (f.rhob_[idx] * f.u_[0][idx] + f.Wmunu_[10][idx]);
        }
    };
    double before[5], after[5];
    content(arena, 1., before);
    evolve.regridTransverse({&arena}, 0., true);
    REQUIRE(arena.nX() == 1);
    content(arena, 2., after);
    for (int m = 0; m < 5; m++) {
        CHECK(after[m] == doctest::Approx(before[m]).epsilon(1e-10));
    }
    // the coarse cell is at rest and its shear is transverse and traceless
    CHECK(std::abs(arena.u_[1][0]) < 1e-10);
    const double uLow[4] = {
        -arena.u_[0][0], arena.u_[1][0], arena.u_[2][0], arena.u_[3][0]};
    const int widx[4][4] = {
        {0, 1, 2, 3}, {1, 4, 5, 6}, {2, 5, 7, 8}, {3, 6, 8, 9}};
    for (int n = 0; n < 4; n++) {
        double uPi = 0.;
        for (int m = 0; m < 4; m++)
            uPi += uLow[m] * arena.Wmunu_[widx[m][n]][0];
        CHECK(std::abs(uPi) < 1e-12);
    }
    const double trace = -arena.Wmunu_[0][0] + arena.Wmunu_[4][0]
                         + arena.Wmunu_[7][0] + arena.Wmunu_[9][0];
    CHECK(std::abs(trace) < 1e-12);
}

TEST_CASE("Evolve::hotCellsAtGridEdge only sees the edge band") {
    InitData DATA = defaultParameters();
    DATA.fac_x = 2;
    DATA.fac_y = 2;
    EOS eos(0);
    Evolve evolve(eos, DATA, nullptr);
    Fields arena(11, 11, 1);
    arena.e_.fill(0.);
    // an interior hot cell that a stride-2 finder can miss: no enlargement
    arena.e_[arena.getFieldIdx(5, 5, 0)] = 1.;
    CHECK(evolve.hotCellsRemain(arena, 0.5));
    CHECK_FALSE(evolve.hotCellsAtGridEdge(arena, 0.5));
    arena.e_[arena.getFieldIdx(10, 4, 0)] = 1.;
    CHECK(evolve.hotCellsAtGridEdge(arena, 0.5));
}

TEST_CASE("Evolve::hotCellsAtGridEdge band matches the finder") {
    // fac = 2: the finder starts cubes at even indices and flags the cube
    // at 0 and those starting at >= n - 4
    InitData DATA = defaultParameters();
    DATA.fac_x = 2;
    DATA.fac_y = 2;
    EOS eos(0);
    Evolve evolve(eos, DATA, nullptr);
    auto onlyHot = [&](int n, int ix, int iy) {
        Fields arena(n, n, 1);
        arena.e_.fill(0.);
        arena.e_[arena.getFieldIdx(ix, iy, 0)] = 1.;
        return (evolve.hotCellsAtGridEdge(arena, 0.5));
    };
    // n = 11: cubes at 0, 2, 4, 6, 8; flagged 0 and 8 (8 >= 7)
    CHECK(onlyHot(11, 1, 5));  // inside the cube at 0, between samples
    CHECK(onlyHot(11, 5, 1));
    CHECK(onlyHot(11, 2, 5));  // far corner of the cube at 0
    CHECK(onlyHot(11, 5, 2));
    CHECK(onlyHot(11, 8, 5));  // near corner of the flagged cube at 8
    CHECK(onlyHot(11, 5, 8));
    CHECK_FALSE(onlyHot(11, 3, 5));  // inside the interior cube at 2
    CHECK_FALSE(onlyHot(11, 5, 3));
    CHECK_FALSE(onlyHot(11, 7, 5));  // inside the interior cube at 6
    CHECK_FALSE(onlyHot(11, 5, 7));
    CHECK(onlyHot(11, 9, 5));  // inside the flagged cube at 8
    CHECK(onlyHot(11, 5, 9));
    // n = 10: cubes at 0, 2, 4, 6; flagged 0 and 6 (6 >= 6); index 9 lies
    // past the last sampled corner (8)
    CHECK_FALSE(onlyHot(10, 5, 5));
    CHECK(onlyHot(10, 7, 5));
    CHECK(onlyHot(10, 9, 5));
    CHECK(onlyHot(10, 5, 9));
}

TEST_CASE("EOS_base: scalar and gradient interpolation agree at a seam") {
    // interpolate2D_with_gradients extrapolated the current table across
    // the seam: 5.4 instead of 5.325 at e = 3.5, n_B = 1.3
    // n_B = 3.2 and 3.8 lie beyond the last row of table 0 (n_B = 3), and
    // 3.8 also beyond table 1 (3.5): the clamped columns are flat in n_B
    SeamEOS s;
    for (double e : {2.5, 3.2, 3.5, 3.9}) {
        for (double nb : {0.2, 1.3, 2.7, 3.2, 3.8}) {
            const int t = s.eos.get_table_idx(e);
            double p, dpde, dpdnb;
            s.eos.interpolate2D_with_gradients(
                e, nb, t, s.eos.pressure_tb, p, dpde, dpdnb);
            auto P = [&](double ee, double nn) {
                return (s.eos.interpolate2D(ee, nn, t, s.eos.pressure_tb));
            };
            CHECK(p == doctest::Approx(P(e, nb)));
            const double h = 1e-6;
            CHECK(
                dpde
                == doctest::Approx((P(e + h, nb) - P(e - h, nb)) / (2. * h))
                       .epsilon(1e-5));
            CHECK(
                dpdnb
                == doctest::Approx((P(e, nb + h) - P(e, nb - h)) / (2. * h))
                       .epsilon(1e-5));
        }
    }
    double p, dpde, dpdnb;
    s.eos.interpolate2D_with_gradients(
        3.5, 1.3, 0, s.eos.pressure_tb, p, dpde, dpdnb);
    CHECK(p == doctest::Approx(5.325));
    // the reviewer's case: only the next table varies at n_B = 3.2
    s.eos.interpolate2D_with_gradients(
        3.5, 3.2, 0, s.eos.pressure_tb, p, dpde, dpdnb);
    CHECK(dpdnb == doctest::Approx(3.25));
}

namespace {
//! integrated T^{tau mu} and J^tau of an (nx, ny, 1) arena
void gridContent(
    const Fields &f, const EOS &eos, double dx, double dy, double out[5]) {
    for (int m = 0; m < 5; m++) out[m] = 0.;
    for (int ix = 0; ix < f.nX(); ix++)
        for (int iy = 0; iy < f.nY(); iy++) {
            const int idx = f.getFieldIdx(ix, iy, 0);
            const double e = f.e_[idx];
            const double P = eos.get_pressure(e, f.rhob_[idx]);
            const double Pi = f.piBulk_[idx] + f.piBulkChem_[idx];
            for (int m = 0; m < 4; m++) {
                const double g = (m == 0 ? -1. : 0.);
                const double uu = f.u_[0][idx] * f.u_[m][idx];
                out[m] +=
                    dx * dy
                    * ((e + P) * uu + P * g + f.Wmunu_[m][idx] + Pi * (g + uu));
            }
            out[4] +=
                dx * dy * (f.rhob_[idx] * f.u_[0][idx] + f.Wmunu_[10][idx]);
        }
}
}  // namespace

TEST_CASE("Evolve: coarsening conserves content on odd and singleton grids") {
    // the clamped 1/4-1/2-1/4 stencil turned 5 uniform cells at dx = 1 into
    // 3 cells at dx = 2 with content 6; a singleton y doubled its dy
    for (int ny : {1, 3}) {
        InitData DATA = defaultParameters("Include_Rhob_Yes_1_No_0 1\n");
        DATA.boost_invariant = true;
        DATA.nx = 5;
        DATA.ny = ny;
        DATA.neta = 1;
        DATA.delta_x = 1.;
        DATA.delta_y = 1.;
        DATA.x_size = 4.;
        DATA.y_size = ny - 1.;
        EOS eos(0);
        Evolve evolve(eos, DATA, nullptr);
        Fields arena(5, ny, 1);
        for (int iy = 0; iy < ny; iy++) {
            for (int ix = 0; ix < 5; ix++) {
                const int idx = arena.getFieldIdx(ix, iy, 0);
                const double ux = 0.3 * (ix - 2) + 0.1 * iy;
                const double uy = 0.2 * (iy - 1);
                arena.e_[idx] = 1. + 0.2 * ix + 0.1 * iy;
                arena.rhob_[idx] = 0.1 + 0.01 * ix;
                arena.u_[1][idx] = ux;
                arena.u_[2][idx] = uy;
                arena.u_[3][idx] = 0.;
                arena.u_[0][idx] = std::sqrt(1. + ux * ux + uy * uy);
                for (int k = 0; k < 14; k++) arena.Wmunu_[k][idx] = 0.;
                arena.piBulk_[idx] = 0.;
                arena.piBulkChem_[idx] = 0.;
            }
        }
        double before[5], after[5];
        gridContent(arena, eos, DATA.delta_x, DATA.delta_y, before);
        evolve.regridTransverse({&arena}, 0., true);
        CHECK(arena.nX() == 3);
        CHECK(arena.nY() == (ny == 1 ? 1 : 2));
        CHECK(DATA.delta_x == doctest::Approx(2.));
        CHECK(DATA.delta_y == doctest::Approx(ny == 1 ? 1. : 2.));
        gridContent(arena, eos, DATA.delta_x, DATA.delta_y, after);
        for (int m = 0; m < 5; m++) {
            CHECK(
                after[m]
                == doctest::Approx(before[m]).epsilon(1e-10).scale(1.));
        }
    }
}

TEST_CASE("Evolve: coarsening a dilute, viscous edge gives physical cells") {
    // the 3x3 parents of a coarse cell at the dilute edge of a beastMode 3
    // run (EoS 20, shear + bulk + diffusion): near-empty cells with
    // u^tau ~ 18 next to e ~ 0.08 fm^-4. Remapping the total T^{tau mu},
    // J^tau returned e = 1.2 fm^-4 at u^tau = 19 (the W projected onto the
    // fast u created energy), and the run blew up after the coarsening.
    // {e, rhob, Pi, u^0..u^3, W[0..13]}
    const double parent[9][21] = {
        {1.2840111748573e-05,     -3.4569408242009398e-05,
         -2.1447662744123178e-09, 1.8039998289256136,
         -1.0223545741864848,     -1.0995665308883349,
         0.012647063438288378,    -5.7200689677964134e-09,
         4.0208345662727588e-09,  5.6471924631265748e-09,
         9.3223021552500641e-11,  -1.3736179365961075e-09,
         -5.3214502126694487e-09, -1.6037403617168184e-10,
         -4.3173176503409831e-09, -4.1685855825102454e-12,
         -2.9133380859322945e-11, -8.6553519487779623e-09,
         9.3932059423236949e-09,  5.4643092400318982e-09,
         -2.1307861119866333e-10},
        {1.7069877265224692e-05,  -7.014875194348459e-05,
         -3.7929400738467811e-09, 1.1438619161413104,
         -0.17541349379179649,    -0.52692302866794782,
         -0.0015169883398727184,  -1.1562928089058486e-09,
         2.5021416513469651e-09,  1.6772579765474527e-09,
         -3.6884837985296445e-11, -1.1272494146917565e-09,
         -5.0565723316556987e-09, 3.6006943658317112e-11,
         -1.9578862685518515e-09, 6.2531003057563841e-11,
         1.928842874337759e-09,   1.3812984625407055e-08,
         -1.2508215338798196e-08, -2.5821435117836917e-08,
         -8.4714379756259541e-11},
        {2.0249137960719431e-07,  -3.0731646081933837e-05,
         -5.3435882411590941e-13, 17.751864345805995,
         16.894296618057695,      5.1890391912429052,
         -1.336151860653443,      -2.5770918499751581e-12,
         -1.4189372035076897e-12, -4.1323969460917965e-12,
         2.4927497521075461e-13,  -3.0570492470780972e-13,
         -3.8210414237286941e-12, 1.4712203228469157e-13,
         -1.6421178978848847e-12, 2.1174889415791515e-13,
         -6.2926902738246338e-13, -2.1427526336489578e-10,
         -1.6907672755064276e-10, -1.8519747410248852e-10,
         -1.0211361578269077e-11},
        {8.1949322457810424e-06,  -1.5096502635249381e-05,
         -8.7327171908005453e-10, 1.2604873230577391,
         -0.64511084915170192,    -0.41505520236303534,
         0.019734774592848232,    -5.7497480978650226e-10,
         5.2887165150141074e-10,  9.2481155689464989e-10,
         1.4199770672895757e-11,  -1.8424143974531409e-10,
         -1.3208942528555086e-09, -2.3520452074216269e-11,
         -7.5502675653500418e-10, 1.0754992900278752e-11,
         3.6429338649381614e-10,  -1.5596055733888732e-09,
         1.2782548349805099e-09,  2.7466873437541157e-09,
         -6.1825935975763402e-11},
        {1.2691560756271292e-07,  -2.0123521160833693e-05,
         -2.0944051321111399e-13, 17.824232279807294,
         10.926194020906056,      13.993927105008254,
         -1.2212881526308041,     -5.0896282972749664e-13,
         -3.2888552895696545e-13, -3.8852299882572439e-13,
         3.3933084934232616e-14,  -5.3509204134859421e-14,
         -3.7690916821293035e-13, 2.4897401378249086e-15,
         -1.9893425610853795e-13, 1.8889889533926446e-14,
         -2.5651936948409921e-13, -3.4037344795889466e-10,
         -2.4283783921281032e-10, -2.4440048577291895e-10,
         -5.3390601385152658e-12},
        {0.0070589725280769712,   -0.0017699681970253769,
         -0.00012626424660334271, 1.3816555211907726,
         -0.66224083532636702,    -0.68586277075652113,
         0.0011455204105877395,   -0.00072879042249979448,
         0.00062596818441582817,  0.00086366444942779597,
         -3.5202702372946641e-05, -8.7139006325394534e-05,
         -0.0011768258318537329,  2.1268950059311201e-05,
         -0.00060345453061741515, 5.0314854516062374e-05,
         -3.8196885556984864e-05, 0.00019154473438996855,
         -0.00017882637206789398, -0.00021320769048625005,
         -7.4232943606509647e-06},
        {5.3741583726050499e-08,  -8.785961146143067e-06,
         -3.7610044361677038e-14, 17.776366151846076,
         10.338290077388271,      14.364785235313363,
         1.3311262049848147,      -9.4857930047591707e-14,
         -6.6853971332981531e-14, -6.9143825282094282e-14,
         -1.3802958233668628e-15, -1.5227010170087745e-14,
         -7.1594504493881577e-14, -1.9229911315097799e-15,
         -3.4397263457888567e-14, 3.8674787459137199e-15,
         -4.523365641961541e-14,  -8.6406239276749051e-11,
         -5.5114329899054732e-11, -6.7280856497887057e-11,
         2.0586836082495403e-13},
        {0.007860044659229069,    -0.0010750417091574644,
         -0.00044494026673855741, 1.5221438742163593,
         -0.80866683639664971,    -0.81422313000031343,
         0.0045402061863042912,   -0.0011449464356661824,
         0.00094097505892206648,  0.0012056289483951188,
         -4.1166268567150532e-05, -0.00025231756987401658,
         -0.0015085155031664924,  -2.1241129882327797e-06,
         -0.00075519686454813509, 7.8301280426543553e-05,
         -0.00013743200124403037, 0.00013084842613087281,
         -0.00011222438556532327, -0.00013316016757724069,
         -9.0238012262817576e-07},
        {0.082312118032461912,    -0.0012546316556891739,
         -0.0040357097940605163,  1.8813238808384398,
         -1.1124574571084593,     -1.140935532132318,
         0.0091682566679006684,   -0.015262009599733216,
         0.011400520124236407,    0.014046871770487005,
         -0.00039422807603059409, -0.0012405741335304685,
         -0.017588427794998075,   7.8814559029658266e-05,
         -0.0060088168575519046,  0.00050882025292087803,
         -0.0080126186086508407,  0.00052589090510674833,
         -0.00039651019372466043, -0.00048064700945221289,
         -1.2745589139941378e-05},
    };
    // P(e, n_B) of the 4D table (EOS/download_Neos4D.sh)
    const char *envPath = std::getenv("HYDROPROGRAMPATH");
    const std::string table =
        std::string(envPath ? envPath : ".") + "/EOS/neos4D/neos4d_p_b.dat";
    if (!std::ifstream(table).good()) {
        MESSAGE("skipped: no 4D EoS table at " << table);
        return;
    }
    InitData DATA = defaultParameters("Include_Rhob_Yes_1_No_0 1\n");
    DATA.boost_invariant = true;
    DATA.nx = 5;
    DATA.ny = 5;
    DATA.neta = 1;
    DATA.delta_x = 0.2;
    DATA.delta_y = 0.2;
    DATA.x_size = 0.8;
    DATA.y_size = 0.8;
    EOS eos(20);
    Evolve evolve(eos, DATA, nullptr);
    Fields arena(5, 5, 1);
    // the coarse cell (1, 1) averages the fine cells 1..3 in x and y; the
    // other fine cells copy their nearest parent
    for (int iy = 0; iy < 5; iy++) {
        for (int ix = 0; ix < 5; ix++) {
            const int a = std::max(0, std::min(2, ix - 1));
            const int b = std::max(0, std::min(2, iy - 1));
            const double *c = parent[3 * a + b];
            const int idx = arena.getFieldIdx(ix, iy, 0);
            arena.e_[idx] = c[0];
            arena.rhob_[idx] = c[1];
            arena.piBulk_[idx] = c[2];
            arena.piBulkChem_[idx] = 0.;
            for (int m = 0; m < 4; m++) arena.u_[m][idx] = c[3 + m];
            for (int k = 0; k < 14; k++) arena.Wmunu_[k][idx] = c[7 + k];
        }
    }
    double u0MaxOld = 1., eMaxOld = 0.;
    for (const auto &c : parent) {
        u0MaxOld = std::max(u0MaxOld, c[3]);
        eMaxOld = std::max(eMaxOld, c[0]);
    }
    evolve.regridTransverse({&arena}, 0., true);
    REQUIRE(arena.nX() == 3);
    const int idx = arena.getFieldIdx(1, 1, 0);
    CHECK(std::isfinite(arena.e_[idx]));
    CHECK(arena.e_[idx] > 0.);
    CHECK(arena.e_[idx] <= eMaxOld);
    CHECK(arena.u_[0][idx] <= u0MaxOld);
    // the dissipative part stays at most of the order of e
    double piMax = 0.;
    for (int k = 4; k < 10; k++)
        piMax = std::max(piMax, std::abs(arena.Wmunu_[k][idx]));
    CHECK(piMax < arena.e_[idx]);
}

TEST_CASE("Evolve: enlarging the grid conserves T^{tau mu}, J^tau") {
    // padding cells were damped copies of the edge cells on top of them,
    // which added their content (3e-4 of the energy per enlargement in a
    // dense beastMode 3 run); the edge cell and its copies now share it
    InitData DATA = defaultParameters("Include_Rhob_Yes_1_No_0 1\n");
    DATA.boost_invariant = true;
    DATA.nx = 3;
    DATA.ny = 3;
    DATA.neta = 1;
    DATA.delta_x = 0.5;
    DATA.delta_y = 0.5;
    DATA.x_size = 1.;
    DATA.y_size = 1.;
    EOS eos(0);
    Evolve evolve(eos, DATA, nullptr);
    Fields arena(3, 3, 1);
    const double a = 0.01;
    for (int iy = 0; iy < 3; iy++) {
        for (int ix = 0; ix < 3; ix++) {
            const int idx = arena.getFieldIdx(ix, iy, 0);
            // a boost along x with rest-frame shear diag(0, 2a, -a, -a)
            const double rap = 0.3 * (ix - 1) + 0.1 * iy;
            const double ch = std::cosh(rap), sh = std::sinh(rap);
            arena.e_[idx] = 1. + 0.2 * ix + 0.1 * iy;
            arena.rhob_[idx] = 0.1 + 0.01 * iy;
            arena.u_[0][idx] = ch;
            arena.u_[1][idx] = sh;
            arena.u_[2][idx] = 0.;
            arena.u_[3][idx] = 0.;
            for (int k = 0; k < 14; k++) arena.Wmunu_[k][idx] = 0.;
            arena.Wmunu_[0][idx] = 2. * a * sh * sh;
            arena.Wmunu_[1][idx] = 2. * a * ch * sh;
            arena.Wmunu_[4][idx] = 2. * a * ch * ch;
            arena.Wmunu_[7][idx] = -a;
            arena.Wmunu_[9][idx] = -a;
            arena.piBulk_[idx] = -0.02;
            arena.piBulkChem_[idx] = 0.;
        }
    }
    double before[5], after[5];
    gridContent(arena, eos, DATA.delta_x, DATA.delta_y, before);
    evolve.regridTransverse({&arena}, 1., false);
    REQUIRE(arena.nX() == 7);
    REQUIRE(arena.nY() == 7);
    gridContent(arena, eos, DATA.delta_x, DATA.delta_y, after);
    for (int m = 0; m < 5; m++) {
        CHECK(after[m] == doctest::Approx(before[m]).epsilon(1e-10).scale(1.));
    }
    // the interior is untouched; along a row the padding decays outwards
    // and keeps the flow of its edge cell
    CHECK(arena.e_[arena.getFieldIdx(3, 3, 0)] == doctest::Approx(1.3));
    for (int ix = 4; ix < 6; ix++) {
        const int idx = arena.getFieldIdx(ix, 3, 0);
        const int idxOut = arena.getFieldIdx(ix + 1, 3, 0);
        CHECK(arena.e_[idxOut] < arena.e_[idx]);
        CHECK(arena.u_[1][idxOut] == doctest::Approx(arena.u_[1][idx]));
    }
    // pi^{mu nu} of every cell is transverse to u and traceless
    const int widx[4][4] = {
        {0, 1, 2, 3}, {1, 4, 5, 6}, {2, 5, 7, 8}, {3, 6, 8, 9}};
    for (int idx = 0; idx < 49; idx++) {
        const int i = arena.getFieldIdx(idx % 7, idx / 7, 0);
        const double uLow[4] = {
            -arena.u_[0][i], arena.u_[1][i], arena.u_[2][i], arena.u_[3][i]};
        for (int n = 0; n < 4; n++) {
            double uPi = 0.;
            for (int m = 0; m < 4; m++)
                uPi += uLow[m] * arena.Wmunu_[widx[m][n]][i];
            CHECK(std::abs(uPi) < 1e-12);
        }
        CHECK(std::isfinite(arena.e_[i]));
        CHECK(arena.e_[i] > 0.);
    }
}

TEST_CASE("the 4D EoS pressure gradients match finite differences of P") {
    // dP/de at fixed n_B, n_Q, n_S needs the mu-tilde term: mu-tilde ~
    // n / sqrt(e) changes with e. Without it dP/de was biased by ~-8% at
    // low Ttilde and high density. Needs the 4D table
    // (EOS/download_Neos4D.sh).
    const char *envPath = std::getenv("HYDROPROGRAMPATH");
    const std::string table =
        std::string(envPath ? envPath : ".") + "/EOS/neos4D/neos4d_p_b.dat";
    if (!std::ifstream(table).good()) {
        MESSAGE("skipped: no 4D EoS table at " << table);
        return;
    }
    EOS eos(20);
    const double alphaNf = (8 / 45.0 + 7 / 60.0 * 3.0) * M_PI * M_PI;
    const double dTtilde = 0.008 / Util::hbarc;  // the table's Ttilde step
    auto eOf = [&](const double Tt) { return 3. * alphaNf * std::pow(Tt, 4); };

    // dP/de against a central difference over +-1 Ttilde cell (which
    // averages over the cell structure), as the mean over a Ttilde window:
    // the T^5-shaped slope scatters by a few % around it within each cell
    for (const auto &window :
         {std::array<double, 2> {0.06, 0.12},
          std::array<double, 2> {0.12, 0.20}}) {
        const int nPts = 40;
        double sum = 0., sum2 = 0.;
        for (int i = 0; i < nPts; i++) {
            const double Tt =
                (window[0] + (window[1] - window[0]) * (i + 0.5) / nPts)
                / Util::hbarc;
            const double nb = 0.3 * Tt * Tt * Tt, nq = 0.4 * nb;
            double p, dpde, dpdnb, dpdnq, dpdns;
            eos.get_pressure_with_gradients(
                eOf(Tt), nb, nq, 0., p, dpde, dpdnb, dpdnq, dpdns);
            const double eL = eOf(Tt - dTtilde), eR = eOf(Tt + dTtilde);
            const double ref = (eos.get_pressure(eR, nb, nq, 0.)
                                - eos.get_pressure(eL, nb, nq, 0.))
                               / (eR - eL);
            sum += dpde / ref - 1.;
            sum2 += (dpde / ref - 1.) * (dpde / ref - 1.);
        }
        CAPTURE(window[0]);
        CHECK(std::abs(sum / nPts) < 0.015);
        CHECK(std::sqrt(sum2 / nPts) < 0.06);
    }

    // dP/dn: linear in mu-tilde within a cell, so a 1% match to a small
    // central difference away from the table nodes (8 MeV apart in Ttilde
    // and mu-tilde, with mu-tilde = 0 on a node)
    int nChecked = 0;
    for (const double TtGeV : {0.123, 0.171, 0.253, 0.405}) {
        const double Tt = TtGeV / Util::hbarc;
        for (const double c : {0.005, 0.03, 0.1}) {
            const std::array<double, 4> q = {
                eOf(Tt), c * Tt * Tt * Tt, 0.4 * c * Tt * Tt * Tt,
                0.1 * c * Tt * Tt * Tt};
            double p, g[4];
            eos.get_pressure_with_gradients(
                q[0], q[1], q[2], q[3], p, g[0], g[1], g[2], g[3]);
            for (int k = 1; k < 4; k++) {
                const double h = 1e-3 * Tt * Tt * Tt;
                std::array<double, 4> a = q, b = q;
                a[k] -= h;
                b[k] += h;
                const double fd = (eos.get_pressure(b[0], b[1], b[2], b[3])
                                   - eos.get_pressure(a[0], a[1], a[2], a[3]))
                                  / (2. * h);
                // small at small n, so at least 1% of P / n_scale
                const double tol =
                    1e-2 * std::max(std::abs(fd), p / (Tt * Tt * Tt));
                CAPTURE(TtGeV);
                CAPTURE(c);
                CAPTURE(k);
                CHECK(std::abs(g[k] - fd) <= tol);
                nChecked++;
            }
        }
    }
    CHECK(nChecked == 36);
}

// ---------------------------------------------------------------------------
// GPU port: the device kernels (gpu_backend 1, the device code on the CPU)
// ---------------------------------------------------------------------------

namespace {
//! a lumpy state with fast transverse and longitudinal flow, n_B, and a few
//! nearly empty cells
void fillLumpyIdeal(Fields &arena) {
    for (int ieta = 0; ieta < arena.nEta(); ieta++)
        for (int iy = 0; iy < arena.nY(); iy++)
            for (int ix = 0; ix < arena.nX(); ix++) {
                const int idx = arena.getFieldIdx(ix, iy, ieta);
                const double s = std::sin(1.3 * ix + 0.7 * iy + 0.9 * ieta);
                const double c = std::cos(0.8 * ix - 1.1 * iy + 0.5 * ieta);
                arena.e_[idx] = 2. + 1.5 * s * c;
                if ((ix + 2 * iy + ieta) % 11 == 0) arena.e_[idx] = 1e-7;
                arena.rhob_[idx] = 0.3 + 0.2 * c;
                arena.rhoq_[idx] = 0.;
                arena.rhos_[idx] = 0.;
                arena.u_[1][idx] = 1.6 * s;
                arena.u_[2][idx] = -1.2 * c;
                arena.u_[3][idx] = (arena.nEta() > 1 ? 0.4 * s * c : 0.);
                arena.u_[0][idx] = std::sqrt(
                    1. + arena.u_[1][idx] * arena.u_[1][idx]
                    + arena.u_[2][idx] * arena.u_[2][idx]
                    + arena.u_[3][idx] * arena.u_[3][idx]);
                arena.piBulk_[idx] = 0.;
                arena.piBulkChem_[idx] = 0.;
                for (int k = 0; k < 14; k++) arena.Wmunu_[k][idx] = 0.;
            }
    arena.updateGhosts();
}

//! steps Heun steps of the ideal update with gpu_backend = backend
std::vector<double> runIdealSteps(
    InitData DATA, int backend, int gpuFloat, int steps) {
    DATA.gpu_backend = backend;
    DATA.gpu_float = gpuFloat;
    EOS eos(0);
    Evolve evolve(eos, DATA, nullptr);
    Fields a(DATA.nx, DATA.ny, DATA.neta), b(a), c(a);
    fillLumpyIdeal(a);
    b = a;
    c = a;
    Fields *fpPrev = &b, *fpCurr = &a, *fpNext = &c;
    double tau = DATA.tau0;
    for (int it = 0; it < steps; it++) {
        evolve.preparePressure(*fpCurr);
        evolve.AdvanceRK(tau, fpPrev, fpCurr, fpNext);
        tau += DATA.delta_tau;
    }
    std::vector<double> out;
    for (int ieta = 0; ieta < DATA.neta; ieta++)
        for (int iy = 0; iy < DATA.ny; iy++)
            for (int ix = 0; ix < DATA.nx; ix++) {
                const int idx = fpCurr->getFieldIdx(ix, iy, ieta);
                out.push_back(fpCurr->e_[idx]);
                out.push_back(fpCurr->rhob_[idx]);
                for (int k = 0; k < 4; k++) out.push_back(fpCurr->u_[k][idx]);
            }
    return out;
}

InitData idealDeviceParameters(bool boostInvariant, int mode) {
    InitData DATA = defaultParameters("Include_Rhob_Yes_1_No_0 1\n");
    DATA.viscosity_flag = 0;
    DATA.turn_on_shear = 0;
    DATA.turn_on_bulk = 0;
    DATA.turn_on_diff = 0;
    DATA.boost_invariant = boostInvariant;
    DATA.face_reconstruction = mode;
    DATA.nx = 13;
    DATA.ny = 11;
    DATA.neta = (boostInvariant ? 1 : 7);
    DATA.delta_x = 0.3;
    DATA.delta_y = 0.3;
    DATA.delta_eta = 0.2;
    DATA.x_size = 0.3 * 12;
    DATA.y_size = 0.3 * 10;
    DATA.eta_size = 0.2 * DATA.neta;
    DATA.tau0 = 1.3;
    DATA.delta_tau = 0.02;
    return DATA;
}
}  // namespace

TEST_CASE("Device kernels reproduce the ideal update of Advance") {
    for (int bi = 0; bi < 2; bi++) {
        for (int mode = 0; mode < 2; mode++) {
            CAPTURE(bi);
            CAPTURE(mode);
            const InitData DATA = idealDeviceParameters(bi == 1, mode);
            const auto host = runIdealSteps(DATA, 0, 0, 4);
            const auto dev = runIdealSteps(DATA, 1, 0, 4);
            REQUIRE(host.size() == dev.size());
            double maxDiff = 0.;
            int nDiff = 0;
            for (size_t i = 0; i < host.size(); i++) {
                const double d = std::abs(host[i] - dev[i]);
                if (d > 0.) nDiff++;
                maxDiff = std::max(maxDiff, d);
                REQUIRE(std::isfinite(dev[i]));
            }
            CAPTURE(nDiff);
            CAPTURE(maxDiff);
            CHECK(maxDiff == 0.);
            // float with tau T^{tau mu} carried on the device: rounding
            // errors of 1e-7 per step, amplified by ~(u^tau)^2 in
            // e = T^{tau tau} - v M (6e-4 in a dilute cell with u^tau = 40)
            const auto devF = runIdealSteps(DATA, 1, 1, 4);
            double l1 = 0., norm = 0., maxRel = 0.;
            for (size_t i = 0; i < host.size(); i += 6) {
                l1 += std::abs(devF[i] - host[i]);
                norm += std::abs(host[i]);
                if (host[i] > 1e-3 && host[i + 2] < 5.)
                    maxRel = std::max(maxRel, std::abs(devF[i] / host[i] - 1.));
            }
            CAPTURE(maxRel);
            CHECK(l1 / norm < 1e-5);
            CHECK(maxRel < 1e-4);
        }
    }
}

