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
        pressure_tb = new double **[1];
        temperature_tb = new double **[1];
        mu_B_tb = new double **[1];
        pressure_tb[0] = Util::mtx_malloc(Nnb, Ne);
        temperature_tb[0] = Util::mtx_malloc(Nnb, Ne);
        mu_B_tb[0] = Util::mtx_malloc(Nnb, Ne);
        for (int j = 0; j < Ne; j++) {
            for (int i = 0; i < Nnb; i++) {
                const double nb = std::min(i, 3 + j) * dnb;
                pressure_tb[0][i][j] =
                    P(1. + j * de, continued ? i * dnb : nb);
                temperature_tb[0][i][j] = 1.;
                mu_B_tb[0][i][j] = mu(nb);
            }
        }
        find_baryon_table_edges();
    }
    ~CappedTableEOS() {
        Util::mtx_free(mu_B_tb[0], Nnb, Ne);
        delete[] mu_B_tb;
    }
    double get_pressure(double e, double rhob) const override {
        return (interpolate2D(e, std::abs(rhob), 0, pressure_tb));
    }
    double get_temperature(double e, double rhob) const override {
        return (interpolate2D(e, std::abs(rhob), 0, temperature_tb));
    }
    double get_muB(double e, double rhob) const override {
        // signed, as in EOS_neos
        return (std::copysign(
            interpolate2D(e, std::abs(rhob), 0, mu_B_tb), rhob));
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
            std::vector<double> thermalVec;
            eos.getThermalVariables(e, nb, thermalVec);
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
        return (a * pow(T, 4)
                + Nc * Nf
                      * (mq * mq * T * T / 6.
                         + pow(mq, 4) / (12. * M_PI * M_PI)));
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
            std::vector<double> thermalVec;
            eos.getThermalVariables(e, n, thermalVec);
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
        CHECK(source.get_hydro_rhob_source(tau, x, 0., 0., {1, 0, 0, 0})
              == 0.);
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
    CHECK(h * std::sqrt(sol.u[1] * sol.u[1] + sol.u[3] * sol.u[3])
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
        U_derivative uder(DATA, eos);
        const double tau = 1.0;
        const int ix = 5;
        const int idx = arena.getFieldIdx(ix, 0, 0);
        uder.MakedU(tau, arena, arena, idx, ix, 0, 0);
        DumuVec a;
        uder.calculate_Du_supmu(tau, arena, idx, a);
        const double e = arena.e_[idx];
        const double P = eos.get_pressure(e, 0.);
        const double dPdx = (eos.get_pressure(e + 0.1, 0.)
                             - eos.get_pressure(e - 0.1, 0.))
                            / (2. * DATA.delta_x);
        const double expected = (config.ideal ? -dPdx / (e + P) : 0.);
        CHECK(a[1] == doctest::Approx(expected).epsilon(1e-6));
        CHECK(a[2] == doctest::Approx(0.));
        CHECK(a[3] == doctest::Approx(0.));
        CHECK(uder.calculate_expansion_rate(tau, arena, idx)
              == doctest::Approx(1. / tau).epsilon(1e-10));
    }
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
        for (int idx = 0; idx < f.nX(); idx++) {
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
    for (auto &e : arena.e_) e = 0.;
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
        for (auto &e : arena.e_) e = 0.;
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
    for (int idx = 0; idx < static_cast<int>(f.e_.size()); idx++) {
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
        out[4] += dx * dy * (f.rhob_[idx] * f.u_[0][idx] + f.Wmunu_[10][idx]);
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
