// Driver of the (2+1)D ideal KT prototype (see kt2d.h).
//
//   kt2d [--backend cpu|metal] [--precision double|float] [--mode 0|1]
//        [--n 201] [--size 10] [--dtau 0.02] [--tau0 1] [--tau1 1.5]
//        [--q 1] [--init file] [--ref file] [--out file] [--threads N]
//        [--write-init file] [--commit K]
//
// The initial state is the ideal Gubser solution with n_B (q, e0 = 1,
// n_B0 = 0.5, the analytic formula of tests/integration), or --init: a file
// in the format MUSIC reads for Initial_profile 0 (x y e n_B u^x u^y).
// --ref compares the final state with MUSIC's Gubser_flow_check_tau_*.dat.
#include <omp.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "kt2d.h"
#ifdef KT2D_WITH_METAL
#include "kt2d_metal.h"
#endif

namespace {

struct Options {
    std::string backend = "cpu", precision = "double";
    int mode = 1, n = 201, threads = 0, commit = 50;
    double absErr = -1., relErr = -1.;  // < 0: the defaults
    double size = 10., dtau = 0.02, tau0 = 1., tau1 = 1.5, q = 1.;
    std::string init, ref, out, writeInit;
    // precision study: round the face fluxes and/or the stored state to
    // float after each stage of a double run (bit 1 fluxes, bit 2 state)
    int roundToFloat = 0;
    int storeQ = 0;
};

void gubser(double x, double y, double tau, double q, double &e, double &nb,
            double &ux, double &uy) {
    const double r = std::sqrt(x * x + y * y);
    const double den = 1. + 2. * q * q * (tau * tau + r * r)
                       + std::pow(q, 4) * std::pow(tau * tau - r * r, 2);
    e = std::pow(tau, -4. / 3.) * std::pow(2. * q, 8. / 3.)
        / std::pow(den, 4. / 3.);
    nb = 0.5 / std::pow(tau, 3) * 4. * std::pow(q * tau, 2) / den;
    const double ur = std::sinh(
        std::atanh(2. * q * q * tau * r / (1. + q * q * tau * tau + q * q * r * r)));
    ux = (r > 0. ? ur * x / r : 0.);
    uy = (r > 0. ? ur * y / r : 0.);
}

double coord(const Options &o, int i) {
    return -o.size / 2. + i * o.size / (o.n - 1);
}

// state: e, n_B, u^tau, u^x, u^y, each n x n with y fastest
std::vector<double> initialState(const Options &o) {
    const int N = o.n * o.n;
    std::vector<double> s(9 * N);
    std::ifstream f;
    if (!o.init.empty()) {
        f.open(o.init);
        if (!f) {
            std::fprintf(stderr, "cannot read %s\n", o.init.c_str());
            std::exit(1);
        }
    }
    for (int ix = 0; ix < o.n; ix++)
        for (int iy = 0; iy < o.n; iy++) {
            const int i = ix * o.n + iy;
            double e, nb, ux, uy;
            if (f.is_open()) {
                double dum;
                f >> dum >> dum >> e >> nb >> ux >> uy;
            } else {
                gubser(coord(o, ix), coord(o, iy), o.tau0, o.q, e, nb, ux, uy);
            }
            s[i] = e;
            s[N + i] = nb;
            s[2 * N + i] = std::sqrt(1. + ux * ux + uy * uy);
            s[3 * N + i] = ux;
            s[4 * N + i] = uy;
            kt::Cell<double> c{e, nb, s[2 * N + i], ux, uy};
            for (int a = 0; a < 4; a++)
                s[(5 + a) * N + i] = o.tau0 * kt::TJ(c, a, 0, e / 3.);
        }
    return s;
}

template <typename Real>
kt::Params<Real> makeParams(const Options &o) {
    kt::Params<Real> P{};
    P.nx = P.ny = o.n;
    P.mode = o.mode;
    P.eosType = 0;
    P.dx = P.dy = Real(o.size / (o.n - 1));
    P.dtau = Real(o.dtau);
    P.tau = Real(o.tau0);
    P.theta = Real(1.8);
    P.storeQ = o.storeQ;
    // MUSIC's tolerances: beastMode 0 in double; in float those of
    // beastMode > 0 (1e-16 is below the float resolution)
    const bool dbl = sizeof(Real) == 8;
    P.absErr = Real(dbl ? 1e-16 : 1e-8);
    P.relErr = Real(dbl ? 1e-14 : 1e-6);
    if (o.absErr >= 0.) P.absErr = Real(o.absErr);
    if (o.relErr >= 0.) P.relErr = Real(o.relErr);
    return P;
}

template <typename Real>
int runCpu(const Options &o, std::vector<double> &sd, int nSteps,
           double &seconds) {
    kt::Params<Real> P = makeParams<Real>(o);
    const int nx = P.nx, ny = P.ny, n = nx * ny;
    std::vector<Real> s0(sd.begin(), sd.end()), s1(9 * n), s2(9 * n);
    std::vector<Real> fx(4 * (nx + 1) * ny), fy(4 * nx * (ny + 1));
    const Real *table = nullptr;
    int flags = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int step = 0; step < nSteps; step++) {
        P.tau = Real(o.tau0 + step * o.dtau);
        for (int rk = 0; rk < 2; rk++) {
            P.rkFlag = rk;
            const Real *sIn = (rk == 0 ? s0.data() : s1.data());
            Real *sOut = (rk == 0 ? s1.data() : s2.data());
#pragma omp parallel
            {
#pragma omp for schedule(static) nowait
                for (int i = 0; i < (nx + 1) * ny; i++)
                    kt::faceFlux(P, table, sIn, 1, i / ny, i % ny, fx.data(),
                                 &flags);
#pragma omp for schedule(static)
                for (int i = 0; i < nx * (ny + 1); i++)
                    kt::faceFlux(P, table, sIn, 2, i / (ny + 1), i % (ny + 1),
                                 fy.data(), &flags);
                if (o.roundToFloat & 1) {
#pragma omp for schedule(static)
                    for (size_t i = 0; i < fx.size(); i++) fx[i] = Real(float(fx[i]));
#pragma omp for schedule(static)
                    for (size_t i = 0; i < fy.size(); i++) fy[i] = Real(float(fy[i]));
                }
#pragma omp for schedule(static)
                for (int i = 0; i < n; i++)
                    kt::cellUpdate(P, table, sIn, s0.data(), fx.data(),
                                   fy.data(), sOut, i / ny, i % ny, &flags);
            }
        }
        std::swap(s0, s2);
        if (o.roundToFloat & 2) {
            for (auto *v : {&s0, &s1})
                for (auto &x : *v) x = Real(float(x));
        }
    }
    seconds = std::chrono::duration<double>(
                  std::chrono::steady_clock::now() - t0)
                  .count();
    for (int i = 0; i < 9 * n; i++) sd[i] = s0[i];
    return flags;
}

void report(const Options &o, const std::vector<double> &s, double tau) {
    const int N = o.n * o.n;
    // relative L1 errors on y = 0 against the analytic solution, as in
    // tests/integration (gubser_face_reconstruction)
    if (o.init.empty() || true) {
        const int iy0 = (o.n - 1) / 2;
        double de = 0, ae = 0, dn = 0, an = 0, du = 0, au = 0;
        for (int ix = 0; ix < o.n; ix++) {
            const int i = ix * o.n + iy0;
            double e, nb, ux, uy;
            gubser(coord(o, ix), coord(o, iy0), tau, o.q, e, nb, ux, uy);
            de += std::fabs(s[i] - e);
            ae += std::fabs(e);
            dn += std::fabs(s[N + i] - nb);
            an += std::fabs(nb);
            du += std::fabs(s[3 * N + i] - ux);
            au += std::fabs(ux);
        }
        std::printf("L1 error vs analytic on y = 0: e %.3e  n_B %.3e  u^x %.3e\n",
                    de / ae, dn / an, du / au);
    }
    if (!o.ref.empty()) {
        // MUSIC: x y e*hbarc n_B T u^x u^y ... (8 digits)
        const double hbarc = 0.19733;
        std::ifstream f(o.ref);
        std::string line;
        double maxRel[4] = {0, 0, 0, 0};
        double sumAbs[4] = {0, 0, 0, 0}, sumRef[4] = {0, 0, 0, 0};
        for (int i = 0; i < N && std::getline(f, line); i++) {
            std::istringstream ls(line);
            double x, y, eh, nb, T, ux, uy;
            ls >> x >> y >> eh >> nb >> T >> ux >> uy;
            const double ref[4] = {eh / hbarc, nb, ux, uy};
            const double mine[4] = {s[i], s[N + i], s[3 * N + i], s[4 * N + i]};
            for (int k = 0; k < 4; k++) {
                const double d = std::fabs(mine[k] - ref[k]);
                const double scale = std::fabs(ref[k]) + 1e-12;
                if (d / scale > maxRel[k] && std::fabs(ref[k]) > 1e-6)
                    maxRel[k] = d / scale;
                sumAbs[k] += d;
                sumRef[k] += std::fabs(ref[k]);
            }
        }
        std::printf("vs MUSIC (%s):\n", o.ref.c_str());
        const char *names[4] = {"e", "n_B", "u^x", "u^y"};
        for (int k = 0; k < 4; k++)
            std::printf("  %-4s L1 rel diff %.3e   max rel diff %.3e\n",
                        names[k], sumAbs[k] / sumRef[k], maxRel[k]);
    }
}

}  // namespace

int main(int argc, char **argv) {
    Options o;
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string k = argv[i], v = argv[i + 1];
        if (k == "--backend") o.backend = v;
        else if (k == "--precision") o.precision = v;
        else if (k == "--mode") o.mode = std::atoi(v.c_str());
        else if (k == "--n") o.n = std::atoi(v.c_str());
        else if (k == "--size") o.size = std::atof(v.c_str());
        else if (k == "--dtau") o.dtau = std::atof(v.c_str());
        else if (k == "--tau0") o.tau0 = std::atof(v.c_str());
        else if (k == "--tau1") o.tau1 = std::atof(v.c_str());
        else if (k == "--q") o.q = std::atof(v.c_str());
        else if (k == "--init") o.init = v;
        else if (k == "--ref") o.ref = v;
        else if (k == "--out") o.out = v;
        else if (k == "--threads") o.threads = std::atoi(v.c_str());
        else if (k == "--commit") o.commit = std::atoi(v.c_str());
        else if (k == "--abserr") o.absErr = std::atof(v.c_str());
        else if (k == "--relerr") o.relErr = std::atof(v.c_str());
        else if (k == "--store-q") o.storeQ = std::atoi(v.c_str());
        else if (k == "--round-to-float") o.roundToFloat = std::atoi(v.c_str());
        else if (k == "--write-init") o.writeInit = v;
        else {
            std::fprintf(stderr, "unknown option %s\n", k.c_str());
            return 1;
        }
    }
    if (o.threads > 0) omp_set_num_threads(o.threads);
    std::vector<double> s = initialState(o);
    if (!o.writeInit.empty()) {
        // the initial state in MUSIC's Initial_profile 0 format
        FILE *f = std::fopen(o.writeInit.c_str(), "w");
        const int N = o.n * o.n;
        for (int ix = 0; ix < o.n; ix++)
            for (int iy = 0; iy < o.n; iy++) {
                const int i = ix * o.n + iy;
                std::fprintf(f, "%.8e %.8e %.8e %.8e %.8e %.8e\n", coord(o, ix),
                             coord(o, iy), s[i], s[N + i], s[3 * N + i],
                             s[4 * N + i]);
            }
        std::fclose(f);
        return 0;
    }
    const int nSteps = int(std::lround((o.tau1 - o.tau0) / o.dtau));
    const double tauEnd = o.tau0 + nSteps * o.dtau;
    double seconds = 0.;
    int flags = 0;
    std::printf("grid %d x %d, dx %.4g fm, %d steps of %.4g fm, mode %d, %s %s\n",
                o.n, o.n, o.size / (o.n - 1), nSteps, o.dtau, o.mode,
                o.backend.c_str(), o.precision.c_str());
    if (o.backend == "cpu") {
        flags = (o.precision == "float" ? runCpu<float>(o, s, nSteps, seconds)
                                        : runCpu<double>(o, s, nSteps, seconds));
        std::printf("cpu: %d threads, %.4f s, %.2f ns per cell and RK stage\n",
                    omp_get_max_threads(), seconds,
                    1e9 * seconds / (2. * nSteps * o.n * o.n));
    } else if (o.backend == "metal") {
#ifdef KT2D_WITH_METAL
        kt::Params<float> P = makeParams<float>(o);
        std::vector<float> sf(s.begin(), s.end());
        KtMetalTiming t;
        char msg[1024] = {0};
        flags = kt2d_metal_run(&P, nullptr, 0, sf.data(), nSteps, o.commit, &t,
                               msg, sizeof(msg));
        if (flags < 0) {
            std::fprintf(stderr, "Metal error: %s\n", msg);
            return 1;
        }
        for (size_t i = 0; i < s.size(); i++) s[i] = sf[i];
        seconds = t.runSeconds;
        std::printf("metal (%s): compile %.3f s, run %.4f s (GPU %.4f s), "
                    "%.2f ns per cell and RK stage\n",
                    t.device, t.compileSeconds, t.runSeconds, t.gpuSeconds,
                    1e9 * t.gpuSeconds / (2. * nSteps * o.n * o.n));
#else
        std::fprintf(stderr, "built without Metal\n");
        return 1;
#endif
    }
    std::printf("flags %d (1 MaxSpeed sqrt<0, 2 speed>1, 4 revert, 8 "
                "regulate, 16 momentum cap)\n",
                flags);
    report(o, s, tauEnd);
    if (!o.out.empty()) {
        FILE *f = std::fopen(o.out.c_str(), "w");
        const int N = o.n * o.n;
        for (int ix = 0; ix < o.n; ix++)
            for (int iy = 0; iy < o.n; iy++) {
                const int i = ix * o.n + iy;
                std::fprintf(f, "%.9e %.9e %.9e %.9e %.9e %.9e\n", coord(o, ix),
                             coord(o, iy), s[i], s[N + i], s[3 * N + i],
                             s[4 * N + i]);
            }
        std::fclose(f);
    }
    return 0;
}
