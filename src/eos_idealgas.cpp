// Copyright 2018 @ Chun Shen
#include "eos_idealgas.h"

#include <algorithm>
#include <cmath>

EOS_idealgas::EOS_idealgas() {
    set_EOS_id(0);
    set_number_of_tables(0);
    set_eps_max(1e5);
    Nc = 3.;
    Nf = 2.5;
    set_flag_muB(false);
    set_flag_muS(false);
    set_flag_muQ(false);
    build_x_table();
}

void EOS_idealgas::initialize_eos() {
    music_message.info("initialze EOS ideal gas ...");
}

double EOS_idealgas::coeff_a() const {
    return (M_PI * M_PI / 90. * (2. * (Nc * Nc - 1.) + 7. / 2. * Nc * Nf));
}

double EOS_idealgas::D_of_x(const double x) const {
    const double x2 = x * x;
    return (coeff_a()
            + coeff_b() * (x2 / 54. + x2 * x2 / (972. * M_PI * M_PI)));
}

double EOS_idealgas::N_of_x(const double x) const {
    return (coeff_b() * (x / 27. + x * x * x / (243. * M_PI * M_PI)));
}

double EOS_idealgas::S_of_x(const double x) const {
    return (4. * coeff_a() + coeff_b() * x * x / 27.);
}

double EOS_idealgas::T_of_e_x(const double e, const double x) const {
    return (sqrt(sqrt(std::max(0., e) / (3. * D_of_x(x)))));
}

//! n_B / e^{3/4} = N(x) / (3 D(x))^{3/4} rises monotonically with x up to
//! the cap; Newton in ln x on ln N - 3/4 ln D
double EOS_idealgas::x_of_rho(const double rho, const double x_guess) const {
    const double target = log(rho * r_cap_) + 0.75 * log(3.);
    const double xMax = muB_over_T_max_;
    double x = std::max(1e-300, std::min(xMax, x_guess));
    for (int iter = 0; iter < 60; iter++) {
        const double N = N_of_x(x);
        const double dN = coeff_b() * (1. / 27. + x * x / (81. * M_PI * M_PI));
        const double h = log(N) - 0.75 * log(D_of_x(x)) - target;
        const double dhdu = x * (dN / N - 0.75 * N / D_of_x(x));
        const double du = -h / dhdu;
        x = std::min(xMax, x * exp(std::max(-2., std::min(2., du))));
        if (std::abs(du) < 1e-13) break;
    }
    return (x);
}

void EOS_idealgas::build_x_table() {
    const double xMax = muB_over_T_max_;
    r_cap_ = N_of_x(xMax) / pow(3. * D_of_x(xMax), 0.75);
    const int Nt = 1025;
    x_table_.assign(Nt, 0.);
    for (int i = 1; i < Nt - 1; i++) {
        const double rho = static_cast<double>(i) / (Nt - 1);
        // the small-x solution, then continuation from the previous point
        const double guess =
            (i == 1) ? 27. * rho * r_cap_ * pow(3. * coeff_a(), 0.75)
                           / coeff_b()
                     : x_table_[i - 1];
        x_table_[i] = x_of_rho(rho, guess);
    }
    x_table_[Nt - 1] = xMax;
}

double EOS_idealgas::get_muB_over_T(const double e, const double rhob) const {
    const double n = std::abs(rhob);
    if (n <= 0. || e <= 0.) return (0.);
    // solve N(x)^4 = C D(x)^3, C = 27 n^4 / e^3, without transcendentals
    const double n2 = n * n;
    const double C = 27. * n2 * n2 / (e * e * e);
    const double rc2 = r_cap_ * r_cap_;
    const double rho = sqrt(sqrt(C / (27. * rc2 * rc2)));
    if (rho >= 1.) return (muB_over_T_max_);
    const double pos = rho * (x_table_.size() - 1);
    const int i = std::min(static_cast<int>(pos), int(x_table_.size()) - 2);
    const double f = pos - i;
    double x = x_table_[i] * (1. - f) + x_table_[i + 1] * f;
    if (x <= 0.) return (x_of_rho(rho, x));
    const double b = coeff_b();
    for (int iter = 0; iter < 8; iter++) {
        const double N = N_of_x(x), D = D_of_x(x);
        const double dN = b * (1. / 27. + x * x / (81. * M_PI * M_PI));
        const double N3 = N * N * N, D2 = D * D;
        const double dx =
            -(N3 * N - C * D2 * D) / (4. * N3 * dN - 3. * C * D2 * N);
        x += dx;
        if (std::abs(dx) < 1e-13 * x) {
            return (std::min(x, double(muB_over_T_max_)));
        }
    }
    // not converged from the tabulated guess: the robust solve
    return (x_of_rho(rho, x_table_[i] * (1. - f) + x_table_[i + 1] * f));
}

double EOS_idealgas::get_temperature(double eps, double rhob) const {
    return (T_of_e_x(eps, get_muB_over_T(eps, rhob)));
}

double EOS_idealgas::get_muB(double e, double rhob) const {
    const double x = get_muB_over_T(e, rhob);
    return (std::copysign(x * T_of_e_x(e, x), rhob));  // [1/fm]
}

//! n_B of the state at (e, rhob): |rhob| up to the mu_B / T cap
double EOS_idealgas::get_rhob_in_table(double e, double rhob) const {
    const double x = get_muB_over_T(e, rhob);
    if (x < muB_over_T_max_) return (rhob);
    const double T = T_of_e_x(e, x);
    const double nCap = T * T * T * N_of_x(x);
    return (std::copysign(std::min(std::abs(rhob), nCap), rhob));
}

//! s = dP/dT = T^3 S(x); beyond the cap this is (e + P - mu_B n_cap) / T
double EOS_idealgas::get_entropy(double e, double rhob) const {
    const double x = get_muB_over_T(e, rhob);
    const double T = T_of_e_x(e, x);
    return (std::max(Util::small_eps, T * T * T * S_of_x(x)));
}

void EOS_idealgas::getThermalVariables(
    double e, double rhob, std::vector<double> &thermalVec) const {
    const double x = get_muB_over_T(e, rhob);
    const double T = T_of_e_x(e, x);
    thermalVec = std::vector<double>(
        {e, rhob, e / 3., 1. / 3., 0., 1. / 3., T, std::copysign(x * T, rhob),
         0., 0., 0., 0., std::max(Util::small_eps, T * T * T * S_of_x(x))});
}

//! e at entropy density s and n_B: n_B / s = N(x) / S(x) rises with x
double EOS_idealgas::get_s2e(double s, double rhob) const {
    if (s <= 0.) return (0.);
    const double r = std::abs(rhob) / s;
    double x = 0.;
    if (r > 0.) {
        const double xMax = muB_over_T_max_;
        if (N_of_x(xMax) / S_of_x(xMax) <= r) {
            x = xMax;
        } else {
            x = 27. * 4. * coeff_a() * r / coeff_b();
            x = std::min(x, xMax);
            for (int iter = 0; iter < 50; iter++) {
                const double N = N_of_x(x), S = S_of_x(x);
                const double dN =
                    coeff_b() * (1. / 27. + x * x / (81. * M_PI * M_PI));
                const double dS = 2. * coeff_b() * x / 27.;
                const double dhdu = x * (dN / N - dS / S);
                const double du = -(log(N / S) - log(r)) / dhdu;
                x = std::min(xMax, x * exp(std::max(-2., std::min(2., du))));
                if (std::abs(du) < 1e-13) break;
            }
        }
    }
    const double T = cbrt(s / S_of_x(x));
    return (3. * T * T * T * T * D_of_x(x));  // in 1/fm^4
}

//! e at temperature T and n_B: n_B / T^3 = N(x) rises with x
double EOS_idealgas::get_T2e(double T_in_GeV, double rhob) const {
    const double T = T_in_GeV / Util::hbarc;
    if (T <= 0.) return (0.);
    const double r = std::abs(rhob) / (T * T * T);
    double x = 0.;
    if (r > 0.) {
        const double xMax = muB_over_T_max_;
        if (N_of_x(xMax) <= r) {
            x = xMax;
        } else {
            x = std::min(xMax, 27. * r / coeff_b());
            for (int iter = 0; iter < 50; iter++) {
                const double dN =
                    coeff_b() * (1. / 27. + x * x / (81. * M_PI * M_PI));
                const double dx = -(N_of_x(x) - r) / dN;
                x = std::max(0., std::min(xMax, x + dx));
                if (std::abs(dx) < 1e-13 * (1. + x)) break;
            }
        }
    }
    return (3. * T * T * T * T * D_of_x(x));
}
