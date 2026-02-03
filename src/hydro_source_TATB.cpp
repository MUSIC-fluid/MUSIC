// Copyright 2019 Chun Shen

#include "hydro_source_TATB.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include "util.h"

using std::string;

HydroSourceTATB::HydroSourceTATB(InitData &DATA_in) : DATA_(DATA_in) {
    double tau_overlap = 2. * 7. / (sinh(DATA_.beam_rapidity));
    tau_source = std::max(DATA_.tau0, tau_overlap);
    set_source_tau_min(tau_source);
    set_source_tau_max(tau_source);
    set_source_tauStart_max(tau_source);

    TA_ = 0.;
    TB_ = 0.;
    if (DATA_.Initial_profile == 113) {
        read_in_participants_and_compute_TATB();
    } else {
        read_in_TATB();
    }

    double N_B = TA_ + TB_;
    double total_energy = DATA_.ecm / 2. * N_B;
    double total_pz = DATA_.ecm / 2. * (TA_ - TB_);
    music_message << "sqrt{s} = " << DATA_.ecm << " GeV, "
                  << "beam rapidity = " << DATA_.beam_rapidity << ", "
                  << "total energy = " << total_energy << " GeV, "
                  << "net Pz = " << total_pz << " GeV, "
                  << "N_B = " << N_B;
    music_message.flush("info");

    music_message << "HydroSourceTATB: tau_min = " << get_source_tau_min()
                  << " fm/c.";
    music_message.flush("info");
    music_message << "HydroSourceTATB: tau_max = " << get_source_tau_max()
                  << " fm/c.";
    music_message.flush("info");

    ybeam_ = DATA_.beam_rapidity;
    tanhYbeam_ = tanh(ybeam_);
    cosh2Ybeam_ = cosh(2 * ybeam_);

    yL_frac_ = DATA_.yL_frac;
    music_message << "Longitudinal velocity fraction yL_frac = " << yL_frac_;
    music_message.flush("info");

    gridXmin_ = -DATA_.x_size / 2.;
    gridYmin_ = -DATA_.y_size / 2.;
    gridDX_ = DATA_.delta_x;
    gridDY_ = DATA_.delta_y;
    gridDtau_ = DATA_.delta_tau;

    eta0_ = DATA_.eta_flat / 2.;
    eta_m_ = DATA_.eta_m;
    sigma_eta_ = DATA_.eta_fall_off;
    C_eta_ =
        (2. * std::sinh(eta0_)
         + std::sqrt(M_PI / 2.0) * sigma_eta_
               * std::exp(sigma_eta_ * sigma_eta_ / 2.0)
               * (std::exp(eta0_) * std::erfc(-sigma_eta_ / sqrt(2))
                  + std::exp(-eta0_) * std::erfc(sigma_eta_ / sqrt(2))));

    beta_ = DATA_.tilted_fraction;
}

HydroSourceTATB::~HydroSourceTATB() {
    profile_TA.clear();
    profile_TB.clear();
}

//! This function reads in the spatal information of the nuclear thickness
//! functions
void HydroSourceTATB::read_in_TATB() {
    music_message << "read in TA and TB from " << DATA_.initName_TA << " and "
                  << DATA_.initName_TB;
    music_message.flush("info");

    std::ifstream TAfile(DATA_.initName_TA.c_str());
    if (!TAfile) {
        music_message << "hydro_source::read_in_TATB: "
                      << "can not open TA file: " << DATA_.initName_TA;
        music_message.flush("error");
        exit(1);
    }
    std::ifstream TBfile(DATA_.initName_TB.c_str());
    if (!TBfile) {
        music_message << "hydro_source::read_in_TATB: "
                      << "can not open TB file: " << DATA_.initName_TB;
        music_message.flush("error");
        exit(1);
    }

    const int nx = DATA_.nx;
    const int ny = DATA_.ny;
    for (int i = 0; i < nx; i++) {
        std::vector<double> TA_temp;
        std::vector<double> TB_temp;
        for (int j = 0; j < ny; j++) {
            double TA, TB;
            TAfile >> TA;
            TBfile >> TB;
            TA_temp.push_back(TA);
            TB_temp.push_back(TB);
            TA_ += TA;
            TB_ += TB;
        }
        profile_TA.push_back(TA_temp);
        profile_TB.push_back(TB_temp);
    }
    TAfile.close();
    TBfile.close();
    TA_ *= DATA_.delta_x * DATA_.delta_y;
    TB_ *= DATA_.delta_x * DATA_.delta_y;
}

//! This function reads in the spatal information of the participants
//! and compute the nuclear thickness functions
void HydroSourceTATB::read_in_participants_and_compute_TATB() {
    music_message << "read in participants from "
                  << DATA_.initName_participants;
    music_message.flush("info");

    std::ifstream partFile(DATA_.initName_participants.c_str());
    if (!partFile) {
        music_message << "hydro_source::read_in_participants_and_compute_TATB: "
                      << "can not open participant file: "
                      << DATA_.initName_participants;
        music_message.flush("error");
        exit(1);
    }

    // read in participants into a list
    std::vector<participant> partList;
    string strDummy;
    double dummy;
    double x_0, y_0;
    double dir, e;
    std::getline(partFile, strDummy);
    partFile >> dummy >> x_0 >> y_0 >> dummy >> dir >> e;
    while (!partFile.eof()) {
        participant part_i;
        part_i.x = x_0;
        part_i.y = y_0;
        part_i.dir = static_cast<int>(dir);
        part_i.e = e;
        if (part_i.dir == 1)
            TA_++;
        else
            TB_++;
        partList.push_back(part_i);
        partFile >> dummy >> x_0 >> y_0 >> dummy >> dir >> e;
    }
    partFile.close();

    // shift the event center to (0, 0)
    double x_CM = 0.;
    double y_CM = 0.;
    for (auto &part_i : partList) {
        x_CM += part_i.x;
        y_CM += part_i.y;
    }
    x_CM /= static_cast<double>(partList.size());
    y_CM /= static_cast<double>(partList.size());
    double x_max = 0.;
    double y_max = 0.;
    for (auto &part_i : partList) {
        part_i.x -= x_CM;
        part_i.y -= y_CM;
        if (x_max < std::abs(part_i.x)) x_max = std::abs(part_i.x);
        if (y_max < std::abs(part_i.y)) y_max = std::abs(part_i.y);
    }

    // adjust transverse grid size
    double gridOffset = std::max(3., 5. * DATA_.nucleonWidth);
    DATA_.x_size = 2. * (x_max + gridOffset);
    DATA_.y_size = 2. * (y_max + gridOffset);
    DATA_.delta_x = DATA_.x_size / (DATA_.nx - 1);
    DATA_.delta_y = DATA_.y_size / (DATA_.ny - 1);
    music_message << "[HydroSourceTATB] Grid info: x_size = " << DATA_.x_size
                  << ", y_size = " << DATA_.y_size << ", dx = " << DATA_.delta_x
                  << " fm, dy = " << DATA_.delta_y << " fm.";
    music_message.flush("info");

    for (int i = 0; i < DATA_.nx; i++) {
        std::vector<double> TA_temp(DATA_.ny, 0.);
        std::vector<double> TB_temp(DATA_.ny, 0.);
        profile_TA.push_back(TA_temp);
        profile_TB.push_back(TB_temp);
    }

    for (const auto &part_i : partList) {
        computeTATB(part_i.x, part_i.y, part_i.dir);
    }
}

void HydroSourceTATB::computeTATB(
    const double x_0, const double y_0, const int dir) {
    const double wsq = DATA_.nucleonWidth * DATA_.nucleonWidth;
    const double preFact = 1. / (2. * M_PI * wsq);
    for (int i = 0; i < DATA_.nx; i++) {
        double x_local = -DATA_.x_size / 2. + i * DATA_.delta_x;
        for (int j = 0; j < DATA_.ny; j++) {
            double y_local = -DATA_.y_size / 2. + j * DATA_.delta_y;
            double dis_sq =
                ((x_local - x_0) * (x_local - x_0)
                 + (y_local - y_0) * (y_local - y_0));
            if (dis_sq < 25. * wsq) {
                double rholocal = preFact * exp(-dis_sq / (2. * wsq));
                if (dir == 1) {
                    profile_TA[i][j] += rholocal;
                } else {
                    profile_TB[i][j] += rholocal;
                }
            }
        }
    }
}

void HydroSourceTATB::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s,
    const FlowVec &u_mu, EnergyFlowVec &j_mu) const {
    j_mu = {0};
    if (std::abs((tau - tau_source)) > 1. / 2. * gridDtau_) return;

    const int ix = static_cast<int>((x - gridXmin_) / gridDX_ + 0.1);
    const int iy = static_cast<int>((y - gridYmin_) / gridDY_ + 0.1);

    const double TA = profile_TA[ix][iy];
    const double TB = profile_TB[ix][iy];
    /*
    if (DATA_.ecm > 2000.) {
        // at LHC energies, we vary eta_flat according to TA + TB
        double total_num_nucleons = 208. * 2.;  // for PbPb runs
        double slope = -0.4;
        eta_flat = eta_flat + ((TA_ + TB_) / total_num_nucleons - 0.5) * slope;
    }
    */
    double y_CM = atanh((TA - TB) / (TA + TB + Util::small_eps) * tanhYbeam_);
    // double y_L = yL_frac_ * y_CM;
    double y_L = compute_yL(TA, TB, y_CM, eta0_, sigma_eta_, eta_m_);

    /*
    double M_inv =
        ((profile_TA[ix][iy] + profile_TB[ix][iy]) * Util::m_N
         * cosh(DATA_.beam_rapidity) / Util::hbarc);  // [1/fm^3]

    double eta0 =
        std::min(eta_flat / 2.0, std::abs(DATA_.beam_rapidity - (y_CM - y_L)));
    */
    double f_plus = eta_envelope_f(eta_s, eta_m_);
    double f_minus = eta_envelope_f(-eta_s, eta_m_);
    /*
    double eta_envelop =
         eta_profile_plateau(eta_s, eta0, DATA_.eta_fall_off);
     */
    double M =
        TA * TA + TB * TB + 2 * TA * TB * Util::m_N * cosh2Ybeam_;  // [1/fm^4]

    /*
    double E_norm =
        tau_source
        * energy_eta_profile_normalisation(y_CM, eta0, DATA_.eta_fall_off);
    */
    // double eta_envelop = eta_profile_plateau_frag(eta_s - (y_CM - y_L), eta0,
    //                                               DATA_.eta_fall_off);
    // double E_norm = tau_source*energy_eta_profile_normalisation_numerical(
    //                                 y_CM, eta0, DATA_.eta_fall_off);
    // double epsilon = M_inv * eta_envelop / E_norm / dtau;  // [1/fm^5]
    double tilted_epsilon = pow(TA, f_plus) * pow(TB, f_minus);
    double shifted_epsilon =
        eta_profile_plateau(eta_s - (y_CM - y_L), eta0_, sigma_eta_);
    double tilted_norm = tau_source
                         * energy_eta_profile_normalisation_tilted(
                             TA, TB, eta0_, eta_m_, sigma_eta_, y_CM, M, y_L);
    double shifted_norm = tau_source * M / C_eta_;

    double epsilon = (beta_ * tilted_epsilon * tilted_norm
                      + (1. - beta_) * shifted_epsilon * shifted_norm)
                     / gridDtau_;  // [1/fm^5]

    j_mu[0] = epsilon * cosh(y_L);  // [1/fm^5]
    j_mu[3] = epsilon * sinh(y_L);  // [1/fm^5]
}

double HydroSourceTATB::get_hydro_rhob_source(
    const double tau, const double x, const double y, const double eta_s,
    const FlowVec &u_mu) const {
    double res = 0.;
    if (std::abs((tau - tau_source)) > 1. / 2. * gridDtau_) return (res);

    const int ix = static_cast<int>((x - gridXmin_) / gridDX_ + 0.1);
    const int iy = static_cast<int>((y - gridYmin_) / gridDY_ + 0.1);
    const double TA = profile_TA[ix][iy];
    const double TB = profile_TB[ix][iy];
    double eta_rhob_plus = eta_rhob_left_factor(eta_s);
    double eta_rhob_minus = eta_rhob_right_factor(eta_s);
    double norm_B = sqrt(2./M_PI)*1/(tau_source*(DATA_.eta_rhob_width_1 + DATA_.eta_rhob_width_2));
    double norm_B_prime = sqrt(2./M_PI)*(TA + TB)/(tau_source*(DATA_.eta_rhob_width_1 + DATA_.eta_rhob_width_2)*(2*TA*TB + Util::small_eps));
    const double omega = DATA_.omega_rhob;
    /*
    res = 0.5*(
          profile_TA[ix][iy]*(  (1. + DATA_.eta_rhob_asym)*eta_rhob_minus
                              + (1. - DATA_.eta_rhob_asym)*eta_rhob_plus)
        + profile_TB[ix][iy]*(  (1. + DATA_.eta_rhob_asym)*eta_rhob_plus
                              + (1. - DATA_.eta_rhob_asym)*eta_rhob_minus)
    );   // [1/fm^3]
    res /= dtau;  // [1/fm^4]
    */
    res = norm_B * (1 - omega) * (TA * eta_rhob_plus + TB * eta_rhob_minus)
          + norm_B_prime * omega * TA * TB * (eta_rhob_plus + eta_rhob_minus);
    return (res);
}
/*
double HydroSourceTATB::eta_rhob_left_factor(const double eta) const {
    double eta_0       = -std::abs(DATA_.eta_rhob_0);
    double delta_eta_1 = DATA_.eta_rhob_width_1;
    double delta_eta_2 = DATA_.eta_rhob_width_2;
    double norm        = 2./(sqrt(M_PI)*tau_source*(delta_eta_1 + delta_eta_2));
    double exp_arg     = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_1;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_2;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}
*/

double HydroSourceTATB::eta_rhob_left_factor(const double eta) const {
    double eta_0_nB = -std::abs(DATA_.eta_rhob_0);
    double sigma_B_plus = DATA_.eta_rhob_width_1;
    double sigma_B_minus = DATA_.eta_rhob_width_2;

    // double norm        = 1./(sqrt(M_PI)*tau_source*sigma_B_plus);
    double exp_arg_1 = (eta - eta_0_nB) / sigma_B_plus;
    double exp_arg_2 = (eta - eta_0_nB) / sigma_B_minus;

    double res =
        Util::theta(eta - eta_0_nB) * std::exp(-exp_arg_1 * exp_arg_1 / 2.)
        + Util::theta(eta_0_nB - eta) * std::exp(-exp_arg_2 * exp_arg_2 / 2.);
    return (res);
}

/*
double HydroSourceTATB::eta_rhob_right_factor(const double eta) const {
    double eta_0       = std::abs(DATA_.eta_rhob_0);
    double delta_eta_1 = DATA_.eta_rhob_width_1;
    double delta_eta_2 = DATA_.eta_rhob_width_2;
    double norm        = 2./(sqrt(M_PI)*tau_source*(delta_eta_1 + delta_eta_2));
    double exp_arg     = 0.0;
    if (eta < eta_0) {
        exp_arg = (eta - eta_0)/delta_eta_2;
    } else {
        exp_arg = (eta - eta_0)/delta_eta_1;
    }
    double res = norm*exp(-exp_arg*exp_arg);
    return(res);
}
*/

double HydroSourceTATB::eta_rhob_right_factor(const double eta) const {
    double eta_0_nB = std::abs(DATA_.eta_rhob_0);
    double sigma_B_plus = DATA_.eta_rhob_width_1;
    double sigma_B_minus = DATA_.eta_rhob_width_2;

    // double norm        = 1./(sqrt(M_PI)*tau_source*sigma_B_plus);
    double exp_arg_1 = (eta + eta_0_nB) / sigma_B_plus;
    double exp_arg_2 = (eta + eta_0_nB) / sigma_B_minus;

    double res =
        Util::theta(eta + eta_0_nB) * std::exp(-exp_arg_1 * exp_arg_1 / 2.)
        + Util::theta(-eta_0_nB - eta) * std::exp(-exp_arg_2 * exp_arg_2 / 2.);
    return (res);
}

double HydroSourceTATB::eta_profile_plateau(
    const double eta, const double eta_0, const double sigma_eta) const {
    // this function return the eta envelope profile for energy density
    // Hirano's plateau + Gaussian fall-off
    double res;
    double exparg1 = (std::abs(eta) - eta_0) / sigma_eta;
    double exparg = exparg1 * exparg1 / 2.0;
    res = exp(-exparg * Util::theta(exparg1));
    return res;
}

double HydroSourceTATB::energy_eta_profile_normalisation(
    const double y_CM, const double eta_0, const double sigma_eta) const {
    // this function returns the normalization of the eta envelope profile
    // for energy density
    //  \int deta eta_profile_plateau(eta - y_CM, eta_0, sigma_eta)*cosh(eta)
    double f1 =
        (exp(eta_0) * erfc(-sqrt(0.5) * sigma_eta)
         + exp(-eta_0) * erfc(sqrt(0.5) * sigma_eta));
    double f2 = sqrt(M_PI / 2.) * sigma_eta * exp(sigma_eta * sigma_eta / 2.);
    double f3 = sinh(eta_0 + y_CM) - sinh(-eta_0 + y_CM);
    double norm = cosh(y_CM) * f2 * f1 + f3;
    return (norm);
}

double HydroSourceTATB::eta_profile_plateau_frag(
    const double eta, const double eta_0, const double sigma_eta) const {
    // this function return the eta envelope profile for energy density
    // Hirano's plateau + Gaussian fall-off
    double res;
    double exparg1 = (std::abs(eta) - eta_0) / sigma_eta;
    double exparg = exparg1 * exparg1 / 2.0;
    double eta_fragmentation = std::max(eta_0, DATA_.beam_rapidity - 2.0);
    double exparg_frag1 = (std::abs(eta) - eta_fragmentation) / 0.5;
    double exparg_frag = exparg_frag1 * exparg_frag1 / 2.0;
    res =
        (exp(-exparg * Util::theta(exparg1))
         * exp(-exparg_frag * Util::theta(exparg_frag1)));
    return res;
}

double HydroSourceTATB::energy_eta_profile_normalisation_numerical(
    const double y_CM, const double eta_0, const double sigma_eta) const {
    // this function returns the normalization of the eta envelope profile
    // for energy density by performing numerical integral
    //  \int deta eta_profile_plateau(eta - y_CM, eta_0, sigma_eta)*cosh(eta)
    const int npoints = 200;
    const double eta_max = DATA_.beam_rapidity + 3.0;
    const double deta = 2. * eta_max / (npoints - 1);
    double f_eta = 0.;
    for (int i = 0; i < npoints; i++) {
        double eta_i = -eta_max + i * deta + y_CM;
        f_eta += eta_profile_plateau_frag(eta_i - y_CM, eta_0, sigma_eta)
                 * cosh(eta_i);
    }
    f_eta *= deta;
    return (f_eta);
}

// the functions defined for the new energy profile
double HydroSourceTATB::compute_yL(
    const double TA, const double TB, const double ycm, const double eta_0,
    const double sigma_eta, const double eta_m) const {
    // this function computes the longitudinal shift yL given TA and TB

    double yL;

    // ------------------------------
    // 1. Basic ratios
    // ------------------------------
    double a = TB / (TA + Util::small_eps);
    // regularize a ≈ 1
    if (std::abs(a - 1.0) < Util::small_eps) a = 1.0 + Util::small_eps;
    double loga = std::log(std::max(a, Util::small_eps));
    double sqrt_a = std::sqrt(a);

    // ------------------------------
    // 2. Denominator base
    // ------------------------------
    double denom_base = 4.0 * eta_m * eta_m - loga * loga;

    // robust regularization
    denom_base = std::copysign(
        std::max(std::abs(denom_base), Util::small_eps), denom_base);

    double denom = sqrt_a * denom_base;

    // ------------------------------
    // 3. denum
    // ------------------------------
    double den_A = (1.0 + a) * (std::sinh(eta_0) - std::sinh(eta_m));

    double den_B =
        (1.0 + a)
        * (std::exp(0.5 * sigma_eta * sigma_eta) * std::sqrt(M_PI / 2.0)
           * sigma_eta
           * (std::cosh(eta_0)
              + std::erf(sigma_eta / std::sqrt(2.0)) * std::sinh(eta_0)));

    double den_C = sqrt_a
                   * (2.0 * eta_m
                      * (-((a - 1.0) * std::cosh(eta_m) * loga)
                         + 2.0 * (1.0 + a) * eta_m * std::sinh(eta_m)))
                   / denom;

    double denum = den_A + den_B + den_C;

    // enforce positivity (by design)
    denum = std::max(denum, Util::small_eps);

    // ------------------------------
    // 4. num
    // ------------------------------
    double num_A = (1.0 - a) * (std::cosh(eta_0) - std::cosh(eta_m));

    double num_B = (1.0 - a)
                   * (std::exp(0.5 * sigma_eta * sigma_eta)
                      * std::sqrt(M_PI / 2.0) * sigma_eta
                      * (std::cosh(eta_0) * std::erf(sigma_eta / std::sqrt(2.0))
                         + std::sinh(eta_0)));

    double num_C = -sqrt_a
                   * (2.0 * eta_m
                      * (2.0 * (a - 1.0) * eta_m * std::cosh(eta_m)
                         - (1.0 + a) * loga * std::sinh(eta_m)))
                   / denom;

    double num = num_A + num_B + num_C;

    // ------------------------------
    // 5. Solve for yL
    // ------------------------------
    double RHS = num / denum;

    // clip for atanh
    RHS = std::min(0.999999, std::max(-0.999999, RHS));
    yL = ycm - std::atanh(RHS);

    return yL;
}

double HydroSourceTATB::eta_envelope_f(
    const double eta, const double eta_m) const {
    // eta_envelope_f, the envelope function for the forward going nucleus
    double f = (eta >= eta_m)    ? 1.0
               : (eta <= -eta_m) ? 0.0
                                 : 0.5 * (1.0 + eta / eta_m);
    return f;
}

double HydroSourceTATB::energy_eta_profile_normalisation_tilted(
    const double TA, const double TB, const double eta_0, const double eta_m,
    const double sigma_eta, const double ycm, const double M,
    const double yL) const {
    // this function returns the normalization of the eta envelope profile
    // for energy density in the tilted source model

    // a = TA / TB
    double a = TA / (TB + Util::small_eps);

    // protect a ~ 1
    if (std::abs(a - 1.0) < Util::small_eps) a = 1.0 + Util::small_eps;

    double loga = std::log(std::max(a, Util::small_eps));
    double sqrt_a = std::sqrt(a);

    // ------------------------------
    // Term 1
    // ------------------------------
    double term1 = (TA + TB) * (std::sinh(eta_0) - std::sinh(eta_m));

    // ------------------------------
    // Term 2
    // ------------------------------
    double term2 =
        (TA + TB)
        * (std::exp(0.5 * sigma_eta * sigma_eta) * std::sqrt(M_PI / 2.0)
           * sigma_eta
           * (std::cosh(eta_0)
              + std::erf(sigma_eta / std::sqrt(2.0)) * std::sinh(eta_0)));

    // ------------------------------
    // Term 3 denominator
    // ------------------------------
    double denom_base = 4.0 * eta_m * eta_m - loga * loga;

    // robust regularization
    denom_base = std::copysign(
        std::max(std::abs(denom_base), Util::small_eps), denom_base);

    double denom = sqrt_a * denom_base;

    // ------------------------------
    // Term 3
    // ------------------------------
    double term3 = std::sqrt(TA * TB)
                   * (2.0 * eta_m
                      * (-((a - 1.0) * std::cosh(eta_m) * loga)
                         + 2.0 * (1.0 + a) * eta_m * std::sinh(eta_m)))
                   / denom;

    // ------------------------------
    // Combine numerator
    // ------------------------------
    double num = term1 + term2 + term3;

    // protect against tiny / negative num
    num = std::max(num, Util::small_eps);

    // ------------------------------
    // epsilon0
    // ------------------------------
    double epsilon0 = M * std::cosh(ycm - yL) / num;
    return epsilon0;
}
