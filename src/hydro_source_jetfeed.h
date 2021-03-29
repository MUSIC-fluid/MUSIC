// Copyright 2019 Chun Shen

#ifndef SRC_HYDRO_SOURCE_JETFEED_H_
#define SRC_HYDRO_SOURCE_JETFEED_H_

#include <vector>
#include <memory>
#include <cmath>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "hydro_source_base.h"
#include "util.h"

class SourceProfile {
 private :

    int nbin_eta_;
    int nbin_x_;
    int nbin_y_;

    double eta_min_;
    double eta_max_;
    double delta_eta_;

    double x_min_;
    double x_max_;
    double delta_x_;

    double y_min_;
    double y_max_;
    double delta_y_;

    double tau_src_;

    std::vector<std::vector<double>> lattice_jvec_;

 public :

    SourceProfile() {
        nbin_eta_ = 0;
        nbin_x_ = 0;
        nbin_y_ = 0;

        eta_min_ = 0.;
        eta_max_ = 0.;
        x_min_ = 0.;
        x_max_ = 0.;
        y_min_ = 0.;
        y_max_ = 0.;

        delta_eta_ = 0.;
        delta_x_ = 0.;
        delta_y_ = 0.;

        tau_src_ = 0.;

        lattice_jvec_.clear();
    }

    SourceProfile(SourceProfile const &src) {
        nbin_eta_ = src.nbin_eta_;
        nbin_x_ = src.nbin_x_;
        nbin_y_ = src.nbin_y_;

        eta_min_ = src.eta_min_;
        eta_max_ = src.eta_max_;
        x_min_ = src.x_min_;
        x_max_ = src.x_max_;
        y_min_ = src.y_min_;
        y_max_ = src.y_max_;

        delta_eta_ = src.delta_eta_;
        delta_x_ = src.delta_x_;
        delta_y_ = src.delta_y_;

        tau_src_ = src.tau_src_;

        int nlat_tot_ = (nbin_eta_ + 1) * (nbin_x_ + 1) * (nbin_y_ + 1);

        lattice_jvec_.clear();
        for (int ilat = 0; ilat < nlat_tot_; ilat++) {
            std::vector<double> new_entry;
            new_entry.clear();
            for (int imu = 0; imu < 5; imu++) {
                new_entry.push_back(src.lattice_jvec_.at(ilat).at(imu));
            }
            lattice_jvec_.push_back(new_entry);
        }
    }

    ~SourceProfile() {
        int nlat_tot_ = (nbin_eta_ + 1) * (nbin_x_ + 1) * (nbin_y_ + 1);
        for (int isrc = 0; isrc < nlat_tot_; isrc++) {
            lattice_jvec_.at(isrc).clear();
        }
        lattice_jvec_.clear();
    }

    SourceProfile &operator=(const SourceProfile &src) {
        nbin_eta_ = src.nbin_eta_;
        nbin_x_ = src.nbin_x_;
        nbin_y_ = src.nbin_y_;

        eta_min_ = src.eta_min_;
        eta_max_ = src.eta_max_;
        x_min_ = src.x_min_;
        x_max_ = src.x_max_;
        y_min_ = src.y_min_;
        y_max_ = src.y_max_;

        delta_eta_ = src.delta_eta_;
        delta_x_ = src.delta_x_;
        delta_y_ = src.delta_y_;

        tau_src_ = src.tau_src_;

        int nlat_tot_ = (nbin_eta_ + 1) * (nbin_x_ + 1) * (nbin_y_ + 1);

        lattice_jvec_.clear();
        for (int ilat = 0; ilat < nlat_tot_; ilat++) {
            std::vector<double> new_entry;
            new_entry.clear();
            for (int imu = 0; imu < 5; imu++) {
                new_entry.push_back(src.lattice_jvec_.at(ilat).at(imu));
            }
            lattice_jvec_.push_back(new_entry);
        }

        return *this;
    }

    void set_nbin_eta(int nbin_in) {nbin_eta_ = nbin_in;}
    void set_nbin_x(int nbin_in) {nbin_x_ = nbin_in;}
    void set_nbin_y(int nbin_in) {nbin_y_ = nbin_in;}
    int get_nbin_eta() const {return nbin_eta_;}
    int get_nbin_x() const {return nbin_x_;}
    int get_nbin_y() const {return nbin_y_;}

    void set_range_eta(double eta_min_in, double eta_max_in) {
        eta_min_ = eta_min_in;
        eta_max_ = eta_max_in;
    }
    double get_eta_min() const {return eta_min_;}
    double get_eta_max() const {return eta_max_;}

    void set_range_x(double x_min_in, double x_max_in) {
        x_min_ = x_min_in;
        x_max_ = x_max_in;
    }
    double get_x_min() const {return x_min_;}
    double get_x_max() const {return x_max_;}

    void set_range_y(double y_min_in, double y_max_in) {
        y_min_ = y_min_in;
        y_max_ = y_max_in;
    }
    double get_y_min() const {return y_min_;}
    double get_y_max() const {return y_max_;}

    double get_delta_eta() const {return delta_eta_;}
    double get_delta_x() const {return delta_x_;}
    double get_delta_y() const {return delta_y_;}

    void set_tau_src(double tau_in) {tau_src_ = tau_in;}
    double get_tau_src() const {return tau_src_;}

    void import_from_file(FILE *fin) {
        char dummy[8][10];
        int nx, ny, neta;
        double tau0;

        fscanf(fin, "%s %s %lf %s %d %s %d %s %d %s %lf %lf %s %lf %lf %s %lf %lf",
            dummy[0], dummy[1], &tau0,
            dummy[2], &neta, dummy[3], &nx, dummy[4], &ny,
            dummy[5], &eta_min_, &eta_max_, dummy[6], &x_min_, &x_max_, dummy[7], &y_min_, &y_max_);

        nbin_eta_ = neta - 1;
        nbin_x_ = nx - 1;
        nbin_y_ = ny - 1;

        delta_eta_ = fabs(eta_max_ - eta_min_) / (double)nbin_eta_;
        delta_x_ = fabs(x_max_ - x_min_) / (double)nbin_x_;
        delta_y_ = fabs(y_max_ - y_min_) / (double)nbin_y_;

        tau_src_ = tau0;

        lattice_jvec_.clear();

        for (int ieta = 0; ieta <= nbin_eta_; ieta++) {
            for (int ix = 0; ix <= nbin_x_; ix++) {
                for (int iy = 0; iy <= nbin_y_; iy++) {
                    double eta_now, x_now, y_now;
                    std::vector<double> new_entry;

                    fscanf(fin, "%lf %lf %lf", &eta_now, &x_now, &y_now);
                    new_entry.clear();
                    for (int imu = 0; imu < 5; imu++) {
                        double src_in;
                        fscanf(fin, "%lf", &src_in);
                        if (imu < 4) {
                            src_in = src_in / Util::hbarc;
                        }
                        new_entry.push_back(src_in);
                    }

                    lattice_jvec_.push_back(new_entry);
                }
            }
        }
    }

    void get_source_jvec(const double eta, const double x, const double y, std::vector<double> &jvec_out) const {
        jvec_out.clear();
        for (int imu = 0; imu < 5; imu++) {
            jvec_out.push_back(0.);
        }

        int ieta = (int)floor((eta - eta_min_) / delta_eta_);
        if (ieta < 0 || ieta >= nbin_eta_) {
            return;
        }

        int ix = (int)floor((x - x_min_) / delta_x_);
        if (ix < 0 || ix >= nbin_x_) {
            return;
        }

        int iy = (int)floor((y - y_min_) / delta_y_);
        if (iy < 0 || iy >= nbin_y_) {
            return;
        }

        double frac_eta[2];
        frac_eta[1] = (eta - eta_min_) / delta_eta_ - (double)ieta;
        frac_eta[0] = 1. - frac_eta[1];

        double frac_x[2];
        frac_x[1] = (x - x_min_) / delta_x_ - (double)ix;
        frac_x[0] = 1. - frac_x[1];

        double frac_y[2];
        frac_y[1] = (y - y_min_) / delta_y_ - (double)iy;
        frac_y[0] = 1. - frac_y[1];

        for (int jeta = 0; jeta < 2; jeta++) {
            for (int jx = 0; jx < 2; jx++) {
                for (int jy = 0; jy < 2; jy++) {
                    int ilat =
                        (ieta + jeta) * (nbin_x_ + 1) * (nbin_y_ + 1) +
                        (ix + jx) * (nbin_y_ + 1) + iy + jy;

                    double frac_overall = frac_eta[jeta] * frac_x[jx] * frac_y[jy];
                    for (int imu = 0; imu < 5; imu++) {
                        jvec_out.at(imu) += frac_overall * lattice_jvec_.at(ilat).at(imu);
                    }
                }
            }
        }
    }
};

class HydroSourceJetFeed : public HydroSourceBase {
 private :
    const InitData &DATA;

    std::vector<SourceProfile> list_src_;

 public :
    HydroSourceJetFeed() = default;
    HydroSourceJetFeed(const InitData &DATA_in);
    ~HydroSourceJetFeed();

    void import_from_file();

    //! this function returns the energy source term J^\mu at a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu, EnergyFlowVec &j_mu) const ;

    //! this function returns the net baryon density source term rho
    //! at a given point (tau, x, y, eta_s)
    double get_hydro_rhob_source(const double tau, const double x,
                                 const double y, const double eta_s,
                                 const FlowVec &u_mu) const ;
};

#endif  // SRC_HYDRO_SOURCE_JETFEED_H_