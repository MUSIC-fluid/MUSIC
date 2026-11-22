// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_4D_H_
#define SRC_EOS_4D_H_

#include <math.h>

#include <array>
#include <vector>

#include "eos_base.h"

class EOS_4D : public EOS_base {
  private:
    // variables for header infos.
    float mubtilde0, muqtilde0, mustilde0, Ttilde0;
    float dmubtilde, dmuqtilde, dmustilde, dTtilde;
    int N_mub, N_muq, N_mus, N_T;

    // max values for the tilde variables
    float T_tilde_max;
    float mub_tilde_max;
    float muq_tilde_max;
    float mus_tilde_max;
    // double Ttilde, mubtilde, muqtilde, mustilde;

    // useful constants
    const double alphaNf = (8 / 45.0 + 7 / 60.0 * 3.0) * M_PI * M_PI;
    const double OneoveralphaNf = 1. / alphaNf;

    bool EoS_file_in_binary_;
    bool file_for_cs_;

    //! Ttilde_{iT+1}^5 - Ttilde_{iT}^5 of each Ttilde cell, for the slope
    //! in FourDLInterp (set in get_eos_max_values)
    std::vector<double> T5diff_;

    // 1D tables.
    std::vector<float> pressure_vec;
    std::vector<float> temp_vec;
    std::vector<float> mub_vec;
    std::vector<float> muq_vec;
    std::vector<float> mus_vec;

    // method to read/mainupalate header info and data
    void read_eos(
        std::string filepath, bool is_cs, std::vector<float> &out,
        int header_size = 2);
    void read_eos_binary(
        std::string filepath, bool is_cs, std::vector<float> &out,
        int header_size = 12);

    void get_eos_max_values();

    void read_header(std::string filepath);
    void read_header_binary(std::string filepath, int header_size = 12);

    // Shift in the index corresponds to the header size.
    int index(int i_T, int i_mub, int i_muq, int i_mus) const;

    void FourDLInterp(
        const std::vector<float> &data, const std::array<float, 4> &TildeVar,
        std::array<float, 5> &ResArr, bool compute_derivatives = false) const;

    //! the table cell and the 16 interpolation weights of a point
    struct FourDCell {
        int iT = 0;
        int base = 0;  // index() of the lower corner
        float w0000, w1111, w1000, w0100, w0010, w0001, w1001, w0101, w0011,
            w1100, w1010, w0110, w0111, w1011, w1101, w1110;
    };
    //! FourDLInterp = FourDSetup + FourDEval: the setup depends only on
    //! TildeVar and is shared by the tables
    void FourDSetup(const std::array<float, 4> &TildeVar, FourDCell &c) const;
    void FourDEval(
        const std::vector<float> &data, const std::array<float, 4> &TildeVar,
        const FourDCell &c, std::array<float, 5> &ResArr,
        bool compute_derivatives = false) const;
    //! get_pressure_with_gradients from the scaled TildeVar and its cell
    void pressureWithGradients(
        const double rhob, const double rhoq, const double rhos,
        const std::array<float, 4> &TildeVar, const float scale,
        const FourDCell &c, double &p, double &dpde, double &dpdrhob,
        double &dpdrhoq, double &dpdrhos) const;
    void get_tilde_variables(
        double e, double nb, double nq, double ns,
        std::array<float, 4> &TildeVar) const;

    //! Below Ttilde_dilute_ the table is not usable: P is interpolated
    //! linearly in Ttilde from P = 0 at Ttilde = 0 (P/e ~ e^{-3/4} near
    //! vacuum), the cs table is zero, and most nodes have the chemical
    //! potentials at the bounds of the table generator. There the EoS is
    //! read at Ttilde_dilute_ with the same mu-tilde and scaled: P and its
    //! e-derivative as e, T and mu as Ttilde (P = e g(mu-tilde)).
    //! Returns Ttilde / Ttilde_dilute_ (1 outside the dilute region) and
    //! moves TildeVar[0] up to Ttilde_dilute_.
    float diluteScale(std::array<float, 4> &TildeVar) const;
    static constexpr double Ttilde_dilute_GeV_ = 0.040;

  public:
    EOS_4D();
    ~EOS_4D();

    void initialize_eos();

    void set_eos_file_in_binary(bool is_binary) {
        EoS_file_in_binary_ = is_binary;
    }
    bool get_eos_file_in_binary() { return EoS_file_in_binary_; }

    double get_temperature(
        double e, double rhob, double rhoq = 0.0, double rhos = 0.0) const;
    double get_muB(
        double e, double rhob, double rhoq = 0.0, double rhos = 0.0) const;
    double get_muS(
        double e, double rhob, double rhoq = 0.0, double rhos = 0.0) const;
    double get_muQ(
        double e, double rhob, double rhoq = 0.0, double rhos = 0.0) const;
    double get_pressure(
        double e, double rhob, double rhoq = 0.0, double rhos = 0.0) const;
    double get_cs2(
        double e, double rhob, double rhoq = 0.0, double rhos = 0.0) const;

    void get_pressure_with_gradients(
        double epsilon, double rhob, double rhoq, double rhos, double &p,
        double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos) const;

    void get_pressure_with_gradients_and_cs2(
        double epsilon, double rhob, double rhoq, double rhos, double &p,
        double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos,
        double &cs2) const;

    void getThermalVariables(
        double epsilon, double rhob, double rhoq, double rhos,
        ThermalVec &thermalVec) const;

    bool has_QS_densities() const { return (true); }

    void check_eos() const { check_4D_eos(); }
};

#endif  // SRC_EOS_4D_H_
