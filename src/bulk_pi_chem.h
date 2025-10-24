// src/bulk_pi_chem.h
#pragma once
#include <algorithm>
#include <cmath>
#include "eos.h"
#include "util.h"

using Util::hbarc;

struct ChemBulkConfig {
    int    chem_bulk_on = 0;     // 0: off (use baseline); 1: on (runtime τ_Π, ζ/s)
    double C_tauPi      = 5.0;   // τ_Π = 1 / [ C_tauPi * T * (1+√Yq)^2 ]
    double tauPi_min    = 1e-3;  // [fm/c] safety floor
    double Yq_eps       = 1e-16; // numerical epsilon inside sqrt/clamp
};

class BulkPiChem {
public:
    explicit BulkPiChem(const ChemBulkConfig& cfg, const EOS* eos_ptr)
        : cfg_(cfg), eos_(eos_ptr) {}

    inline bool enabled() const { return cfg_.chem_bulk_on != 0; }

    // ---- EOS hooks --------------------------------------------------------
    inline double I_QCD(double e, double rhob) const {
    #ifdef HAVE_EOS_IQCD
        return eos_->get_IQCD(e);
    #else
        return e - 3.0*eos_->get_pressure(e, rhob);
    #endif
    }
    inline double I0(double /*e*/, double /*rhob*/) const {
    #ifdef HAVE_EOS_I0
        return eos_->get_I0(e);
    #else
        return 0.0;
    #endif
    }
    inline double T(double e, double rhob) const { return eos_->get_temperature(e, rhob); }

    // convenience
    inline double dI(double e, double rhob) const { return I_QCD(e, rhob) - I0(e, rhob); }

    // ---- Utilities --------------------------------------------------------
    static inline double clamp01(double x) {
        return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
    }

    // Π → Y_q :  Y_q = clamp01( [ 1 - 3 Π / (I_QCD - I0) ]^2 )
    inline double Yq_from_Pi(double e, double rhob, double Pi_total) const {
        if (!enabled()) return 1.0;
        const double dI_val = dI(e, rhob);
        if (std::abs(dI_val) < 1e-30) return 1.0;
        const double y = 1.0 - 3.0 * Pi_total / dI_val;
        return clamp01(y * y);
    }

    // Y_q → Π :  Π = (1/3) (I_QCD - I0) (1 - √Y_q)
    inline double Pi_from_Yq(double e, double rhob, double Yq) const {
        const double y = std::sqrt(clamp01(Yq));
        return (1.0 / 3.0) * dI(e, rhob) * (1.0 - y);
    }

    // ---- Transport coefficients ------------------------------------------
    // τ_Π(e,Π) = 1 / [ C * T(e) * (1 + √Y_q)^2 ]  (GeV⁻¹ → fm/c via ħc)
    inline double tauPi_runtime(double e, double rhob, double Pi_total) const {
        if (!cfg_.chem_bulk_on) return cfg_.tauPi_min;
        const double yq   = std::max(cfg_.Yq_eps, Yq_from_Pi(e, rhob, Pi_total));
        const double temp = std::max(1e-16, T(e, rhob));
        const double denom = cfg_.C_tauPi * temp * std::pow(1.0 + std::sqrt(yq), 2.0);
        const double tau   = hbarc / denom;  // convert GeV⁻¹ → fm/c
        return std::max(cfg_.tauPi_min, tau);
    }

    // (ζ/s)_runtime(e,Π): (ζ/s)_eq(e) × (1 - √Y_q)^2  (→0 as Y_q→1)
    inline double zeta_over_s_runtime(double e, double rhob, double Pi_total,
                                      double zetas_eq) const {
        if (!cfg_.chem_bulk_on) return std::max(0.0, zetas_eq);
        const double yq    = std::max(cfg_.Yq_eps, Yq_from_Pi(e, rhob, Pi_total));
        const double shape = std::pow(1.0 - std::sqrt(yq), 2.0);
        return std::max(0.0, zetas_eq * shape);
    }

    // Optional: placeholder if you later split Π = Π_kin + Π_chem
    inline double PiChem(double /*e*/, double /*rhob*/, double /*Pi_kin*/) const {
        if (!cfg_.chem_bulk_on) return 0.0;
        return 0.0;
    }

private:
    ChemBulkConfig cfg_;
    const EOS* eos_;
};
