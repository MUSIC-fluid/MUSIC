// Copyright @ Chun Shen
#ifndef SRC_CRITICAL_MODES_H_
#define SRC_CRITICAL_MODES_H_

#include <array>
#include <vector>
#include "data.h"
#include "grid.h"
#include "eos.h"
#include "pretty_ostream.h"
#include "data_struct.h"

class CriticalSlowModes {
 private:
    InitData &DATA;
    const EOS &eos;
    std::vector<double> Qvec;
    pretty_ostream music_message;

 public:
    CriticalSlowModes(const EOS &eos_in, InitData &DATA_in);
    ~CriticalSlowModes();

    void InitializeFields(const int nQ, SCGrid &arena_current);

    int get_Qvec_size() const {return(Qvec.size());}

    double phiQbar_f2(const double x) const;
    double phiQbar_0(const double e, const double rho_b) const;
    //! This function gets the local correlation length
    double get_xi(const double e, const double rho_b) const;

    //! This function computes the equilibrium value of phi_Q
    double compute_phiQ_equilibrium(const double Qxi,
                                    const double e, const double rho_b) const;

    //! This function computes the relaxation rate for the phiQ fields
    double get_GammaQ(const double Q, const double xi,
                      const double T, const double eta) const;

    //! This function evolves the phi_Q field at ix, iy, ieta
    //! by one RK step in time
    void evolve_phiQfields(const double tau, const int ix, const int iy,
                           const int ieta, const int rk_flag,
                           const double T_local, FlowVec umu);
};

#endif  // SRC_CRITICAL_MODES_H_
