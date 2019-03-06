// Copyright @ Chun Shen
#ifndef SRC_CRITICAL_MODES_H_
#define SRC_CRITICAL_MODES_H_

#include <vector>
#include "data.h"
#include "grid.h"
#include "eos.h"
#include "pretty_ostream.h"

class CriticalSlowModes {
 private:
    InitData &DATA;
    const EOS &eos;
    std::vector<double> Qvec;
    std::vector<GridphiQ> phiQfields;
    pretty_ostream music_message;

 public:
    CriticalSlowModes(const EOS &eos_in, InitData &DATA_in);
    ~CriticalSlowModes() {}

    void InitializeFields(const int nQ, const int nX,
                          const int nY, const int nEta);
    int get_Qvec_size() const {return(Qvec.size());}
    int get_nphiQfields() const {return(phiQfields.size());}

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
};

#endif  // SRC_CRITICAL_MODES_H_
