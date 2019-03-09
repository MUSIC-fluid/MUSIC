// Copyright @ Chun Shen
#ifndef SRC_CRITICAL_MODES_H_
#define SRC_CRITICAL_MODES_H_

#include <array>
#include <vector>
#include "data.h"
#include "grid.h"
#include "eos.h"
#include "minmod.h"
#include "pretty_ostream.h"
#include "data_struct.h"

class CriticalSlowModes {
 private:
    const InitData &DATA;
    const EOS &eos;
    std::vector<double> Qvec;
    Minmod minmod;
    pretty_ostream music_message;

 public:
    CriticalSlowModes(const EOS &eos_in, const InitData &DATA_in);
    ~CriticalSlowModes();

    void InitializeFields(const int nQ, SCGrid &arena_current);

    int get_Qvec_size() const {return(Qvec.size());}
    double get_Qi(int iQ) const {return(Qvec[iQ]);}
    std::vector<double> get_Qvec() const {return(Qvec);}

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
    void evolve_phiQfields(const double tau, SCGrid &arena_prev,
                           SCGrid &arena_current, SCGrid &arena_future,
                           const double theta_local,
                           const int ix, const int iy, const int ieta,
                           const int rk_flag);

    void compute_KTflux(
        const double tau, SCGrid &arena, const double theta_local,
        const int ix, const int iy, const int ieta, const int iQ,
        const DeltaXVec delta, double &flux_term) const;

    double compute_relaxation_source_term(
        const double tau, Cell_small *grid_pt, Cell_small *grid_pt_prev,
        const int iQ, const int rk_flag) const;
};

#endif  // SRC_CRITICAL_MODES_H_
