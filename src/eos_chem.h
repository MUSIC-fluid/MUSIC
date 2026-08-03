#ifndef SRC_EOS_CHEM_H_
#define SRC_EOS_CHEM_H_

#include "eos_base.h"

// Equation of state for the chemically non-equilibrated QGP.
// P, s, and T are tabulated on a two-dimensional (e, Y_q) grid, where
// Y_q in [0,1] is the chemical-composition variable of the chemical
// equilibration paper: Y_q = 1 is full chemical equilibrium (matches the
// lattice QCD EOS), Y_q = 0 is the pure-glue conformal limit (P = e/3).
// The table is built from
//     I(e,Y_q) = (1 - sqrt(Y_q)) * I_0(e) + sqrt(Y_q) * I_QCD(e)
//     P(e,Y_q) = (e - I(e,Y_q)) / 3
// with s(e,Y_q) obtained from integrating the entropy relation at fixed
// Y_q starting from a common reference entropy s(e0,Y_q) = s_QCD(e0) for
// every Y_q, and T(e,Y_q) = (e + P(e,Y_q)) / s(e,Y_q).
//
// Y_q enters through the third argument of the extended EOS_base
// interface (see eos_base.h), whose default value of 1 corresponds to
// full chemical equilibrium. The rhob argument is accepted for interface
// compatibility and ignored: this table has no baryon density
// dependence. Internally the generic 2D table machinery inherited from
// EOS_base stores Y_q on its "rhob" axis; that is a private table-layout
// detail and no longer visible in any function signature.
//
// This EOS is still not selectable through EOS_to_use (see the NOTE in
// eos.cpp): the generic call sites in the evolution code omit Y_q and
// would therefore always evaluate the Y_q = 1 slice, which is identical
// to the lattice EOS they already use. The physically relevant Y_q of a
// fluid cell is derived from the evolved pi_b_chem field (see the Y_q
// calculation in grid_info.cpp) and must be passed explicitly by code
// that has computed it.
//
// Note: get_entropy(double, double, double) below hides (does not
// override, since EOS_base::get_entropy is not virtual) the base class
// version. Calling get_entropy on an EOS_chem object or pointer uses the
// tabulated s(e,Y_q) directly. Calling it through a generic EOS_base
// pointer (for example via the EOS wrapper class) falls back to the base
// class's (e + P - ...) / T formula built from the (correctly virtual)
// get_pressure/get_temperature overrides below; that fallback is
// numerically consistent with the table to high accuracy since the table
// already satisfies T*s = e + P by construction, but it is not a direct
// table lookup.
class EOS_chem : public EOS_base {
  private:
    const int eos_id;
    double ***entropy_tb;

  public:
    EOS_chem(const int eos_id_in);
    ~EOS_chem();

    void initialize_eos();
    double p_e_func(double e, double rhob, double Y_q) const;
    double get_temperature(double e, double rhob, double Y_q) const;
    double get_pressure(double e, double rhob, double Y_q) const;
    double get_entropy(double e, double rhob, double Y_q) const;
    double get_s2e(double s, double rhob, double Y_q) const;
    double get_T2e(double T_in_GeV, double rhob, double Y_q) const;
    void get_pressure_with_gradients(
        double e, double rhob, double &p, double &dpde, double &dpdYq,
        double &cs2, double Y_q) const;

    void check_eos() const { check_eos_no_muB(); }
};

#endif  // SRC_EOS_CHEM_H_
