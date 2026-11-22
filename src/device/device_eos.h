// Copyright 2026 Chun Shen
//
// The equation of state on the device: a POD descriptor of the tables and
// the pressure functions of the EoSs the device path supports, written
// operation by operation as their host versions (EOS_base::interpolate1D,
// interpolate2D and their _with_gradients variants, EOS_idealgas,
// EOS_hotQCD, EOS_neos), so that in double the device reproduces the host.
//   type 0: conformal ideal gas (EoS 0), P = e/3
//   type 1: one (e) table, hotQCD (EoS 9, 91)
//   type 2: (e, n_B) tables, NEOS (EoS 10 - 14)

#ifndef SRC_DEVICE_DEVICE_EOS_H_
#define SRC_DEVICE_DEVICE_EOS_H_

#ifndef MUSIC_DEVICE_BUNDLED
#include "device_defs.h"
#endif

namespace music_device {

enum { kMaxEosTables = 16 };

template <typename Real>
struct EosDesc {
    int type;
    int nTables;
    int eLength[kMaxEosTables];
    int nbLength[kMaxEosTables];
    int offset[kMaxEosTables];  // of each table in the pressure buffer
    Real eBounds[kMaxEosTables];
    Real eSpacing[kMaxEosTables];
    Real nbBounds[kMaxEosTables];
    Real nbSpacing[kMaxEosTables];
    // the thermal variables of the viscous update (eosThermal)
    int tOffset;          // type 1: T^5 on the e grid of table 0
    int xOffset, xLength;  // type 0: mu_B / T on a grid in rho
    Real coeffA, coeffB;  // type 0: P = T^4 D(x), D = a + b (...)
    Real rCap, xMax;      // type 0: n_B / e^{3/4} at the cap, the cap
};

template <typename Real>
struct EosGrad {
    Real p, dpde, dpdrhob, cs2;
};

//! table[t][inb][ie]
template <typename Real>
MD_FN Real eosTab(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, int t,
    int inb, int ie) {
    return tab[d.offset[t] + inb * d.eLength[t] + ie];
}

template <typename Real>
MD_FN int eosTableIdx(MD_CONST const EosDesc<Real> &d, Real e) {
    for (int itable = 1; itable < d.nTables; itable++) {
        if (e < d.eBounds[itable]) return itable - 1;
    }
    return dmax(0, d.nTables - 1);
}

template <typename Real>
MD_FN Real eosInterp1D(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e) {
    const Real e0 = d.eBounds[0];
    const Real de = d.eSpacing[0];
    const int Ne = d.eLength[0];
    int idx = int((e - e0) / de);
    idx = dmax(0, idx);
    idx = dmin(Ne - 2, idx);
    if (e < e0) return eosTab(d, tab, 0, 0, 0) * e / e0;
    const Real frac = (e - (Real(idx) * de + e0)) / de;
    const Real t1 = eosTab(d, tab, 0, 0, idx);
    const Real t2 = eosTab(d, tab, 0, 0, idx + 1);
    return t1 * (Real(1) - frac) + t2 * frac;
}

template <typename Real>
MD_FN void eosInterp1DGrad(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    MD_THREAD Real &p, MD_THREAD Real &dpde) {
    const Real e0 = d.eBounds[0];
    const Real de = d.eSpacing[0];
    const int Ne = d.eLength[0];
    int idx = int((e - e0) / de);
    idx = dmin(Ne - 2, dmax(0, idx));
    if (e < e0) {
        dpde = eosTab(d, tab, 0, 0, 0) / e0;
        p = dpde * e;
    } else {
        const Real t1 = eosTab(d, tab, 0, 0, idx);
        const Real t2 = eosTab(d, tab, 0, 0, idx + 1);
        dpde = (t2 - t1) / de;
        const Real frac = e - (Real(idx) * de + e0);
        p = t1 + dpde * frac;
    }
}

// the next table's first e column at n_B (the seam of interpolate2D)
template <typename Real>
MD_FN Real eosNextTableAtNb(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, int next,
    Real nb, MD_THREAD Real &fRaw) {
    const Real nb0 = d.nbBounds[next];
    const Real dnb = d.nbSpacing[next];
    const int Nnb = d.nbLength[next];
    int i = int((nb - nb0) / dnb);
    i = dmax(0, dmin(Nnb - 2, i));
    fRaw = (nb - (Real(i) * dnb + nb0)) / dnb;
    const Real f = dmax(Real(0), dmin(Real(1), fRaw));
    return eosTab(d, tab, next, i, 0) * (Real(1) - f)
           + eosTab(d, tab, next, i + 1, 0) * f;
}

template <typename Real>
MD_FN Real eosInterp2D(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real nb, int t) {
    const Real e0 = d.eBounds[t];
    const Real nb0 = d.nbBounds[t];
    const Real de = d.eSpacing[t];
    const Real dnb = d.nbSpacing[t];
    const int Ne = d.eLength[t];
    const int Nnb = d.nbLength[t];
    int idx_e = int((e - e0) / de);
    int idx_nb = int((nb - nb0) / dnb);
    idx_e = dmin(Ne - 1, idx_e);
    if (t == d.nTables - 1) idx_e = dmin(Ne - 2, idx_e);
    idx_nb = dmin(Nnb - 2, idx_nb);
    idx_e = dmax(0, idx_e);
    idx_nb = dmax(0, idx_nb);
    const Real frac_e = (e - (Real(idx_e) * de + e0)) / de;
    Real frac_rhob = (nb - (Real(idx_nb) * dnb + nb0)) / dnb;
    frac_rhob = dmin(Real(1), frac_rhob);
    const Real temp1 = eosTab(d, tab, t, idx_nb, idx_e);
    const Real temp4 = eosTab(d, tab, t, idx_nb + 1, idx_e);
    if (idx_e == Ne - 1) {
        Real fRaw;
        const Real curr = temp1 * (Real(1) - frac_rhob) + temp4 * frac_rhob;
        const Real nextVal = eosNextTableAtNb(d, tab, t + 1, nb, fRaw);
        return curr * (Real(1) - frac_e) + nextVal * frac_e;
    }
    const Real temp2 = eosTab(d, tab, t, idx_nb, idx_e + 1);
    const Real temp3 = eosTab(d, tab, t, idx_nb + 1, idx_e + 1);
    return ((temp1 * (Real(1) - frac_e) + temp2 * frac_e) * (Real(1) - frac_rhob)
            + (temp3 * frac_e + temp4 * (Real(1) - frac_e)) * frac_rhob);
}

template <typename Real>
MD_FN void eosInterp2DGrad(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real nb, int t, MD_THREAD Real &p, MD_THREAD Real &dpde,
    MD_THREAD Real &dpdrhob) {
    const Real e0 = d.eBounds[t];
    const Real nb0 = d.nbBounds[t];
    const Real de = d.eSpacing[t];
    const Real dnb = d.nbSpacing[t];
    const int Ne = d.eLength[t];
    const int Nnb = d.nbLength[t];
    int idx_e = int((e - e0) / de);
    int idx_nb = int((nb - nb0) / dnb);
    idx_e = dmin(Ne - 1, idx_e);
    if (t == d.nTables - 1) idx_e = dmin(Ne - 2, idx_e);
    idx_nb = dmin(Nnb - 2, idx_nb);
    idx_e = dmax(0, idx_e);
    idx_nb = dmax(0, idx_nb);
    const Real frac_e = (e - (Real(idx_e) * de + e0)) / de;
    Real frac_rhob = (nb - (Real(idx_nb) * dnb + nb0)) / dnb;
    const bool nbClamped = (frac_rhob > Real(1));
    frac_rhob = dmin(Real(1), frac_rhob);
    const Real temp1 = eosTab(d, tab, t, idx_nb, idx_e);
    const Real temp4 = eosTab(d, tab, t, idx_nb + 1, idx_e);
    if (idx_e == Ne - 1) {
        const int next = t + 1;
        Real fRaw;
        const Real nextVal = eosNextTableAtNb(d, tab, next, nb, fRaw);
        const Real dnb_next = d.nbSpacing[next];
        int i = int((nb - d.nbBounds[next]) / dnb_next);
        i = dmax(0, dmin(d.nbLength[next] - 2, i));
        const Real dnextdnb = (fRaw < Real(0) || fRaw > Real(1))
                                  ? Real(0)
                                  : (eosTab(d, tab, next, i + 1, 0)
                                     - eosTab(d, tab, next, i, 0))
                                        / dnb_next;
        const Real curr = temp1 * (Real(1) - frac_rhob) + temp4 * frac_rhob;
        const Real dcurrdnb = nbClamped ? Real(0) : (temp4 - temp1) / dnb;
        p = curr * (Real(1) - frac_e) + nextVal * frac_e;
        dpde = (nextVal - curr) / de;
        dpdrhob = dcurrdnb * (Real(1) - frac_e) + dnextdnb * frac_e;
        return;
    }
    const Real temp2 = eosTab(d, tab, t, idx_nb, idx_e + 1);
    const Real temp3 = eosTab(d, tab, t, idx_nb + 1, idx_e + 1);
    Real p1 = temp1 * (Real(1) - frac_rhob) + temp4 * frac_rhob;
    Real p2 = temp2 * (Real(1) - frac_rhob) + temp3 * frac_rhob;
    dpde = (p2 - p1) / de;
    p1 = temp1 * (Real(1) - frac_e) + temp2 * frac_e;
    p2 = temp4 * (Real(1) - frac_e) + temp3 * frac_e;
    dpdrhob = nbClamped ? Real(0) : (p2 - p1) / dnb;
    p = p1 * (Real(1) - frac_rhob) + p2 * frac_rhob;
}

//! EOS_base::interpolate1D of the table at offset off on the e grid of
//! table 0
template <typename Real>
MD_FN Real eosInterp1DAt(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, int off,
    Real e) {
    const Real e0 = d.eBounds[0];
    const Real de = d.eSpacing[0];
    const int Ne = d.eLength[0];
    int idx = int((e - e0) / de);
    idx = dmax(0, idx);
    idx = dmin(Ne - 2, idx);
    if (e < e0) return tab[off] * e / e0;
    const Real frac = (e - (Real(idx) * de + e0)) / de;
    const Real t1 = tab[off + idx];
    const Real t2 = tab[off + idx + 1];
    return t1 * (Real(1) - frac) + t2 * frac;
}

// ----- EOS_idealgas: x = mu_B / T at (e, n_B) -------------------------------

template <typename Real>
MD_FN Real idealD(MD_CONST const EosDesc<Real> &d, Real x) {
    const Real x2 = x * x;
    return (d.coeffA
            + d.coeffB
                  * (x2 / Real(54.)
                     + x2 * x2
                           / (Real(972.) * Real(3.14159265358979323846)
                              * Real(3.14159265358979323846))));
}

template <typename Real>
MD_FN Real idealN(MD_CONST const EosDesc<Real> &d, Real x) {
    return (d.coeffB
            * (x / Real(27.)
               + x * x * x
                     / (Real(243.) * Real(3.14159265358979323846)
                        * Real(3.14159265358979323846))));
}

template <typename Real>
MD_FN Real idealS(MD_CONST const EosDesc<Real> &d, Real x) {
    return (Real(4.) * d.coeffA + d.coeffB * x * x / Real(27.));
}

template <typename Real>
MD_FN Real idealTofEX(MD_CONST const EosDesc<Real> &d, Real e, Real x) {
    return (sqrt(sqrt(dmax(Real(0.), e) / (Real(3.) * idealD(d, x)))));
}

template <typename Real>
MD_FN Real idealDN(MD_CONST const EosDesc<Real> &d, Real x) {
    return (d.coeffB
            * (Real(1.) / Real(27.)
               + x * x
                     / (Real(81.) * Real(3.14159265358979323846)
                        * Real(3.14159265358979323846))));
}

//! EOS_idealgas::x_of_rho
template <typename Real>
MD_FN Real idealXofRho(MD_CONST const EosDesc<Real> &d, Real rho, Real xGuess) {
    const Real target = log(rho * d.rCap) + Real(0.75) * log(Real(3.));
    const Real xMax = d.xMax;
    // 1e-300 as on the host; the smallest normal float in float
#if defined(__METAL_VERSION__)
    const Real xMin = Real(1.2e-38);
#else
    const Real xMin = (sizeof(Real) == 8 ? Real(1e-300) : Real(1.2e-38));
#endif
    Real x = dmax(xMin, dmin(xMax, xGuess));
    for (int iter = 0; iter < 60; iter++) {
        const Real N = idealN(d, x);
        const Real dN = idealDN(d, x);
        const Real h = log(N) - Real(0.75) * log(idealD(d, x)) - target;
        const Real dhdu = x * (dN / N - Real(0.75) * N / idealD(d, x));
        const Real du = -h / dhdu;
        x = dmin(xMax, x * exp(dmax(Real(-2.), dmin(Real(2.), du))));
        if (fabs(du) < Real(1e-13)) break;
    }
    return (x);
}

//! EOS_idealgas::get_muB_over_T
template <typename Real>
MD_FN Real idealMuBOverT(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob) {
    const Real n = fabs(rhob);
    if (n <= Real(0.) || e <= Real(0.)) return (Real(0.));
    const Real n2 = n * n;
    const Real C = Real(27.) * n2 * n2 / (e * e * e);
    const Real rc2 = d.rCap * d.rCap;
    const Real rho = sqrt(sqrt(C / (Real(27.) * rc2 * rc2)));
    if (rho >= Real(1.)) return (d.xMax);
    const Real pos = rho * Real(d.xLength - 1);
    const int i = dmin(int(pos), d.xLength - 2);
    const Real f = pos - Real(i);
    MD_GLOBAL const Real *xt = tab + d.xOffset;
    Real x = xt[i] * (Real(1.) - f) + xt[i + 1] * f;
    if (x <= Real(0.)) return (idealXofRho(d, rho, x));
    const Real b = d.coeffB;
    for (int iter = 0; iter < 8; iter++) {
        const Real N = idealN(d, x), D = idealD(d, x);
        const Real dN =
            b
            * (Real(1.) / Real(27.)
               + x * x
                     / (Real(81.) * Real(3.14159265358979323846)
                        * Real(3.14159265358979323846)));
        const Real N3 = N * N * N, D2 = D * D;
        const Real dx =
            -(N3 * N - C * D2 * D) / (Real(4.) * N3 * dN - Real(3.) * C * D2 * N);
        x += dx;
        if (fabs(dx) < Real(1e-13) * x) {
            return (dmin(x, d.xMax));
        }
    }
    return (idealXofRho(d, rho, xt[i] * (Real(1.) - f) + xt[i + 1] * f));
}

//! EOS::getThermalVariables, the entries the viscous update reads
template <typename Real>
struct Thermal {
    Real e, rhob, p, dpde, dpdrhob, cs2, T, muB, s;
};

template <typename Real>
MD_FN Thermal<Real> eosThermal(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob, Real rhoq, Real rhos) {
    const Real smallEps = Real(1e-16);
    Thermal<Real> t;
    t.e = e;
    t.rhob = rhob;
    if (d.type == 0) {
        // EOS_idealgas::getThermalVariables
        const Real x = idealMuBOverT(d, tab, e, rhob);
        const Real T = idealTofEX(d, e, x);
        t.p = e / Real(3.);
        t.dpde = Real(1.) / Real(3.);
        t.dpdrhob = Real(0.);
        t.cs2 = Real(1.) / Real(3.);
        t.T = T;
        t.muB = copysign(x * T, rhob);
        t.s = dmax(smallEps, T * T * T * idealS(d, x));
        return (t);
    }
    // EOS_base::getThermalVariables with EOS_hotQCD: mu = 0
    const EosGrad<Real> g = eosPressureGrad(d, tab, e, rhob, rhoq, rhos);
    t.p = g.p;
    t.dpde = g.dpde;
    t.dpdrhob = g.dpdrhob;
    t.cs2 = g.cs2;
    const Real T5 = eosInterp1DAt(d, tab, d.tOffset, e);
    t.T = dmax(smallEps, Real(pow(T5, Real(0.2))));
    t.muB = Real(0.);
    const Real muS = Real(0.), muQ = Real(0.);
    const Real rhobT = rhob, rhoS = Real(0.), rhoQ = Real(0.4) * rhob;
    Real entropy = ((e + t.p - t.muB * rhobT - muS * rhoS - muQ * rhoQ)
                    / (t.T + smallEps));
    t.s = dmax(smallEps, entropy);
    return (t);
}

//! EOS::get_pressure
template <typename Real>
MD_FN Real eosPressure(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob) {
    const Real smallEps = Real(1e-16);
    if (d.type == 0) return Real(1) / Real(3) * e;
    if (d.type == 1) return dmax(smallEps, eosInterp1D(d, tab, e));
    const int t = eosTableIdx(d, e);
    return dmax(smallEps, eosInterp2D(d, tab, e, fabs(rhob), t));
}

//! EOS::get_pressure_with_gradients_and_cs2 (n_Q, n_S do not enter these
//! EoSs: dP/dn_Q = dP/dn_S = 0)
template <typename Real>
MD_FN EosGrad<Real> eosPressureGrad(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob, Real rhoq, Real rhos) {
    const Real smallEps = Real(1e-16);
    EosGrad<Real> g;
    if (d.type == 0) {
        g.dpde = Real(1) / Real(3);
        g.p = g.dpde * e;
        g.dpdrhob = Real(0);
        g.cs2 = g.dpde;
        return g;
    }
    if (d.type == 1) {
        eosInterp1DGrad(d, tab, e, g.p, g.dpde);
        g.p = dmax(smallEps, g.p);
        g.dpdrhob = Real(0);
    } else {
        const int t = eosTableIdx(d, e);
        eosInterp2DGrad(d, tab, e, fabs(rhob), t, g.p, g.dpde, g.dpdrhob);
        g.p = dmax(smallEps, g.p);
    }
    Real cs2 = g.dpde
               + (g.dpdrhob * rhob + Real(0) * rhos + Real(0) * rhoq)
                     / (e + g.p);
    g.cs2 = dmax(Real(0.01), dmin(Real(1) / Real(3), cs2));
    return g;
}

}  // namespace music_device

#endif  // SRC_DEVICE_DEVICE_EOS_H_
