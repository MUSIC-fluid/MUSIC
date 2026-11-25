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
    // type 2 (NEOS), thermal variables: the T^5, mu_B, mu_S, mu_Q tables
    // (laid out as the pressure table, -1 if absent) and the n_B edge of
    // each e column (EOS_base::nb_edge, the tables one after the other)
    int t2Offset, muBOffset, muSOffset, muQOffset, edgeOffset;
    // type 3: EOS_4D. The tables are float (as on the host) and start at
    // index f4Offset of the Real table buffer, f4Size floats each: P, then
    // (with the thermal variables) T, mu_B, mu_Q, mu_S. T5diff_ (Real) is
    // at index t5Offset.
    int f4Offset, f4Size, f4NTables, t5Offset;
    int nT4, nMub4, nMuq4, nMus4;
    float T04, mub04, muq04, mus04, dT4, dmub4, dmuq4, dmus4;
    float Tmax4, mubMax4, muqMax4, musMax4, TrefDilute4;
    Real alphaNf4, oneOverAlphaNf4;
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
        // EOS_idealgas::get_muB_over_T: a slope <= 0, the robust solve
        const Real slope = Real(4.) * N3 * dN - Real(3.) * C * D2 * N;
        if (!(slope > Real(0.))) break;
        const Real dx = -(N3 * N - C * D2 * D) / slope;
        x += dx;
        if (fabs(dx) < Real(1e-13) * x) {
            return (dmin(x, d.xMax));
        }
    }
    return (idealXofRho(d, rho, xt[i] * (Real(1.) - f) + xt[i + 1] * f));
}


// ----- EOS_4D: quadrilinear tables in (Ttilde, mu-tilde_{B,Q,S}) ------------
// Every operation as in eos_4D.cpp: float where the host computes in float,
// Real where it computes in double (float on Metal).

enum { kEos4DP = 0, kEos4DT = 1, kEos4DMuB = 2, kEos4DMuQ = 3, kEos4DMuS = 4 };

template <typename Real>
MD_FN MD_GLOBAL const float *eos4DTable(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, int which) {
    return ((MD_GLOBAL const float *)(tab + d.f4Offset)) + which * d.f4Size;
}

//! EOS_4D::get_tilde_variables
template <typename Real>
MD_FN void eos4DTilde(
    MD_CONST const EosDesc<Real> &d, Real e, Real rhob, Real rhoq, Real rhos,
    MD_THREAD float *TildeVar) {
    Real Ttilde_sq = sqrt(mdDiv(e, Real(3.0)) * d.oneOverAlphaNf4);
    Real Ttilde = sqrt(Ttilde_sq);
    // EOS_4D::get_tilde_variables: Ttilde^2 floored at the smallest normal
    const Real Tsq_div = dmax(Ttilde_sq, tinyNormal<Real>());
    Real mubtilde, muqtilde, mustilde;
    // the divisions (device_defs.h): mdDivF for finite densities over a
    // normal Ttilde^2, mdDiv otherwise
    if (Tsq_div >= Real(1.1754943508222875e-38)
        && Tsq_div < Real(8.507059173023462e37)
        && isfinite(rhob) && isfinite(rhoq) && isfinite(rhos)) {
        mubtilde = mdDivF(Real(5) * rhob - rhoq + Real(2) * rhos, Tsq_div);
        muqtilde = mdDivF(-rhob + Real(2) * rhoq - rhos, Tsq_div);
        mustilde = mdDivF(Real(2) * rhob - rhoq + Real(2) * rhos, Tsq_div);
    } else {
        mubtilde = mdDiv(Real(5) * rhob - rhoq + Real(2) * rhos, Tsq_div);
        muqtilde = mdDiv(-rhob + Real(2) * rhoq - rhos, Tsq_div);
        mustilde = mdDiv(Real(2) * rhob - rhoq + Real(2) * rhos, Tsq_div);
    }
    TildeVar[0] = float(Ttilde);
    TildeVar[1] = float(mubtilde);
    TildeVar[2] = float(muqtilde);
    TildeVar[3] = float(mustilde);
}

//! EOS_4D::diluteScale
template <typename Real>
MD_FN float eos4DDilute(
    MD_CONST const EosDesc<Real> &d, MD_THREAD float *TildeVar) {
    const float Tref = d.TrefDilute4;
    if (TildeVar[0] >= Tref) return 1.f;
    // TildeVar[0] < Tref (or NaN); Tref checked at the upload
    const float scale = mdDivF(TildeVar[0], Tref);
    TildeVar[0] = Tref;
    return scale;
}

//! EOS_4D::FourDCell
struct Eos4DCell {
    int iT;
    int base;
    float w0000, w1111, w1000, w0100, w0010, w0001, w1001, w0101, w0011,
        w1100, w1010, w0110, w0111, w1011, w1101, w1110;
};

//! EOS_4D::FourDSetup
template <typename Real>
MD_FN Eos4DCell eos4DSetup(
    MD_CONST const EosDesc<Real> &d, MD_THREAD const float *TildeVar) {
    Eos4DCell c;
    // the grid spacings are checked at the upload; mdDivF needs finite
    // tilde variables
    Real Tfrac, mubfrac, muqfrac, musfrac;
    if (isfinite(TildeVar[0]) && isfinite(TildeVar[1])
        && isfinite(TildeVar[2]) && isfinite(TildeVar[3])) {
        Tfrac = Real(mdDivF(TildeVar[0] - d.T04, d.dT4));
        mubfrac = Real(mdDivF(TildeVar[1] - d.mub04, d.dmub4));
        muqfrac = Real(mdDivF(TildeVar[2] - d.muq04, d.dmuq4));
        musfrac = Real(mdDivF(TildeVar[3] - d.mus04, d.dmus4));
    } else {
        Tfrac = Real(mdDiv(TildeVar[0] - d.T04, d.dT4));
        mubfrac = Real(mdDiv(TildeVar[1] - d.mub04, d.dmub4));
        muqfrac = Real(mdDiv(TildeVar[2] - d.muq04, d.dmuq4));
        musfrac = Real(mdDiv(TildeVar[3] - d.mus04, d.dmus4));
    }

    int iT = int(Tfrac);
    iT = dmax(0, dmin(d.nT4 - 2, iT));
    int ib = int(mubfrac);
    ib = dmax(0, dmin(d.nMub4 - 2, ib));
    int iq = int(muqfrac);
    iq = dmax(0, dmin(d.nMuq4 - 2, iq));
    int is = int(musfrac);
    is = dmax(0, dmin(d.nMus4 - 2, is));

    float dx = float(dmax(Real(0.), dmin(Real(1.), Tfrac - Real(iT))));
    float dy = float(dmax(Real(0.), dmin(Real(1.), mubfrac - Real(ib))));
    float dz = float(dmax(Real(0.), dmin(Real(1.), muqfrac - Real(iq))));
    float dt = float(dmax(Real(0.), dmin(Real(1.), musfrac - Real(is))));

    c.w0000 = (1 - dx) * (1 - dy) * (1 - dz) * (1 - dt);
    c.w1111 = dx * dy * dz * dt;

    c.w1000 = dx * (1 - dy) * (1 - dz) * (1 - dt);
    c.w0100 = (1 - dx) * dy * (1 - dz) * (1 - dt);
    c.w0010 = (1 - dx) * (1 - dy) * dz * (1 - dt);
    c.w0001 = (1 - dx) * (1 - dy) * (1 - dz) * dt;

    c.w1001 = dx * (1 - dy) * (1 - dz) * dt;
    c.w0101 = (1 - dx) * dy * (1 - dz) * dt;
    c.w0011 = (1 - dx) * (1 - dy) * dz * dt;
    c.w1100 = dx * dy * (1 - dz) * (1 - dt);
    c.w1010 = dx * (1 - dy) * dz * (1 - dt);
    c.w0110 = (1 - dx) * dy * dz * (1 - dt);

    c.w0111 = (1 - dx) * dy * dz * dt;
    c.w1011 = dx * (1 - dy) * dz * dt;
    c.w1101 = dx * dy * (1 - dz) * dt;
    c.w1110 = dx * dy * dz * (1 - dt);

    c.iT = iT;
    // EOS_4D::index
    c.base = ((iT * d.nMus4 + is) * d.nMuq4 + iq) * d.nMub4 + ib;
    return c;
}

//! EOS_4D::FourDEval: res[0] the value, res[1..4] its derivatives in
//! (e, n_B, n_Q, n_S) if derivatives
template <typename Real>
MD_FN void eos4DEval(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab,
    MD_GLOBAL const float *data, MD_THREAD const float *TildeVar,
    MD_THREAD const Eos4DCell &c, MD_THREAD float *ResArr, bool derivatives) {
    const float w0000 = c.w0000;
    const float w1111 = c.w1111;
    const float w1000 = c.w1000;
    const float w0100 = c.w0100;
    const float w0010 = c.w0010;
    const float w0001 = c.w0001;
    const float w1001 = c.w1001;
    const float w0101 = c.w0101;
    const float w0011 = c.w0011;
    const float w1100 = c.w1100;
    const float w1010 = c.w1010;
    const float w0110 = c.w0110;
    const float w0111 = c.w0111;
    const float w1011 = c.w1011;
    const float w1101 = c.w1101;
    const float w1110 = c.w1110;
    const int iT = c.iT;
    const int oB = 1;
    const int oQ = d.nMub4;
    const int oS = d.nMub4 * d.nMuq4;
    const int oT = d.nMub4 * d.nMuq4 * d.nMus4;
    MD_GLOBAL const float *d0 = data + c.base;
    float data_0000 = d0[0];
    float data_1111 = d0[oT + oB + oQ + oS];

    float data_1000 = d0[oT];
    float data_0100 = d0[oB];
    float data_0010 = d0[oQ];
    float data_0001 = d0[oS];

    float data_1001 = d0[oT + oS];
    float data_0101 = d0[oB + oS];
    float data_0011 = d0[oQ + oS];
    float data_1100 = d0[oT + oB];
    float data_1010 = d0[oT + oQ];
    float data_0110 = d0[oB + oQ];

    float data_0111 = d0[oB + oQ + oS];
    float data_1011 = d0[oT + oQ + oS];
    float data_1101 = d0[oT + oB + oS];
    float data_1110 = d0[oT + oB + oQ];

    float interpolated_value =
        (w0000 * data_0000 + w1111 * data_1111 + w1000 * data_1000
         + w0100 * data_0100 + w0010 * data_0010 + w0001 * data_0001
         + w1001 * data_1001 + w0101 * data_0101 + w0011 * data_0011
         + w1100 * data_1100 + w1010 * data_1010 + w0110 * data_0110
         + w0111 * data_0111 + w1011 * data_1011 + w1101 * data_1101
         + w1110 * data_1110);

    ResArr[0] = interpolated_value;

    if (!derivatives) return;

    float wT000 = w0000 + w1000;
    float wT111 = w0111 + w1111;
    float wT100 = w0100 + w1100;
    float wT110 = w0110 + w1110;
    float wT010 = w0010 + w1010;
    float wT011 = w0011 + w1011;
    float wT001 = w0001 + w1001;
    float wT101 = w0101 + w1101;

    float tempT1 =
        (wT000 * data_0000 + wT100 * data_0100 + wT010 * data_0010
         + wT001 * data_0001 + wT101 * data_0101 + wT011 * data_0011
         + wT110 * data_0110 + wT111 * data_0111);

    float tempT2 =
        (wT111 * data_1111 + wT000 * data_1000 + wT001 * data_1001
         + wT100 * data_1100 + wT010 * data_1010 + wT011 * data_1011
         + wT101 * data_1101 + wT110 * data_1110);

    // after eos4DDilute TildeVar[0] >= Tref: T lies in
    // [max(T04, min(Tmax4, Tref)), max(T04, Tmax4)], checked at the upload
    // for the divisions by T (mdDivF), as are T5diff and the spacings; the
    // tables are finite (checked), so are the numerators
    float T = dmax(d.T04, dmin(d.Tmax4, TildeVar[0]));
    const Real Tsq = Real(T) * Real(T);
    MD_GLOBAL const Real *T5diff = tab + d.t5Offset;
    float dXdTtilde =
        float(mdDivF(Real(5 * (tempT2 - tempT1)), T5diff[iT]) * (Tsq * Tsq));

    float wb000 = w0000 + w0100;
    float wb111 = w1011 + w1111;
    float wb100 = w1000 + w1100;
    float wb110 = w1010 + w1110;
    float wb010 = w0010 + w0110;
    float wb011 = w0011 + w0111;
    float wb001 = w0001 + w0101;
    float wb101 = w1001 + w1101;

    float tempb1 =
        (wb000 * data_0000 + wb100 * data_1000 + wb010 * data_0010
         + wb001 * data_0001 + wb101 * data_1001 + wb011 * data_0011
         + wb110 * data_1010 + wb111 * data_1011);

    float tempb2 =
        (wb111 * data_1111 + wb000 * data_0100 + wb001 * data_0101
         + wb100 * data_1100 + wb010 * data_0110 + wb011 * data_0111
         + wb101 * data_1101 + wb110 * data_1110);

    float dXdmubtilde = mdDivF(tempb2 - tempb1, d.dmub4);

    float wq000 = w0000 + w0010;
    float wq111 = w1101 + w1111;
    float wq100 = w1000 + w1010;
    float wq110 = w1100 + w1110;
    float wq010 = w0100 + w0110;
    float wq011 = w0101 + w0111;
    float wq001 = w0001 + w0011;
    float wq101 = w1001 + w1011;

    float tempq1 =
        (wq000 * data_0000 + wq100 * data_1000 + wq010 * data_0100
         + wq001 * data_0001 + wq101 * data_1001 + wq011 * data_0101
         + wq110 * data_1100 + wq111 * data_1101);

    float tempq2 =
        (wq111 * data_1111 + wq000 * data_0010 + wq001 * data_0011
         + wq100 * data_1010 + wq010 * data_0110 + wq011 * data_0111
         + wq101 * data_1011 + wq110 * data_1110);

    float dXdmuqtilde = mdDivF(tempq2 - tempq1, d.dmuq4);

    float ws000 = w0000 + w0001;
    float ws111 = w1110 + w1111;
    float ws100 = w1000 + w1001;
    float ws110 = w1100 + w1101;
    float ws010 = w0100 + w0101;
    float ws011 = w0110 + w0111;
    float ws001 = w0010 + w0011;
    float ws101 = w1010 + w1011;

    float temps1 =
        (ws000 * data_0000 + ws100 * data_1000 + ws010 * data_0100
         + ws001 * data_0010 + ws110 * data_1100 + ws101 * data_1010
         + ws011 * data_0110 + ws111 * data_1110);

    float temps2 =
        (ws111 * data_1111 + ws000 * data_0001 + ws100 * data_1001
         + ws010 * data_0101 + ws001 * data_0011 + ws011 * data_0111
         + ws101 * data_1011 + ws110 * data_1101);

    float dXdmustilde = mdDivF(temps2 - temps1, d.dmus4);

    if (TildeVar[0] > d.Tmax4) dXdTtilde = 0.f;
    if (TildeVar[1] < d.mub04 || TildeVar[1] > d.mubMax4) dXdmubtilde = 0.f;
    if (TildeVar[2] < d.muq04 || TildeVar[2] > d.muqMax4) dXdmuqtilde = 0.f;
    if (TildeVar[3] < d.mus04 || TildeVar[3] > d.musMax4) dXdmustilde = 0.f;

    const Real pi = Real(3.14159265358979323846);
    const Real T0sq = Real(TildeVar[0]) * Real(TildeVar[0]);
    const Real eLocal = Real(3.) * d.alphaNf4 * (T0sq * T0sq);
    ResArr[1] = float(
        mdDivF(Real(3.0),
               Real(19.0) * pi * pi * Real(T) * Real(T) * Real(T))
            * Real(dXdTtilde)
        - mdDiv(Real(0.5)
                    * Real(TildeVar[1] * dXdmubtilde
                           + TildeVar[2] * dXdmuqtilde
                           + TildeVar[3] * dXdmustilde),
                dmax(eLocal, tinyNormal<Real>())));
    ResArr[2] = float(mdDivF(
        Real(5.0) * Real(dXdmubtilde) - Real(dXdmuqtilde)
         + Real(2.0) * Real(dXdmustilde),
        Real(T * T)));
    ResArr[3] = float(mdDivF(
        Real(-1.0) * Real(dXdmubtilde) + Real(2.0) * Real(dXdmuqtilde)
            - Real(dXdmustilde),
        Real(T * T)));
    ResArr[4] = float(mdDivF(
        Real(2.0) * Real(dXdmubtilde) - Real(dXdmuqtilde)
            + Real(2.0) * Real(dXdmustilde),
        Real(T * T)));
}

//! EOS_4D::pressureWithGradients
template <typename Real>
MD_FN void eos4DPressureGrad(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real rhob,
    Real rhoq, Real rhos, MD_THREAD const float *TildeVar, float scale,
    MD_THREAD const Eos4DCell &c, MD_THREAD Real &p, MD_THREAD Real &dpde,
    MD_THREAD Real &dpdrhob, MD_THREAD Real &dpdrhoq,
    MD_THREAD Real &dpdrhos) {
    const Real smallEps = Real(1e-16);
    float ResArr[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
    eos4DEval(d, tab, eos4DTable(d, tab, kEos4DP), TildeVar, c, ResArr, true);
    if (scale < 1.f) {
        const Real s2 = dmax(Real(scale) * Real(scale), tinyNormal<Real>());
        const Real Tsq = Real(TildeVar[0]) * Real(TildeVar[0]);
        const Real eref = Real(3.) * d.alphaNf4 * (Tsq * Tsq);
        const Real pref = ResArr[0];
        const Real ndpdn = mdDiv(
            rhob * Real(ResArr[2]) + rhoq * Real(ResArr[3])
             + rhos * Real(ResArr[4]),
            s2);
        p = dmax(smallEps, pref * s2 * s2);
        dpde = mdDiv(pref - Real(0.5) * ndpdn, eref);
        dpdrhob = s2 * Real(ResArr[2]);
        dpdrhoq = s2 * Real(ResArr[3]);
        dpdrhos = s2 * Real(ResArr[4]);
        return;
    }
    p = dmax(smallEps, Real(ResArr[0]));
    dpde = Real(ResArr[1]);
    dpdrhob = Real(ResArr[2]);
    dpdrhoq = Real(ResArr[3]);
    dpdrhos = Real(ResArr[4]);
}

//! EOS_4D::get_pressure_with_gradients_and_cs2
template <typename Real>
MD_FN EosGrad<Real> eos4DGrad(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob, Real rhoq, Real rhos) {
    const Real smallEps = Real(1e-16);
    float TildeVar[4];
    eos4DTilde(d, e, rhob, rhoq, rhos, TildeVar);
    const float scale = eos4DDilute(d, TildeVar);
    const Eos4DCell c = eos4DSetup(d, TildeVar);
    EosGrad<Real> g;
    Real dpdrhoq, dpdrhos;
    eos4DPressureGrad(
        d, tab, rhob, rhoq, rhos, TildeVar, scale, c, g.p, g.dpde, g.dpdrhob,
        dpdrhoq, dpdrhos);
    const Real enthalpy = dmax(e + g.p, Real(0.)) + smallEps;
    Real cs2 =
        g.dpde
        + mdDiv(rhob * g.dpdrhob + rhoq * dpdrhoq + rhos * dpdrhos, enthalpy);
    g.cs2 = dmax(Real(0.01), dmin(Real(1.) / Real(3), cs2));
    return g;
}

//! EOS_4D::get_pressure
template <typename Real>
MD_FN Real eos4DPressure(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob, Real rhoq, Real rhos) {
    const Real smallEps = Real(1e-16);
    float TildeVar[4];
    eos4DTilde(d, e, rhob, rhoq, rhos, TildeVar);
    const float scale = eos4DDilute(d, TildeVar);
    const Eos4DCell c = eos4DSetup(d, TildeVar);
    float ResArr[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
    eos4DEval(d, tab, eos4DTable(d, tab, kEos4DP), TildeVar, c, ResArr, false);
    ResArr[0] *= scale * scale * scale * scale;
    return (dmax(smallEps, Real(ResArr[0])));
}

//! EOS_base::get_rhob_edge
template <typename Real>
MD_FN Real eosRhobEdge(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e) {
    if (d.edgeOffset < 0) return (realMax<Real>());
    const int t = eosTableIdx(d, e);
    int base = d.edgeOffset;
    for (int k = 0; k < t; k++) base += d.eLength[k];
    MD_GLOBAL const Real *edge = tab + base;
    const int N_e = d.eLength[t];
    const Real x = (e - d.eBounds[t]) / d.eSpacing[t];
    if (x <= Real(0.)) return (edge[0]);
    const int j = int(x);
    if (j >= N_e - 1) {
        if (t + 1 >= d.nTables) return (edge[N_e - 1]);
        const Real f = dmin(Real(1.), x - Real(N_e - 1));
        return (edge[N_e - 1] * (Real(1.) - f) + edge[N_e] * f);
    }
    const Real f = x - Real(j);
    return (edge[j] * (Real(1.) - f) + edge[j + 1] * f);
}

//! EOS::getThermalVariables, the entries the viscous update reads
template <typename Real>
struct Thermal {
    Real e, rhob, p, dpde, dpdrhob, cs2, T, muB, s;
};

template <typename Real, int EosType = -1>
MD_FN EosGrad<Real> eosPressureGrad(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob, Real rhoq, Real rhos);

template <typename Real, int EosType = -1>
MD_FN Thermal<Real> eosThermal(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob, Real rhoq, Real rhos) {
    const int type = (EosType < 0 ? MD_EOS_TYPE(d) : EosType);
    const Real smallEps = Real(1e-16);
    Thermal<Real> t;
    t.e = e;
    t.rhob = rhob;
    if (type == 0) {
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
    if (type == 2) {
        // EOS_base::getThermalVariables with EOS_neos
        const EosGrad<Real> g =
            eosPressureGrad<Real, EosType>(d, tab, e, rhob, rhoq, rhos);
        t.p = g.p;
        t.dpde = g.dpde;
        t.dpdrhob = g.dpdrhob;
        t.cs2 = g.cs2;
        const Real nb = fabs(rhob);
        const int it = eosTableIdx(d, e);
        // EOS_neos::get_temperature
        Real T5 = eosInterp2D(d, tab + d.t2Offset, e, nb, it);
        T5 = dmax(smallEps, T5);
        t.T = Real(pow(T5, Real(0.2)));
        // EOS_neos::get_muB, get_muS, get_muQ
        const Real sign = rhob / (fabs(rhob) + smallEps);
        t.muB = sign * eosInterp2D(d, tab + d.muBOffset, e, nb, it);
        const Real muS =
            (d.muSOffset < 0
                 ? Real(0.)
                 : sign * eosInterp2D(d, tab + d.muSOffset, e, nb, it));
        const Real muQ =
            (d.muQOffset < 0
                 ? Real(0.)
                 : sign * eosInterp2D(d, tab + d.muQOffset, e, nb, it));
        // EOS_base::get_densities_in_table: n_B up to the table edge,
        // n_Q = 0.4 n_B, n_S = 0
        const Real nEdge = eosRhobEdge(d, tab, e);
        const Real rhobT = (fabs(rhob) <= nEdge ? rhob : copysign(nEdge, rhob));
        const Real rhoQ = Real(0.4) * rhobT, rhoS = Real(0.);
        Real entropy = ((e + t.p - t.muB * rhobT - muS * rhoS - muQ * rhoQ)
                        / (dmax(t.T, Real(0.)) + smallEps));
        t.s = dmax(smallEps, entropy);
        return (t);
    }
    if (type == 3) {
        // EOS_4D::getThermalVariables
        float TildeVar[4];
        eos4DTilde(d, e, rhob, rhoq, rhos, TildeVar);
        const float scale = eos4DDilute(d, TildeVar);
        const Eos4DCell c = eos4DSetup(d, TildeVar);
        Real dpdrhoq, dpdrhos;
        eos4DPressureGrad(
            d, tab, rhob, rhoq, rhos, TildeVar, scale, c, t.p, t.dpde,
            t.dpdrhob, dpdrhoq, dpdrhos);
        const Real enthalpy = dmax(e + t.p, Real(0.)) + smallEps;
        Real cs2 = t.dpde
                   + (rhob * t.dpdrhob + rhoq * dpdrhoq + rhos * dpdrhos)
                         / enthalpy;
        t.cs2 = dmax(Real(0.01), dmin(Real(1.) / Real(3), cs2));
        float ResArr[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
        eos4DEval(
            d, tab, eos4DTable(d, tab, kEos4DT), TildeVar, c, ResArr, false);
        ResArr[0] *= scale;
        t.T = dmax(smallEps, Real(ResArr[0]));
        eos4DEval(
            d, tab, eos4DTable(d, tab, kEos4DMuB), TildeVar, c, ResArr, false);
        ResArr[0] *= scale;
        t.muB = Real(ResArr[0]);
        eos4DEval(
            d, tab, eos4DTable(d, tab, kEos4DMuS), TildeVar, c, ResArr, false);
        ResArr[0] *= scale;
        const Real muS = Real(ResArr[0]);
        eos4DEval(
            d, tab, eos4DTable(d, tab, kEos4DMuQ), TildeVar, c, ResArr, false);
        ResArr[0] *= scale;
        const Real muQ = Real(ResArr[0]);
        Real entropy =
            ((e + t.p - t.muB * rhob - muS * rhos - muQ * rhoq)
             / (dmax(t.T, Real(0.)) + smallEps));
        t.s = dmax(smallEps, entropy);
        return (t);
    }
    // EOS_base::getThermalVariables with EOS_hotQCD: mu = 0
    const EosGrad<Real> g =
        eosPressureGrad<Real, EosType>(d, tab, e, rhob, rhoq, rhos);
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
                    / (dmax(t.T, Real(0.)) + smallEps));
    t.s = dmax(smallEps, entropy);
    return (t);
}

//! EOS::get_pressure
//! A negative EosType uses the descriptor or the backend's run specialization.
template <typename Real, int EosType = -1>
MD_FN Real eosPressure(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob, Real rhoq, Real rhos) {
    const int type = (EosType < 0 ? MD_EOS_TYPE(d) : EosType);
    const Real smallEps = Real(1e-16);
    if (type == 0) return Real(1) / Real(3) * e;
    if (type == 3) return eos4DPressure(d, tab, e, rhob, rhoq, rhos);
    if (type == 1) return dmax(smallEps, eosInterp1D(d, tab, e));
    const int t = eosTableIdx(d, e);
    return dmax(smallEps, eosInterp2D(d, tab, e, fabs(rhob), t));
}

//! EOS::get_pressure_with_gradients_and_cs2 (n_Q, n_S enter only the 4D
//! EoS; for the others dP/dn_Q = dP/dn_S = 0)
template <typename Real, int EosType>
MD_FN EosGrad<Real> eosPressureGrad(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob, Real rhoq, Real rhos) {
    const int type = (EosType < 0 ? MD_EOS_TYPE(d) : EosType);
    const Real smallEps = Real(1e-16);
    EosGrad<Real> g;
    if (type == 0) {
        g.dpde = Real(1) / Real(3);
        g.p = g.dpde * e;
        g.dpdrhob = Real(0);
        g.cs2 = g.dpde;
        return g;
    }
    if (type == 3) return eos4DGrad(d, tab, e, rhob, rhoq, rhos);
    if (type == 1) {
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
                     / dmax(e + g.p, tinyNormal<Real>());
    g.cs2 = dmax(Real(0.01), dmin(Real(1) / Real(3), cs2));
    return g;
}

//! EOS::get_temperature and EOS::get_muB (the mu_B / T gradients of the
//! baryon diffusion)
template <typename Real, int EosType = -1>
MD_FN void eosTMuB(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob, Real rhoq, Real rhos, MD_THREAD Real &T, MD_THREAD Real &muB) {
    const int type = (EosType < 0 ? MD_EOS_TYPE(d) : EosType);
    const Real smallEps = Real(1e-16);
    if (type == 0) {
        // EOS_idealgas
        const Real x = idealMuBOverT(d, tab, e, rhob);
        T = idealTofEX(d, e, x);
        muB = copysign(x * idealTofEX(d, e, x), rhob);
        return;
    }
    if (type == 1) {
        // EOS_hotQCD: mu_B = 0
        const Real T5 = eosInterp1DAt(d, tab, d.tOffset, e);
        T = dmax(smallEps, Real(pow(T5, Real(0.2))));
        muB = Real(0.);
        return;
    }
    if (type == 2) {
        // EOS_neos
        const int it = eosTableIdx(d, e);
        const Real nb = fabs(rhob);
        Real T5 = eosInterp2D(d, tab + d.t2Offset, e, nb, it);
        T5 = dmax(smallEps, T5);
        T = Real(pow(T5, Real(0.2)));
        const Real sign = rhob / (fabs(rhob) + smallEps);
        muB = sign * eosInterp2D(d, tab + d.muBOffset, e, nb, it);
        return;
    }
    // EOS_4D
    float TildeVar[4];
    eos4DTilde(d, e, rhob, rhoq, rhos, TildeVar);
    const float scale = eos4DDilute(d, TildeVar);
    const Eos4DCell c = eos4DSetup(d, TildeVar);
    float ResArr[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
    eos4DEval(d, tab, eos4DTable(d, tab, kEos4DT), TildeVar, c, ResArr, false);
    ResArr[0] *= scale;
    T = dmax(smallEps, Real(ResArr[0]));
    eos4DEval(d, tab, eos4DTable(d, tab, kEos4DMuB), TildeVar, c, ResArr, false);
    ResArr[0] *= scale;
    muB = Real(ResArr[0]);
}

//! EOS::get_cs2, which Advance::nCausalityConstraints reads: for the
//! tables EOS_base::calculate_velocity_of_sound_sq with the finite
//! differences of get_dpOverde3 and get_dpOverdrhob2 (NEOS; hotQCD has
//! dP/dn_B = 0), clamped to [0.01, 1]
template <typename Real, int EosType = -1>
MD_FN Real eosCs2(
    MD_CONST const EosDesc<Real> &d, MD_GLOBAL const Real *tab, Real e,
    Real rhob, Real rhoq, Real rhos) {
    const int type = (EosType < 0 ? MD_EOS_TYPE(d) : EosType);
    const Real smallEps = Real(1e-16);
    if (type == 0) return (Real(1.) / Real(3.));
    if (type == 3) return (eos4DGrad(d, tab, e, rhob, rhoq, rhos).cs2);
    // EOS_base::get_dpOverde3
    const Real de = dmax(Real(0.01), Real(0.01) * e);
    const Real eLeft = dmax(Real(1e-16), e - de);
    const Real eRight = e + de;
    const Real pL = eosPressure<Real, EosType>(
        d, tab, eLeft, rhob, Real(0.), Real(0.));
    const Real pR = eosPressure<Real, EosType>(
        d, tab, eRight, rhob, Real(0.), Real(0.));
    const Real dpde = (pR - pL) / (eRight - eLeft);
    Real dpdrho = Real(0.);
    if (type == 2) {
        // EOS_base::get_dpOverdrhob2
        const int t = eosTableIdx(d, e);
        const Real deltaRhob = d.nbSpacing[t];
        const Real rhobLeft = rhob - deltaRhob * Real(0.5);
        const Real rhobRight = rhob + deltaRhob * Real(0.5);
        const Real qL = eosPressure<Real, EosType>(
            d, tab, e, rhobLeft, Real(0.), Real(0.));
        const Real qR = eosPressure<Real, EosType>(
            d, tab, e, rhobRight, Real(0.), Real(0.));
        dpdrho = (qR - qL) / (rhobRight - rhobLeft);
    }
    const Real pressure = eosPressure<Real, EosType>(
        d, tab, e, rhob, Real(0.), Real(0.));
    Real v_sound = dpde + rhob / (dmax(e + pressure, Real(0.)) + smallEps) * dpdrho;
    v_sound = dmax(Real(0.01), dmin(Real(1.0), v_sound));
    return (v_sound);
}

}  // namespace music_device

#endif  // SRC_DEVICE_DEVICE_EOS_H_
