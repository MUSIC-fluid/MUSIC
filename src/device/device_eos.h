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
        const Real dx =
            -(N3 * N - C * D2 * D) / (Real(4.) * N3 * dN - Real(3.) * C * D2 * N);
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
    Real Ttilde_sq = sqrt(e / Real(3.0) * d.oneOverAlphaNf4);
    Real Ttilde = sqrt(Ttilde_sq);
    Real mubtilde = (Real(5) * rhob - rhoq + Real(2) * rhos) / Ttilde_sq;
    Real muqtilde = (-rhob + Real(2) * rhoq - rhos) / Ttilde_sq;
    Real mustilde = (Real(2) * rhob - rhoq + Real(2) * rhos) / Ttilde_sq;
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
    const float scale = TildeVar[0] / Tref;
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
    Real Tfrac = Real((TildeVar[0] - d.T04) / d.dT4);
    Real mubfrac = Real((TildeVar[1] - d.mub04) / d.dmub4);
    Real muqfrac = Real((TildeVar[2] - d.muq04) / d.dmuq4);
    Real musfrac = Real((TildeVar[3] - d.mus04) / d.dmus4);

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

    float T = dmax(d.T04, dmin(d.Tmax4, TildeVar[0]));
    const Real Tsq = Real(T) * Real(T);
    MD_GLOBAL const Real *T5diff = tab + d.t5Offset;
    float dXdTtilde =
        float(Real(5 * (tempT2 - tempT1)) / T5diff[iT] * (Tsq * Tsq));

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

    float dXdmubtilde = (tempb2 - tempb1) / d.dmub4;

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

    float dXdmuqtilde = (tempq2 - tempq1) / d.dmuq4;

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

    float dXdmustilde = (temps2 - temps1) / d.dmus4;

    if (TildeVar[0] > d.Tmax4) dXdTtilde = 0.f;
    if (TildeVar[1] < d.mub04 || TildeVar[1] > d.mubMax4) dXdmubtilde = 0.f;
    if (TildeVar[2] < d.muq04 || TildeVar[2] > d.muqMax4) dXdmuqtilde = 0.f;
    if (TildeVar[3] < d.mus04 || TildeVar[3] > d.musMax4) dXdmustilde = 0.f;

    const Real pi = Real(3.14159265358979323846);
    const Real T0sq = Real(TildeVar[0]) * Real(TildeVar[0]);
    const Real eLocal = Real(3.) * d.alphaNf4 * (T0sq * T0sq);
    ResArr[1] = float(
        Real(3.0) / (Real(19.0) * pi * pi * Real(T) * Real(T) * Real(T))
            * Real(dXdTtilde)
        - Real(0.5)
              * Real(TildeVar[1] * dXdmubtilde + TildeVar[2] * dXdmuqtilde
                     + TildeVar[3] * dXdmustilde)
              / eLocal);
    ResArr[2] = float(
        (Real(5.0) * Real(dXdmubtilde) - Real(dXdmuqtilde)
         + Real(2.0) * Real(dXdmustilde))
        / Real(T * T));
    ResArr[3] = float(
        (Real(-1.0) * Real(dXdmubtilde) + Real(2.0) * Real(dXdmuqtilde)
         - Real(dXdmustilde))
        / Real(T * T));
    ResArr[4] = float(
        (Real(2.0) * Real(dXdmubtilde) - Real(dXdmuqtilde)
         + Real(2.0) * Real(dXdmustilde))
        / Real(T * T));
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
        const Real s2 = Real(scale) * Real(scale);
        const Real Tsq = Real(TildeVar[0]) * Real(TildeVar[0]);
        const Real eref = Real(3.) * d.alphaNf4 * (Tsq * Tsq);
        const Real pref = ResArr[0];
        const Real ndpdn =
            (rhob * Real(ResArr[2]) + rhoq * Real(ResArr[3])
             + rhos * Real(ResArr[4]))
            / s2;
        p = dmax(smallEps, pref * s2 * s2);
        dpde = (pref - Real(0.5) * ndpdn) / eref;
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
    const Real enthalpy = e + g.p + smallEps;
    Real cs2 =
        g.dpde + (rhob * g.dpdrhob + rhoq * dpdrhoq + rhos * dpdrhos) / enthalpy;
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
    if (d.type == 3) {
        // EOS_4D::getThermalVariables
        float TildeVar[4];
        eos4DTilde(d, e, rhob, rhoq, rhos, TildeVar);
        const float scale = eos4DDilute(d, TildeVar);
        const Eos4DCell c = eos4DSetup(d, TildeVar);
        Real dpdrhoq, dpdrhos;
        eos4DPressureGrad(
            d, tab, rhob, rhoq, rhos, TildeVar, scale, c, t.p, t.dpde,
            t.dpdrhob, dpdrhoq, dpdrhos);
        const Real enthalpy = e + t.p + smallEps;
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
             / (t.T + smallEps));
        t.s = dmax(smallEps, entropy);
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
    Real rhob, Real rhoq, Real rhos) {
    const Real smallEps = Real(1e-16);
    if (d.type == 0) return Real(1) / Real(3) * e;
    if (d.type == 3) return eos4DPressure(d, tab, e, rhob, rhoq, rhos);
    if (d.type == 1) return dmax(smallEps, eosInterp1D(d, tab, e));
    const int t = eosTableIdx(d, e);
    return dmax(smallEps, eosInterp2D(d, tab, e, fabs(rhob), t));
}

//! EOS::get_pressure_with_gradients_and_cs2 (n_Q, n_S enter only the 4D
//! EoS; for the others dP/dn_Q = dP/dn_S = 0)
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
    if (d.type == 3) return eos4DGrad(d, tab, e, rhob, rhoq, rhos);
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
