// Copyright 2026 Chun Shen
//
// Portability layer of the device kernels (src/device/*.h). The same
// headers compile as
//   C++11 (host, OpenMP; the reference backend),
//   CUDA (nvcc, or clang -x cuda),
//   Metal Shading Language (bundled into one source string at build time
//   and compiled at run time by the Metal backend).
// The kernels are templates on the floating type Real, use only POD data and
// fixed-size arrays, no virtual functions, no recursion, no standard
// library. Errors are reported as flag bits (see DeviceFlags).

#ifndef SRC_DEVICE_DEVICE_DEFS_H_
#define SRC_DEVICE_DEVICE_DEFS_H_

#if defined(__METAL_VERSION__)
#include <metal_stdlib>
#define MD_FN inline
#define MD_GLOBAL device
#define MD_CONST constant
#define MD_THREAD thread
#define MD_FLAGS device metal::atomic_int *
#elif defined(__CUDACC__)
#define MD_FN __host__ __device__ inline
#define MD_GLOBAL
#define MD_CONST
#define MD_THREAD
#define MD_FLAGS int *
#else
#include <cmath>
#define MD_FN inline
#define MD_GLOBAL
#define MD_CONST
#define MD_THREAD
#define MD_FLAGS int *
#endif

// The switches of a run (StepParams and EosDesc fields that are fixed for
// the run). The Metal backend compiles the kernels for one run with
// MUSIC_SPEC and MUSIC_SPEC_<name> defined (DeviceHydro::specKey), so the
// code of the features that are off is removed; elsewhere they are read at
// run time.
#if defined(__METAL_VERSION__) && defined(MUSIC_SPEC)
#define MD_RUN(P, name) (MUSIC_SPEC_##name)
#define MD_EOS_TYPE(d) (MUSIC_SPEC_eosType)
#else
#define MD_RUN(P, name) ((P).name)
#define MD_EOS_TYPE(d) ((d).type)
#endif

namespace music_device {

#if defined(__METAL_VERSION__)
using metal::acos;
using metal::acosh;
using metal::asinh;
using metal::atan2;
using metal::cosh;
using metal::sin;
using metal::sinh;
using metal::atanh;
using metal::copysign;
using metal::cos;
using metal::exp;
using metal::fabs;
using metal::isfinite;
using metal::isnan;
using metal::log;
using metal::pow;
using metal::sqrt;
using metal::tanh;
MD_FN void flagOr(MD_FLAGS f, int v) {
    metal::atomic_fetch_or_explicit(f, v, metal::memory_order_relaxed);
}
#elif defined(__CUDACC__)
// host and device, as the MD_FN kernels that call it (nvcc rejects a call
// of a __device__ function from a __host__ __device__ one)
__host__ __device__ inline void flagOr(MD_FLAGS f, int v) {
#ifdef __CUDA_ARCH__
    atomicOr(f, v);
#else
    __atomic_fetch_or(f, v, __ATOMIC_RELAXED);
#endif
}
#else
using std::acos;
using std::acosh;
using std::asinh;
using std::atan2;
using std::cosh;
using std::sin;
using std::sinh;
using std::atanh;
using std::copysign;
using std::cos;
using std::exp;
using std::fabs;
using std::isfinite;
using std::isnan;
using std::log;
using std::pow;
using std::sqrt;
using std::tanh;
inline void flagOr(MD_FLAGS f, int v) {
    __atomic_fetch_or(f, v, __ATOMIC_RELAXED);
}
#endif

//! mdDivF(a, b) and mdDiv(a, b) = a / b. On Metal the division of
//! MTLMathModeSafe is the correctly rounded IEEE division; these give the
//! same bits at about a quarter of its cost (13.5 against 3.1 ns per
//! dependent operation on the M4): the fast reciprocal, one Newton step on
//! 1/b and Markstein's correction of the quotient, with a numerator of
//! exponent field 1 ... 25 scaled by 2^64 (its remainder would underflow).
//! mdDiv also gives +-inf, NaN and the divisors 0 and +-inf as the division
//! does. Checked on the M4 against a / b: all 2^46 pairs of float mantissas
//! (exhaustive) and every pair of exponent fields (4096 random mantissas
//! and signs each), with no difference inside the preconditions:
//!   mdDivF: a not +-inf, b normal with |b| < 2^126, |a / b| < 2^127;
//!   mdDiv: |b| < 2^126 or b = +-inf, |a / b| < 2^127.
//! Elsewhere (CPU, CUDA) both are the division.
#if defined(__METAL_VERSION__)
MD_FN float mdDivF(float a, float b) {
    const uint ua = as_type<uint>(a) & 0x7fffffffu;
    const bool tiny = (ua - 0x00800000u) < 0x0C800000u;
    const float as = a * (tiny ? 0x1p64f : 1.0f);
    float y = metal::fast::divide(1.0f, b);
    const float e = metal::fma(-b, y, 1.0f);
    y = metal::fma(e, y, y);
    const float q = as * y;
    const float r = metal::fma(-b, q, as);
    // a = 0: q has the sign of a / b, the correction would give +0 for -0
    return ((a == 0.0f ? q : metal::fma(r, y, q))
            * (tiny ? 0x1p-64f : 1.0f));
}
MD_FN float mdDiv(float a, float b) {
    const uint ua = as_type<uint>(a) & 0x7fffffffu;
    const uint ub = as_type<uint>(b) & 0x7fffffffu;
    const uint sgn = (as_type<uint>(a) ^ as_type<uint>(b)) & 0x80000000u;
    float l = mdDivF(a, b);
    if (ua >= 0x7f800000u) l = a * metal::copysign(1.0f, b);
    if (b == 0.0f) {
        l = (a == 0.0f || ua > 0x7f800000u ? NAN
                                           : as_type<float>(sgn | 0x7f800000u));
    }
    if (ub == 0x7f800000u)
        l = (ua >= 0x7f800000u ? NAN : as_type<float>(sgn));
    if (ub > 0x7f800000u) l = b;
    return (l);
}
#else
template <typename T>
MD_FN T mdDivF(T a, T b) {
    return (a / b);
}
template <typename T>
MD_FN T mdDiv(T a, T b) {
    return (a / b);
}
#endif

//! erf: MSL has none. Float there: erfc of Numerical Recipes (erfcc, a
//! Chebyshev fit with relative error below 1.2e-7 everywhere), so erf has an
//! absolute error below 1.2e-7, the rounding of float near 1
#if defined(__METAL_VERSION__)
MD_FN float derfcPositive(float z) {
    const float t = 1.f / (1.f + 0.5f * z);
    const float ans =
        t
        * exp(-z * z - 1.26551223f
              + t * (1.00002368f
                     + t * (0.37409196f
                            + t * (0.09678418f
                                   + t * (-0.18628806f
                                          + t * (0.27886807f
                                                 + t * (-1.13520398f
                                                        + t * (1.48851587f
                                                               + t * (-0.82215223f
                                                                      + t * 0.17087277f)))))))));
    return ans;
}
MD_FN float derf(float x) {
    const float ans = derfcPositive(fabs(x));
    return (x >= 0.f ? 1.f - ans : ans - 1.f);
}
#elif defined(__CUDACC__)
MD_FN double derfcPositive(double x) { return erfc(x); }
MD_FN float derfcPositive(float x) { return erfc(x); }
MD_FN double derf(double x) { return erf(x); }
MD_FN float derf(float x) { return erf(x); }
#else
inline double derfcPositive(double x) { return std::erfc(x); }
inline float derfcPositive(float x) { return std::erfc(x); }
inline double derf(double x) { return std::erf(x); }
inline float derf(float x) { return std::erf(x); }
#endif

//! std::numeric_limits<Real>::max(); MSL has no double literals
template <typename Real>
MD_FN Real realMax() {
#if defined(__METAL_VERSION__)
    return (FLT_MAX);
#else
    return (sizeof(Real) == 8 ? Real(1.7976931348623157e308)
                              : Real(3.40282346638528859812e38f));
#endif
}

//! the smallest normal Real (Util::tiny_normal in double): a floor for a
//! denominator that is >= 0 by construction, which changes no result unless
//! it is 0 or subnormal
template <typename Real>
MD_FN Real tinyNormal() {
#if defined(__METAL_VERSION__)
    return (FLT_MIN);
#else
    return (sizeof(Real) == 8 ? Real(2.2250738585072014e-308)
                              : Real(1.17549435082228751e-38f));
#endif
}

//! exp(y) - 1 without the cancellation at small y (MSL has no expm1):
//! Kahan's (u - 1) y / log(u), u = exp(y)
template <typename Real>
MD_FN Real dexpm1(Real y) {
    const Real u = exp(y);
    if (u == Real(1)) return (y);
    const Real um1 = u - Real(1);
    if (um1 == Real(-1)) return (Real(-1));
    return (um1 * y / log(u));
}

// std::max / std::min semantics (the second argument wins only if it is
// strictly larger / smaller), so that NaNs and ties propagate as on the host
template <typename T>
MD_FN T dmax(T a, T b) {
    return (a < b) ? b : a;
}
template <typename T>
MD_FN T dmin(T a, T b) {
    return (b < a) ? b : a;
}

//! error and event flags of a step, ORed over all cells and faces
enum DeviceFlags {
    kFlagMaxSpeedFatal = 1,  // MaxSpeed: a case where MUSIC exits
    kFlagRevert = 2,         // a velocity solve failed, cell reverted
    kFlagRegulate = 4,       // T^{tau tau} < abs_err, cell set to rest
    kFlagMomentumCap = 8,    // superluminal T^{tau mu} capped
    kFlagHybridFail = 16,    // the bracketing velocity solver found no root
    kFlagSourceNan = 32,     // qi_source is nan (strings source on the device)
};
enum { kFatalFlags = kFlagMaxSpeedFatal | kFlagSourceNan };

//! number of components of a device state slice: e, n_B, n_Q, n_S,
//! u^tau, u^x, u^y, u^eta, and (with storeQ) tau T^{tau mu}, tau J^tau_{B,Q,S}
// enums: Metal allows no const variables at program scope
enum { kNPrim = 8, kNComp = 15, kNTJ = 7 };

//! viscous runs append W^{mu nu} (14 components, as Fields::Wmunu_), Pi and
//! Pi_chem to the slice: kNCompVisc components
enum { kOffW = 15, kOffPi = 29, kOffPiChem = 30, kNCompVisc = 31 };

//! values per face of the W kernel (Diss::computeWDivergences): the KT
//! fluxes of the up to 10 evolved W components, the face averages of
//! tau W^{d alpha} (alpha < 5) and of the bulk part (alpha < 4)
enum { kNWFlux = 10, kNWFaceVals = 19 };

}  // namespace music_device

#endif  // SRC_DEVICE_DEVICE_DEFS_H_
