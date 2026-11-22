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

namespace music_device {

#if defined(__METAL_VERSION__)
using metal::acos;
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
__device__ inline void flagOr(MD_FLAGS f, int v) { atomicOr(f, v); }
#else
using std::acos;
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
};
enum { kFatalFlags = kFlagMaxSpeedFatal };

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
