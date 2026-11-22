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
using metal::fabs;
using metal::isfinite;
using metal::sqrt;
MD_FN void flagOr(MD_FLAGS f, int v) {
    metal::atomic_fetch_or_explicit(f, v, metal::memory_order_relaxed);
}
#elif defined(__CUDACC__)
__device__ inline void flagOr(MD_FLAGS f, int v) { atomicOr(f, v); }
#else
using std::fabs;
using std::isfinite;
using std::sqrt;
inline void flagOr(MD_FLAGS f, int v) {
    __atomic_fetch_or(f, v, __ATOMIC_RELAXED);
}
#endif

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

}  // namespace music_device

#endif  // SRC_DEVICE_DEVICE_DEFS_H_
