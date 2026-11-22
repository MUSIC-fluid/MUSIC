// Copyright 2026 Chun Shen
//
// The C interface of the GPU backends. Metal (metal_backend.mm, Apple clang)
// and CUDA (cuda_backend.cu, nvcc) are compiled by other compilers than the
// rest of MUSIC and may use another C++ standard library, so only C types
// cross this boundary. The structures behind the void pointers are
// music_device::StepParams<Real> and music_device::EosDesc<Real>, plain
// structs with the same layout in every compiler, with Real = float
// (realSize 4) or double (8).
//
// A backend holds three state slices (the current one, the predictor and
// the next, which at the start of a step holds the previous slice), the
// pressure of the first two, the face fluxes and the EoS table. A slice has
// nComp components of nx ny neta cells: kNComp (device_ideal.h), or
// kNCompVisc with the viscous update (device_visc.h), which also needs the
// W face values. state(0) is a host-visible copy of the current slice,
// state(1) of the previous one (read by the viscous update only):
// upload() makes them the device state, step() advances it by one Heun
// step (two RK stages), download() makes state(0) the new current slice. On Apple's unified
// memory state() is the device buffer itself and upload/download copy
// nothing. Functions returning int return 0 on success; on failure they
// write a message to err.

#ifndef SRC_DEVICE_DEVICE_API_H_
#define SRC_DEVICE_DEVICE_API_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

//! byte offset of StepParams<Real>::rkFlag (after five ints, for both Real;
//! checked in device_hydro.cpp): the backends set it for each stage
#define MUSIC_DEVICE_RKFLAG_OFFSET 16

typedef struct MusicDeviceApi {
    void *(*create)(int realSize, char *err, int errLen);
    void (*destroy)(void *dev);
    const char *(*name)(void *dev);
    int (*resize)(
        void *dev, int nx, int ny, int neta, int nComp, char *err, int errLen);
    int (*setEos)(
        void *dev, const void *desc, size_t descBytes, const void *table,
        size_t tableBytes, char *err, int errLen);
    void *(*state)(void *dev, int which);
    int (*upload)(void *dev, char *err, int errLen);
    // params: StepParams<Real> of the first stage (rkFlag is set inside;
    // nComp selects the ideal or the viscous update); flags: the OR of the
    // DeviceFlags of both stages
    int (*step)(
        void *dev, const void *params, size_t paramBytes, int *flags,
        char *err, int errLen);
    int (*download)(void *dev, char *err, int errLen);
} MusicDeviceApi;

//! the backends compiled into this build (NULL if not)
const MusicDeviceApi *music_metal_api(void);
const MusicDeviceApi *music_cuda_api(void);

#ifdef __cplusplus
}
#endif

#endif  // SRC_DEVICE_DEVICE_API_H_
