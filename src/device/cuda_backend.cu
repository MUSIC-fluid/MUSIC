// Copyright 2026 Chun Shen
//
// gpu_backend = 3: the device kernels on an NVIDIA GPU through CUDA, in
// double (default) or float (gpu_float = 1), behind the C interface of
// device_api.h (built by nvcc with CMake option cuda=ON).
//
// The state stays in device memory; state(0), state(1) and state(2) are
// pinned host copies of the current, the previous (viscous update) and the
// predictor slice. upload() sends the first two to the device; download()
// fills state(0) after a step, fetch() components of a slice between the
// stages. One RK stage is five kernel launches on the default stream
// (pressure, three face directions and the cells), eight with the viscous
// update (also the W faces of the three directions).
//
// nvcc contracts multiply-adds by default (--fmad=true), as g++ does on
// the host: in double the CUDA backend agrees with MUSIC to rounding
// (1e-15 relative per step), not bit for bit; nvcc --fmad=false comes
// closer.
//
// Without an NVIDIA toolchain the source is type-checked by clang's CUDA
// front end (MUSIC_CUDA_SYNTAX_ONLY stubs the runtime API):
//   clang++ -x cuda -nocudainc -nocudalib --cuda-gpu-arch=sm_80 \
//       -std=c++11 -DMUSIC_CUDA_SYNTAX_ONLY -fsyntax-only \
//       [--cuda-host-only | --cuda-device-only] cuda_backend.cu

#ifdef MUSIC_CUDA_SYNTAX_ONLY
#define __host__ __attribute__((host))
#define __device__ __attribute__((device))
#define __global__ __attribute__((global))
#define __shared__ __attribute__((shared))
#define __constant__ __attribute__((constant))
#include <__clang_cuda_builtin_vars.h>
#include <stddef.h>
struct dim3 {
    unsigned x, y, z;
    __host__ __device__ dim3(unsigned a = 1, unsigned b = 1, unsigned c = 1)
        : x(a), y(b), z(c) {}
};
typedef struct CUstream_st *cudaStream_t;
extern "C" int cudaConfigureCall(dim3, dim3, size_t = 0, cudaStream_t = 0);
extern "C" unsigned __cudaPushCallConfiguration(
    dim3, dim3, size_t = 0, void * = 0);
typedef int cudaError_t;
enum { cudaSuccess = 0 };
enum cudaMemcpyKind { cudaMemcpyHostToDevice = 1, cudaMemcpyDeviceToHost = 2 };
struct cudaDeviceProp {
    char name[256];
};
cudaError_t cudaMalloc(void **, size_t);
cudaError_t cudaMallocHost(void **, size_t);
cudaError_t cudaMemcpy(void *, const void *, size_t, cudaMemcpyKind);
cudaError_t cudaMemset(void *, int, size_t);
cudaError_t cudaFree(void *);
cudaError_t cudaFreeHost(void *);
cudaError_t cudaGetLastError();
cudaError_t cudaGetDevice(int *);
cudaError_t cudaGetDeviceProperties(cudaDeviceProp *, int);
cudaError_t cudaDeviceSynchronize();
const char *cudaGetErrorString(cudaError_t);
struct cudaFuncAttributes {
    int maxThreadsPerBlock, numRegs;
};
template <class T>
cudaError_t cudaFuncGetAttributes(cudaFuncAttributes *, T *);
typedef struct CUevent_st *cudaEvent_t;
cudaError_t cudaEventCreate(cudaEvent_t *);
cudaError_t cudaEventDestroy(cudaEvent_t);
cudaError_t cudaEventRecord(cudaEvent_t, cudaStream_t = 0);
cudaError_t cudaEventSynchronize(cudaEvent_t);
cudaError_t cudaEventElapsedTime(float *, cudaEvent_t, cudaEvent_t);
enum { cudaStreamNonBlocking = 1 };
cudaError_t cudaStreamCreateWithFlags(cudaStream_t *, unsigned);
cudaError_t cudaStreamDestroy(cudaStream_t);
cudaError_t cudaStreamSynchronize(cudaStream_t);
cudaError_t cudaMemcpyAsync(
    void *, const void *, size_t, cudaMemcpyKind, cudaStream_t = 0);
__device__ int atomicOr(int *, int);
#define MUSIC_CUDA_MATH(f)                \
    __host__ __device__ double f(double); \
    __host__ __device__ float f(float);
MUSIC_CUDA_MATH(sqrt)
MUSIC_CUDA_MATH(fabs)
MUSIC_CUDA_MATH(acos)
MUSIC_CUDA_MATH(atanh)
MUSIC_CUDA_MATH(cos)
MUSIC_CUDA_MATH(tanh)
MUSIC_CUDA_MATH(acosh)
MUSIC_CUDA_MATH(asinh)
MUSIC_CUDA_MATH(cosh)
MUSIC_CUDA_MATH(sinh)
MUSIC_CUDA_MATH(sin)
MUSIC_CUDA_MATH(erf)
__host__ __device__ double atan2(double, double);
__host__ __device__ float atan2(float, float);
MUSIC_CUDA_MATH(exp)
MUSIC_CUDA_MATH(log)
__host__ __device__ double pow(double, double);
__host__ __device__ float pow(float, float);
__host__ __device__ double copysign(double, double);
__host__ __device__ float copysign(float, float);
__host__ __device__ bool isfinite(double);
__host__ __device__ bool isfinite(float);
__host__ __device__ bool isnan(double);
__host__ __device__ bool isnan(float);
#define __CUDACC__ 1
#else
#include <cuda_runtime.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "device_api.h"
#include "device_diag.h"
#include "device_strings.h"

namespace {

using music_device::EosDesc;
using music_device::StepParams;
using music_device::faceCount;
using music_device::kNComp;
using music_device::kNCompVisc;
using music_device::kNTJ;
using music_device::kNWFaceVals;

template <typename Real>
__global__ void cudaPressure(
    const StepParams<Real> P, const EosDesc<Real> D, const Real *tab,
    const Real *s, Real *p) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= P.nx * P.ny * P.neta) return;
    music_device::kernelPressure(P, D, tab, s, p, i);
}

template <typename Real, int EosType>
__global__ void cudaFace(
    const StepParams<Real> P, const EosDesc<Real> D, const Real *tab,
    const Real *s, const Real *p, int direction, Real *flux, int *flags) {
    const int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= faceCount(P.nx, P.ny, P.neta, direction)) return;
    music_device::kernelFace<Real, EosType>(
        P, D, tab, s, p, direction, f, flux, flags);
}

template <typename Real>
__global__ void cudaCell(
    const StepParams<Real> P, const EosDesc<Real> D, const Real *tab,
    const Real *sCurr, const Real *sPrev, const Real *pCurr,
    const Real *pPrev, const Real *fluxX, const Real *fluxY,
    const Real *fluxEta, const Real *src, Real *sNext, int *flags) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= P.nx * P.ny * P.neta) return;
    music_device::kernelCell(
        P, D, tab, sCurr, sPrev, pCurr, pPrev, fluxX, fluxY, fluxEta, src,
        sNext, i, flags);
}

template <typename Real>
__global__ void cudaWFace(
    const StepParams<Real> P, const Real *s, int direction, Real *out) {
    const int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= faceCount(P.nx, P.ny, P.neta, direction)) return;
    music_device::kernelWFace(P, s, direction, f, out);
}

template <typename Real>
__global__ void cudaMuBOverT(
    const StepParams<Real> P, const EosDesc<Real> D, const Real *tab,
    const Real *s, Real *a) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= P.nx * P.ny * P.neta) return;
    music_device::kernelMuBOverT(P, D, tab, s, a, i);
}

template <typename Real, int EosType>
__global__ void cudaCellVisc(
    const StepParams<Real> P, const EosDesc<Real> D, const Real *tab,
    const Real *sCurr, const Real *sPrev, const Real *pCurr,
    const Real *pPrev, const Real *aCurr, const Real *aPrev,
    const Real *fluxX, const Real *fluxY, const Real *fluxEta,
    const Real *wFaceX, const Real *wFaceY, const Real *wFaceEta,
    const Real *src, Real *sNext, int *flags) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= P.nx * P.ny * P.neta) return;
    music_device::kernelCellVisc<Real, -1, EosType>(
        P, D, tab, sCurr, sPrev, pCurr, pPrev, aCurr, aPrev, fluxX, fluxY,
        fluxEta, wFaceX, wFaceY, wFaceEta, src, sNext, i, flags);
}

template <typename Real>
__global__ void cudaStringSource(
    const StepParams<Real> P, const Real *S, const Real *sIn, Real *src,
    int *flags) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= P.nx * P.ny * P.neta) return;
    music_device::kernelStringSource(P, S, sIn, src, i, flags);
}

//! all nComp components of the listed cells, component-major
template <typename Real>
__global__ void cudaGather(
    const Real *s, size_t n, const int *cells, int nCells, int nComp,
    Real *out) {
    const size_t t = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= (size_t)nCells * nComp) return;
    const size_t c = t / nCells;
    out[t] = s[c * n + cells[t - c * nCells]];
}

template <typename Real>
__global__ void cudaConservation(
    const StepParams<Real> P, const EosDesc<Real> D, const Real *tab,
    const Real *sCurr, const Real *sPrev, Real *out) {
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= P.nx * P.neta) return;
    music_device::kernelConservationRow(P, D, tab, sCurr, sPrev, out, row);
}

template <typename Real>
__global__ void cudaCoarsening(
    int nx, int ny, int neta, Real eThreshold, const Real *s, Real *out) {
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= nx * neta) return;
    music_device::kernelCoarseningRow(nx, ny, neta, eThreshold, s, out, row);
}

template <typename Real>
__global__ void cudaEdgeMax(
    int nx, int ny, int neta, int bandX, int bandY, const Real *s,
    Real *out) {
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= nx * neta) return;
    music_device::kernelEdgeMaxRow(nx, ny, neta, bandX, bandY, s, out, row);
}

//! the kernels, for their block sizes and the profile
enum {
    kKStrings,
    kKPressure,
    kKFace,
    kKWFace,
    kKCell,
    kKCellVisc,
    kKConservation,
    kKCoarsening,
    kKEdgeMax,
    kKMuBOverT,
    kNKernels
};
const char *const kKernelNames[kNKernels] = {
    "strings", "pressure", "face", "wface", "cell", "cellVisc",
    "conservation", "coarsening", "edgeMax", "muBOverT"};

//! threads per block: MUSIC_CUDA_BLOCK (default 128; a multiple of 32),
//! capped by what each kernel's register use allows
int requestedBlock() {
    const char *s = getenv("MUSIC_CUDA_BLOCK");
    int b = (s != NULL ? atoi(s) : 128);
    if (b < 32) b = 32;
    if (b > 1024) b = 1024;
    return b - b % 32;
}

template <typename F>
int blockFor(F *kernel, int want) {
    cudaFuncAttributes a;
    if (cudaFuncGetAttributes(&a, kernel) != cudaSuccess) return 128;
    const int most = a.maxThreadsPerBlock - a.maxThreadsPerBlock % 32;
    return (want < most ? want : most);
}

bool ok(cudaError_t status, const char *what, char *err, int len) {
    if (status == cudaSuccess) return true;
    if (err != NULL && len > 0)
        snprintf(err, len, "CUDA: %s: %s", what, cudaGetErrorString(status));
    return false;
}

struct DeviceBase {
    virtual ~DeviceBase() {}
    char name[300];
    int nx = 0, ny = 0, neta = 0;
    virtual int resize(
        int nx, int ny, int neta, int nComp, char *err, int len) = 0;
    virtual int setEos(
        const void *desc, size_t descBytes, const void *table,
        size_t tableBytes, char *err, int len) = 0;
    virtual void *state(int which) = 0;
    virtual int upload(char *err, int len) = 0;
    virtual int stage(
        const void *params, size_t bytes, int rk, int *flags, char *err,
        int len) = 0;
    virtual void *source() = 0;
    virtual int uploadSource(char *err, int len) = 0;
    virtual int fetch(
        int which, int firstComp, int nComp, char *err, int len) = 0;
    virtual int download(char *err, int len) = 0;
    virtual int setStrings(const void *data, size_t bytes, char *err, int len) = 0;
    virtual int conservation(
        const void *params, size_t bytes, void *out, size_t outBytes,
        char *err, int len) = 0;
    virtual int coarsening(
        double eThreshold, void *out, size_t outBytes, char *err,
        int len) = 0;
    virtual int edgeMax(
        int bandX, int bandY, void *out, size_t outBytes, char *err,
        int len) = 0;
    virtual int finish(int *flags, char *err, int len) = 0;
    virtual int fetchCells(
        int which, const int *cells, int nCells, void **out, char *err,
        int len) = 0;
    virtual int fetchUnwritten(
        int which, int firstComp, int nComp, char *err, int len) = 0;
    virtual int conservationBegin(
        const void *params, size_t bytes, size_t outBytes, char *err,
        int len) = 0;
    virtual int conservationEnd(
        void *out, size_t outBytes, char *err, int len) = 0;
    virtual int init(char *err, int len) = 0;
};

template <typename Real>
struct Device : public DeviceBase {
    // The hotQCD viscous specialization is enabled for single precision.
    static constexpr int kHotQCDViscousEosType =
        (sizeof(Real) == sizeof(float) ? 1 : -1);

    Real *slice[3] = {NULL, NULL, NULL};
    Real *pCur = NULL, *pPred = NULL;
    // mu_B / T (baryon diffusion only, allocated on first use): aCur of
    // the first stage's input slice, aOld of the slice before it in the
    // first stage and of the second stage's input slice in the second
    Real *aCur = NULL, *aOld = NULL;
    Real *flux[3] = {NULL, NULL, NULL};
    Real *wFace[3] = {NULL, NULL, NULL};  // viscous update only
    Real *src = NULL, *srcHost = NULL;   // the source terms
    Real *table = NULL;
    Real *strings = NULL;  // the packed strings (device_strings.h)
    size_t stringsBytes = 0;
    Real *consOut = NULL;  // the conservation rows
    size_t consBytes = 0;
    Real *coarsenOut = NULL;  // the coarsening rows
    size_t coarsenBytes = 0;
    Real *edgeOut = NULL;  // the grid-edge rows
    size_t edgeBytes = 0;
    // pinned copies of the current, the previous and the predictor slice
    Real *host[3] = {NULL, NULL, NULL};
    int *flags = NULL;
    EosDesc<Real> desc;
    // next holds tau_{n-1} at the start of a step
    int cur = 0, pred = 1, next = 2;
    int nComp = kNComp;
    size_t sliceBytes = 0;
    // the second stage of a step is in flight: finish() reads its flags
    bool pending = false;
    int block[kNKernels];
    // MUSIC_CUDA_PROFILE=1: the GPU time of each kernel (each launch is
    // timed with events and waited for, so the stages do not overlap the
    // host), printed at the end of the run
    bool profile = false;
    cudaEvent_t ev[2] = {NULL, NULL};
    double profMs[kNKernels];
    long profCalls[kNKernels];

    int init(char *err, int len) {
        const int want = requestedBlock();
        block[kKStrings] = blockFor(cudaStringSource<Real>, want);
        block[kKPressure] = blockFor(cudaPressure<Real>, want);
        block[kKMuBOverT] = blockFor(cudaMuBOverT<Real>, want);
        block[kKWFace] = blockFor(cudaWFace<Real>, want);
        block[kKCell] = blockFor(cudaCell<Real>, want);
        // the generic kernels until setEos picks the specialized ones
        block[kKFace] = blockFor(cudaFace<Real, -1>, want);
        block[kKCellVisc] = blockFor(cudaCellVisc<Real, -1>, want);
        // one thread per (eta, x) row: few threads, small blocks
        block[kKConservation] = 32;
        block[kKCoarsening] = 32;
        block[kKEdgeMax] = 32;
        const char *p = getenv("MUSIC_CUDA_PROFILE");
        profile = (p != NULL && atoi(p) != 0);
        for (int k = 0; k < kNKernels; k++) {
            profMs[k] = 0.;
            profCalls[k] = 0;
        }
        if (!ok(cudaStreamCreateWithFlags(&side, cudaStreamNonBlocking),
                "cudaStreamCreateWithFlags", err, len))
            return 1;
        if (profile
            && (!ok(cudaEventCreate(&ev[0]), "cudaEventCreate", err, len)
                || !ok(cudaEventCreate(&ev[1]), "cudaEventCreate", err, len)))
            return 1;
        return 0;
    }
    int grid(int n, int k) const { return (n + block[k] - 1) / block[k]; }
    void tic() {
        if (profile) cudaEventRecord(ev[0]);
    }
    void toc(int k) {
        if (!profile) return;
        cudaEventRecord(ev[1]);
        cudaEventSynchronize(ev[1]);
        float ms = 0.f;
        cudaEventElapsedTime(&ms, ev[0], ev[1]);
        profMs[k] += ms;
        profCalls[k]++;
    }
    void printProfile() const {
        if (!profile) return;
        double total = 0.;
        for (int k = 0; k < kNKernels; k++) total += profMs[k];
        fprintf(stderr, "[CUDA profile] %s, %d-byte reals\n", name,
                (int)sizeof(Real));
        for (int k = 0; k < kNKernels; k++) {
            if (profCalls[k] == 0) continue;
            fprintf(stderr,
                    "[CUDA profile] %-12s %4d threads/block %8ld calls "
                    "%9.3f s %5.1f%%\n",
                    kKernelNames[k], block[k], profCalls[k],
                    profMs[k] * 1e-3, 100. * profMs[k] / total);
        }
        fprintf(stderr, "[CUDA profile] total        %9.3f s\n",
                total * 1e-3);
    }

    // fetchCells: the cell list on the device (re-sent when it changes),
    // the gathered values on the device and in pinned host memory
    int *dCells = NULL;
    int *hCells = NULL;  // the list dCells holds
    int nCellsCap = 0, nCellsHeld = -1;
    Real *dGather = NULL, *hGather = NULL;
    size_t gatherCap = 0;

    // the stream of a read of slice which: during a step (after stage()
    // with rk = 1, before finish()) slice 1 is tau_n, which the stages in
    // flight read and do not write, so its reads go to a stream of their
    // own (non-blocking: no implicit order with the default stream) and do
    // not wait for the step. Every other read is on the default stream,
    // behind the stages
    cudaStream_t side = NULL;
    cudaStream_t streamFor(int which) const {
        return (pending && which == 1 ? side : (cudaStream_t)0);
    }

    int fetchCells(int which, const int *cells, int nCells, void **out,
                   char *err, int len) {
        if (nCells <= 0) {
            *out = hGather;
            return 0;
        }
        const cudaStream_t st = streamFor(which);
        // the buffers grow rarely (cudaFree waits for the device)
        if (nCells > nCellsCap) {
            cudaFree(dCells);
            cudaFreeHost(hCells);
            dCells = NULL;
            hCells = NULL;
            nCellsCap = nCellsHeld = 0;
            if (!ok(cudaMalloc((void **)&dCells, sizeof(int) * nCells),
                    "cudaMalloc", err, len)
                || !ok(cudaMallocHost((void **)&hCells, sizeof(int) * nCells),
                       "cudaMallocHost", err, len))
                return 1;
            nCellsCap = nCells;
        }
        if (nCells != nCellsHeld
            || memcmp(hCells, cells, sizeof(int) * nCells) != 0) {
            memcpy(hCells, cells, sizeof(int) * nCells);
            if (!ok(cudaMemcpyAsync(dCells, hCells, sizeof(int) * nCells,
                                    cudaMemcpyHostToDevice, st),
                    "fetchCells", err, len))
                return 1;
            nCellsHeld = nCells;
        }
        const size_t nv = (size_t)nCells * nComp;
        if (nv > gatherCap) {
            cudaFree(dGather);
            cudaFreeHost(hGather);
            dGather = hGather = NULL;
            gatherCap = 0;
            if (!ok(cudaMalloc((void **)&dGather, sizeof(Real) * nv),
                    "cudaMalloc", err, len)
                || !ok(cudaMallocHost((void **)&hGather, sizeof(Real) * nv),
                       "cudaMallocHost", err, len))
                return 1;
            gatherCap = nv;
        }
        const Real *from = slice[which == 0 ? cur : (which == 1 ? next : pred)];
        const int b = 256;
        cudaGather<Real><<<(unsigned)((nv + b - 1) / b), b, 0, st>>>(
            from, (size_t)nx * ny * neta, dCells, nCells, nComp, dGather);
        if (!ok(cudaGetLastError(), "kernel launch", err, len)) return 1;
        if (!ok(cudaMemcpyAsync(hGather, dGather, sizeof(Real) * nv,
                                cudaMemcpyDeviceToHost, st),
                "fetchCells", err, len)
            || !ok(cudaStreamSynchronize(st), "fetchCells", err, len))
            return 1;
        *out = hGather;
        return 0;
    }

    // fetch() of slice 1 (tau_n) while the step may run: on the side
    // stream, without waiting for the stages
    int fetchUnwritten(int which, int firstComp, int nComp_, char *err,
                       int len) {
        if (which != 1) {
            if (err) snprintf(err, len, "CUDA: fetchUnwritten: only the "
                                        "slice of tau_n");
            return 1;
        }
        const cudaStream_t st = streamFor(1);
        const size_t n = (size_t)nx * ny * neta;
        return ok(cudaMemcpyAsync(host[1] + firstComp * n,
                                  slice[next] + firstComp * n,
                                  sizeof(Real) * nComp_ * n,
                                  cudaMemcpyDeviceToHost, st),
                  "fetchUnwritten", err, len)
                       && ok(cudaStreamSynchronize(st), "fetchUnwritten",
                             err, len)
                   ? 0
                   : 1;
    }

    // the conservation rows of the new slice, queued behind the step (on
    // the default stream) into a buffer of their own, so that a
    // conservation() call in between does not overwrite them
    Real *consAhead = NULL;
    size_t consAheadCap = 0, consAheadBytes = 0;
    int conservationBegin(const void *params, size_t bytes, size_t outBytes,
                          char *err, int len) {
        if (bytes != sizeof(StepParams<Real>)) {
            if (err) snprintf(err, len, "CUDA: StepParams size mismatch");
            return 1;
        }
        StepParams<Real> P;
        memcpy(&P, params, sizeof(P));
        if (outBytes > consAheadCap) {
            cudaFree(consAhead);
            consAhead = NULL;
            consAheadCap = 0;
            if (!ok(cudaMalloc((void **)&consAhead, outBytes), "cudaMalloc",
                    err, len))
                return 1;
            consAheadCap = outBytes;
        }
        tic();
        cudaConservation<Real>
            <<<grid(nx * neta, kKConservation), block[kKConservation]>>>(
                P, desc, table, slice[cur], slice[next], consAhead);
        if (!ok(cudaGetLastError(), "kernel launch", err, len)) return 1;
        toc(kKConservation);
        consAheadBytes = outBytes;
        return 0;
    }
    int conservationEnd(void *out, size_t outBytes, char *err, int len) {
        if (outBytes != consAheadBytes || consAhead == NULL) {
            if (err) snprintf(err, len, "CUDA: conservationEnd: not the size "
                                        "of conservationBegin");
            return 1;
        }
        // behind the kernel on the default stream
        return ok(cudaMemcpy(out, consAhead, outBytes,
                             cudaMemcpyDeviceToHost),
                  "conservationEnd", err, len)
                   ? 0
                   : 1;
    }

    // waits for the step in flight (every read of the device state is a
    // blocking copy on the default stream and waits as well)
    int finish(int *flagsOut, char *err, int len) {
        int f = 0;
        if (pending) {
            pending = false;
            if (!ok(cudaMemcpy(&f, flags, sizeof(int),
                               cudaMemcpyDeviceToHost),
                    "the step", err, len))
                return 1;
        }
        *flagsOut = f;
        return 0;
    }

    void release() {
        for (int k = 0; k < 3; k++) {
            cudaFree(slice[k]);
            cudaFree(flux[k]);
            cudaFree(wFace[k]);
            slice[k] = flux[k] = wFace[k] = NULL;
        }
        cudaFree(pCur);
        cudaFree(pPred);
        cudaFree(aCur);
        cudaFree(aOld);
        aCur = aOld = NULL;
        cudaFree(src);
        cudaFreeHost(srcHost);
        for (int k = 0; k < 3; k++) {
            cudaFreeHost(host[k]);
            host[k] = NULL;
        }
        pCur = pPred = src = srcHost = NULL;
    }
    ~Device() {
        cudaDeviceSynchronize();
        printProfile();
        for (int k = 0; k < 2; k++)
            if (ev[k] != NULL) cudaEventDestroy(ev[k]);
        cudaFree(dCells);
        cudaFreeHost(hCells);
        cudaFree(dGather);
        cudaFreeHost(hGather);
        cudaFree(consAhead);
        if (side != NULL) cudaStreamDestroy(side);
        release();
        cudaFree(table);
        cudaFree(flags);
        cudaFree(strings);
        cudaFree(consOut);
        cudaFree(coarsenOut);
        cudaFree(edgeOut);
    }

    int setStrings(const void *data, size_t bytes, char *err, int len) {
        if (bytes > stringsBytes) {
            cudaFree(strings);
            strings = NULL;
            stringsBytes = 0;
            if (!ok(cudaMalloc((void **)&strings, bytes), "cudaMalloc", err,
                    len))
                return 1;
            stringsBytes = bytes;
        }
        return ok(cudaMemcpy(strings, data, bytes, cudaMemcpyHostToDevice),
                  "setStrings", err, len)
                   ? 0
                   : 1;
    }

    int conservation(
        const void *params, size_t bytes, void *out, size_t outBytes,
        char *err, int len) {
        if (bytes != sizeof(StepParams<Real>)) {
            if (err) snprintf(err, len, "CUDA: StepParams size mismatch");
            return 1;
        }
        StepParams<Real> P;
        memcpy(&P, params, sizeof(P));
        if (outBytes > consBytes) {
            cudaFree(consOut);
            consOut = NULL;
            consBytes = 0;
            if (!ok(cudaMalloc((void **)&consOut, outBytes), "cudaMalloc",
                    err, len))
                return 1;
            consBytes = outBytes;
        }
        tic();
        cudaConservation<Real>
            <<<grid(nx * neta, kKConservation), block[kKConservation]>>>(
                P, desc, table, slice[cur], slice[next], consOut);
        if (!ok(cudaGetLastError(), "kernel launch", err, len)) return 1;
        toc(kKConservation);
        return ok(cudaMemcpy(out, consOut, outBytes, cudaMemcpyDeviceToHost),
                  "conservation", err, len)
                   ? 0
                   : 1;
    }

    int coarsening(
        double eThreshold, void *out, size_t outBytes, char *err, int len) {
        if (outBytes < sizeof(Real) * music_device::kCoarsenVals * nx * neta) {
            if (err) snprintf(err, len, "CUDA: the coarsening buffer is too "
                                        "small");
            return 1;
        }
        if (outBytes > coarsenBytes) {
            cudaFree(coarsenOut);
            coarsenOut = NULL;
            coarsenBytes = 0;
            if (!ok(cudaMalloc((void **)&coarsenOut, outBytes), "cudaMalloc",
                    err, len))
                return 1;
            coarsenBytes = outBytes;
        }
        tic();
        cudaCoarsening<Real>
            <<<grid(nx * neta, kKCoarsening), block[kKCoarsening]>>>(
                nx, ny, neta, static_cast<Real>(eThreshold), slice[cur],
                coarsenOut);
        if (!ok(cudaGetLastError(), "kernel launch", err, len)) return 1;
        toc(kKCoarsening);
        return ok(cudaMemcpy(out, coarsenOut, outBytes,
                             cudaMemcpyDeviceToHost),
                  "coarsening", err, len)
                   ? 0
                   : 1;
    }

    int edgeMax(
        int bandX, int bandY, void *out, size_t outBytes, char *err,
        int len) {
        if (outBytes < sizeof(Real) * nx * neta) {
            if (err) snprintf(err, len, "CUDA: the edge buffer is too small");
            return 1;
        }
        if (outBytes > edgeBytes) {
            cudaFree(edgeOut);
            edgeOut = NULL;
            edgeBytes = 0;
            if (!ok(cudaMalloc((void **)&edgeOut, outBytes), "cudaMalloc",
                    err, len))
                return 1;
            edgeBytes = outBytes;
        }
        tic();
        cudaEdgeMax<Real><<<grid(nx * neta, kKEdgeMax), block[kKEdgeMax]>>>(
            nx, ny, neta, bandX, bandY, slice[cur], edgeOut);
        if (!ok(cudaGetLastError(), "kernel launch", err, len)) return 1;
        toc(kKEdgeMax);
        return ok(cudaMemcpy(out, edgeOut, outBytes, cudaMemcpyDeviceToHost),
                  "edgeMax", err, len)
                   ? 0
                   : 1;
    }

    int resize(int nx_, int ny_, int neta_, int nComp_, char *err, int len) {
        release();
        nx = nx_;
        ny = ny_;
        neta = neta_;
        nComp = nComp_;
        const size_t n = (size_t)nx * ny * neta;
        sliceBytes = sizeof(Real) * nComp * n;
        for (int k = 0; k < 3; k++) {
            if (!ok(cudaMalloc((void **)&slice[k], sliceBytes), "cudaMalloc",
                    err, len)
                || !ok(cudaMemset(slice[k], 0, sliceBytes), "cudaMemset",
                       err, len)
                || !ok(cudaMalloc(
                           (void **)&flux[k],
                           sizeof(Real) * kNTJ
                               * faceCount(nx, ny, neta, k + 1)),
                       "cudaMalloc", err, len))
                return 1;
            if (nComp == kNCompVisc
                && !ok(cudaMalloc(
                           (void **)&wFace[k],
                           sizeof(Real) * kNWFaceVals
                               * faceCount(nx, ny, neta, k + 1)),
                       "cudaMalloc", err, len))
                return 1;
        }
        if (!ok(cudaMalloc((void **)&pCur, sizeof(Real) * n), "cudaMalloc",
                err, len)
            || !ok(cudaMalloc((void **)&pPred, sizeof(Real) * n),
                   "cudaMalloc", err, len)
            || !ok(cudaMalloc((void **)&src, sizeof(Real) * kNTJ * n),
                   "cudaMalloc", err, len)
            || !ok(cudaMemset(src, 0, sizeof(Real) * kNTJ * n),
                   "cudaMemset", err, len)
            || !ok(cudaMallocHost((void **)&srcHost, sizeof(Real) * kNTJ * n),
                   "cudaMallocHost", err, len))
            return 1;
        for (int k = 0; k < 3; k++) {
            if (!ok(cudaMallocHost((void **)&host[k], sliceBytes),
                    "cudaMallocHost", err, len))
                return 1;
            memset(host[k], 0, sliceBytes);
        }
        cur = 0;
        pred = 1;
        next = 2;
        return 0;
    }

    int setEos(
        const void *d, size_t descBytes, const void *tab, size_t tableBytes,
        char *err, int len) {
        if (descBytes != sizeof(desc)) {
            if (err) snprintf(err, len, "CUDA: EoS descriptor size mismatch");
            return 1;
        }
        memcpy(&desc, d, sizeof(desc));
        // Use the register and block-size limits of the selected face kernel.
        const int want = requestedBlock();
        switch (desc.type) {
            case 0:
                block[kKFace] = blockFor(cudaFace<Real, 0>, want);
                break;
            case 1:
                block[kKFace] = blockFor(cudaFace<Real, 1>, want);
                break;
            case 2:
                block[kKFace] = blockFor(cudaFace<Real, 2>, want);
                break;
            case 3:
                block[kKFace] = blockFor(cudaFace<Real, 3>, want);
                break;
            default:
                block[kKFace] = blockFor(cudaFace<Real, -1>, want);
                break;
        }
        block[kKCellVisc] =
            (desc.type == 1
                 ? blockFor(cudaCellVisc<Real, kHotQCDViscousEosType>, want)
                 : blockFor(cudaCellVisc<Real, -1>, want));
        cudaFree(table);
        table = NULL;
        const size_t bytes = (tableBytes > 0 ? tableBytes : sizeof(Real));
        if (!ok(cudaMalloc((void **)&table, bytes), "cudaMalloc", err, len))
            return 1;
        if (tableBytes > 0
            && !ok(cudaMemcpy(table, tab, tableBytes, cudaMemcpyHostToDevice),
                   "cudaMemcpy", err, len))
            return 1;
        return 0;
    }

    void *state(int which) { return host[which < 0 || which > 2 ? 0 : which]; }

    void *source() { return srcHost; }

    int uploadSource(char *err, int len) {
        const size_t bytes = sizeof(Real) * kNTJ * nx * ny * neta;
        return ok(cudaMemcpy(src, srcHost, bytes, cudaMemcpyHostToDevice),
                  "uploadSource", err, len)
                   ? 0
                   : 1;
    }

    int fetch(int which, int firstComp, int nComp_, char *err, int len) {
        const size_t n = (size_t)nx * ny * neta;
        const Real *from = slice[which == 0 ? cur : (which == 1 ? next : pred)];
        Real *to = host[which < 0 || which > 2 ? 0 : which];
        return ok(cudaMemcpy(to + firstComp * n, from + firstComp * n,
                             sizeof(Real) * nComp_ * n,
                             cudaMemcpyDeviceToHost),
                  "fetch", err, len)
                   ? 0
                   : 1;
    }

    int upload(char *err, int len) {
        if (!ok(cudaMemcpy(slice[cur], host[0], sliceBytes,
                           cudaMemcpyHostToDevice),
                "upload", err, len))
            return 1;
        if (nComp == kNCompVisc
            && !ok(cudaMemcpy(slice[next], host[1], sliceBytes,
                              cudaMemcpyHostToDevice),
                   "upload", err, len))
            return 1;
        return 0;
    }

    int download(char *err, int len) {
        return ok(cudaMemcpy(host[0], slice[cur], sliceBytes,
                             cudaMemcpyDeviceToHost),
                  "download", err, len)
                   ? 0
                   : 1;
    }

    template <int EosType>
    void launchFace(const StepParams<Real> &P, const Real *s, const Real *p,
                    int direction, int nFaces) {
        cudaFace<Real, EosType><<<grid(nFaces, kKFace), block[kKFace]>>>(
            P, desc, table, s, p, direction, flux[direction - 1], flags);
    }

    template <int EosType>
    void launchViscousCell(
        int n, const StepParams<Real> &P, const Real *sIn, const Real *sPrev,
        const Real *pIn, const Real *aIn, const Real *aPrev, Real *sOut) {
        cudaCellVisc<Real, EosType>
            <<<grid(n, kKCellVisc), block[kKCellVisc]>>>(
                P, desc, table, sIn, sPrev, pIn, pCur, aIn, aPrev, flux[0],
                flux[1], flux[2], wFace[0], wFace[1], wFace[2], src, sOut,
                flags);
    }

    int stage(const void *params, size_t bytes, int rk, int *flagsOut,
              char *err, int len) {
        if (bytes != sizeof(StepParams<Real>)) {
            if (err) snprintf(err, len, "CUDA: StepParams size mismatch");
            return 1;
        }
        StepParams<Real> P;
        memcpy(&P, params, sizeof(P));
        if (P.nComp != nComp) {
            if (err)
                snprintf(err, len, "CUDA: the slices have the wrong size");
            return 1;
        }
        const bool viscous = (nComp == kNCompVisc);
        // the flags of both stages add up; finish() reads them. The kernels
        // are queued on the default stream: everything later on it (the
        // copies of fetch, uploadSource, setStrings, download) runs after
        // them, so neither stage waits here, and the host goes on while
        // the step runs (DeviceHydro::beginStep / advanceStep)
        if (rk == 0) {
            pending = false;
            if (!ok(cudaMemset(flags, 0, sizeof(int)), "cudaMemset", err,
                    len))
                return 1;
        }
        const int n = nx * ny * neta;
        {
            P.rkFlag = rk;
            const Real *sIn = slice[rk == 0 ? cur : pred];
            const Real *sPrev = slice[rk == 0 && viscous ? next : cur];
            Real *sOut = slice[rk == 0 ? pred : next];
            Real *pIn = (rk == 0 ? pCur : pPred);
            if (P.stringSource) {
                if (strings == NULL) {
                    if (err) snprintf(err, len, "CUDA: the strings were not set");
                    return 1;
                }
                tic();
                cudaStringSource<Real>
                    <<<grid(n, kKStrings), block[kKStrings]>>>(
                        P, strings, sIn, src, flags);
                toc(kKStrings);
            }
            tic();
            cudaPressure<Real><<<grid(n, kKPressure), block[kKPressure]>>>(
                P, desc, table, sIn, pIn);
            toc(kKPressure);
            // mu_B / T of sIn and sPrev for the diffusion derivatives; in
            // the second stage sPrev is the first stage's sIn (aCur)
            const bool withMuBOverT = (viscous && P.diff == 1);
            if (withMuBOverT && aCur == NULL) {
                // both or neither: a failed second allocation must not
                // leave aCur set, or a later stage would skip this and
                // write through aOld = NULL
                Real *a = NULL, *b = NULL;
                const bool okA = ok(cudaMalloc((void **)&a, sizeof(Real) * n),
                                    "cudaMalloc", err, len);
                const bool okB =
                    okA
                    && ok(cudaMalloc((void **)&b, sizeof(Real) * n),
                          "cudaMalloc", err, len);
                if (!okB) {
                    if (okA) cudaFree(a);
                    return 1;
                }
                aCur = a;
                aOld = b;
            }
            // without diffusion the cell kernel does not read them
            const Real *aIn =
                (withMuBOverT ? (rk == 0 ? aCur : aOld) : pCur);
            const Real *aPrev =
                (withMuBOverT ? (rk == 0 ? aOld : aCur) : pCur);
            if (withMuBOverT) {
                tic();
                cudaMuBOverT<Real><<<grid(n, kKMuBOverT), block[kKMuBOverT]>>>(
                    P, desc, table, sIn, rk == 0 ? aCur : aOld);
                if (rk == 0) {
                    cudaMuBOverT<Real>
                        <<<grid(n, kKMuBOverT), block[kKMuBOverT]>>>(
                            P, desc, table, sPrev, aOld);
                }
                toc(kKMuBOverT);
            }
            for (int d = 1; d <= 3; d++) {
                const int nf = faceCount(nx, ny, neta, d);
                tic();
                // Dispatch by EoS before launch so the primitive solver
                // compiles only the selected pressure and gradient paths.
                switch (desc.type) {
                    case 0:
                        launchFace<0>(P, sIn, pIn, d, nf);
                        break;
                    case 1:
                        launchFace<1>(P, sIn, pIn, d, nf);
                        break;
                    case 2:
                        launchFace<2>(P, sIn, pIn, d, nf);
                        break;
                    case 3:
                        launchFace<3>(P, sIn, pIn, d, nf);
                        break;
                    default:
                        launchFace<-1>(P, sIn, pIn, d, nf);
                        break;
                }
                toc(kKFace);
                if (viscous) {
                    tic();
                    cudaWFace<Real><<<grid(nf, kKWFace), block[kKWFace]>>>(
                        P, sIn, d, wFace[d - 1]);
                    toc(kKWFace);
                }
            }
            tic();
            if (viscous) {
                if (desc.type == 1) {
                    launchViscousCell<kHotQCDViscousEosType>(
                        n, P, sIn, sPrev, pIn, aIn, aPrev, sOut);
                } else {
                    launchViscousCell<-1>(
                        n, P, sIn, sPrev, pIn, aIn, aPrev, sOut);
                }
                toc(kKCellVisc);
            } else {
                cudaCell<Real><<<grid(n, kKCell), block[kKCell]>>>(
                    P, desc, table, sIn, sPrev, pIn, pCur, flux[0], flux[1],
                    flux[2], src, sOut, flags);
                toc(kKCell);
            }
        }
        if (!ok(cudaGetLastError(), "kernel launch", err, len)) return 1;
        if (rk == 1) {
            pending = true;
            const int t = cur;
            cur = next;
            next = t;
        }
        *flagsOut = 0;
        return 0;
    }
};

void *cudaCreate(int realSize, char *err, int len) {
    int id = 0;
    cudaDeviceProp prop;
    if (!ok(cudaGetDevice(&id), "no CUDA device", err, len)
        || !ok(cudaGetDeviceProperties(&prop, id), "cudaGetDeviceProperties",
               err, len))
        return NULL;
    DeviceBase *d = NULL;
    if (realSize == 8) {
        d = new Device<double>();
    } else if (realSize == 4) {
        d = new Device<float>();
    } else {
        if (err) snprintf(err, len, "CUDA: realSize must be 4 or 8");
        return NULL;
    }
    snprintf(d->name, sizeof(d->name), "CUDA (%s)", prop.name);
    int *flags = NULL;
    if (!ok(cudaMalloc((void **)&flags, sizeof(int)), "cudaMalloc", err,
            len)) {
        delete d;
        return NULL;
    }
    if (realSize == 8) {
        static_cast<Device<double> *>(d)->flags = flags;
    } else {
        static_cast<Device<float> *>(d)->flags = flags;
    }
    if (d->init(err, len)) {
        delete d;
        return NULL;
    }
    return d;
}

DeviceBase *get(void *dev) { return static_cast<DeviceBase *>(dev); }

void cudaDestroy(void *dev) { delete get(dev); }
const char *cudaName(void *dev) { return get(dev)->name; }
int cudaResize(
    void *dev, int nx, int ny, int neta, int nComp, char *err, int len) {
    return get(dev)->resize(nx, ny, neta, nComp, err, len);
}
int cudaSetEos(
    void *dev, const void *desc, size_t descBytes, const void *table,
    size_t tableBytes, char *err, int len) {
    return get(dev)->setEos(desc, descBytes, table, tableBytes, err, len);
}
void *cudaState(void *dev, int which) { return get(dev)->state(which); }
int cudaUpload(void *dev, char *err, int len) {
    return get(dev)->upload(err, len);
}
int cudaStage(
    void *dev, const void *params, size_t bytes, int rk, int *flags,
    char *err, int len) {
    return get(dev)->stage(params, bytes, rk, flags, err, len);
}
void *cudaSource(void *dev) { return get(dev)->source(); }
int cudaUploadSource(void *dev, char *err, int len) {
    return get(dev)->uploadSource(err, len);
}
int cudaFetch(
    void *dev, int which, int firstComp, int nComp, char *err, int len) {
    return get(dev)->fetch(which, firstComp, nComp, err, len);
}
// fetch() of tau_n without waiting for the step in flight
int cudaFetchUnwritten(
    void *dev, int which, int firstComp, int nComp, char *err, int len) {
    return get(dev)->fetchUnwritten(which, firstComp, nComp, err, len);
}
int cudaConservationBegin(
    void *dev, const void *params, size_t bytes, size_t outBytes, char *err,
    int len) {
    return get(dev)->conservationBegin(params, bytes, outBytes, err, len);
}
int cudaConservationEnd(
    void *dev, void *out, size_t outBytes, char *err, int len) {
    return get(dev)->conservationEnd(out, outBytes, err, len);
}
int cudaDownload(void *dev, char *err, int len) {
    return get(dev)->download(err, len);
}
// the kernels read the switches at run time
int cudaSpecialize(void *, const char *, char *, int) { return 0; }
// waits for the second stage and gives the flags of the step
int cudaFinish(void *dev, int *flags, char *err, int len) {
    return get(dev)->finish(flags, err, len);
}
int cudaSetStrings(
    void *dev, const void *data, size_t bytes, char *err, int len) {
    return get(dev)->setStrings(data, bytes, err, len);
}
int cudaFetchCells(
    void *dev, int which, const int *cells, int nCells, void **out,
    char *err, int len) {
    return get(dev)->fetchCells(which, cells, nCells, out, err, len);
}
int cudaCoarseningRows(
    void *dev, double eThreshold, void *out, size_t outBytes, char *err,
    int len) {
    return get(dev)->coarsening(eThreshold, out, outBytes, err, len);
}
int cudaEdgeMaxRows(
    void *dev, int bandX, int bandY, void *out, size_t outBytes, char *err,
    int len) {
    return get(dev)->edgeMax(bandX, bandY, out, outBytes, err, len);
}
int cudaConservationRows(
    void *dev, const void *params, size_t bytes, void *out, size_t outBytes,
    char *err, int len) {
    return get(dev)->conservation(params, bytes, out, outBytes, err, len);
}

const MusicDeviceApi kCudaApi = {
    cudaCreate, cudaDestroy, cudaName,   cudaResize,       cudaSetEos,
    cudaState,  cudaUpload,  cudaStage,  cudaSource,       cudaUploadSource,
    cudaFetch,  cudaDownload, cudaSetStrings, cudaConservationRows,
    cudaSpecialize, cudaFinish, cudaFetchUnwritten,
    cudaConservationBegin, cudaConservationEnd, cudaFetchCells,
    cudaCoarseningRows, cudaEdgeMaxRows,
};

}  // namespace

extern "C" const MusicDeviceApi *music_cuda_api(void) { return &kCudaApi; }
