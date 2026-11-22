// Copyright 2026 Chun Shen
//
// gpu_backend = 3: the device kernels on an NVIDIA GPU through CUDA, in
// double (default) or float (gpu_float = 1), behind the C interface of
// device_api.h (built by nvcc with CMake option cuda=ON).
//
// The state stays in device memory; state() is a pinned host copy of the
// current slice, which upload() sends to the device and download() fills
// after a step. One Heun step is ten kernel launches on the default stream
// (pressure, three face directions and the cells, for both stages).
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
__device__ int atomicOr(int *, int);
#define MUSIC_CUDA_MATH(f)                \
    __host__ __device__ double f(double); \
    __host__ __device__ float f(float);
MUSIC_CUDA_MATH(sqrt)
MUSIC_CUDA_MATH(fabs)
__host__ __device__ bool isfinite(double);
__host__ __device__ bool isfinite(float);
#define __CUDACC__ 1
#else
#include <cuda_runtime.h>
#endif

#include <stdio.h>
#include <string.h>

#include "device_api.h"
#include "device_ideal.h"

namespace {

using music_device::EosDesc;
using music_device::StepParams;
using music_device::faceCount;
using music_device::kNComp;
using music_device::kNTJ;

template <typename Real>
__global__ void cudaPressure(
    const StepParams<Real> P, const EosDesc<Real> D, const Real *tab,
    const Real *s, Real *p) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= P.nx * P.ny * P.neta) return;
    music_device::kernelPressure(P, D, tab, s, p, i);
}

template <typename Real>
__global__ void cudaFace(
    const StepParams<Real> P, const EosDesc<Real> D, const Real *tab,
    const Real *s, const Real *p, int direction, Real *flux, int *flags) {
    const int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= faceCount(P.nx, P.ny, P.neta, direction)) return;
    music_device::kernelFace(P, D, tab, s, p, direction, f, flux, flags);
}

template <typename Real>
__global__ void cudaCell(
    const StepParams<Real> P, const EosDesc<Real> D, const Real *tab,
    const Real *sCurr, const Real *sPrev, const Real *pCurr,
    const Real *pPrev, const Real *fluxX, const Real *fluxY,
    const Real *fluxEta, Real *sNext, int *flags) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= P.nx * P.ny * P.neta) return;
    music_device::kernelCell(
        P, D, tab, sCurr, sPrev, pCurr, pPrev, fluxX, fluxY, fluxEta, sNext,
        i, flags);
}

const int kBlock = 128;

int blocks(int n) { return (n + kBlock - 1) / kBlock; }

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
    virtual int resize(int nx, int ny, int neta, char *err, int len) = 0;
    virtual int setEos(
        const void *desc, size_t descBytes, const void *table,
        size_t tableBytes, char *err, int len) = 0;
    virtual void *state() = 0;
    virtual int upload(char *err, int len) = 0;
    virtual int step(
        const void *params, size_t bytes, int *flags, char *err, int len) = 0;
    virtual int download(char *err, int len) = 0;
};

template <typename Real>
struct Device : public DeviceBase {
    Real *slice[3] = {NULL, NULL, NULL};
    Real *pCur = NULL, *pPred = NULL;
    Real *flux[3] = {NULL, NULL, NULL};
    Real *table = NULL;
    Real *host = NULL;  // pinned copy of the current slice
    int *flags = NULL;
    EosDesc<Real> desc;
    int cur = 0, pred = 1, next = 2;
    size_t sliceBytes = 0;

    void release() {
        for (int k = 0; k < 3; k++) {
            cudaFree(slice[k]);
            cudaFree(flux[k]);
            slice[k] = flux[k] = NULL;
        }
        cudaFree(pCur);
        cudaFree(pPred);
        cudaFreeHost(host);
        pCur = pPred = host = NULL;
    }
    ~Device() {
        release();
        cudaFree(table);
        cudaFree(flags);
    }

    int resize(int nx_, int ny_, int neta_, char *err, int len) {
        release();
        nx = nx_;
        ny = ny_;
        neta = neta_;
        const size_t n = (size_t)nx * ny * neta;
        sliceBytes = sizeof(Real) * kNComp * n;
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
        }
        if (!ok(cudaMalloc((void **)&pCur, sizeof(Real) * n), "cudaMalloc",
                err, len)
            || !ok(cudaMalloc((void **)&pPred, sizeof(Real) * n),
                   "cudaMalloc", err, len)
            || !ok(cudaMallocHost((void **)&host, sliceBytes),
                   "cudaMallocHost", err, len))
            return 1;
        memset(host, 0, sliceBytes);
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

    void *state() { return host; }

    int upload(char *err, int len) {
        return ok(cudaMemcpy(slice[cur], host, sliceBytes,
                             cudaMemcpyHostToDevice),
                  "upload", err, len)
                   ? 0
                   : 1;
    }

    int download(char *err, int len) {
        return ok(cudaMemcpy(host, slice[cur], sliceBytes,
                             cudaMemcpyDeviceToHost),
                  "download", err, len)
                   ? 0
                   : 1;
    }

    int step(const void *params, size_t bytes, int *flagsOut, char *err,
             int len) {
        if (bytes != sizeof(StepParams<Real>)) {
            if (err) snprintf(err, len, "CUDA: StepParams size mismatch");
            return 1;
        }
        StepParams<Real> P;
        memcpy(&P, params, sizeof(P));
        if (!ok(cudaMemset(flags, 0, sizeof(int)), "cudaMemset", err, len))
            return 1;
        const int n = nx * ny * neta;
        for (int rk = 0; rk < 2; rk++) {
            P.rkFlag = rk;
            const Real *sIn = slice[rk == 0 ? cur : pred];
            Real *sOut = slice[rk == 0 ? pred : next];
            Real *pIn = (rk == 0 ? pCur : pPred);
            cudaPressure<Real><<<blocks(n), kBlock>>>(P, desc, table, sIn, pIn);
            for (int d = 1; d <= 3; d++) {
                cudaFace<Real><<<blocks(faceCount(nx, ny, neta, d)), kBlock>>>(
                    P, desc, table, sIn, pIn, d, flux[d - 1], flags);
            }
            cudaCell<Real><<<blocks(n), kBlock>>>(
                P, desc, table, sIn, slice[cur], pIn, pCur, flux[0], flux[1],
                flux[2], sOut, flags);
        }
        if (!ok(cudaGetLastError(), "kernel launch", err, len)) return 1;
        int f = 0;
        if (!ok(cudaMemcpy(&f, flags, sizeof(int), cudaMemcpyDeviceToHost),
                "the step", err, len))
            return 1;
        const int t = cur;
        cur = next;
        next = t;
        *flagsOut = f;
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
    return d;
}

DeviceBase *get(void *dev) { return static_cast<DeviceBase *>(dev); }

void cudaDestroy(void *dev) { delete get(dev); }
const char *cudaName(void *dev) { return get(dev)->name; }
int cudaResize(void *dev, int nx, int ny, int neta, char *err, int len) {
    return get(dev)->resize(nx, ny, neta, err, len);
}
int cudaSetEos(
    void *dev, const void *desc, size_t descBytes, const void *table,
    size_t tableBytes, char *err, int len) {
    return get(dev)->setEos(desc, descBytes, table, tableBytes, err, len);
}
void *cudaState(void *dev) { return get(dev)->state(); }
int cudaUpload(void *dev, char *err, int len) {
    return get(dev)->upload(err, len);
}
int cudaStep(
    void *dev, const void *params, size_t bytes, int *flags, char *err,
    int len) {
    return get(dev)->step(params, bytes, flags, err, len);
}
int cudaDownload(void *dev, char *err, int len) {
    return get(dev)->download(err, len);
}

const MusicDeviceApi kCudaApi = {
    cudaCreate, cudaDestroy, cudaName, cudaResize, cudaSetEos,
    cudaState,  cudaUpload,  cudaStep, cudaDownload,
};

}  // namespace

extern "C" const MusicDeviceApi *music_cuda_api(void) { return &kCudaApi; }
