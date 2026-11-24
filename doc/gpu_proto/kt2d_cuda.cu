// CUDA backend of the (2+1)D ideal KT prototype, the same kernels as the
// CPU and Metal backends (kt2d.h), here in double or float.
//
//   nvcc -O3 -std=c++17 -arch=native kt2d_cuda.cu -o kt2d_cuda
//   ./kt2d_cuda [n] [double|float] [mode] [steps]
//
// Without an NVIDIA toolchain (macOS) the source is only type-checked by
// clang's CUDA front end: make cuda-syntax (KT2D_NO_CUDA_SDK stubs the
// runtime API).
#ifdef KT2D_NO_CUDA_SDK
#define __host__ __attribute__((host))
#define __device__ __attribute__((device))
#define __global__ __attribute__((global))
#define __shared__ __attribute__((shared))
#define __constant__ __attribute__((constant))
#include <__clang_cuda_builtin_vars.h>
#include <cstddef>
struct dim3 {
    unsigned x, y, z;
    __host__ __device__ dim3(unsigned a = 1, unsigned b = 1, unsigned c = 1)
        : x(a), y(b), z(c) {}
};
typedef struct CUstream_st *cudaStream_t;
extern "C" int cudaConfigureCall(dim3, dim3, size_t = 0, cudaStream_t = 0);
extern "C" unsigned __cudaPushCallConfiguration(dim3, dim3, size_t = 0,
                                                void * = 0);
typedef int cudaError_t;
typedef struct CUevent_st *cudaEvent_t;
enum { cudaMemcpyHostToDevice = 1, cudaMemcpyDeviceToHost = 2 };
cudaError_t cudaMalloc(void **, size_t);
cudaError_t cudaMemcpy(void *, const void *, size_t, int);
cudaError_t cudaMemset(void *, int, size_t);
cudaError_t cudaFree(void *);
cudaError_t cudaEventCreate(cudaEvent_t *);
cudaError_t cudaEventRecord(cudaEvent_t, void * = 0);
cudaError_t cudaEventSynchronize(cudaEvent_t);
cudaError_t cudaEventElapsedTime(float *, cudaEvent_t, cudaEvent_t);
cudaError_t cudaGetLastError();
const char *cudaGetErrorString(cudaError_t);
__device__ int atomicOr(int *, int);
// the overloads of CUDA's math headers
#define KT_MATH(f)                 \
    __device__ double f(double);   \
    __device__ float f(float);
KT_MATH(sqrt)
KT_MATH(fabs)
KT_MATH(log)
__device__ double fmax(double, double);
__device__ float fmax(float, float);
__device__ double fmin(double, double);
__device__ float fmin(float, float);
__device__ bool isfinite(double);
__device__ bool isfinite(float);
#define __CUDACC__ 1
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "kt2d.h"

template <typename Real>
__global__ void faceX(const kt::Params<Real> P, const Real *table,
                      const Real *s, Real *flux, int *flags) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= (P.nx + 1) * P.ny) return;
    kt::faceFlux(P, table, s, 1, i / P.ny, i % P.ny, flux, flags);
}

template <typename Real>
__global__ void faceY(const kt::Params<Real> P, const Real *table,
                      const Real *s, Real *flux, int *flags) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= P.nx * (P.ny + 1)) return;
    kt::faceFlux(P, table, s, 2, i / (P.ny + 1), i % (P.ny + 1), flux, flags);
}

template <typename Real>
__global__ void cells(const kt::Params<Real> P, const Real *table,
                      const Real *sCurr, const Real *sPrev, const Real *fluxX,
                      const Real *fluxY, Real *sNext, int *flags) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= P.nx * P.ny) return;
    kt::cellUpdate(P, table, sCurr, sPrev, fluxX, fluxY, sNext, i / P.ny,
                   i % P.ny, flags);
}

#define CK(x)                                                              \
    do {                                                                   \
        cudaError_t e_ = (x);                                              \
        if (e_ != 0) {                                                     \
            std::fprintf(stderr, "%s: %s\n", #x, cudaGetErrorString(e_));  \
            std::exit(1);                                                  \
        }                                                                  \
    } while (0)

template <typename Real>
int run(int n, int mode, int nSteps) {
    const double size = (n - 1) * 0.05, q = 10. / size, dtau = 0.005;
    kt::Params<Real> P{};
    P.nx = P.ny = n;
    P.mode = mode;
    P.dx = P.dy = Real(0.05);
    P.dtau = Real(dtau);
    P.theta = Real(1.8);
    P.absErr = Real(sizeof(Real) == 8 ? 1e-16 : 1e-8);
    P.relErr = Real(sizeof(Real) == 8 ? 1e-14 : 1e-6);
    const int N = n * n;
    // the Gubser initial state (e0 = 1, n_B0 = 0.5) at tau = 1
    std::vector<Real> h(9 * N);  // storeQ = 0: components 5-8 unused
    for (int ix = 0; ix < n; ix++)
        for (int iy = 0; iy < n; iy++) {
            const double x = -size / 2 + ix * 0.05, y = -size / 2 + iy * 0.05;
            const double r = std::sqrt(x * x + y * y), tau = 1.;
            const double den = 1. + 2. * q * q * (tau * tau + r * r)
                               + std::pow(q, 4) * std::pow(tau * tau - r * r, 2);
            const double ur = std::sinh(std::atanh(
                2. * q * q * tau * r / (1. + q * q * tau * tau + q * q * r * r)));
            const int i = ix * n + iy;
            const double ux = (r > 0 ? ur * x / r : 0.), uy = (r > 0 ? ur * y / r : 0.);
            h[i] = Real(std::pow(2. * q, 8. / 3.) / std::pow(den, 4. / 3.));
            h[N + i] = Real(2. * q * q / den);
            h[2 * N + i] = Real(std::sqrt(1. + ux * ux + uy * uy));
            h[3 * N + i] = Real(ux);
            h[4 * N + i] = Real(uy);
        }
    Real *s[3], *fx, *fy;
    int *flags;
    for (auto &p : s) CK(cudaMalloc((void **)&p, sizeof(Real) * 9 * N));
    CK(cudaMalloc((void **)&fx, sizeof(Real) * 4 * (n + 1) * n));
    CK(cudaMalloc((void **)&fy, sizeof(Real) * 4 * (n + 1) * n));
    CK(cudaMalloc((void **)&flags, sizeof(int)));
    CK(cudaMemset(flags, 0, sizeof(int)));
    CK(cudaMemcpy(s[0], h.data(), sizeof(Real) * 9 * N, cudaMemcpyHostToDevice));
    const int B = 256;
    cudaEvent_t t0, t1;
    CK(cudaEventCreate(&t0));
    CK(cudaEventCreate(&t1));
    CK(cudaEventRecord(t0));
    Real *cur = s[0], *pred = s[1], *next = s[2];
    for (int step = 0; step < nSteps; step++) {
        P.tau = Real(1. + step * dtau);
        for (int rk = 0; rk < 2; rk++) {
            P.rkFlag = rk;
            Real *sIn = (rk == 0 ? cur : pred), *sOut = (rk == 0 ? pred : next);
            faceX<Real><<<((n + 1) * n + B - 1) / B, B>>>(P, nullptr, sIn, fx, flags);
            faceY<Real><<<((n + 1) * n + B - 1) / B, B>>>(P, nullptr, sIn, fy, flags);
            cells<Real><<<(N + B - 1) / B, B>>>(P, nullptr, sIn, cur, fx, fy, sOut, flags);
        }
        Real *t = cur;
        cur = next;
        next = t;
    }
    CK(cudaEventRecord(t1));
    CK(cudaEventSynchronize(t1));
    CK(cudaGetLastError());
    float ms = 0.f;
    CK(cudaEventElapsedTime(&ms, t0, t1));
    int hflags = 0;
    CK(cudaMemcpy(h.data(), cur, sizeof(Real) * 9 * N, cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(&hflags, flags, sizeof(int), cudaMemcpyDeviceToHost));
    std::printf("cuda %s: %d x %d, mode %d, %d steps: %.4f s, %.2f ns per cell "
                "and RK stage, flags %d, e(0,0) = %.9e\n",
                sizeof(Real) == 8 ? "double" : "float", n, n, mode, nSteps,
                1e-3 * ms, 1e6 * ms / (2. * nSteps * N), hflags,
                double(h[(n / 2) * n + n / 2]));
    for (auto &p : s) cudaFree(p);
    cudaFree(fx);
    cudaFree(fy);
    cudaFree(flags);
    return hflags;
}

int main(int argc, char **argv) {
    const int n = (argc > 1 ? std::atoi(argv[1]) : 1001);
    const bool dbl = (argc > 2 ? std::strcmp(argv[2], "float") != 0 : true);
    const int mode = (argc > 3 ? std::atoi(argv[3]) : 1);
    const int steps = (argc > 4 ? std::atoi(argv[4]) : 50);
    return dbl ? run<double>(n, mode, steps) : run<float>(n, mode, steps);
}
