// Metal backend of the (2+1)D ideal KT prototype: the kernels of kt2d.h are
// compiled at run time (newLibraryWithSource, no offline metal compiler
// needed) in float, the only floating type Metal has.
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

#include "kt2d.h"
#include "kt2d_metal.h"

static const char *kEntryPoints = R"MSL(
using namespace kt;
kernel void faceX(constant Params<float> &P [[buffer(0)]],
                  device const float *table [[buffer(1)]],
                  device const float *s [[buffer(2)]],
                  device float *flux [[buffer(3)]],
                  device atomic_int *flags [[buffer(4)]],
                  uint gid [[thread_position_in_grid]]) {
    const int nf = (P.nx + 1) * P.ny;
    if (int(gid) >= nf) return;
    faceFlux(P, table, s, 1, int(gid) / P.ny, int(gid) % P.ny, flux, flags);
}
kernel void faceY(constant Params<float> &P [[buffer(0)]],
                  device const float *table [[buffer(1)]],
                  device const float *s [[buffer(2)]],
                  device float *flux [[buffer(3)]],
                  device atomic_int *flags [[buffer(4)]],
                  uint gid [[thread_position_in_grid]]) {
    const int nf = P.nx * (P.ny + 1);
    if (int(gid) >= nf) return;
    faceFlux(P, table, s, 2, int(gid) / (P.ny + 1), int(gid) % (P.ny + 1),
             flux, flags);
}
kernel void cells(constant Params<float> &P [[buffer(0)]],
                  device const float *table [[buffer(1)]],
                  device const float *sCurr [[buffer(2)]],
                  device const float *sPrev [[buffer(3)]],
                  device const float *fluxX [[buffer(4)]],
                  device const float *fluxY [[buffer(5)]],
                  device float *sNext [[buffer(6)]],
                  device atomic_int *flags [[buffer(7)]],
                  uint gid [[thread_position_in_grid]]) {
    if (int(gid) >= P.nx * P.ny) return;
    cellUpdate(P, table, sCurr, sPrev, fluxX, fluxY, sNext,
               int(gid) / P.ny, int(gid) % P.ny, flags);
}
)MSL";

int kt2d_metal_run(const kt::Params<float> *p0, const float *table,
                   int tableLen, float *state, int nSteps, int stepsPerCommit,
                   KtMetalTiming *timing, char *errMsg, int errLen) {
    @autoreleasepool {
        auto fail = [&](const std::string &m) {
            std::snprintf(errMsg, errLen, "%s", m.c_str());
            return -1;
        };
        id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
        if (!dev) return fail("no Metal device");
        std::ifstream hdr(KT2D_HEADER_PATH);
        if (!hdr) return fail(std::string("cannot read ") + KT2D_HEADER_PATH);
        std::stringstream ss;
        ss << hdr.rdbuf() << kEntryPoints;
        const std::string src = ss.str();

        auto t0 = std::chrono::steady_clock::now();
        MTLCompileOptions *opt = [MTLCompileOptions new];
        // IEEE division/sqrt and no isfinite() folding
        if (@available(macOS 15.0, *)) {
            opt.mathMode = MTLMathModeSafe;
        } else {
            opt.fastMathEnabled = NO;
        }
        NSError *err = nil;
        id<MTLLibrary> lib = [dev
            newLibraryWithSource:[NSString stringWithUTF8String:src.c_str()]
                         options:opt
                           error:&err];
        if (!lib) return fail([[err localizedDescription] UTF8String]);
        id<MTLComputePipelineState> psX, psY, psC;
        for (int k = 0; k < 3; k++) {
            const char *name = (k == 0 ? "faceX" : (k == 1 ? "faceY" : "cells"));
            id<MTLFunction> f =
                [lib newFunctionWithName:[NSString stringWithUTF8String:name]];
            id<MTLComputePipelineState> ps =
                [dev newComputePipelineStateWithFunction:f error:&err];
            if (!ps) return fail([[err localizedDescription] UTF8String]);
            (k == 0 ? psX : (k == 1 ? psY : psC)) = ps;
        }
        auto t1 = std::chrono::steady_clock::now();

        const int nx = p0->nx, ny = p0->ny, n = nx * ny;
        const size_t stBytes = sizeof(float) * 9 * n;
        const MTLResourceOptions shared = MTLResourceStorageModeShared;
        // unified memory: the state buffers are shared with the CPU, so the
        // host reads them (freeze-out, output) without a copy
        id<MTLBuffer> bS0 = [dev newBufferWithBytes:state length:stBytes options:shared];
        id<MTLBuffer> bS1 = [dev newBufferWithLength:stBytes options:shared];
        id<MTLBuffer> bS2 = [dev newBufferWithLength:stBytes options:shared];
        id<MTLBuffer> bFx = [dev newBufferWithLength:sizeof(float) * 4 * (nx + 1) * ny
                                             options:MTLResourceStorageModePrivate];
        id<MTLBuffer> bFy = [dev newBufferWithLength:sizeof(float) * 4 * nx * (ny + 1)
                                             options:MTLResourceStorageModePrivate];
        id<MTLBuffer> bTab = [dev newBufferWithBytes:(tableLen > 0 ? table : state)
                                              length:sizeof(float) * (tableLen > 0 ? tableLen : 1)
                                             options:shared];
        id<MTLBuffer> bFlag = [dev newBufferWithLength:sizeof(int) options:shared];
        std::memset([bFlag contents], 0, sizeof(int));
        id<MTLCommandQueue> q = [dev newCommandQueue];

        auto dispatch = [&](id<MTLComputeCommandEncoder> enc,
                            id<MTLComputePipelineState> ps, int count) {
            [enc setComputePipelineState:ps];
            NSUInteger tg = std::min<NSUInteger>(256, ps.maxTotalThreadsPerThreadgroup);
            [enc dispatchThreads:MTLSizeMake(count, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
        };

        // s0: tau_n, s1: predictor, s2: tau_{n+1}; rotate s0 <-> s2
        id<MTLBuffer> cur = bS0, pred = bS1, next = bS2;
        kt::Params<float> P = *p0;
        int step = 0;
        while (step < nSteps) {
            id<MTLCommandBuffer> cb = [q commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
            [enc setBuffer:bTab offset:0 atIndex:1];
            const int last = std::min(nSteps, step + stepsPerCommit);
            for (; step < last; step++) {
                P.tau = p0->tau + step * p0->dtau;
                for (int rk = 0; rk < 2; rk++) {
                    P.rkFlag = rk;
                    id<MTLBuffer> sIn = (rk == 0 ? cur : pred);
                    id<MTLBuffer> sOut = (rk == 0 ? pred : next);
                    // the faces
                    [enc setBytes:&P length:sizeof(P) atIndex:0];
                    [enc setBuffer:sIn offset:0 atIndex:2];
                    [enc setBuffer:bFx offset:0 atIndex:3];
                    [enc setBuffer:bFlag offset:0 atIndex:4];
                    dispatch(enc, psX, (nx + 1) * ny);
                    [enc setBuffer:bFy offset:0 atIndex:3];
                    dispatch(enc, psY, nx * (ny + 1));
                    // the cells
                    [enc setBuffer:sIn offset:0 atIndex:2];
                    [enc setBuffer:cur offset:0 atIndex:3];
                    [enc setBuffer:bFx offset:0 atIndex:4];
                    [enc setBuffer:bFy offset:0 atIndex:5];
                    [enc setBuffer:sOut offset:0 atIndex:6];
                    [enc setBuffer:bFlag offset:0 atIndex:7];
                    dispatch(enc, psC, n);
                }
                std::swap(cur, next);
            }
            [enc endEncoding];
            [cb commit];
            [cb waitUntilCompleted];
            if (cb.error) return fail([[cb.error localizedDescription] UTF8String]);
            if (timing) timing->gpuSeconds += cb.GPUEndTime - cb.GPUStartTime;
        }
        auto t2 = std::chrono::steady_clock::now();
        std::memcpy(state, [cur contents], stBytes);
        if (timing) {
            timing->compileSeconds = std::chrono::duration<double>(t1 - t0).count();
            timing->runSeconds = std::chrono::duration<double>(t2 - t1).count();
            std::snprintf(timing->device, sizeof(timing->device), "%s",
                          [[dev name] UTF8String]);
        }
        return *static_cast<int *>([bFlag contents]);
    }
}
