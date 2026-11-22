// Copyright 2026 Chun Shen
//
// gpu_backend = 2: the device kernels on an Apple GPU through Metal, in
// float (Metal has no double). The kernel source (device_defs.h,
// device_eos.h, device_ideal.h, device_visc.h and metal_kernels.metal,
// bundled by CMake
// into metal_source.inc) is compiled when the backend starts, so no offline
// Metal compiler is needed. Plain Objective-C (ARC) behind the C interface
// of device_api.h: this file is compiled by clang and linked into a MUSIC
// built by any compiler.
//
// The buffers are in shared storage: on Apple's unified memory state(0)
// and state(1) are the device buffers of the current and the previous
// slice, and upload/download copy nothing. One Heun step (two stages of
// pressure, three face (and with the viscous update three W face) and one
// cell dispatch) is one command buffer; a serial compute encoder orders
// the dispatches.

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <stdio.h>
#include <string.h>

#include "device_api.h"
#include "metal_source.inc"

// device_defs.h: kNComp, kNTJ, kNCompVisc, kNWFaceVals
enum {
    kMetalNComp = 15,
    kMetalNTJ = 7,
    kMetalNCompVisc = 31,
    kMetalNWFaceVals = 19
};

@interface MusicMetalDevice : NSObject {
  @public
    id<MTLDevice> dev;
    id<MTLCommandQueue> queue;
    id<MTLComputePipelineState> psPressure, psFace, psCell;
    id<MTLComputePipelineState> psWFace, psCellVisc;
    // next holds tau_{n-1} at the start of a step
    id<MTLBuffer> slice[3];
    id<MTLBuffer> pCur, pPred;
    id<MTLBuffer> flux[3];
    id<MTLBuffer> wFace[3];  // viscous update only
    id<MTLBuffer> eosDesc, eosTable, flags;
    int cur, pred, next;
    int nx, ny, neta, nComp;
    char name[256];
}
@end

@implementation MusicMetalDevice
@end

static void setErr(char *err, int len, const char *what, NSError *e) {
    if (err == NULL || len <= 0) return;
    if (e != nil) {
        snprintf(err, len, "Metal: %s: %s", what,
                 [[e localizedDescription] UTF8String]);
    } else {
        snprintf(err, len, "Metal: %s", what);
    }
}

static MusicMetalDevice *get(void *dev) {
    return (__bridge MusicMetalDevice *)dev;
}

static int faceCount(int nx, int ny, int neta, int d) {
    return (nx + (d == 1)) * (ny + (d == 2)) * (neta + (d == 3));
}

static void *metalCreate(int realSize, char *err, int errLen) {
    @autoreleasepool {
        if (realSize != 4) {
            setErr(err, errLen, "Metal GPUs have no double: use float", nil);
            return NULL;
        }
        MusicMetalDevice *m = [MusicMetalDevice new];
        m->dev = MTLCreateSystemDefaultDevice();
        if (m->dev == nil) {
            setErr(err, errLen, "no Metal device", nil);
            return NULL;
        }
        MTLCompileOptions *opt = [MTLCompileOptions new];
        // IEEE division and sqrt, no folding of isfinite()
        if (@available(macOS 15.0, *)) {
            opt.mathMode = MTLMathModeSafe;
        } else {
            opt.fastMathEnabled = NO;
        }
        NSError *e = nil;
        id<MTLLibrary> lib = [m->dev
            newLibraryWithSource:[NSString
                                     stringWithUTF8String:kMusicMetalSource]
                         options:opt
                           error:&e];
        if (lib == nil) {
            setErr(err, errLen, "kernel compilation failed", e);
            return NULL;
        }
        const char *names[5] = {"musicPressure", "musicFace", "musicCell",
                                "musicWFace", "musicCellVisc"};
        id<MTLComputePipelineState> ps[5];
        for (int k = 0; k < 5; k++) {
            id<MTLFunction> f = [lib
                newFunctionWithName:[NSString stringWithUTF8String:names[k]]];
            if (f == nil) {
                setErr(err, errLen, "a kernel is missing", nil);
                return NULL;
            }
            ps[k] = [m->dev newComputePipelineStateWithFunction:f error:&e];
            if (ps[k] == nil) {
                setErr(err, errLen, "pipeline creation failed", e);
                return NULL;
            }
        }
        m->psPressure = ps[0];
        m->psFace = ps[1];
        m->psCell = ps[2];
        m->psWFace = ps[3];
        m->psCellVisc = ps[4];
        m->queue = [m->dev newCommandQueue];
        m->flags = [m->dev newBufferWithLength:sizeof(int)
                                       options:MTLResourceStorageModeShared];
        m->cur = 0;
        m->pred = 1;
        m->next = 2;
        snprintf(m->name, sizeof(m->name), "Metal (%s)",
                 [[m->dev name] UTF8String]);
        return (__bridge_retained void *)m;
    }
}

static void metalDestroy(void *dev) {
    if (dev == NULL) return;
    MusicMetalDevice *m = (__bridge_transfer MusicMetalDevice *)dev;
    (void)m;  // released at the end of the scope
}

static const char *metalName(void *dev) { return get(dev)->name; }

static int metalResize(
    void *dev, int nx, int ny, int neta, int nComp, char *err, int errLen) {
    @autoreleasepool {
        MusicMetalDevice *m = get(dev);
        if (nComp != kMetalNComp && nComp != kMetalNCompVisc) {
            setErr(err, errLen, "unknown number of slice components", nil);
            return 1;
        }
        const size_t n = (size_t)nx * ny * neta;
        const MTLResourceOptions shared = MTLResourceStorageModeShared;
        for (int k = 0; k < 3; k++) {
            m->slice[k] = [m->dev
                newBufferWithLength:sizeof(float) * nComp * n
                            options:shared];
            m->flux[k] = [m->dev
                newBufferWithLength:sizeof(float) * kMetalNTJ
                                    * faceCount(nx, ny, neta, k + 1)
                            options:MTLResourceStorageModePrivate];
            m->wFace[k] = nil;
            if (nComp == kMetalNCompVisc) {
                m->wFace[k] = [m->dev
                    newBufferWithLength:sizeof(float) * kMetalNWFaceVals
                                        * faceCount(nx, ny, neta, k + 1)
                                options:MTLResourceStorageModePrivate];
                if (m->wFace[k] == nil) {
                    setErr(err, errLen, "out of device memory", nil);
                    return 1;
                }
            }
            if (m->slice[k] == nil || m->flux[k] == nil) {
                setErr(err, errLen, "out of device memory", nil);
                return 1;
            }
            memset([m->slice[k] contents], 0, sizeof(float) * nComp * n);
        }
        m->pCur = [m->dev newBufferWithLength:sizeof(float) * n options:shared];
        m->pPred = [m->dev newBufferWithLength:sizeof(float) * n
                                       options:shared];
        if (m->pCur == nil || m->pPred == nil) {
            setErr(err, errLen, "out of device memory", nil);
            return 1;
        }
        m->nx = nx;
        m->ny = ny;
        m->neta = neta;
        m->nComp = nComp;
        m->cur = 0;
        m->pred = 1;
        m->next = 2;
        return 0;
    }
}

static int metalSetEos(
    void *dev, const void *desc, size_t descBytes, const void *table,
    size_t tableBytes, char *err, int errLen) {
    @autoreleasepool {
        MusicMetalDevice *m = get(dev);
        const MTLResourceOptions shared = MTLResourceStorageModeShared;
        m->eosDesc = [m->dev newBufferWithBytes:desc
                                         length:descBytes
                                        options:shared];
        if (tableBytes > 0) {
            m->eosTable = [m->dev newBufferWithBytes:table
                                              length:tableBytes
                                             options:shared];
        } else {
            m->eosTable = [m->dev newBufferWithLength:sizeof(float)
                                              options:shared];
        }
        if (m->eosDesc == nil || m->eosTable == nil) {
            setErr(err, errLen, "cannot allocate the EoS table", nil);
            return 1;
        }
        return 0;
    }
}

static void *metalState(void *dev, int which) {
    MusicMetalDevice *m = get(dev);
    return [m->slice[which == 0 ? m->cur : m->next] contents];
}

static int metalNoCopy(void *dev, char *err, int errLen) {
    (void)dev;
    (void)err;
    (void)errLen;
    return 0;
}

static void dispatch(
    id<MTLComputeCommandEncoder> enc, id<MTLComputePipelineState> ps,
    int count) {
    [enc setComputePipelineState:ps];
    NSUInteger tg = ps.maxTotalThreadsPerThreadgroup;
    if (tg > 256) tg = 256;
    [enc dispatchThreads:MTLSizeMake(count, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
}

static int metalStep(
    void *dev, const void *params, size_t paramBytes, int *flags, char *err,
    int errLen) {
    @autoreleasepool {
        MusicMetalDevice *m = get(dev);
        if (paramBytes > 1024) {
            setErr(err, errLen, "StepParams too large", nil);
            return 1;
        }
        unsigned char P[2][1024];
        for (int rk = 0; rk < 2; rk++) {
            memcpy(P[rk], params, paramBytes);
            memcpy(P[rk] + MUSIC_DEVICE_RKFLAG_OFFSET, &rk, sizeof(int));
        }
        const bool viscous = (m->nComp == kMetalNCompVisc);
        memset([m->flags contents], 0, sizeof(int));
        const int n = m->nx * m->ny * m->neta;
        id<MTLCommandBuffer> cb = [m->queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setBuffer:m->eosDesc offset:0 atIndex:1];
        [enc setBuffer:m->eosTable offset:0 atIndex:2];
        for (int rk = 0; rk < 2; rk++) {
            id<MTLBuffer> sIn = m->slice[rk == 0 ? m->cur : m->pred];
            // the second stage averages with tau_n; the viscous update of
            // the first stage takes its time derivatives against tau_{n-1}
            id<MTLBuffer> sPrev =
                m->slice[rk == 0 && viscous ? m->next : m->cur];
            id<MTLBuffer> sOut = m->slice[rk == 0 ? m->pred : m->next];
            id<MTLBuffer> pIn = (rk == 0 ? m->pCur : m->pPred);
            [enc setBytes:P[rk] length:paramBytes atIndex:0];
            // the pressure of the stage's input slice
            [enc setBuffer:sIn offset:0 atIndex:3];
            [enc setBuffer:pIn offset:0 atIndex:4];
            dispatch(enc, m->psPressure, n);
            // the faces
            [enc setBuffer:m->flags offset:0 atIndex:6];
            for (int d = 1; d <= 3; d++) {
                const int nf = faceCount(m->nx, m->ny, m->neta, d);
                [enc setBuffer:m->flux[d - 1] offset:0 atIndex:5];
                [enc setBytes:&d length:sizeof(int) atIndex:7];
                dispatch(enc, m->psFace, nf);
                if (viscous) {
                    [enc setBuffer:m->wFace[d - 1] offset:0 atIndex:5];
                    dispatch(enc, m->psWFace, nf);
                }
            }
            // the cells
            [enc setBuffer:sPrev offset:0 atIndex:4];
            [enc setBuffer:pIn offset:0 atIndex:5];
            [enc setBuffer:m->pCur offset:0 atIndex:6];
            for (int d = 0; d < 3; d++)
                [enc setBuffer:m->flux[d] offset:0 atIndex:7 + d];
            [enc setBuffer:sOut offset:0 atIndex:10];
            [enc setBuffer:m->flags offset:0 atIndex:11];
            if (viscous) {
                for (int d = 0; d < 3; d++)
                    [enc setBuffer:m->wFace[d] offset:0 atIndex:12 + d];
                dispatch(enc, m->psCellVisc, n);
            } else {
                dispatch(enc, m->psCell, n);
            }
        }
        [enc endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
        if (cb.error != nil) {
            setErr(err, errLen, "the step failed", cb.error);
            return 1;
        }
        const int t = m->cur;
        m->cur = m->next;
        m->next = t;
        *flags = *(int *)[m->flags contents];
        return 0;
    }
}

static const MusicDeviceApi kMetalApi = {
    metalCreate, metalDestroy, metalName,   metalResize, metalSetEos,
    metalState,  metalNoCopy,  metalStep,   metalNoCopy,
};

const MusicDeviceApi *music_metal_api(void) { return &kMetalApi; }
