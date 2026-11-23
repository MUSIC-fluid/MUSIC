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
// slice, and upload/download copy nothing. One RK stage (pressure, three
// face (and with the viscous update three W face) and one cell dispatch)
// is one command buffer; a serial compute encoder orders the dispatches.

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <stdio.h>
#include <stdlib.h>
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
    // psCellVisc[rk]: the viscous cell update of each RK stage
    id<MTLComputePipelineState> psWFace, psCellVisc[2];
    id<MTLComputePipelineState> psStrings, psConservation;
    // next holds tau_{n-1} at the start of a step
    id<MTLBuffer> slice[3];
    id<MTLBuffer> pCur, pPred;
    id<MTLBuffer> flux[3];
    id<MTLBuffer> wFace[3];  // viscous update only
    id<MTLBuffer> src;       // the source terms, kNTJ components
    id<MTLBuffer> strings;   // the packed strings (device_strings.h)
    id<MTLBuffer> consOut;   // the conservation rows
    id<MTLCommandBuffer> pending;  // the first stage, in flight
    id<MTLBuffer> eosDesc, eosTable, flags;
    NSString *specKey;  // the switches the kernels were compiled for
    // MUSIC_METAL_PROFILE set: the GPU time of each kernel (s)
    BOOL profile;
    NSMutableDictionary *profileSeconds;
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

//! compiles the kernels; macros (may be nil): the switches of the run
//! (MUSIC_SPEC, MUSIC_SPEC_<name>, see MD_RUN in device_defs.h)
static int buildPipelines(
    MusicMetalDevice *m, NSDictionary *macros, char *err, int errLen) {
    MTLCompileOptions *opt = [MTLCompileOptions new];
    // IEEE division and sqrt, no folding of isfinite()
    if (@available(macOS 15.0, *)) {
        opt.mathMode = MTLMathModeSafe;
    } else {
        opt.fastMathEnabled = NO;
    }
    if (macros != nil) opt.preprocessorMacros = macros;
    NSError *e = nil;
    id<MTLLibrary> lib = [m->dev
        newLibraryWithSource:[NSString stringWithUTF8String:kMusicMetalSource]
                     options:opt
                       error:&e];
    if (lib == nil) {
        setErr(err, errLen, "kernel compilation failed", e);
        return 1;
    }
    const char *names[8] = {"musicPressure",     "musicFace",
                            "musicCell",         "musicWFace",
                            "musicCellVisc0",    "musicStringSource",
                            "musicConservation", "musicCellVisc1"};
    id<MTLComputePipelineState> ps[8];
    for (int k = 0; k < 8; k++) {
        id<MTLFunction> f =
            [lib newFunctionWithName:[NSString stringWithUTF8String:names[k]]];
        if (f == nil) {
            setErr(err, errLen, "a kernel is missing", nil);
            return 1;
        }
        ps[k] = [m->dev newComputePipelineStateWithFunction:f error:&e];
        if (ps[k] == nil) {
            setErr(err, errLen, "pipeline creation failed", e);
            return 1;
        }
    }
    m->psPressure = ps[0];
    m->psFace = ps[1];
    m->psCell = ps[2];
    m->psWFace = ps[3];
    m->psCellVisc[0] = ps[4];
    m->psCellVisc[1] = ps[7];
    m->psStrings = ps[5];
    m->psConservation = ps[6];
    return 0;
}

//! the generic kernels, if no specialize() came first
static int ensurePipelines(MusicMetalDevice *m, char *err, int errLen) {
    if (m->psCell != nil) return 0;
    return buildPipelines(m, nil, err, errLen);
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
        m->queue = [m->dev newCommandQueue];
        m->profile = (getenv("MUSIC_METAL_PROFILE") != NULL);
        m->profileSeconds = [NSMutableDictionary new];
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
    if (m->profile) {
        double total = 0.;
        for (NSString *k in m->profileSeconds)
            total += [m->profileSeconds[k] doubleValue];
        fprintf(stderr,
                "Metal profile (GPU time per kernel); max threads per "
                "threadgroup: pressure %lu, face %lu, wface %lu, cell %lu, "
                "cellVisc %lu / %lu\n",
                (unsigned long)m->psPressure.maxTotalThreadsPerThreadgroup,
                (unsigned long)m->psFace.maxTotalThreadsPerThreadgroup,
                (unsigned long)m->psWFace.maxTotalThreadsPerThreadgroup,
                (unsigned long)m->psCell.maxTotalThreadsPerThreadgroup,
                (unsigned long)m->psCellVisc[0].maxTotalThreadsPerThreadgroup,
                (unsigned long)m->psCellVisc[1].maxTotalThreadsPerThreadgroup);
        for (NSString *k in [m->profileSeconds keysSortedByValueUsingSelector:
                                                   @selector(compare:)]) {
            const double t = [m->profileSeconds[k] doubleValue];
            fprintf(stderr, "  %-13s %8.3f s  %5.1f%%\n", [k UTF8String], t,
                    100. * t / (total > 0. ? total : 1.));
        }
    }
    // released at the end of the scope
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
        m->pending = nil;
        m->src = [m->dev newBufferWithLength:sizeof(float) * kMetalNTJ * n
                                     options:shared];
        if (m->src == nil) {
            setErr(err, errLen, "out of device memory", nil);
            return 1;
        }
        memset([m->src contents], 0, sizeof(float) * kMetalNTJ * n);
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
    const int k = (which == 0 ? m->cur : (which == 1 ? m->next : m->pred));
    return [m->slice[k] contents];
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
    // the largest threadgroup the pipeline allows: on the M4 1024 runs the
    // viscous update 11% faster than 256 (tiles in x, y, eta instead of
    // whole x rows were slower). MUSIC_METAL_TG overrides it (tuning)
    static int tgEnv = -1;
    if (tgEnv < 0) {
        const char *v = getenv("MUSIC_METAL_TG");
        tgEnv = (v != NULL ? atoi(v) : 0);
    }
    if (tgEnv > 0 && (NSUInteger)tgEnv < tg) tg = (NSUInteger)tgEnv;
    [enc dispatchThreads:MTLSizeMake(count, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
}

//! waits for the stage in flight (if any)
static int waitPending(MusicMetalDevice *m, char *err, int errLen) {
    if (m->pending == nil) return 0;
    id<MTLCommandBuffer> cb = m->pending;
    m->pending = nil;
    [cb waitUntilCompleted];
    if (cb.error != nil) {
        setErr(err, errLen, "the step failed", cb.error);
        return 1;
    }
    return 0;
}

//! one RK stage as one command buffer. The first stage is committed
//! without waiting (the queue runs the second after it); fetch() and the
//! second stage wait.
static int metalStage(
    void *dev, const void *params, size_t paramBytes, int rk, int *flags,
    char *err, int errLen) {
    @autoreleasepool {
        MusicMetalDevice *m = get(dev);
        if (ensurePipelines(m, err, errLen)) return 1;
        if (paramBytes > 1024) {
            setErr(err, errLen, "StepParams too large", nil);
            return 1;
        }
        unsigned char P[1024];
        memcpy(P, params, paramBytes);
        memcpy(P + MUSIC_DEVICE_RKFLAG_OFFSET, &rk, sizeof(int));
        const bool viscous = (m->nComp == kMetalNCompVisc);
        if (rk == 0) {
            if (waitPending(m, err, errLen)) return 1;
            memset([m->flags contents], 0, sizeof(int));
        }
        const int n = m->nx * m->ny * m->neta;
        id<MTLBuffer> sIn = m->slice[rk == 0 ? m->cur : m->pred];
        // the second stage averages with tau_n; the viscous update of the
        // first stage takes its time derivatives against tau_{n-1}
        id<MTLBuffer> sPrev = m->slice[rk == 0 && viscous ? m->next : m->cur];
        id<MTLBuffer> sOut = m->slice[rk == 0 ? m->pred : m->next];
        id<MTLBuffer> pIn = (rk == 0 ? m->pCur : m->pPred);
        const size_t pb = paramBytes;
        // each kernel binds everything it reads, so that the profile mode
        // can give each one a command buffer of its own
        typedef void (^Encode)(id<MTLComputeCommandEncoder>);
        NSMutableArray *work = [NSMutableArray new];
        NSMutableArray *names = [NSMutableArray new];
        const unsigned char *Pp = P;
        // the strings source from the u of the stage's input slice
        int stringSource = 0;
        memcpy(&stringSource, P + MUSIC_DEVICE_STRINGS_OFFSET, sizeof(int));
        if (stringSource) {
            if (m->strings == nil) {
                setErr(err, errLen, "the strings were not set", nil);
                return 1;
            }
            [names addObject:@"strings"];
            [work addObject:[^(id<MTLComputeCommandEncoder> enc) {
                [enc setBytes:Pp length:pb atIndex:0];
                [enc setBuffer:m->strings offset:0 atIndex:16];
                [enc setBuffer:sIn offset:0 atIndex:3];
                [enc setBuffer:m->src offset:0 atIndex:15];
                [enc setBuffer:m->flags offset:0 atIndex:11];
                dispatch(enc, m->psStrings, n);
            } copy]];
        }
        // the pressure of the stage's input slice
        [names addObject:@"pressure"];
        [work addObject:[^(id<MTLComputeCommandEncoder> enc) {
            [enc setBytes:Pp length:pb atIndex:0];
            [enc setBuffer:m->eosDesc offset:0 atIndex:1];
            [enc setBuffer:m->eosTable offset:0 atIndex:2];
            [enc setBuffer:sIn offset:0 atIndex:3];
            [enc setBuffer:pIn offset:0 atIndex:4];
            dispatch(enc, m->psPressure, n);
        } copy]];
        // the faces
        for (int d = 1; d <= 3; d++) {
            const int nf = faceCount(m->nx, m->ny, m->neta, d);
            [names addObject:@"face"];
            [work addObject:[^(id<MTLComputeCommandEncoder> enc) {
                int dir = d;
                [enc setBytes:Pp length:pb atIndex:0];
                [enc setBuffer:m->eosDesc offset:0 atIndex:1];
                [enc setBuffer:m->eosTable offset:0 atIndex:2];
                [enc setBuffer:sIn offset:0 atIndex:3];
                [enc setBuffer:pIn offset:0 atIndex:4];
                [enc setBuffer:m->flux[d - 1] offset:0 atIndex:5];
                [enc setBuffer:m->flags offset:0 atIndex:6];
                [enc setBytes:&dir length:sizeof(int) atIndex:7];
                dispatch(enc, m->psFace, nf);
            } copy]];
            if (viscous) {
                [names addObject:@"wface"];
                [work addObject:[^(id<MTLComputeCommandEncoder> enc) {
                    int dir = d;
                    [enc setBytes:Pp length:pb atIndex:0];
                    [enc setBuffer:sIn offset:0 atIndex:3];
                    [enc setBuffer:m->wFace[d - 1] offset:0 atIndex:5];
                    [enc setBytes:&dir length:sizeof(int) atIndex:7];
                    dispatch(enc, m->psWFace, nf);
                } copy]];
            }
        }
        // the cells
        [names addObject:(viscous ? (rk == 0 ? @"cellVisc_rk0" : @"cellVisc_rk1")
                                  : @"cell")];
        [work addObject:[^(id<MTLComputeCommandEncoder> enc) {
            [enc setBytes:Pp length:pb atIndex:0];
            [enc setBuffer:m->eosDesc offset:0 atIndex:1];
            [enc setBuffer:m->eosTable offset:0 atIndex:2];
            [enc setBuffer:sIn offset:0 atIndex:3];
            [enc setBuffer:sPrev offset:0 atIndex:4];
            [enc setBuffer:pIn offset:0 atIndex:5];
            [enc setBuffer:m->pCur offset:0 atIndex:6];
            for (int d = 0; d < 3; d++)
                [enc setBuffer:m->flux[d] offset:0 atIndex:7 + d];
            [enc setBuffer:sOut offset:0 atIndex:10];
            [enc setBuffer:m->flags offset:0 atIndex:11];
            [enc setBuffer:m->src offset:0 atIndex:15];
            if (viscous) {
                for (int d = 0; d < 3; d++)
                    [enc setBuffer:m->wFace[d] offset:0 atIndex:12 + d];
                dispatch(enc, m->psCellVisc[rk], n);
            } else {
                dispatch(enc, m->psCell, n);
            }
        } copy]];
        id<MTLCommandBuffer> cb = nil;
        if (m->profile) {
            // one command buffer per kernel, timed on the GPU
            for (NSUInteger k = 0; k < [work count]; k++) {
                id<MTLCommandBuffer> one = [m->queue commandBuffer];
                id<MTLComputeCommandEncoder> enc = [one computeCommandEncoder];
                ((Encode)work[k])(enc);
                [enc endEncoding];
                [one commit];
                [one waitUntilCompleted];
                if (one.error != nil) {
                    setErr(err, errLen, "a kernel failed", one.error);
                    return 1;
                }
                NSString *key = names[k];
                NSNumber *old = m->profileSeconds[key];
                m->profileSeconds[key] =
                    @((old ? [old doubleValue] : 0.)
                      + (one.GPUEndTime - one.GPUStartTime));
                cb = one;
            }
        } else {
            cb = [m->queue commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
            for (Encode e in work) e(enc);
            [enc endEncoding];
            [cb commit];
        }
        *flags = 0;
        if (rk == 0) {
            m->pending = cb;
            return 0;
        }
        // the second stage runs while the host goes on (metalFinish waits
        // for it); every read of the slices waits for it first
        m->pending = cb;
        const int t = m->cur;
        m->cur = m->next;
        m->next = t;
        return 0;
    }
}

static void *metalSource(void *dev) { return [get(dev)->src contents]; }

static int metalSetStrings(
    void *dev, const void *data, size_t bytes, char *err, int errLen) {
    @autoreleasepool {
        MusicMetalDevice *m = get(dev);
        // the stages in flight may still read the old buffer
        if (waitPending(m, err, errLen)) return 1;
        if (m->strings == nil || [m->strings length] < bytes) {
            m->strings = [m->dev newBufferWithLength:(bytes > 0 ? bytes : 4)
                                             options:MTLResourceStorageModeShared];
            if (m->strings == nil) {
                setErr(err, errLen, "out of device memory", nil);
                return 1;
            }
        }
        memcpy([m->strings contents], data, bytes);
        return 0;
    }
}

static int metalFinish(void *dev, int *flags, char *err, int errLen) {
    @autoreleasepool {
        MusicMetalDevice *m = get(dev);
        if (waitPending(m, err, errLen)) return 1;
        *flags = *(int *)[m->flags contents];
        return 0;
    }
}

static int metalSpecialize(void *dev, const char *key, char *err, int errLen) {
    @autoreleasepool {
        MusicMetalDevice *m = get(dev);
        // MUSIC_METAL_NOSPEC: the generic kernels (debugging)
        if (getenv("MUSIC_METAL_NOSPEC") != NULL) return 0;
        NSString *k = [NSString stringWithUTF8String:key];
        if (m->specKey != nil && [m->specKey isEqualToString:k]) return 0;
        if (waitPending(m, err, errLen)) return 1;
        NSMutableDictionary *macros = [NSMutableDictionary new];
        macros[@"MUSIC_SPEC"] = @1;
        for (NSString *item in [k componentsSeparatedByString:@";"]) {
            NSArray<NSString *> *kv = [item componentsSeparatedByString:@"="];
            if ([kv count] != 2) {
                setErr(err, errLen, "a malformed specialisation key", nil);
                return 1;
            }
            macros[[@"MUSIC_SPEC_" stringByAppendingString:kv[0]]] =
                @([kv[1] intValue]);
        }
        // MUSIC_METAL_SPEC_SKIP=name,...: these switches stay run-time
        // fields (debugging; the kernels name their parameters P and d)
        const char *skip = getenv("MUSIC_METAL_SPEC_SKIP");
        for (NSString *name in (skip != NULL
                                    ? [[NSString stringWithUTF8String:skip]
                                          componentsSeparatedByString:@","]
                                    : @[])) {
            if ([name length] == 0) continue;
            macros[[@"MUSIC_SPEC_" stringByAppendingString:name]] =
                ([name isEqualToString:@"eosType"]
                     ? @"((d).type)"
                     : [NSString stringWithFormat:@"((P).%@)", name]);
        }
        if (buildPipelines(m, macros, err, errLen)) return 1;
        m->specKey = k;
        return 0;
    }
}

static int metalConservation(
    void *dev, const void *params, size_t paramBytes, void *out,
    size_t outBytes, char *err, int errLen) {
    @autoreleasepool {
        MusicMetalDevice *m = get(dev);
        if (waitPending(m, err, errLen)) return 1;
        if (ensurePipelines(m, err, errLen)) return 1;
        if (m->consOut == nil || [m->consOut length] < outBytes) {
            m->consOut = [m->dev newBufferWithLength:outBytes
                                             options:MTLResourceStorageModeShared];
            if (m->consOut == nil) {
                setErr(err, errLen, "out of device memory", nil);
                return 1;
            }
        }
        id<MTLCommandBuffer> cb = [m->queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setBytes:params length:paramBytes atIndex:0];
        [enc setBuffer:m->eosDesc offset:0 atIndex:1];
        [enc setBuffer:m->eosTable offset:0 atIndex:2];
        [enc setBuffer:m->slice[m->cur] offset:0 atIndex:3];
        [enc setBuffer:m->slice[m->next] offset:0 atIndex:4];
        [enc setBuffer:m->consOut offset:0 atIndex:5];
        // one thread per (eta, x) row: few threads, so small threadgroups
        // spread them over the GPU cores
        [enc setComputePipelineState:m->psConservation];
        [enc dispatchThreads:MTLSizeMake(m->nx * m->neta, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
        [enc endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
        if (cb.error != nil) {
            setErr(err, errLen, "the conservation kernel failed", cb.error);
            return 1;
        }
        if (m->profile) {
            NSNumber *old = m->profileSeconds[@"conservation"];
            m->profileSeconds[@"conservation"] =
                @((old ? [old doubleValue] : 0.)
                  + (cb.GPUEndTime - cb.GPUStartTime));
        }
        memcpy(out, [m->consOut contents], outBytes);
        return 0;
    }
}

static int metalFetch(
    void *dev, int which, int firstComp, int nComp, char *err, int errLen) {
    (void)which;
    (void)firstComp;
    (void)nComp;
    // shared memory: only wait for the stage in flight
    return waitPending(get(dev), err, errLen);
}

static int metalFetchUnwritten(
    void *dev, int which, int firstComp, int nComp, char *err, int errLen) {
    (void)dev;
    (void)firstComp;
    (void)nComp;
    // shared memory, and no stage in flight writes the slice
    if (which != 1) {
        setErr(err, errLen, "fetchUnwritten: only the slice of tau_n", nil);
        return 1;
    }
    return 0;
}

static const MusicDeviceApi kMetalApi = {
    metalCreate, metalDestroy, metalName,   metalResize, metalSetEos,
    metalState,  metalNoCopy,  metalStage,  metalSource, metalNoCopy,
    metalFetch,  metalNoCopy,  metalSetStrings, metalConservation,
    metalSpecialize, metalFinish, metalFetchUnwritten,
};

const MusicDeviceApi *music_metal_api(void) { return &kMetalApi; }
