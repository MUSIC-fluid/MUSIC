// Copyright 2026 Chun Shen
//
// The Metal entry points of the device kernels. CMake bundles this file
// after device_defs.h, device_eos.h, device_ideal.h and device_visc.h into
// one source string (metal_source.inc) that metal_backend.m compiles at run
// time (with device_strings.h and device_diag.h).

using namespace music_device;
typedef StepParams<float> MSP;
typedef EosDesc<float> MED;

kernel void musicPressure(
    constant MSP &P [[buffer(0)]], constant MED &D [[buffer(1)]],
    device const float *tab [[buffer(2)]],
    device const float *s [[buffer(3)]], device float *p [[buffer(4)]],
    uint gid [[thread_position_in_grid]]) {
    if (int(gid) >= P.nx * P.ny * P.neta) return;
    kernelPressure(P, D, tab, s, p, int(gid));
}

kernel void musicFace(
    constant MSP &P [[buffer(0)]], constant MED &D [[buffer(1)]],
    device const float *tab [[buffer(2)]],
    device const float *s [[buffer(3)]],
    device const float *p [[buffer(4)]], device float *flux [[buffer(5)]],
    device metal::atomic_int *flags [[buffer(6)]],
    constant int &direction [[buffer(7)]],
    uint gid [[thread_position_in_grid]]) {
    if (int(gid) >= faceCount(P.nx, P.ny, P.neta, direction)) return;
    kernelFace(P, D, tab, s, p, direction, int(gid), flux, flags);
}

kernel void musicCell(
    constant MSP &P [[buffer(0)]], constant MED &D [[buffer(1)]],
    device const float *tab [[buffer(2)]],
    device const float *sCurr [[buffer(3)]],
    device const float *sPrev [[buffer(4)]],
    device const float *pCurr [[buffer(5)]],
    device const float *pPrev [[buffer(6)]],
    device const float *fluxX [[buffer(7)]],
    device const float *fluxY [[buffer(8)]],
    device const float *fluxEta [[buffer(9)]],
    device float *sNext [[buffer(10)]],
    device metal::atomic_int *flags [[buffer(11)]],
    device const float *src [[buffer(15)]],
    uint gid [[thread_position_in_grid]]) {
    if (int(gid) >= P.nx * P.ny * P.neta) return;
    kernelCell(P, D, tab, sCurr, sPrev, pCurr, pPrev, fluxX, fluxY, fluxEta,
               src, sNext, int(gid), flags);
}

kernel void musicWFace(
    constant MSP &P [[buffer(0)]], device const float *s [[buffer(3)]],
    device float *out [[buffer(5)]], constant int &direction [[buffer(7)]],
    uint gid [[thread_position_in_grid]]) {
    if (int(gid) >= faceCount(P.nx, P.ny, P.neta, direction)) return;
    kernelWFace(P, s, direction, int(gid), out);
}

kernel void musicCellVisc(
    constant MSP &P [[buffer(0)]], constant MED &D [[buffer(1)]],
    device const float *tab [[buffer(2)]],
    device const float *sCurr [[buffer(3)]],
    device const float *sPrev [[buffer(4)]],
    device const float *pCurr [[buffer(5)]],
    device const float *pPrev [[buffer(6)]],
    device const float *fluxX [[buffer(7)]],
    device const float *fluxY [[buffer(8)]],
    device const float *fluxEta [[buffer(9)]],
    device float *sNext [[buffer(10)]],
    device metal::atomic_int *flags [[buffer(11)]],
    device const float *wFaceX [[buffer(12)]],
    device const float *wFaceY [[buffer(13)]],
    device const float *wFaceEta [[buffer(14)]],
    device const float *src [[buffer(15)]],
    uint gid [[thread_position_in_grid]]) {
    if (int(gid) >= P.nx * P.ny * P.neta) return;
    kernelCellVisc(P, D, tab, sCurr, sPrev, pCurr, pPrev, fluxX, fluxY,
                   fluxEta, wFaceX, wFaceY, wFaceEta, src, sNext, int(gid),
                   flags);
}

kernel void musicStringSource(
    constant MSP &P [[buffer(0)]], device const float *S [[buffer(16)]],
    device const float *sIn [[buffer(3)]], device float *src [[buffer(15)]],
    device metal::atomic_int *flags [[buffer(11)]],
    uint gid [[thread_position_in_grid]]) {
    if (int(gid) >= P.nx * P.ny * P.neta) return;
    kernelStringSource(P, S, sIn, src, int(gid), flags);
}

kernel void musicConservation(
    constant MSP &P [[buffer(0)]], constant MED &D [[buffer(1)]],
    device const float *tab [[buffer(2)]],
    device const float *sCurr [[buffer(3)]],
    device const float *sPrev [[buffer(4)]], device float *out [[buffer(5)]],
    uint gid [[thread_position_in_grid]]) {
    if (int(gid) >= P.nx * P.neta) return;
    kernelConservationRow(P, D, tab, sCurr, sPrev, out, int(gid));
}
