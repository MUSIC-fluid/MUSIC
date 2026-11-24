#ifndef GPU_PROTO_KT2D_METAL_H_
#define GPU_PROTO_KT2D_METAL_H_

#include "kt2d.h"

struct KtMetalTiming {
    double compileSeconds = 0.;  // run-time MSL compile + pipeline creation
    double runSeconds = 0.;      // wall time of the time steps
    double gpuSeconds = 0.;      // GPU execution time of the command buffers
    char device[128] = {0};
};

// advances state (9 n floats: e, n_B, u^tau, u^x, u^y, then tau T^{tau mu},
// tau J^tau; y fastest) by nSteps
// Heun steps. Returns the OR of the error flags, or -1 on a Metal error
int kt2d_metal_run(const kt::Params<float> *p, const float *table,
                   int tableLen, float *state, int nSteps, int stepsPerCommit,
                   KtMetalTiming *timing, char *errMsg, int errLen);

#endif  // GPU_PROTO_KT2D_METAL_H_
