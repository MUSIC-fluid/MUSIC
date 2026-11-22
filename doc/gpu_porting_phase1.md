# GPU port, phase 1: ideal hydro on the device

Branch `gpu`, from a7fd784. 2026-11-22. Machine: Apple M4 (10 CPU cores, 10-core GPU, 16 GB unified memory). There is no NVIDIA hardware or CUDA toolkit here. Follows the plan of [gpu_porting_investigation.md](gpu_porting_investigation.md).

## Summary

- **New input parameter `gpu_backend`.** When it is non-zero, the ideal hydro update of `Evolve::AdvanceRK` runs on a device:
  - `0`: off, the default. The CPU code runs unchanged.
  - `1`: the device kernels on the CPU (OpenMP), in double. This is the reference for the GPU backends.
  - `2`: Metal (Apple GPUs), in float.
  - `3`: CUDA (NVIDIA GPUs), in double by default.
- **Supported cases:** ideal hydro (`Viscosity_Flag_Yes_1_No_0 = 0`), 3D or boost-invariant, without source terms. Both `face_reconstruction` modes. EoS 0 (conformal), 9 and 91 (hotQCD) and 10–14 (NEOS).
  - Anything else prints a warning and runs on the CPU code.
  - Freeze-out, output, the conservation diagnostic, beastMode regridding and Δτ changes stay on the host, unchanged.
- **Validation:**
  - Backend 1 reproduces `Advance` **bit for bit** in double when compiled with clang. The unit test checks this for four configurations, and the integration test finds 0 difference for EoS 0, 9 and 14.
  - Metal (float) agrees with the CPU code to **~1e-6 of E** on a lumpy 3D event with n_B, for all three EoSs. On the Gubser test it has the same L1 errors against the analytic solution, to 1e-4 relative.
  - CUDA is type-checked only (no toolkit here).
- **Speed:** lumpy 3D event, 161×161×49 cells, 100 steps, against the CPU code on 10 threads.
  - **Metal is 10× faster in the hydro kernels.** Wall time is 5.6× faster for EoS 0 (6.0 s against 33.6 s) and 7× faster for NEOS (7.7 s against 53.9 s).
- **Results:** with `gpu_backend = 0` (the default) nothing changes. The Gubser output is byte-identical to the build before the port in both modes, and all integration tests pass.

## 1. Code

| file | what it does |
|---|---|
| `src/device/device_defs.h` | portability macros (C++ / CUDA / MSL), `dmax`/`dmin` with `std::max`/`std::min` semantics, error flags, sizes |
| `src/device/device_eos.h` | EoS as data: a POD descriptor plus one pressure table; `eosPressure`, `eosPressureGrad` written operation by operation like `EOS_base::interpolate1D/2D` |
| `src/device/device_ideal.h` | kernels `kernelPressure`, `kernelFace` (`computeFaceFlux[Primitive]`, `MaxSpeed`), `kernelCell` (`FirstRKStepT` + `MakeDeltaQI` + `ReconstIt_shell`, with the velocity Newton, hybrid and u⁰ solvers) |
| `src/device/device_api.h` | the C interface of the GPU backends (built by other compilers) |
| `src/device/device_hydro.h/.cpp` | `DeviceHydro`: backend choice, EoS tables, upload/download, the pointer rotation of `AdvanceRK`, flags, timing; the OpenMP backend |
| `src/device/metal_kernels.metal`, `metal_backend.m` | Metal entry points; Objective-C backend that compiles the bundled MSL at start-up |
| `src/device/cuda_backend.cu` | CUDA backend (double or float) |
| `cmake/bundle_msl.cmake` | bundles the kernel headers and the Metal entry points into one C string |

Changes to existing files:

- `evolve.cpp/h`: calls `DeviceHydro::advanceStep` in place of `AdvanceRK`. It invalidates the device copy on regrid and Δτ change, and adds `preparePressure` for the unit test.
- `read_in_parameters.cpp` and `data.h`: `gpu_backend`, `gpu_float`, `gpu_store_q`.
- Small accessors in `eos.h`, `reconst.h` and `advance.h`.
- CMake options `metal` and `cuda`.

**The state stays on the device.** The host Fields are uploaded only:

- at the first step,
- when the slice pointer or the grid changes,
- after an invalidate (regrid, Δτ change).

After each step the new slice is copied into the Fields that `AdvanceRK` would have filled. This costs 3.5 ms per step at 1.3M cells, about 10% of a Metal step. Freeze-out and output therefore see exactly what they saw before.

**Error handling.** The kernels set flags where MUSIC prints and exits. The host reports the MaxSpeed failure and exits like MUSIC. The other flags (revert, regulate, momentum cap, failed bracket) are printed at `echo_level > 5`.

## 2. Bit-for-bit reproduction and compiler contraction

At τ₀ = 1.3 the first comparison found differences of 1e-15. They were not in the fluxes, the pressure or the root solve. They came from the two-statement average in stage 2:

- `FirstRKStepT` writes `qi += rk_flag * q0; qi *= 1/(1 + rk_flag)` in one loop body.
- The device had `qi += q0; qi *= 0.5`.

clang fuses `a*b + c` within one statement, so `rk_flag * q0 + qi` becomes one FMA on the host and not on the device. With the host's form copied exactly, the device matches with clang:

- unit test (`Device kernels reproduce the ideal update of Advance`): maxDiff = 0 for 3D and boost-invariant, modes 0 and 1, at τ₀ = 1.3, over 4 steps;
- integration test: 0 for EoS 0, 9 and 14.

**g++ also fuses across statements**, differently in each translation unit. With g++-15 -O3, backend 1 differs from the CPU code by:

- 1e-16 of E for EoS 0 and EoS 9;
- for NEOS, 1e-19 noise that grows to 4e-8 by τ = 2.3 fm/c, as soon as one capped or reverted cell takes the other branch.

This is rounding noise, not a bug. Building with clang, or with `-ffp-contract=off` on both paths, gives bit equality. The integration test's limit for backend 1 is 1e-6, and the exact check stays in the unit test.

## 3. Metal (float)

- **Build:** `cmake -Dmetal=ON` (macOS). `metal_backend.m` is compiled by the OBJC language (clang) whatever compiles the C++, and links with `-framework Metal -framework Foundation`.
- **Kernel compile:** the MSL source is compiled at start-up with `MTLMathModeSafe` (IEEE division and sqrt). It takes ~0.5 s and needs no offline Metal compiler.
- **Porting fixes:** MSL needs `metal::`-qualified atomics, and it allows no `const` variables at program scope, so the constants are now enums.
- **Precision:**
  - Float runs with stored τT^{τμ} (`gpu_store_q`, default on in float).
  - The tolerances are raised to the beastMode > 0 values (1e-8, 1e-6), because MUSIC's 1e-16 and 1e-14 are below float resolution.
- **Memory:** the buffers are shared memory, so `upload` and `download` copy nothing. One Heun step (10 dispatches) is one command buffer.

Accuracy against the CPU code, from `tests/integration/music_integration_tests.py gpu_backend`. The lumpy 3D event has 41×41×17 cells, with n_B, and runs 1 fm/c:

| EoS | max \|Δ(E, P, N_B)\| / E, backend 1 (g++) | backend 2 (Metal) |
|---|---|---|
| 0 | 3.9e-16 | 1.1e-6 |
| 9 | 1.3e-16 | 5.6e-7 |
| 14 | 4.2e-8 | 7.5e-7 |

Gubser flow (201², mode 1, τ = 1.5 fm/c), L1 errors against the analytic solution:

| | e | n_B | u^x |
|---|---|---|---|
| CPU code | 5.186e-4 | 2.635e-4 | 6.485e-4 |
| Metal | 5.187e-4 | 2.635e-4 | 6.485e-4 |

The energy drift over the 3D event is 1.0e-6 on Metal and 1.2e-6 on the CPU. The N_B drift is 1.37e-6 on Metal and 1.32e-6 on the CPU.

## 4. Speed

Lumpy 3D event with n_B, 161×161×49 = 1.27M cells, Δτ = 0.01, 100 steps, mode 1, no freeze-out. CPU code on 10 threads (g++-15 -O3). Two repetitions agree to 1%.

| | EoS 0 wall | EoS 0 hydro kernels | EoS 14 wall | EoS 14 hydro kernels |
|---|---|---|---|---|
| CPU code (`gpu_backend 0`) | 33.6 s | — | 53.9 s | — |
| device kernels on CPU (1) | 33.9 s | 30.7 s | 45.9 s | 42.2 s |
| Metal M4 (2) | **6.0 s** | **3.1 s** | **7.7 s** | **4.2 s** |

- **Per cell and RK stage:** Metal takes 12 ns (EoS 0) and 16.5 ns (NEOS). The device kernels on the CPU take 120 ns and 165 ns.
- **Where the rest of Metal's wall time goes:** about 2.6 s of the 6 s is outside the hydro kernels. That is initialisation (the TATB profile), the host copy (0.35 s) and the per-step diagnostics on the host.
- **One unexplained slow run:** the first Metal run of this benchmark took 441 s, and it did not reproduce in three later runs (6.0 s each). That build had no timing output, so the cause is not known; a one-off first-time shader compile or GPU power state is the likely explanation. Watch for it.
- **Further gains:**
  - fuse the pressure into the face kernel, as the prototype did;
  - copy to the host only on steps that freeze-out or output need;
  - move the conservation diagnostic to the device.

## 5. CUDA

- **Build:** `cmake -Dcuda=ON` uses `enable_language(CUDA)`, `CUDA::cudart` and `CMAKE_CUDA_ARCHITECTURES=native`.
- **Memory:** the state is in device memory. `state()` is a pinned host copy, and `upload`/`download` are a `cudaMemcpy` of one slice each.
- **Check here:** both host and device passes type-check with clang's CUDA front end (`-DMUSIC_CUDA_SYNTAX_ONLY`, see the file header). Apple clang has no NVPTX target, so no PTX was generated and nothing was run.
- **Precision:** nvcc contracts FMAs by default. In double, expect rounding-level agreement with the CPU code, not bit equality; `--fmad=false` comes closer.
- **To do on an NVIDIA machine:** build, run `gpu_backend` with `gpu_backend=3` (the integration test can be extended the same way as for Metal), and benchmark double against float.

## 6. Tests

- `unittest_bugfix`, `Device kernels reproduce the ideal update of Advance`:
  - backend 1 in double against `Advance` over 4 steps (3D and boost-invariant, modes 0 and 1, τ₀ = 1.3): exact;
  - backend 1 in float with stored q: L1 < 1e-5, and max relative error < 1e-4 in cells with u^τ < 5. The largest error, 6e-4, is in a dilute cell with u^τ = 40, where e = T^{ττ} − vM loses (u^τ)² digits.
- Integration test `gpu_backend`:
  - backends 1 and 2 (when built) against 0 for EoS 0, 9 and 14;
  - Gubser accuracy on Metal.
- All 19 integration tests pass with the Metal build, and all unit tests pass.
- **Not yet checked:** a freeze-out surface produced with the device on. Freeze-out reads the downloaded Fields, so with backend 1 under clang it should be identical.

## 7. Next steps

1. Viscous hydro on the device: W fluxes, the `u_derivative` gradients, `Make_uWRHS`, the relaxation and the stage-2 fixed point. Test float precision there; it may need double for the W variables.
2. Source terms (TATB, strings) as device input per step.
3. Fewer host copies, and the conservation diagnostic on the device.
4. CUDA validation and timing on NVIDIA hardware.
