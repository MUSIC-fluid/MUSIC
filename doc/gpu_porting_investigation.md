# Porting MUSIC to GPUs (CUDA and Metal): investigation

Branch `gpu`, from `optimization` at 3c2beac. 2026-11-22. Machine: Apple M4 (10 CPU cores, 4 performance and 6 efficiency; 10-core GPU; Metal 4; 16 GB unified memory). No NVIDIA hardware or CUDA toolkit is available here.

## Summary

- **A port is feasible, and the hydro core suits GPUs well.** Every hot loop in a time step is a loop over cells or faces with no dependencies inside it. That covers the KT fluxes, the velocity root solves, the cell update, the W relaxation and the stage-2 fixed point. Freeze-out, output and diagnostics happen between steps and can stay on the host.
- **Prototype:** `gpu/proto/`. It is the ideal, boost-invariant KT update (both `face_reconstruction` modes, Heun RK2, MUSIC's velocity solve, the conformal EoS with n_B). It is written once, in a portable header, and compiled as C++/OpenMP, CUDA and Metal Shading Language.
  - In double on the CPU it reproduces MUSIC's own output to 5e-9, the 8 digits MUSIC writes, on the Gubser integration test.
  - On the M4 GPU, in float, it agrees with MUSIC to 1e-5 and gives the same L1 errors against the analytic solution.
- **Speed (ideal (2+1)D, M4 GPU against MUSIC on all 10 CPU cores):**
  - 1.7 ns per cell and RK stage in mode 1 and 3.1 ns in mode 0, from 1001² cells up. MUSIC takes 86–98 ns and 135–145 ns.
  - That is **40–60× faster than MUSIC**, and **9–15× faster than the same specialised kernels on the 10 CPU cores**. The second number is the fair, like-for-like gain from the GPU.
- **Precision is the main issue for Metal, which has no fp64.** As written, the scheme rebuilds τT^{τμ} from the primitive variables in every stage. In float that accumulates roundoff: 2e-4 relative in e after 400 steps, at fast-moving cells.
  - Storing the conserved τT^{τμ}, τJ^τ next to the primitives cuts this **50×, to 4e-6**. In double it changes nothing (1e-13). It costs 15–27% speed.
  - With that change, float is good enough for the ideal part. The viscous part still has to be tested.
  - CUDA on data-center GPUs (A100, H100) should run in double. Consumer NVIDIA cards have 1/32–1/64-rate fp64, so they need the float path too.
- **Recommended approach:** a single-source kernel layer, a header of plain functions with three thin backends (OpenMP, CUDA/HIP, Metal), as in the prototype.
  - Metal accepts only MSL, which has no virtual functions, no recursion, no standard library and no double. That rules out Kokkos, SYCL, OpenMP offload and std::par as a common layer: none of them has a Metal backend.
  - The main refactoring work is the EoS: a virtual class hierarchy has to become data (tables plus a type tag). The rest is mechanical: logging and `exit()` in the hot code, Eigen in the causality check, and the host-side source terms.

## 1. Where the time goes and how it maps to kernels

One RK stage of `Advance::AdvanceIt` (src/advance.cpp):

| step | MUSIC code | parallel over | GPU kernel |
|---|---|---|---|
| ghost cells | `Fields::updateGhosts` | ghost cells | copy kernel, or clamped indices as in the prototype |
| pressure | `computePressure` | cells | fused into its consumers (prototype) |
| KT face fluxes | `computeFaceFluxes` → `computeFaceFlux` / `computeFaceFluxPrimitive`, `MaxSpeed`, `ReconstIt_shell` | faces, 3 directions | one kernel per direction (prototype: `faceFlux`) |
| W fluxes | `diss_helper.computeWDivergences` | lines and faces | face kernel plus difference |
| cell update | `FirstRKStepT` (ideal) | cells | `cellUpdate` (prototype) |
| viscous update | `MakedU`, gradients, `Make_uWRHS`, `getThermalVariables`, `FirstRKStepW`, `stage2ViscousFixedPoint` | cells | 1–3 cell kernels |
| between steps | freeze-out (Cornelius), output, conservation diagnostic, beastMode regrid and Δτ changes | — | host |

Profile of a 3D EoS 9 shear+bulk run (from [performance_optimizations_vs_chun_dev.md](performance_optimizations_vs_chun_dev.md), 10 threads):

- **mode 1:** faces 1.1 s, cell loop 4.5 s, W divergences 2.0 s, rest 0.6 s.
- **mode 0:** faces 4.8 s, cell loop 4.4 s.

So a viscous run needs the dissipative part (src/dissipative.cpp, 1740 lines) on the GPU as well. With only the ideal part ported, the speed-up is capped at about 1.5×.

## 2. The prototype

Files in `gpu/proto/`:

- `kt2d.h`: the kernels.
- `kt2d_main.cpp`: the CPU backend and the driver.
- `kt2d_metal.mm`: the Metal backend. The MSL is compiled at run time with `newLibraryWithSource`, so no offline `metal` compiler is needed (the Command Line Tools do not include one).
- `kt2d_cuda.cu`: the CUDA backend.

**Portability layer.** A dozen macros cover the differences between the three compilers:

- Address spaces: `device`, `constant` and `thread` in MSL.
- `__device__` for CUDA.
- Atomic flag updates.
- Math functions.

The kernels are templates on the floating type. Errors come back as flag bits ORed atomically (revert, regulate, momentum cap, MaxSpeed anomalies) instead of messages and `exit()`.

**What the prototype ports.** `computeFaceFlux`, `computeFaceFluxPrimitive`, `MakeDeltaQI` for a boost-invariant grid, `FirstRKStepT`, `MaxSpeed`, `ReconstIt_shell`, `ReconstIt_velocity_Newton`, `solve_v_from_guess` and `solve_v_Hybrid`. Each keeps MUSIC's branches, tolerances and fallbacks (including the transverse momentum cap and `regulate_grid` keeping n_B).

It leaves out the u^τ solve for v > v_critical, which MUSIC uses as a refinement; the double results agree with MUSIC without it. Data layout: structure of arrays with y fastest, zero-gradient edges through clamped indices (MUSIC's ghost cells hold the same values).

**CUDA.** `kt2d_cuda.cu` cannot be run here. `make cuda-syntax` type-checks it with clang's CUDA front end, on both the host side and the sm_80 device side, with all template instantiations, stubbing the runtime API. It builds unchanged with `nvcc` (`make cuda`).

### 2.1 Validation against MUSIC

Gubser test of `tests/integration` (`gubser_face_reconstruction`): 201² cells, dx 0.05 fm, τ 1 → 1.5, conformal EoS with n_B. MUSIC resets Δτ to 0.1·dx = 0.005 fm (100 steps). The initial state is MUSIC's own input file. The table gives the relative L1 difference from MUSIC's `Gubser_flow_check_tau_1.5.dat`, with the largest single-cell difference in brackets.

| | mode | e | n_B | u^x |
|---|---|---|---|---|
| CPU double | 0 | 8.8e-10 (4.9e-9) | 1.3e-9 (5.0e-9) | 8.8e-10 (5.0e-9) |
| CPU double | 1 | 8.8e-10 (4.9e-9) | 1.3e-9 (4.8e-9) | 9.0e-10 (4.9e-9) |
| CPU float | 1 | 2.1e-5 | 7.6e-6 | 8.4e-6 |
| Metal float | 0 | 1.1e-5 | 3.9e-6 | 4.7e-6 |
| Metal float | 1 | 1.1e-5 | 4.1e-6 | 4.7e-6 |

The double kernel agrees with MUSIC to the precision MUSIC writes, so it is the same algorithm. Its L1 errors against the analytic solution on y = 0 equal MUSIC's:

- **mode 1:** 5.186e-4 (e), 2.635e-4 (n_B), 6.485e-4 (u^x).
- **mode 0:** 1.92e-3, 1.36e-3, 1.33e-3.

In float the errors against the analytic solution are 5.09e-4, 2.59e-4 and 6.49e-4 (Metal, mode 1). The roundoff is far below the discretization error.

### 2.2 Precision in float

The same Gubser run continued to τ = 3 (400 steps). The table gives the relative L1 difference from the double run.

| | e | n_B | u^x |
|---|---|---|---|
| CPU float, as MUSIC (τ = 1.5 / 2 / 3) | 2.1e-5 / 6.1e-5 / 2.2e-4 | 7.6e-6 / 2.5e-5 / 9.8e-5 | 8.4e-6 / 2.3e-5 / 8.0e-5 |
| Metal float, as MUSIC (τ = 3) | 1.0e-4 | 5.5e-5 | 4.8e-5 |
| CPU float, stored τT^{τμ} (τ = 3) | 3.9e-6 | 2.3e-6 | 1.9e-6 |
| Metal float, stored τT^{τμ} (τ = 3) | 3.8e-6 | 2.3e-6 | 1.8e-6 |
| CPU double, stored τT^{τμ} (τ = 3) | 1.6e-13 | 0 | 9e-18 |

For comparison, the discretization error against the analytic solution at τ = 3 is 1.8e-3 in e. Grid noise in the result, not the solver, set the limit:

- **Error growth:** about steps^1.6, which is systematic, not a random walk.
- **Location:** the error is largest (up to 5e-4) where |u| ≈ 2.4–2.8. There e = T^{ττ} − v·M cancels by about (4/3)u₀² ≈ 10.
- **Ruled out:**
  - The velocity-solve tolerance: tightening the float tolerances from (1e-8, 1e-6) to (0, 0) gives 2.0e-4.
  - The renormalisation of u^μ: written as (v u₀)² instead of u₀² − 1, no change.
  - The storage of the face fluxes: rounded to float in a double run, 1e-7.
  - The storage of the state: rounded to float, 2e-6.
- **Cause:** the round trip from primitives to τT^{τμ} and back in every stage. MUSIC starts each stage from `get_TJb(cell)·τ`.
- **Fix:** store the conserved variables and update those (`--store-q 1`). Then only the fluxes see the reconstruction roundoff, and the error stays at a few 1e-6.

On the M4 the fix costs 15–27% speed: 8 more floats per cell, and the kernels are memory-bound.

Open: the viscous part takes differences of u^μ between neighbouring cells (gradients, θ, σ^{μν}). In float those carry about 1e-7·u/(dx·∂u) relative error, which may be visible where gradients are small. It has to be measured once the dissipative kernels exist. The fallback is double-float (float2) arithmetic for the accumulations only.

### 2.3 Speed

Ideal (2+1)D Gubser flow (q scaled so the flow fills the grid), dx 0.05 fm, Δτ 0.005, 50 steps. The table gives ns per cell and RK stage (wall time ÷ (2 × steps × cells)). The machine load average was about 1.3.

| grid | mode | MUSIC, 10 threads | prototype CPU double, 10 thr | prototype CPU float, 10 thr | Metal float | Metal float, stored q |
|---|---|---|---|---|---|---|
| 1001² | 0 | 135 | 53.2 | 43.3 | 3.35 | 3.70 |
| 1001² | 1 | 85.8 | 20.5 | 17.6 | 1.66 | 2.12 |
| 2001² | 0 | 145 | 47.1 | 36.0 | 3.07 | 3.54 |
| 2001² | 1 | 97.8 | 18.7 | 15.6 | 1.67 | 2.12 |
| 4001² | 0 | — | 44.9 | 32.1 | 3.03 | 3.48 |
| 4001² | 1 | — | 18.1 | 14.8 | 1.65 | 2.09 |

- **MUSIC's timing:** the median wall time of one time step, from time-stamped log lines. That excludes the file reading and output (their full runs took 19–102 s).
- **Why the prototype CPU is 3–5× faster than MUSIC:** it is specialised to (2+1)D with 4 components. MUSIC runs its general 3D path: 7 components, an η direction of one cell, ghost updates and per-step diagnostics.
- **The GPU gain:** about 9–15× over the 10 CPU cores running the same kernels, and 40–60× over MUSIC.
- **Run-time MSL compile:** 0.7 s the first time, then cached by the OS.

**The kernels are memory-bound.** Mode 1 moves about 150 B per cell and stage, which at 1.66 ns is about 90 GB/s against the M4's 120 GB/s. The extra 32 B of the stored q costs 0.45 ns, the same rate. The next gain is therefore to fuse the face fluxes into the cell kernel through a tile in threadgroup or shared memory. That removes the 64 B of flux written and read back per cell, worth an estimated 1.5×.

On a data-center NVIDIA GPU (1.5–3 TB/s), bandwidth suggests 10–25× the M4's float throughput, and about half of that in double. This is an estimate, not a measurement.

## 3. What has to change in MUSIC

| issue | where | remedy |
|---|---|---|
| EoS is a virtual hierarchy, 22 virtual methods in `EOS_base`, 10 derived classes (EOS_BEST, 1DGenerator, 4D, neos, eosQ, hotQCD, idealgas, WB, s95p, UH) | `EOS::eos_ptr`, called in every kernel | A device EoS: a POD descriptor (type tag, table sizes, offsets into the flat `EosTable` buffer that 7609fb7 already introduced) and `switch`-dispatched interpolation functions in the portable header. The host classes keep loading and checking the tables. The analytic EoSs (ideal gas with the x(ρ) Newton solve, eosQ) become small device functions. |
| Logging and `exit()` in hot code (72 sites in advance.cpp, reconst.cpp, dissipative.cpp) | velocity solve, MaxSpeed, W checks | error flag bits plus the first failing cell index (atomic), checked by the host after each step, which then reports or aborts as now |
| Eigen (`EigenSolver` on the 4×4 π^μ_ν in the causality constraints) | `solveEigenvaluesWmunu` | a closed-form or Jacobi solver: π^{μν} is symmetric after lowering one index with g^{1/2}, so a 4×4 symmetric Jacobi with a fixed number of sweeps fits a GPU |
| `std::vector`, `std::array` (`TJbVec`, `ThermalVec`, `FlowVec`) | everywhere | fixed-size POD arrays in the header; MSL has no standard library |
| Hydro source terms (TATB, strings) | `hydro_source_terms_ptr->get_hydro_*_source`, per cell and stage | evaluate on the host into a per-cell array for each stage, or port the smearing kernel (it is a sum over sources: parallel over cells) |
| Iterative solvers (velocity solve, u^τ solve, stage-2 fixed point) | reconst.cpp, advance.cpp | fine on a GPU: bounded loops without recursion. Divergence between threads costs little at about 3 secant steps per solve (prototype). |
| OpenMP pragmas, dynamic scheduling | all loops | the loops become kernel launches; the CPU backend keeps OpenMP |
| No fp64 on Metal | everything | stored conserved variables (§2.2); float for the ideal part; the viscous part still to be measured |
| beastMode regrid, Δτ changes, freeze-out, output | evolve.cpp, freeze-out | stay on the host between steps. On Apple's unified memory the host reads shared buffers without copies. On CUDA, copy the slices the surface finder needs, or port the Cornelius cube scan (parallel over cells). |

**Memory.** `Fields` holds 24 doubles per cell. Three slices, the face fluxes (3 × 7), the W buffers and the pressure come to about 1 KB per cell in double.

- A 261² × 64 3D grid (4.4 M cells) needs about 4.4 GB in double and 2.2 GB in float.
- That fits on any current data-center GPU (40–80 GB) and in the M4's 16 GB in float.
- (2+1)D grids are small in comparison.

## 4. Backend options

| option | CUDA | AMD (HIP) | Metal | notes |
|---|---|---|---|---|
| **portable header plus thin backends (prototype)** | yes | yes (HIP is CUDA-compatible) | yes | The only option that covers Metal; the CPU code path stays the same source; no new dependency. |
| Kokkos | yes | yes | no | mature, but no Metal backend |
| SYCL (DPC++, AdaptiveCpp) | yes | yes | no | no production Metal backend |
| OpenMP target offload | yes (nvc++, clang) | yes | no | |
| C++ std::par (nvc++) | yes | no | no | NVIDIA only |
| native CUDA only | yes | via hipify | no | |

**Recommendation.** Build out the portable kernel header that the prototype introduces:

- The kernels are plain templated functions on POD data.
- The CPU build calls them in OpenMP loops, so the CPU path remains the reference and is tested by the existing integration tests.
- CUDA and HIP use it in double by default.
- Metal uses it in float, compiled at run time from the same header.

The cost is that MSL limits the language used in kernels (no virtual functions, no recursion, no standard library, explicit address spaces). The prototype shows these limits are workable for MUSIC's kernels.

## 5. Suggested plan

Effort is given in rough working weeks for one developer familiar with the code. These are estimates.

1. **Ideal 3D on device, conformal and tabulated EoS (2–3 weeks).**
   - The prototype extended to 7 components and the η direction with its geometric terms.
   - The u^τ solve for v > v_critical.
   - A device EoS for the table EoSs (EoS 9, 14 first) and the ideal gas.
   - Wired into `Evolve` behind a run-time switch (`gpu_backend = none/cuda/metal`), with host freeze-out from device slices.
   - Validation: the ideal integration tests in double on CUDA (expected agreement about 1e-12 with the CPU) and in float on Metal.
2. **Dissipative part (3–5 weeks):**
   - W fluxes and divergences, gradients, `Make_uWRHS`, `FirstRKStepW`, the stage-2 fixed point and the causality limiter.
   - Measure float precision on Metal against the viscous Gubser and lumpy tests, with the conservation diagnostic.
3. **Sources and the remaining EoSs (1–2 weeks):** the TATB/string sources and the 4D EoS (EoS 20).
4. **Performance (1–2 weeks):** fuse the face fluxes into the cell kernel, keep slices resident between steps, overlap freeze-out on the host with the next steps.

beastMode's regrid and Δτ doubling stay on the host and copy the state back and forth; they happen rarely.

## Reproducing

```
cd gpu/proto && make                      # CPU + Metal (Apple clang, Homebrew libomp)
make cuda-syntax                          # type-check the CUDA backend (host and sm_80)
python3 tests/integration/music_integration_tests.py gubser_face_reconstruction \
    --music <build>/src/MUSIChydro --workdir <dir>      # MUSIC reference
./kt2d --mode 1 --dtau 0.005 --threads 10 \
    --init <dir>/gubser_face_reconstruction_1/tests/Gubser_flow/y=0_tau=1.00_ideal.dat \
    --ref  <dir>/gubser_face_reconstruction_1/Gubser_flow_check_tau_1.5.dat \
    [--backend metal] [--precision float] [--store-q 1]
./kt2d --n 2001 --size 100 --q 0.1 --mode 1 --dtau 0.005 --tau1 1.25 --backend metal
```

MUSIC's large-grid runs use the same input as the integration test, with the grid, `face_reconstruction`, Δτ 0.005 and 0.25 fm/c of evolution changed. Their initial files come from `kt2d --write-init`.
