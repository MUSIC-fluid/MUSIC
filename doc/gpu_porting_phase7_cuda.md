# GPU port, phase 7: the CUDA backend on NVIDIA hardware

Branch `gpu`, from fe606d7 (phase 6, [gpu_porting_phase6_performance.md](gpu_porting_phase6_performance.md)). 2026-11-22.

**Machine:**

- Intel i7-9750H (6 cores, 12 threads), 6 GB RAM.
- NVIDIA GeForce GTX 1650 (Turing, sm_75): 4 GB, PCIe 3 ×16, driver 620.70.06.
- CUDA 13.5, g++ 15.2.

This is the first time the CUDA backend (`gpu_backend = 3`) has been compiled by nvcc and run. Phases 1–6 could only type-check it with clang on the Mac.

## Summary

- **It needed one fix to compile.** `flagOr` was `__device__`-only but is called from `__host__ __device__` kernels; nvcc rejects that and clang's syntax check had not.
- **Once built, the backend was correct without further changes.** All six GPU integration tests passed in both precisions.
  - **CUDA double** agrees with the CPU code to 1e-16 of E. W, Π, Π_chem and μ_B agree exactly for EoS 0 and 9, and to 1e-7 or better (L1) for NEOS and EoS 20.
  - **CUDA float** agrees to 2–5e-6 of E, inside the limits of Metal float.
- **One latent bug, CUDA only, fixed:** with host-evaluated sources (e.g. TATB profile 112), the first RK stage read u from a host copy that is only refreshed when a diagnostic syncs (section 2).
- **Three performance changes.** All 156 output files of the benchmark are byte-identical before and after.
  1. **Steps run while the host works:** the asynchronous stage and finish that Metal got in phase 6.
  2. **Partial copies gather on the device:** the diagnostics' η windows and edge shells copy only their cells over PCIe, not the whole slice.
  3. **Block size and profiling knobs:** `MUSIC_CUDA_BLOCK` and `MUSIC_CUDA_PROFILE=1`.
- **Benchmark** (the phase 1–6 event: 161×161×49 = 1.27M cells, 101 steps; mean of 2 runs):

  | case | CPU code, 12 threads | CUDA before | **CUDA now** | speed-up over the CPU code |
  |---|---|---|---|---|
  | ideal, EoS 9, float | 100.2 s | 13.6 s | **9.9 s** (−27%) | **10.1×** |
  | viscous, EoS 9, float | 159.7 s | 31.3 s | **24.8 s** (−21%) | **6.4×** |
  | viscous, EoS 20, float | 260.1 s | 42.1 s | **35.9 s** (−15%) | **7.2×** |
  | ideal, EoS 9, double | 100.2 s | 59.3 s | **52.5 s** (−12%) | 1.9× |
  | viscous, EoS 9, double | 159.7 s | 105.3 s | **93.9 s** (−11%) | 1.7× |
  | viscous, EoS 20, double | 260.1 s | 160.0 s | **144.1 s** (−10%) | 1.8× |

  - **Use `gpu_float = 1` on GeForce cards.** The GTX 1650 runs double at 1/32 of its float rate, so double is only ~1.8× faster than the CPU.
  - On a data-centre card (A100, H100: 1/2 rate) double would be the natural choice; that was not tested here.
- **Tests:** all integration tests pass (25/25, including `memory_history` with the unit-test driver). The unit tests pass except one pre-existing float tolerance in `unittest_bugfix` under g++ (section 6).

## 1. Building

```bash
export PATH=/usr/local/cuda/bin:$PATH
cmake -S . -B build_cuda -DCMAKE_BUILD_TYPE=Release -Dcuda=ON
cmake --build build_cuda -j
```

`CMAKE_CUDA_ARCHITECTURES` defaults to `native` (sm_75 here).

**No external libraries are needed** (GSL stays optional). Eigen was removed after this phase: its only use was the eigenvalue solve of W^μ_ν in the causality limiter (`Advance::solveEigenvaluesWmunu`, `FlagResetCausality = 1`). That now uses the closed form the device code already had (`ViscSolver::wEigenvalues`):
- W is transverse to u, so its eigenvalues are 0 and the three real roots of a cubic in tr A, tr A², tr A³.
- **Host and device:** they now agree bit for bit in the unit test (before: up to 8e-12 after 4 steps).
- **Against Eigen,** on the lumpy viscous event:
  - conservation sums within 1.3e-16;
  - W/(e+P) within 3e-16, or 3e-11 with η/s = 0.5, where the limiter acts on many cells;
  - T and μ identical.
- All 25 integration tests pass without Eigen.

**`ptxas` resource use**, sm_75, from `nvcc -Xptxas -v`:

| kernel | float: registers / stack | double: registers / stack (spills) |
|---|---|---|
| face | 121 / 64 B | 199 / 1824 B |
| cellVisc | 255 / 1696 B | 255 / 3072 B (384 B of spill stores) |
| wface | 55 / 0 | 74 / 0 |
| cell (ideal) | — / 0 | — / 0 |

The large stack frames of cellVisc are its local arrays (the gradients, W, the transport coefficients), held in local memory. This is where the viscous kernel spends its time (section 5).

## 2. Fixes

| | file | what |
|---|---|---|
| compile | `device_defs.h` | `flagOr` is now `__host__ __device__`: `atomicOr` under `__CUDA_ARCH__`, `__atomic_fetch_or` on the host pass |
| correctness | `device_hydro.cpp` (`beginStep`) | host-evaluated source terms fetch u of the current slice in **both** RK stages (before: only the second) |
| build | `CMakeLists.txt`, `advance.cpp` | Eigen removed: the causality limiter's eigenvalues in closed form (section 1) |
| warnings | `device_hydro.cpp` | three `-Wenum-compare` (`kNCompVisc` and `kNComp` are in different anonymous enums) |

**The source-term bug.**
- **What goes wrong:** `DeviceHydro::beginStep` evaluates `HydroSourceBase` on the host from u in `state(which)`. For stage 0 it skipped the `fetch`. On Metal and backend 1 `state(0)` is the device slice itself, so that was fine. On CUDA `state(0)` is a pinned copy, refreshed only when a sync happens to run that step.
- **When it matters:** runs with beastMode > 0 skip the diagnostics, and so do boost-invariant runs (no shell). There, and between output steps, the source saw a stale u.
- **Why the tests passed anyway:** the integration tests run beastMode 0 with diagnostics every step, so they did not catch it.
- **The fix:** stage 0 now fetches the 4 components of u as well. That is a no-op on Metal and backend 1, and on CUDA a 20 MB (float) or 40 MB (double) copy per step at this grid size.

## 3. Performance changes

### 3.1 The step runs while the host works

- **Before:** CUDA `stage()` ended each RK stage with a blocking `cudaMemcpy` of the flags, and `finish()` was a stub that returned 0. So the GPU sat idle while the host did the output, the diagnostics and the freeze-out. The phase-6 overlap in `Evolve::EvolveIt` did nothing on CUDA.
- **Now, as on Metal:**
  - The flags are reset at stage 0 only and are OR'd over both stages.
  - Both stages are queued on the default stream without waiting.
  - `finish()` copies the flags back, which waits for the step.
- **Why this is safe:** every read of device state (`fetch`, `fetchCells`, `download`, `conservation`) is a blocking copy on the same stream, so it runs after the step. `DeviceHydro::requireIdle` still guards every read in Evolve.
- **The fatal flags still work:** `finish()` returns them, and `DeviceHydro::advanceStep` exits on MaxSpeed and NaN, as before.

### 3.2 Partial host copies gather on the device

**The problem:** each step, Evolve reads the η window −0.5 < η < 0.5 (5 of 49 slices, about 10% of the cells) and the 1-cell edge shell of the current and previous slices (about 6% each). `downloadCells` called `fetch()`, which on CUDA copies the **whole** slice: 31 components × 1.27M cells = 157 MB in float. That happened three times per step over PCIe. On Metal's unified memory it was free.

**New optional API call** `fetchCells(dev, which, cells, nCells, &out)` (`device_api.h`):

- The backend gathers all components of the listed cells into a compact component-major buffer and copies just that.
- The cell list stays on the device and is re-sent only when it changes.
- `DeviceHydro::downloadCells` builds the list from the same `keep` ranges as before. It writes the gathered values straight into the Fields, through the per-cell writer that `downloadT` now shares (`storeCell`).
- The full syncs (output and freeze-out steps) keep the old path.
- **Other backends:** a backend that leaves `fetchCells` NULL keeps `fetch()`. Metal sets it to NULL explicitly and backend 1 uses the base-class default, so neither changes.

| copy time per run (196 partial + 4 full syncs) | before | now |
|---|---|---|
| ideal, EoS 9, float | 2.70 s | 0.96 s |
| viscous, EoS 9, float | 6.05 s | 2.46 s |
| viscous, EoS 9, double | 10.40 s | 3.13 s |

**What is left:** at 11.5 GB/s (pinned, measured) the gathered data takes only about 0.4 s per run. Most of the remaining ~2 s is the host-side scatter into the Fields: ~5 ms for 82k cells × 31 components. Its cost does not change with the OpenMP thread count, so it is not a parallelism problem. It was not pursued further: kernels take 85% of the run.

### 3.3 Launch configuration and profiling

- **Block size:** `MUSIC_CUDA_BLOCK=<n>` sets the threads per block (default 128, a multiple of 32). It is capped per kernel at the `maxThreadsPerBlock` its registers allow, so cellVisc (255 registers) can never fail to launch with "too many resources".
- **The sweep** (viscous EoS 9, float, kernel time from the profile):

  | block | face | cellVisc | total |
  |---|---|---|---|
  | 64 | 4.53 s | 12.72 s | 18.68 s |
  | **128** | **3.94 s** | **12.35 s** | **17.71 s** |
  | 256 | 4.57 s | 13.51 s | 19.53 s |
  | 512 (cellVisc capped at 256) | 3.85 s | 13.51 s | 18.83 s |

  128 stays the default. Unlike Metal (phase 6: 1024 was best), the size barely matters here, because the viscous kernels are bound by latency and local memory.
- **`MUSIC_CUDA_PROFILE=1`:** times each kernel with CUDA events and prints a table at the end of the run, like `MUSIC_METAL_PROFILE`. Each launch is waited for, so the steps no longer overlap the host and the wall time rises a little.

## 4. Benchmark details

**Event and setup:**
- Lumpy TATB with n_B, generated by `_write_lumpy_TATB` of the integration tests, as in phases 1–6.
- Grid: 161×161×49, 16 fm × 16 fm × 10 units (dx = 0.1 fm), Δτ = 0.01, τ = 1 → 2.01 fm/c, 101 steps plus the initial one.
- No freeze-out, all diagnostics on.
- Viscous: η/s = 0.12, bulk, second-order terms.
- With the 12 fm box that I tried first (dx = 0.075), MUSIC shrinks Δτ to 0.0075 and runs 135 steps instead.

**How it was run:**
- The CPU code used all 12 threads.
- "CUDA before" is HEAD (fe606d7) plus only the `flagOr` compile fix; "CUDA now" is this phase.
- Each GPU case ran twice, before and after interleaved.

**Repeatability:**
- The "before" builds agreed with their repeats to ≤ 2%.
- The "now" builds varied by up to 12% (ideal double: 55.7 / 49.2 s). This is a laptop GPU (clocks and thermals), and the faster pipeline leaves the GPU less idle time to cool.
- All runs, wall time in seconds:

| case | CUDA before (run 1 / run 2) | CUDA now (run 1 / run 2) |
|---|---|---|
| ideal, EoS 9, float | 13.60 / 13.51 | 10.52 / 9.34 |
| viscous, EoS 9, float | 31.32 / 31.24 | 24.99 / 24.59 |
| viscous, EoS 20, float | 42.25 / 41.96 | 37.21 / 34.58 |
| ideal, EoS 9, double | 59.59 / 59.08 | 55.68 / 49.24 |
| viscous, EoS 9, double | 105.57 / 104.97 | 97.83 / 90.00 |
| viscous, EoS 20, double | 161.67 / 158.26 | 145.32 / 142.77 |

**Results:**
- All 156 output files of the six cases are byte-identical between "before" and "now":
  - `global_conservation_laws.dat`;
  - the eccentricities, momentum anisotropy, ⟨p_T⟩ and inverse Reynolds number files.
- CUDA double against the CPU code, max |Δ| / E over the run: 2e-17 (viscous EoS 9), 4e-17 (ideal EoS 9), 1.4e-10 (EoS 20, whose tables are float).

**Kernel time** (`MUSIC_CUDA_PROFILE=1`, 204 RK stages):

| kernel | ideal EoS 9 float | viscous EoS 9 float | viscous EoS 20 float | viscous EoS 9 double |
|---|---|---|---|---|
| pressure | 0.03 s | 0.03 s | 0.06 s | 0.11 s |
| face | 4.39 s (83%) | 5.19 s (25%) | 9.88 s (35%) | 36.87 s (44%) |
| wface | — | 1.26 s (6%) | 1.24 s (4%) | 13.14 s (16%) |
| cell / cellVisc | 0.78 s (15%) | 13.92 s (68%) | 16.51 s (59%) | 32.39 s (39%) |
| conservation | 0.06 s | 0.16 s | 0.19 s | 0.50 s |
| **total** | **5.26 s** | **20.55 s** | **27.87 s** | **83.00 s** |

**Per cell and RK stage:**

| | ideal | viscous |
|---|---|---|
| GTX 1650, float | 20 ns | 79 ns (EoS 9) |
| Apple M4, float (phases 1 and 6) | 12 ns | ~25 ns |

**Notes:**
- The ideal and viscous EoS 9 float columns are from separate profiled runs. Kernel times vary by up to 13% between them (cellVisc: 12.35 s in the block sweep, 13.92 s here), which is the thermal drift noted above.
- The ideal update is host-bound: 9.9 s of wall time against 5.3 s of kernels. The rest is the TATB initialisation and the diagnostics, as in phase 1.

## 5. What is left

**The viscous update:**
- **cellVisc** takes 60–70% of GPU time in float, at 255 registers with a 1.7 KB local-memory stack. Ways to shrink its live arrays:
  - split it into a gradients/coefficients kernel and an update kernel;
  - keep the 4×4 tensors symmetric (10 values, not 16);
  - use `__launch_bounds__` with a register budget.

  Each changes the order of float operations, so the result changes at rounding level.
- **Specialisation:** the switches of a run could be template parameters on CUDA, as `MD_RUN` does with macros on Metal. On Metal that removed 25% of the viscous kernel time (phase 6). On CUDA `cudaSpecialize` is still a no-op.

**The other kernels:**
- **face with EoS 20** (35%): the 4D table lookups in the conserved-to-primitive solve.
- **face in ideal hydro** (83% of its GPU time): at 121 registers it runs at low occupancy.

**Host and memory:**
- **Host scatter of the partial syncs:** ~2 s per viscous run (section 3.2).
- **Memory:** at 1.27M cells, viscous double uses ~2 GB of device memory. Peak host memory, as MUSIC logs it, is 2.55 GB (EoS 9) and 2.9 GB (EoS 20); viscous float uses 2.05 GB. Part of the host memory is the three pinned slice copies. That fits a 4 GB card but not much more: 161×161×97 would not fit in double.

## 6. Tests

**Integration** (`tests/integration`, CUDA build, `--memory-driver build_unittest/src/integration_hydro_memory.e`): **25/25 pass**.

The GPU tests now run every backend in the build:
- 0: the CPU code;
- 1: the device kernels on the CPU;
- 2: Metal;
- 3: CUDA double;
- 4: CUDA float (`gpu_backend 3`, `gpu_float 1`).

How the harness handles them:
- `GPU_KEYS`, `GPU_LIMITS` and `gpu()` in `music_integration_tests.py` set the backends and their limits. A backend not in the build ("built without") is skipped.
- **Limits:** double 1e-6, as backend 1; float 1e-4, as Metal.
- **Gubser tests:** both now check every GPU backend against the analytic solution, not just Metal.

CUDA results on this machine (worst case over the cases of each test):

| test | CUDA double | CUDA float |
|---|---|---|
| gpu_backend, conservation / E (EoS 0, 9, 14, 20) | ≤ 5.6e-8 | ≤ 3.6e-6 |
| gpu_backend_viscous, conservation / E | ≤ 1.4e-8 | ≤ 3.6e-6 |
| gpu_backend_viscous, L1 W, Π, μ_B (hot cells) | ≤ 1.0e-6 | ≤ 7.3e-6 |
| gpu_backend_sources (TATB, TATB viscous, strings) | ≤ 8.4e-12 | ≤ 8.3e-6 |
| gpu_backend_pi_chem, L1 Π_bulk, Π_chem | 0 | ≤ 4.1e-6 |
| Gubser ideal and viscous, error / CPU error − 1 | ≤ 7e-15 | ≤ 2.6e-4 |
| gpu_backend_freezeout, surface E, S, u_⊥ | 0 | ≤ 4.4e-6 |

**Unit tests** (`-Dunittest=ON`, g++ 15):
- cell, grid, minmod and reconst pass.
- **unittest_bugfix:** 41 of 42 cases pass. The failure:
  - **Where:** "Device kernels reproduce the viscous update of Advance", EoS 9, config 2: `CHECK(relE < 1e-5)` with relE = 1.32e-5.
  - **Not caused by this phase:** that check compares backend 1 **in float** with the CPU code, and unmodified fe606d7 gives the same 1.31559e-05.
  - **Cause:** phase 6 ran the unit tests with clang (bit for bit). g++ contracts float multiply-adds differently, so the float result moves slightly.
  - **Left as it is:** the tolerance or the compiler setting is the author's call.

**Not checked here:**
- **Metal:** this phase touched it only for the explicit `NULL` of `fetchCells` in `kMetalApi`.
- **The `MUSIC_CUDA_SYNTAX_ONLY` stubs** for the new CUDA calls (`cudaFuncGetAttributes`, `cudaEvent*`): there is no clang here. nvcc compiled the real code.

## 7. Merge with origin/gpu (7afc004, 8aa8b51, 6880b71)

The three upstream commits add three backend calls for Metal: `fetchUnwritten`, `conservationBegin` and `conservationEnd`. Their CUDA stubs assumed synchronous stages. With this phase's asynchronous stages, those stubs stayed correct, but `conservationBegin = NULL` made `beginStep` compute the rows on the host's call, which waits for the whole step and removes the overlap. The merge implements them on CUDA:

- **`conservationBegin` / `conservationEnd`:** the row kernel is queued behind the step on the default stream, into a buffer of its own, and `End` copies it back. All 98 row sets of the benchmark are computed this way.
- **`fetchUnwritten`, and `fetchCells` of slice 1 during a step:** they run on a non-blocking side stream and do not wait for the stages. Slice 1 is τ_n, which the stages in flight read and never write.
- **`outflowFaces`:** it now gathers only the cells of the boundary faces through `fetchCells`, instead of fetching the whole slice every step. `gatherOutflow` reads through an accessor, with the same arithmetic.
- **Struct order:** `fetchCells` is the last member of `MusicDeviceApi`, after upstream's three.

**Tests after the merge:**
- **Integration:** 25/25 pass.
- **Unit tests:** 44/45 pass; the failure is the pre-existing g++ float tolerance (section 6).
- **Upstream's three new unit tests** (the reads during a step, the copy during a step, the outflow faces) are written for backend 1. Pointed at CUDA (`gpu_backend = 3`) in a local run, all 256 assertions pass.
- **Outputs:** byte-identical to the pre-merge CUDA build (78 output files: ideal and viscous EoS 9, viscous EoS 20, float).

| float, 1.27M cells, 101 steps (mean of 2) | before the merge | merged |
|---|---|---|
| ideal, EoS 9 | 9.9 s | **8.9 s** |
| viscous, EoS 9 | 24.8 s | **23.9 s** |
| viscous, EoS 20 | 35.9 s | **35.7 s** |

## Files

| file | change |
|---|---|
| `src/device/cuda_backend.cu` | async stages and `finish`; `fetchCells` with the `cudaGather` kernel; per-kernel block sizes (`MUSIC_CUDA_BLOCK`); `MUSIC_CUDA_PROFILE`; stubs for the syntax check |
| `src/device/device_api.h` | `fetchCells` (optional, last member) |
| `src/device/device_hydro.{h,cpp}` | `Backend::fetchCells`; gathered partial downloads (`storeCell`, `downloadGathered`); the source-term fetch in both stages; enum warnings |
| `src/device/device_defs.h` | `flagOr` for nvcc |
| `src/device/metal_backend.m` | `fetchCells = NULL` written out |
| `CMakeLists.txt`, `src/CMakeLists.txt`, `.github/workflows/BuildTest.yml` | no Eigen |
| `src/advance.cpp` | closed-form eigenvalues of W^μ_ν instead of `Eigen::EigenSolver` |
| `tests/integration/music_integration_tests.py` | CUDA double and float in the six GPU tests |
