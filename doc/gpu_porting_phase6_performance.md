# GPU port, phase 6: performance (item 4)

Branch `gpu`, from f69f3c1 (phase 5, [gpu_porting_phase5_strings_lazy_sync.md](gpu_porting_phase5_strings_lazy_sync.md)). 2026-11-22. Machine: Apple M4 (10 GPU cores), as before.

## Summary

- **Four changes, none of which changes a result:** every output file of the three benchmarks is byte-identical to phase 5, and backend 1 stays bit for bit with the CPU code.
  1. **Kernels compiled for the run (Metal).** The switches of a run are compile-time constants. This removes the 25% that phase 4 added to the viscous kernels: they are now 3% faster than before phase 4.
  2. **Threadgroups of 1024 threads instead of 256 (Metal):** 11% faster viscous kernels and 7% faster ideal ones.
  3. **Partial host copies by x-ranges** instead of a test of every cell; the conservation kernel in small threadgroups.
  4. **The step runs while the host works.** Evolve takes what it reads from the device, starts the next step, and does the output, diagnostics and freeze-out on its copies while the GPU computes.
- **Metal against phase 5, 1.27M cells, 101 steps:**

  | case | phase 5 | phase 6 | speed-up over the CPU code |
  |---|---|---|---|
  | viscous, EoS 9 | 11.37 s | **7.97 s** (−30%) | 7.6× (CPU 60.9 s) |
  | viscous, EoS 20 | 16.6 s | **13.05 s** (−21%) | 8.8× (CPU 114.8 s) |
  | ideal, EoS 9 | 5.15 s | **3.47 s** (−33%) | |

  Both builds ran back to back, twice each, and the repeats agreed to 1%.
- **A Metal compiler problem, found by the tests and worked around:** with one loop bound compiled in, the causality-limiter run went wrong on Metal (section 1).

- **New tooling:**
  - `MUSIC_METAL_PROFILE=1` prints the GPU time of each kernel at the end of a run.
  - `MUSIC_METAL_TG=<n>` caps the threadgroup size, for tuning.
  - `MUSIC_METAL_NOSPEC=1` uses the generic kernels.
  - `MUSIC_METAL_SPEC_SKIP=name,...` leaves the named switches as run-time fields.
- **Tests:** all unit tests and all integration tests pass (section 6).

## 1. Kernels compiled for the run

- **Cause of the phase 4 slowdown:** phase 4 left the diffusion, Π_chem and causality code in the viscous kernels behind run-time flags. On the benchmark the kernels went from 7.36 s (phase 3) to 9.22 s (phase 4) with all three off.
- **The fix:**
  - The device code reads the switches of a run through `MD_RUN(P, name)` and `MD_EOS_TYPE(d)` (`device_defs.h`). These are the StepParams and EosDesc fields that stay fixed for a run:
    - the face scheme and stored Q;
    - shear, bulk, second order, vorticity and Du;
    - resummation, causality, Π_chem, diffusion and the 1+1D profile;
    - μ_B-dependent shear, the bulk relaxation type, QuestRevert and nWComp;
    - the EoS type.
  - `DeviceHydro::specKey` writes them as `name=value;...`. The new API call `specialize` hands them to the backend once, at start-up.
  - Metal compiles the library with `MUSIC_SPEC` and `MUSIC_SPEC_<name>` as preprocessor macros, so the code of the features that are off is removed. The library is now compiled lazily, generic if no key came.
- **The other backends are unchanged:** on the CPU and in CUDA `MD_RUN(P, name)` is `(P).name`, so their code is the same. `hasSource` and `stringSource` stay run-time flags, because they change from call to call.

**The causality case.** The integration test `gpu_backend_viscous` ("causality": EoS 0, `FlagResetCausality = 1`) failed with specialisation: W differed from the CPU code by 2.6e-2 (L1, hot cells) instead of 2.0e-6.
- **What is wrong:** the specialised build differs from the generic one by 5e-5 per step from the first step on. That is far above float rounding (CPU float: 7.8e-6 after 50 steps).
- **Bisection:**
  - `MUSIC_METAL_SPEC_SKIP` narrowed it to `nWComp`.
  - Making each of its three loops a run-time bound in turn narrowed it to one loop: the add of the W flux divergence to the right-hand side in `Gradients` (`device_visc.h`).
- **Not a bug in the source, as far as I can tell:** the loop is in bounds, its arrays are fully set, and the key at the first step equals the key at start-up.
- **It needs the limiter:** without it (all other viscous cases), specialised and generic builds give byte-identical output.
- **Conclusion and workaround:** this is treated as a Metal compiler problem (M4, macOS 26.7.3). That one loop reads `P.nWComp` at run time, with a comment. The causality case now agrees with the generic build exactly, and the loop costs 0.2 s of cellVisc on the benchmark.
- **Risk:** other compile-time switches could expose similar problems in configurations the tests do not cover. `MUSIC_METAL_NOSPEC=1` and `MUSIC_METAL_SPEC_SKIP` make that easy to check.

| Metal kernels (s) | phase 3 | phase 4/5 | specialised |
|---|---|---|---|
| viscous, EoS 9 | 7.36 | 9.23 | **7.13** |
| viscous, EoS 20 | 11.0 | 12.95 | **10.62** |

## 2. The threadgroup size

**Where the GPU time goes** (`MUSIC_METAL_PROFILE=1`, viscous EoS 9, 256-thread groups):

| kernel | GPU time | share |
|---|---|---|
| cellVisc (the cell update) | 3.88 s | 55% |
| face (KT fluxes; the conserved-to-primitive solve of both face states) | 2.09 s | 30% |
| wface (the W fluxes) | 1.01 s | 14% |
| pressure | 0.04 s | 1% |

For EoS 20, face rises to 4.4 s (42%) because of the 4D table lookups in the face solve.

**Threadgroup sizes from 32 to 1024** (viscous EoS 9):
- 1024 (the pipeline maximum on the M4) is the only one that helps: cellVisc goes from 3.88 to 3.22 s and face from 2.09 to 1.93 s.
- Tiles of the grid instead of x-rows (32×8×4, 16×16×4, 8×8×8, 64×4×4, …) were all slower, by 18 to 120%. The gain is therefore not cache reuse of neighbours: long runs in x, which keep the loads contiguous, matter most.
- The backend now takes the largest threadgroup the pipeline allows; `MUSIC_METAL_TG` caps it.

| Metal wall / kernels | 256 | 1024 |
|---|---|---|
| viscous, EoS 9 | 9.67 / 7.12 s | 8.45 / 6.30 s |
| viscous, EoS 20 | 14.22 / 10.62 s | 13.99 / 10.36 s |
| ideal, EoS 9 | 5.15 / 3.01 s | 4.47 / 2.81 s |

## 3. Host copies and the conservation kernel

- **Partial copies:** `syncShell` and `syncEtaWindow` gave `downloadT` a test per cell, run over all 1.27M cells. Now each (y, η) row gives at most two x-ranges to copy. The copy time per run falls from 0.52 to 0.40 s.
- **Conservation kernel:** it has one thread per (η, x) row, 7,889 threads, which is only 8 groups of 1024. It now runs in groups of 32. Its GPU time is 0.19 s per run (1.8 ms per call), and it is bound by the serial sum over y of each row. Splitting the rows would change the float summation order on Metal; for 2% of the run it is not done.

## 4. The step overlaps the host work

**Before:** a step took 79 ms, of which 62 ms was kernels. Most of the remaining 16 ms was host work while the GPU waited (`sample` of the process):
- `output_momentum_anisotropy_vs_tau` and the eccentricities (cexp, pow, atan2);
- the boundary outflow;
- the copies.

In runs with freeze-out the host share is larger: in the fireball test (216 steps), Metal ran for 4.0 s with 2.6 s of kernels.

**Now:**
- `DeviceHydro::beginStep` starts the step and returns while it runs. On Metal, the second stage is committed without a wait. `advanceStep` waits for it (the new API call `finish`), checks the flags and rotates the slices, as before.
- In each iteration, `Evolve::EvolveIt` first takes everything the host reads of the step: the full sync on output and freeze-out steps, the η windows, the conservation rows and the shell. Then it calls `beginStep`. The output, the diagnostics and the freeze-out run on the host copies while the GPU computes the next step, and `AdvanceRK` finishes it.

**Why this is safe:**
- During the step the device writes only its own buffers: the new slice, and the slice of τ_{n−1}, which has already been read.
- All reads of the device state (`syncHost`, `syncEtaWindow`, `syncShell`, `conservationRows`) fail while a step runs (`requireIdle`), so a misplaced read in Evolve cannot go unnoticed.
- The returns to re-run hydro leave a step running. Nothing reads the Fields after them, and the DeviceHydro destructor waits for the step.
- The MaxSpeed and NaN checks now come at the end of the step (in `advanceStep`), before anything reads the new state.

**Gain** (both runs of each case, outputs byte-identical):

| Metal wall | without overlap | with overlap |
|---|---|---|
| viscous, EoS 9 | 8.35 s | 7.74 / 8.12 s |
| viscous, EoS 20 | 13.99 s | 12.86 / 12.89 s |
| ideal, EoS 9 | 4.47 s | 3.45 / 3.48 s |

The log line `device: …` now reports the time the host spends launching and waiting for the device. That is no longer the kernel time, because the kernels overlap the host work; `MUSIC_METAL_PROFILE=1` gives the kernel time.

## 5. What is left

- **cellVisc** (half of the GPU time) is bound by compute and latency, not by bandwidth. Its memory traffic would take about 6 ms per stage at full bandwidth; it takes 16 ms.
- **face** (30%, 42% with EoS 20) is the iterative solve for the primitive variables of both face states, with EoS lookups.
- **wface** (15%) is memory-bound: it writes 11 values per face that cellVisc reads back. Fusing it into cellVisc would compute each face twice to save about 0.5 s.
- **Rounding:** each of these changes the order of operations, so the float results on Metal would change at rounding level. Keeping backend 1 bit for bit with the CPU code would need a separate device path.
- **Copies on the critical path:** with the overlap, the per-step copies take 0.82 s per run instead of 0.40 s. Without the overlap (profile mode) they take 0.40 s. They run while the GPU waits, right after the host waited for the step. The cause is not established; OpenMP threads waking from idle, or caches after the GPU's writes, are the likely candidates. Copying on a thread of its own, or fewer parallel regions, would be the next thing to try.
- **CUDA:** this phase does not change it except for the two new no-op API calls. `MD_RUN` could become template parameters there.

## 6. Tests

- **Unit tests (bgunit, clang):** all pass:
  - unittest_bugfix: 42 cases, 1,519,832 assertions;
  - cell, grid, minmod and reconst: all pass.
- **CUDA:** syntax check (host and device passes) passes.
- **Integration (bmetal, with the bgunit memory driver):** 25/25 pass. The causality case is back at 2.0e-6.
