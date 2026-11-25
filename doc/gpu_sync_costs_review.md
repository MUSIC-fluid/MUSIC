# Remaining host costs of GPU runs (branch `gpu`, after 5456c02)

## Summary

- **Copies from the device are no longer the main cost.** In the (3+1)D strings event on Metal they take 2.3 s of a 47 s run, and most of that overlaps the GPU.
- **The host work is.** About 13 s of the run is host work that runs while the GPU is idle, or that takes longer than the GPU step it overlaps:
  1. **The beastMode regrid: 5.1–6.4 s for one coarsening, 11–13 % of the run.** This is the largest item, and the GPU is idle the whole time. The remap runs on one thread, because Apple clang in the Metal build has no OpenMP. One of the five Fields it remaps, `fpNext`, is overwritten by the next step anyway.
  2. **Freeze-out steps are bound by the host.** On those steps the wait for the device is 0. The host spends about 100 ms per step (39 ms copying the current slice, 43 ms in the finder, 20 ms copying Fields) against about 66 ms for the GPU step: about 1.7 s exposed over 48 steps.
     - Half of the Fields copying, the copy of `fpPrev` into `freezeoutFieldPrev`, copies host data that is stale on the device path and is never read.
  3. **The momentum-anisotropy output: 1.2 s for 4 calls,** host computation after a full sync.
- **Everything else is small:**
  - Full syncs: 0.36 s for 6.
  - Grid-edge checks: 0.23 s for 566.
  - The Δτ extrapolation: 0.23 s.
  - Coarsening checks: 0.02 s.
  - Ghost-cell updates: 0.1 s.
- **Not measured:**
  - CUDA, where every copy also crosses PCIe, so the freeze-out copy of the current slice costs more.
  - The per-step diagnostics of beastMode 0: conservation rows, η windows and outflow faces.

## How this was measured

- **The run:** `strings_event_0.dat` with `example_inputfiles/3D_dynamical/music_input_mode_2`, beastMode 3 with `beastModeCoarsenTolerance` 0.01, `gpu_backend = 2` (Metal, Apple M4, float, viscous). The grid starts at 120 × 120 × 107 and is coarsened once, at τ = 5.70 fm/c. The run has 566 steps and freezes out at τ = 6.0 fm/c.
- **The instrumentation:** a scratch copy of `src/evolve.cpp` with wall-clock timers around each part of the main loop. The wait for the device (`AdvanceRK`) is timed separately on freeze-out steps, check steps and other steps.
- **Results unchanged:** the observables of the instrumented build are bit-identical to those of the committed build.

| part of the loop | calls | time [s] | GPU idle? |
|---|---|---|---|
| wait for the device, ordinary steps | 494 | 32.9 | — (the GPU steps) |
| wait for the device, coarsening-check steps | 24 | 1.76 | — |
| wait for the device, freeze-out steps | 48 | **0.0003** | host-bound |
| regrid (`regridTransverse`, one coarsening) | 1 | **5.06–6.39** | yes |
| freeze-out finder + Fields copies | 48 | 3.02 | partly |
| of which `store_previous_step_for_freezeout` (Fields copies) | 103 | 0.85–1.06 | partly |
| copy of the current slice during the step (freeze-out) | 48 | 1.87 | partly |
| momentum-anisotropy output (`itSinceFreezeStart` 0, 10, 30, 50) | 4 | 1.17 | partly |
| full syncs (those 4 output steps) | 4 | 0.24 | partly |
| grid-edge checks (device reduction) | 566 | 0.23 | yes |
| Δτ extrapolation of `fpPrev` | 2 | 0.23 | yes |
| syncs before the regrid | 2 | 0.11 | yes |
| ghost-cell updates after the grid changes | 3 | 0.10 | yes |
| coarsening checks (device reduction) | 24 | 0.02 | yes |
| whole loop | | 47.3 | |

The device's own counters give the same picture: 2.28 s copying in total (7 full syncs and 48 overlapped copies of the current slice), with no partial syncs left.

## Done: skip the fpNext remap and the unused fpPrev copy

- **The regrid remaps only the Fields that are read afterwards:** `fpPrev`, `fpCurr` and `freezeoutFieldCurr`, plus `freezeoutFieldPrev` when vorticity output is on. `fpNext` is only resized, and so is `freezeoutFieldPrev` without vorticity output. This covers both the coarsening and the grid enlargement (`src/evolve.cpp`, the `regrid` lambda in `EvolveIt`).
  - **Why `fpNext` can be resized:** the CPU update writes every cell of it before reading it. `FirstRKStepT` reconstructs e, n and u into it, and `FirstRKStepW` then reads those and writes W and Π. `resizeFields` leaves a fluid at rest with W = Π = 0, which is what the ideal update keeps. The device path uploads only `fpCurr` and `fpPrev`.
- **`fpPrev` is copied into `freezeoutFieldPrev` only when the finder reads it,** that is with `output_vorticity = 1` in (3+1)D, at the start of freeze-out and on each freeze-out step. Otherwise the copy, and the Δτ extrapolation of the stored slice, are skipped, on the CPU and the device path alike. On the device path the copy was of stale data, since `fpPrev` is not synced on those steps.
- **Results are bit-identical.** I compared the surface, the evolution and the FO files of the previous and the new executable (md5) on the 2D fireball of the `beastmode3` integration test:
  - CPU and Metal;
  - ideal and viscous (shear + bulk);
  - grid enlargement only (the default tolerance), and enlargement + coarsening (tolerance 0).

  The (3+1)D strings event on Metal also gives identical observables. All 48 unit tests pass.
- **Saved, in the (3+1)D strings event on Metal:**
  - Each remapped Fields costs 1.3–1.4 s on the single host thread, so two fewer save about 2.7 s per coarsening or enlargement. The regrid's measured total is noisy, though: 5.1–6.4 s before, 4.1–5.2 s after, and identical work varied by 1.3 s.
  - The Fields copies on freeze-out steps fell from 0.97 s to 0.42 s, and the Δτ extrapolation from 0.23 s to 0.11 s.
  - The loop took 47.3 s before and 46.5–46.8 s after, in two runs with only one coarsening.
- **Not run:** a (3+1)D run with vorticity output. Its path is unchanged: it still remaps, copies and extrapolates `freezeoutFieldPrev`.

## Done: OpenMP on the host of GPU runs

**Before,** CMake skipped OpenMP for Apple clang, so every Metal build ran its host work on one thread.

**The build change:**
- With Apple clang, CMake now looks for libomp: `LIBOMP_ROOT`, or Homebrew's `/opt/homebrew/opt/libomp` or `/usr/local/opt/libomp`.
- When it finds it, it compiles with `-Xpreprocessor -fopenmp` and links `libomp`, into the library, the executable and the unit tests.
- Other compilers use `find_package(OpenMP)` as before.
- `-Dopenmp=OFF` builds without OpenMP. Without libomp, the build goes on single-threaded and says so.

**The remap loop is now parallel.** `coarseGrainAndEnlargeGrid` had no OpenMP pragma. Its new cells are independent: each reads the old copy and writes only itself, and the shared `Reconst` is already used this way in `Advance::AdvanceIt`. So the loop over (η, y) is now `omp parallel for`, with a reduction over the count of cells remapped with the ideal part only. This also speeds up CPU runs.

**Results:**
- The 2D fireball (beastMode 3 with enlargement and coarsening; CPU and Metal; ideal and viscous) gives bit-identical outputs (md5) to the build without OpenMP.
- The (3+1)D strings event on Metal differs only in the last digit of S, ⟨u_⊥⟩ and the v2 proxy. That comes from the order in which the parallel freeze-out finder writes the surface; the parallel remap alone changes nothing.
- All 48 unit tests pass with OpenMP.
- With `-Dopenmp=OFF` the build has no libomp.

**3D strings event on Metal, 10 threads:**

| part | 1 thread | OpenMP | OpenMP + parallel remap |
|---|---|---|---|
| regrid (`regridTransverse`, 3 Fields) | 4.0–5.2 s | 4.0 s | **0.77 s** |
| freeze-out finder + Fields copies (48 steps) | 2.7 s | 1.36 s | 1.34–1.44 s |
| copy of the current slice on freeze-out steps (48) | 1.86 s | 0.65 s | 0.67–0.71 s |
| momentum-anisotropy output (4) | 1.14 s | 0.28 s | 0.28 s |
| wait for the device on freeze-out steps (48) | 0.0003 s | 1.5 s | 1.5–1.6 s |
| whole loop | 46.5–46.8 s | 42.4–42.6 s | 38.8–40.7 s |
| run wall time | 47.7 s | 43.5–43.8 s | 39.9–42.5 s |

- **Freeze-out steps are no longer bound by the host.** The wait for the device on those steps is 1.5 s again, which is GPU time the host work now fits inside.
- **The run is 12–16 % faster than with one host thread,** and the regrid 5× faster.
- **Against beastMode 1** on the same GPU (65.8–70.3 s), beastMode 3 is now 1.55–1.76× faster, close to the 1.74× fewer cell updates.
- **What is left of the host work** is the grid-edge check (0.2 s), the Fields copies of the freeze-out steps (0.8–0.9 s, serial `std::copy`), and the syncs around the regrid.

## Details and options

### 1. The regrid (largest item)

**Why it is slow.** `regridTransverse` remaps five Fields:
- `fpPrev`, `fpCurr` and `fpNext`;
- the two stored freeze-out slices.

For each new cell, `coarseGrainAndEnlargeGrid`:
- sums the conserved densities of up to 9 parent cells, with one EoS pressure call per parent;
- reconstructs e, n and u with `Reconst`, a root solve that calls the EoS.

At 60 × 60 × 107 new cells per Field, on one thread, that takes 5–6 s. The GPU waits all that time. A grid enlargement goes through the same function and costs about as much.

**Options, cheapest first:**
- **Don't remap `fpNext`; just resize it.** (Done; see above.)
  - On the device path, the host `fpNext` is never uploaded: only `fpCurr` and `fpPrev` are, after `invalidate()`.
  - On the CPU path, the next RK step writes every cell of `fpNext`. Before relying on this, check that no stage reads `fpNext` before writing it, for example as the starting guess of the velocity solve (`Advance` passes it near `advance.cpp:206`).
  - This saves a fifth, about 1.1–1.3 s.
- **(Done; see above.) Build the Metal executable with OpenMP** (Homebrew libomp with Apple clang, or g++ for the host files). The remap loop is independent per new cell. This also speeds up the freeze-out finder and the anisotropy output, which are single-threaded in the Metal build too.
- **Remap on the device.** It needs a kernel for the filter, plus the device versions of the reconstruction, which already exist for the RK update. It also saves the two syncs and the upload around the regrid. This is the most work.

### 2. Freeze-out steps

**What a freeze-out step does on the host,** every `facTau` steps:
1. copies all of `fpCurr` from the device while the step runs (39 ms);
2. runs the Cornelius finder (about 43 ms);
3. copies `fpPrev` and `fpCurr` into the freeze-out Fields (two full Fields copies, about 20 ms).

The GPU step takes about 66 ms, so roughly 36 ms per freeze-out step is exposed.

**Options:**
- **(Done; see above.) Skip the copy of `fpPrev` into `freezeoutFieldPrev` on the device path** when vorticity output is off. `fpPrev` is not synced on those steps (the code says so where it chooses `syncCurrDuringStep`), so the copy moves stale data that is never read. This saves about half of the 1 s of Fields copying.
- **Swap `freezeoutFieldCurr` with the host `fpCurr` instead of copying it,** on the device path. `fpCurr`'s host Fields are filled again by the next sync before they are read. Check this against `DeviceHydro::hostCurrent` and the pointer checks in `DeviceHydro` before doing it.
- **Copy less.** The finder only needs e (and n_B, μ) to place the surface, and the other components only for the cells next to it. The device could send e first, and then only the cells of the cubes the surface crosses. That replaces most of the 39 ms copy, and most of the PCIe traffic on CUDA.

### 3. The momentum-anisotropy output

**What it costs.** `output_momentum_anisotropy_vs_etas` runs at `itSinceFreezeStart` = 0, 10, 30 and 50, each time after a full sync. It takes about 290 ms of single-threaded host work per call, against a 66 ms GPU step.

**Options:**
- OpenMP (as above);
- a row reduction on the device, as for the conservation sums;
- an input switch to turn it off in production runs.

### 4. Small items

- **The grid-edge check, every step:** 0.4 ms of kernel launch and wait. It could be queued behind the step, as the conservation rows of beastMode 0 are (`conservationAhead`), so that the host never waits for it. That saves at most 0.2 s here.
- **The Δτ extrapolation and the syncs around the regrid:** 0.3 s. They would go away with a device remap.

## What was not covered

- **CUDA.** There was no CUDA machine. There a full sync is a PCIe transfer: about 0.4 GB per slice in double for this grid. So the freeze-out copy of the current slice and the full syncs weigh more than on Metal, while the two device reductions (coarsening and grid edge) are cheap on both.
- **beastMode 0.** It reads the conservation rows, two η windows and the outflow faces every step. These were optimised in earlier phases (`conservationAhead`, `fetchCells`) and were not measured again here.
