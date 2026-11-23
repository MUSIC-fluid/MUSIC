# GPU against CPU: benchmark of fe606d7, and the next optimisation

Branch `gpu` at fe606d7 (phase 6, [gpu_porting_phase6_performance.md](gpu_porting_phase6_performance.md)). 2026-11-23.

- **Machine:** Apple M4 (4 performance and 6 efficiency CPU cores, 10 GPU cores, 16 GB).
- **CPU code:** `gpu_backend 0`, OpenMP with 10 threads.
- **GPU:** `gpu_backend 2`, Metal, float.
- **Build:** the same binary for both (clang, Release).

## Summary

- **Speed-up:** Metal is **7.5–11× faster** than the CPU code on production-size grids (0.4–1.3M cells). It is 2.7× on a small 41×41×17 grid, where launches and host work dominate.
- **Agreement:** Metal agrees with the CPU code at float rounding level in every case:
  - conservation series to about 1e-6 of E;
  - freeze-out surface: the same 370,636 cells, and energy and entropy fluxes to 7e-8.
- **Where the time goes:**
  - The GPU computes 84% of the wall time of a large viscous run, so the kernels decide the speed.
  - cellVisc takes half of the GPU time, and the KT face kernel 30–40%.
- **The next step: the floating-point mode of the Metal compiler.**
  - The kernels are compiled with `MTLMathModeSafe` (exact IEEE arithmetic). `MTLMathModeRelaxed` makes the face kernel 2.6× faster (3.2× with EoS 20), and the cell kernel 12–20% faster.
  - Measured wall times drop by **22–37%** (viscous EoS 9: 8.0 → 6.2 s; EoS 20: 13.0 → 8.8 s). The accuracy against the CPU code stays the same.
  - **It is not a drop-in change:** one integration case (an ideal run with a TATB source) aborts under relaxed math with an out-of-range KT signal speed. Some consistency checks of the solver rely on exact rounding and must be made robust first (section 4).
- **Second step, done (section 6):** the per-step boundary-shell copy cost 0.55 s of an 8 s run. The outflow now reads its cells from a compact face array gathered on the device side, and that array is carried over for the previous slice. The copies fall from 0.91 to 0.56 s and the wall time by 5%, with every output file bit-identical.
- **Hypotheses ruled out:**
  - the solver tolerances in float;
  - the number of fixed-point passes;
  - the EoS lookup cost;
  - OpenMP thread placement;
  - start-up and shader compilation.

## 1. Benchmark

- **Events:** lumpy 3D initial conditions (`Initial_profile 111`, the phase 2–6 benchmark), viscous with shear, bulk, second order and n_B. Δτ = 0.01 fm/c and 101 steps, unless stated otherwise.
- **Timing:** wall time of the whole run, start-up included.
- **Throughput:** cells × steps / wall time, in millions of cell-steps per second (Mc/s).

| case | grid (cells) | steps | CPU code | Metal | speed-up | Metal Mc/s | CPU Mc/s |
|---|---|---|---|---|---|---|---|
| viscous, EoS 9 | 41×41×17 (28.6k) | 101 | 1.54 s | 0.56 s | 2.7× | 5.1 | 1.9 |
| viscous, EoS 9 | 81×81×25 (164k) | 101 | 7.83 s | 1.27 s | 6.2× | 13.1 | 2.1 |
| viscous, EoS 9 | 121×121×37 (542k) | 101 | 25.27 s | 3.58 s | 7.1× | 15.3 | 2.2 |
| viscous, EoS 9 | 161×161×49 (1.27M) | 101 | 59.26 s | 7.88 s | 7.5× | 16.3 | 2.2 |
| ideal, EoS 9 | 161×161×49 | 101 | 39.67 s | 3.49 s | **11.4×** | 36.8 | 3.2 |
| viscous, EoS 20 (4D, n_B, n_Q, n_S) | 161×161×49 | 101 | 115.97 s | 13.02 s | 8.9× | 9.8 | 1.1 |
| viscous + diffusion, EoS 14 | 161×161×49 | 101 | 99.86 s | 12.65 s | 7.9× | 10.1 | 1.3 |
| full run with freeze-out, EoS 9 (e_fo = 0.26 GeV/fm³, Cornelius every 5 steps, Δτ = 0.02, until all cells froze out) | 101×101×41 (418k) | 381 | 99.45 s | 11.69 s | 8.5× | 13.6 | 1.6 |

- **Scaling:**
  - The CPU code has a flat throughput, about 2.2 Mc/s for viscous EoS 9.
  - Metal saturates at about 0.5M cells.
  - Below 0.2M cells, the per-step host work and the launch latency dominate. A 41×41×17 event runs at a third of Metal's peak throughput.
- **Agreement with the CPU code** (Metal float against CPU double):
  - **conservation series:** max |Δ| / E is 0.8–2.6e-6 in every case;
  - **freeze-out run:**

    | quantity | difference |
    |---|---|
    | surface cells | the same 370,636 |
    | energy flux | 8.3e-8 |
    | entropy flux | 6.7e-8 |
    | P_z | 6.8e-7 |
    | ⟨u_⊥⟩ | 5.1e-7 |

## 2. Where the time goes

**GPU time per kernel** (`MUSIC_METAL_PROFILE=1`, whole run):

| case | cellVisc / cell | face | wface | conservation | pressure | GPU total | wall |
|---|---|---|---|---|---|---|---|
| viscous, EoS 9 | 3.45 s (52%) | 1.91 s (29%) | 1.02 s (16%) | 0.19 s | 0.03 s | 6.60 s | 7.88 s |
| ideal, EoS 9 | 0.70 s (25%) | 2.01 s (72%) | — | 0.08 s | 0.02 s | 2.80 s | 3.49 s |
| viscous, EoS 20 | 5.16 s (49%) | 4.14 s (39%) | 1.04 s (10%) | 0.19 s | 0.07 s | 10.6 s | 13.02 s |
| diffusion, EoS 14 | 6.08 s (59%) | 2.71 s (26%) | 1.34 s (13%) | 0.19 s | 0.05 s | 10.4 s | 12.65 s |
| freeze-out run | 5.23 s (53%) | 3.21 s (33%) | 1.22 s (12%) | 0.19 s | 0.03 s | 9.9 s | 11.69 s |

**Host time on the critical path** (viscous EoS 9, 1.27M cells, 101 steps):
- **Start-up:** 0.3 s on Metal, including the specialised compile of the kernels (a 2-step run takes 0.39 s).
- **Copies to the host:** 0.82 s, split as follows:

  | copy | time | per step |
  |---|---|---|
  | the boundary shell (both slices, every step) | **0.55 s** | 5.7 ms |
  | the η window of the diagnostics | 0.15 s | |
  | 4 full syncs | 0.12 s | |

  - **The shell:** only about 80k cells per slice, but they are mostly 2 cells at the two ends of each x-row. In the SoA layout each of its 31 components then costs separate cache-line misses, so the shell moves more memory than a full copy of a slice.
  - **With the overlap** the copies run 2× slower than in profile mode (0.26 s for the shell there). That fits host cores clocked down while they wait for the GPU. The OpenMP thread count and `OMP_WAIT_POLICY` make no difference.

## 3. Hypotheses tested

Tested with temporary probes, all removed again: CPU-side counters in the device kernels (backend 1, float), and Metal compile variants.

| hypothesis | test | result |
|---|---|---|
| The solver tolerances are below float resolution, so the solves run to their caps | read `DeviceHydro` | no: in float they are raised to 1e-6 relative and 1e-8 absolute |
| The stage-2 fixed point (tolerance 1e-8) never converges in float | count of passes per cell, viscous EoS 9, 11.4M cell-stages | no: mean 0.59 passes, the same as in double. 57% of cells need 0 passes, 28% 1, 15% 2, 0.4% 3. In a 32-lane SIMD group someone usually needs 2, so the stage costs about 3 solves and W updates. |
| The face solves iterate a lot | count of EoS calls per solve: 138M face solves, 30M cell solves | mean 1.95 per face solve, with a tail to 25. 47% need none (dilute cells). Cell solves average 1.26. |
| The solve is what makes the face kernel expensive | face kernel with the solve removed (timing only) | face 1.91 → 0.73 s (EoS 9) and 4.17 → 1.46 s (EoS 20). The solve is 60–65% of the face kernel. |
| **The exact IEEE arithmetic of `MTLMathModeSafe`** | the same kernels compiled with relaxed or fast math | **face 1.92 → 0.74 s (EoS 9) and 4.17 → 1.31 s (EoS 20); cellVisc 3.43 → 3.04 s and 5.16 → 4.14 s.** With relaxed math, the solve costs almost nothing. |
| The cost is in the built-in functions (sqrt, cos, …) | relaxed arithmetic with precise functions | no: face 0.79 s, nearly the full gain. It is the arithmetic mode. |
| The cost is the correctly rounded division in the solver and the EoS lookup | `fast::divide` in the 15 divisions of `velocityF`, the secant step and the 2D table lookup, in safe mode | only 9% of the face kernel (1.91 → 1.74 s). The gain of relaxed mode is spread over the kernel (FMA contraction, reassociation and the remaining divisions). |
| OpenMP placement slows the copies | `OMP_NUM_THREADS` 4 or 10, `OMP_WAIT_POLICY=active` | no effect (0.80–0.85 s) |

**Relaxed math, end to end** (wall, two runs each, the same binary):

| case | safe | relaxed | change |
|---|---|---|---|
| viscous, EoS 9 | 8.29 / 8.02 s | 6.44 / 6.21 s | −22% |
| viscous, EoS 20 | 13.38 / 13.01 s | 8.86 / 8.79 s | **−33%** |
| ideal, EoS 9 | 3.85 / 3.47 s | 2.50 / 2.31 s | −34% |
| diffusion, EoS 14 | 13.20 / 12.20 s | 9.36 / 8.63 s | −29% |

**Accuracy against the CPU code (double):** unchanged.
- **The 1.27M-cell runs:**

  | case | safe: E / N_B | relaxed: E / N_B |
  |---|---|---|
  | viscous, EoS 9 | 1.0e-6 / 8.0e-8 | 1.2e-6 / 1.1e-7 |
  | viscous, EoS 20 | 1.1e-6 / 8.5e-8 | 1.3e-6 / 1.1e-7 |

- **The integration suite under relaxed math:** 24 of 25 pass, every Metal metric at the safe-mode level. For example, the causality case gives W 2.1e-6 (safe 2.0e-6), and EoS 20 gives μ_B 2.2e-6 (safe 1.9e-6).
- **The failure:** `gpu_backend_sources` (TATB source, ideal, EoS 0) aborts at step 31 with "MaxSpeed: the signal speed of a KT face is out of range". It passes in safe mode.
  - Compiling MaxSpeed itself in safe mode (`#pragma METAL fp math_mode(safe)`) does not help. The face state already arrives inconsistent.
  - The likely cause is the normalisation of u after the solve: `|u0² − u² − 1| > absErr` and the rescaling. Relaxed math reassociates it, so u_x can exceed u^τ by a few ulp in dilute, fast cells.

## 4. The next steps, in order of gain

1. **Relaxed math on Metal, with robust consistency checks.** Expected −22% to −37% wall on the benchmarks, from the measurements above.
   - **What it needs:**
     - **Find the failing operation.** Run `gpu_backend_sources` under relaxed math with backend-1 style flags at the face, recording which check fires and the state that caused it.
     - **Make the state physical after each solve,** whatever the rounding: a final clamp or renormalisation that enforces u^τ ≥ |u^i| and f ≤ 1. Alternatively, compile those few lines in safe mode (the pragma compiles; the scope has to be the normalisation, not MaxSpeed).
     - **Keep exact rounding where the CPU code compares with exact thresholds:** `isfinite` (relaxed mode keeps INF/NaN semantics; fast mode would not) and the MaxSpeed range checks.
   - **Validation:**
     - the whole integration suite on Metal;
     - the 1.27M-cell runs against the CPU code;
     - a run per EoS type;
     - a deliberately broken input, to confirm that NaN detection still fires.
   - **What it does not change:** backend 1 and CUDA. The math mode is a Metal compile option, so the CPU code and backend 1 stay bit-exact.
   - **Risk:** the compiler problem with `nWComp` in phase 6 showed that this compiler can change results in ways the tests catch only when a case covers them. Run the integration suite after every change to the math mode or specialisation.
2. **The boundary-shell copy: done (section 6).** The copy went from 0.55 to 0.18 s per run, with bit-identical output. Packing the shell on the GPU would not have paid off: the scattered writes into the host Fields cost more than the device reads.
3. **Divergence in the cell kernel's second stage.** 15% of cells need 2 fixed-point passes, so most SIMD groups run 3 solves and W updates.
   - **Measure first:** time stage 2 against stage 1 in the profile. With relaxed math the solve gets cheaper, so this matters less.
   - **Possible fix:** a compacted second pass over only the cells that did not converge after one pass. It changes no arithmetic, so it stays bit-exact.
4. **wface into cellVisc** (about 0.5 s, from phase 6): it computes each face twice to save the write and read of 11 values per face. It changes the order of operations at rounding level.
5. **Small grids:** below 0.2M cells, the per-step host diagnostics and the two command buffers per step dominate. To do: batch several steps per command buffer when no output is due, and make the per-step diagnostics optional.

**Not worth doing:** the start-up (0.3 s), the OpenMP settings, and the EoS lookup code (9% at most, as a micro-optimisation of the face kernel).

## 5. Reproduction

- **Driver:** `run_bench.py` (kept in the session scratchpad, not in the repository) writes the inputs from the `visc9_0` benchmark input, with the integration suite's lumpy T_A/T_B. It runs both backends of each case and stores wall times, the `device:` lines and the conservation differences. Cases: `visc9_{41,81,121,161}`, `ideal9_161`, `visc20_161`, `diff14_161`, `fullFO9_101`.
- **Kernel times:** `MUSIC_METAL_PROFILE=1`.
- **The probes** (math mode, solve counters, copy split) were temporary and are not in the tree. The source is unchanged at fe606d7.

## 6. Done: the boundary outflow without the shell copy

**What the shell was for:** in each step, Evolve copied the outermost cells of both slices (`syncShell(curr, prev, 1, 1, 1)`) into the host Fields. Only `Cell_info::boundary_outflow_rate` read them.

**Why packing the shell would not help** (a stand-alone copy benchmark, 82k shell cells × 31 components, caches flushed):

| copy | time per slice |
|---|---|
| before: scattered reads from the device and scattered writes into the Fields | 2.1–2.6 ms |
| packed on the GPU, then scattered writes into the Fields | 1.75 ms |
| packed reads only | 0.3–0.5 ms |

- **The Fields are per component too (SoA),** so writing two cells per x-row into 31 arrays costs as much as reading them.
- **So the fix removes the writes,** not the scattered reads.

**The change** (no new kernel; the backend API and CUDA are unchanged):
- **`boundary_outflow_rate` is one loop over the faces in a fixed order.** `Cell_info::outflowFace` maps a face number to its cell, direction and side, in the order of the former nested loops. The arithmetic is the same statement by statement, and it reads its 24 values per cell (e, n_B, n_Q, n_S, Π, Π_chem, u^μ, W^0…W^13) either from the Fields or from a face array (a new overload).
- **`DeviceHydro::outflowFaces(arena, which)`** gathers that array from the device slice in parallel, converting as `downloadT` would. It runs in Evolve's gather before the step starts. `check_conservation_law` then reads it while the GPU computes, through two new optional arguments.
- **The previous slice is not gathered again.** It is the device buffer that was current in the step before, and its host Fields rotate with it. So `advanceStep` swaps the two arrays; an upload clears them, and so does a step whose current faces were not taken.
- **Ideal runs:** the device has no Π or W, so the array takes them from the host Fields, as the old path did. The ideal update never changes them, so they are taken only after an upload.
- **Unchanged:** `syncShell` stays for beastMode 3. Output steps (full sync) and boost-invariant runs use the Fields as before.

**Validation:**
- **Byte-identical to fe606d7, every output file:**
  - the CPU code (`gpu_backend 0`);
  - backend 1 in double and float, viscous and ideal;
  - Metal viscous, ideal, and the freeze-out run, surfaces included.
- **New unit test "Device outflow faces reproduce boundary_outflow_rate"** (backend 1; ideal and viscous × EoS 0 and 9 × double and float): the rates from the device face values, current slice and carried-over previous one, equal those from the synced Fields bit for bit.
  - **The test can fail:** with a perturbed carried-over array, 16 of its 152 assertions fail.
- **All unit tests pass** (43 cases, 1,519,984 assertions), and the CUDA syntax check passes.
- **Integration suite (Metal):** 25/25 pass, with the metrics of fe606d7 (the causality case still gives W 2.015e-6).

**Timing:** interleaved with the fe606d7 binary, three runs each. The machine was more loaded than for section 1 (load average 4–7), so the absolute times are higher.

| Metal | copies, fe606d7 | copies, now | wall, fe606d7 | wall, now |
|---|---|---|---|---|
| viscous, EoS 9, 1.27M cells | 0.87–0.94 s | 0.55–0.56 s | 8.67–8.87 s | 8.26–8.35 s (−5%) |
| ideal, EoS 9, 1.27M cells | 0.14–0.16 s | 0.11 s | 3.68–3.90 s | 3.64–4.04 s |
| freeze-out run, 418k cells | 1.12–1.22 s | 1.14–1.26 s | 12.4–13.0 s | 12.3–12.8 s |

- **The gather:** 1.8 ms per step for the current slice. The previous slice is gathered only after the 4 full syncs.
- **The freeze-out run gains nothing:** its copies are the 77 full syncs of the freeze-out steps, not the shell. A partial sync for Cornelius would be the next host-copy target there.
