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
- **The floating-point mode of the Metal compiler: tried, not adopted (section 7).**
  - The kernels are compiled with `MTLMathModeSafe` (exact IEEE arithmetic). `MTLMathModeRelaxed` makes the face kernel 2.6× faster (3.2× with EoS 20), the cell kernel 12–20% faster, and the runs 22–37% faster, at the same accuracy against the CPU code on most metrics.
  - **It is not reliable for MUSIC.** Relaxed rounding occasionally flips a cell-level decision of the primitive solve (solve, cap, regulate or revert), and each flip is a discrete, non-conservative event: +2 GeV in one step of a long freeze-out run, +0.1 GeV in the strings test.
  - **The cause:** two concrete holes were found and closed (a NaN test folded away, and the bracketing solver failing). But which case breaks depends on the rounding pattern, so each variant failed a different test. Exact float never did.
  - **The code stays in safe mode.** The attempt is kept as `doc/gpu_relaxed_math_attempt.patch`.
- **Second step, done (section 6):** the per-step boundary-shell copy cost 0.55 s of an 8 s run. The outflow now reads its cells from a compact face array gathered on the device side, and that array is carried over for the previous slice. The copies fall from 0.91 to 0.56 s and the wall time by 5%, with every output file bit-identical.
- **Items 3 and 6, done (section 8):** −4% to −11% of the wall time, every output byte-identical except the freeze-out run's conservation file (5e-8 of E).
  - **The cell kernel per RK stage.** The stage-2 fixed-point loop costs registers in the stage-1 kernel too, where it never runs. Compiling one pipeline per stage makes the stage-1 cellVisc 15% faster (EoS 9) and 26% faster (EoS 20).
  - **The stage-2 compaction is not worth building.** Divergence doubles the loop's cost, but the whole loop is 0.37 s of a 7.3 s run, and a compaction would gain under 1%.
  - **Freeze-out steps copy the current slice while the next step runs,** and skip the previous slice, which only the vorticity output reads.
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

1. **Relaxed math on Metal: tried, not adopted (section 7).** It would gain 22–37% of the wall time, but it makes the cell-level decisions of the primitive solve depend on rounding. The integration suite fails in a different case for each variant tried.
   - **A future attempt would need:** a primitive solve whose decisions (cap, regulate, revert, the u⁰ sanity check) are insensitive to rounding by construction, for example with margins around each threshold, applied identically on the host. That changes the CPU code's results and is a physics decision.
2. **The boundary-shell copy: done (section 6).** The copy went from 0.55 to 0.18 s per run, with bit-identical output. Packing the shell on the GPU would not have paid off: the scattered writes into the host Fields cost more than the device reads.
3. **The cell kernel's second stage: done (section 8),** as one pipeline per RK stage. The compaction of the fixed-point passes was measured and not built.
4. **wface into cellVisc** (about 0.5 s, from phase 6): it computes each face twice to save the write and read of 11 values per face. It changes the order of operations at rounding level.
5. **Small grids:** below 0.2M cells, the per-step host diagnostics and the two command buffers per step dominate. To do: batch several steps per command buffer when no output is due, and make the per-step diagnostics optional.

6. **The freeze-out steps' full syncs: done (section 8).** The current slice is copied while the next step runs.

**Not worth doing:** the start-up (0.3 s), the OpenMP settings, and the EoS lookup code (9% at most, as a micro-optimisation of the face kernel).

## 5. Reproduction

- **Driver:** `run_bench.py` (kept in the session scratchpad, not in the repository) writes the inputs from the `visc9_0` benchmark input, with the integration suite's lumpy T_A/T_B. It runs both backends of each case and stores wall times, the `device:` lines and the conservation differences. Cases: `visc9_{41,81,121,161}`, `ideal9_161`, `visc20_161`, `diff14_161`, `fullFO9_101`.
- **Kernel times:** `MUSIC_METAL_PROFILE=1`.
- **The probes** (math mode, solve counters, pass histogram, loop caps, copy split, freeze-out timers) were temporary and are not in the tree.

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

## 7. Tried: the relaxed floating-point mode

**First failure: a NaN folded away** (`gpu_backend_sources`, TATB source, ideal, EoS 0, aborted at step 31 with "MaxSpeed: the signal speed of a KT face is out of range").
- **Capture:** a temporary probe wrote the first failing face into a debug buffer.
- **The face:** the reconstructed conserved state is nearly lightlike: τT^{ττ} = 0.54253 and |τT^{τi}| = 0.54253, at the edge of the source region between e = 1.25 and e = 0.11.
- **Exact rounding:** replayed on the CPU with the device solver, the secant solve fails, the bracketing solve fails, and the face reverts to the cell state (flag 2), in float and in double. A scan of ±4e-6 around this state gives no NaN in exact arithmetic.
- **Relaxed rounding:**
  - The secant solve converges to v = 0.9999938 (u⁰ = 284), and `solveU0Hybrid` then returns success with u⁰ = NaN. Its NaN test `root != root` had been folded to false by the compiler.
  - The NaN then passes the sanity check |u⁰ − u⁰_guess| / u⁰_guess > 100, because comparisons with NaN are false.
- **Fix:** NaN and finiteness tests on the bit pattern of the float (integer operations, which no math mode folds). With it the case runs to the end. An injected NaN in e, u^x or τT^{ττ} is still caught at the first step, exactly as in safe mode and on backend 1.

**Second failure: energy conservation** (the 101×101×41 run to freeze-out, 381 steps).

| Metal build | max \|E/E0 − 1\| | steps | steps where the bracketing solver failed | wall |
|---|---|---|---|---|
| CPU code (double) | 7.4e-7 | 381 | — | 99.5 s |
| backend 1, float (exact, on the CPU) | 4.0e-6 | 381 | — | 74.6 s |
| safe (default) | 1.1e-6 | 381 | 0 | 11.7–12.6 s |
| relaxed | **2.7e-3** (+2.04 GeV in the step at τ = 7.33) | 386 | 25 | 8.5 s |
| relaxed, solver exact (`#pragma METAL fp math_mode(safe)` in `Solver`) | 1.9e-5 | 381 | 0 | 10.0 s |
| relaxed, solver and EoS exact | 1.7e-4 | 381 | 0 | 10.5 s |
| relaxed, EoS exact | 3.9e-4 | 381 | 0 | 8.9 s |

- **Every relaxed variant violates conservation in discrete steps,** not by drift. The jumps range from 1 MeV to 2 GeV, at the times when a cell takes a different decision.
- **The size of the violation is not monotonic** in how much is kept exact. It depends on which threshold the rounding happens to flip.
- **Exact float, on Metal or on the CPU, stays at the 1e-6 level.**
- **The freeze-out surface** of the fully relaxed run still agrees with the CPU code: E to 1e-6 and S to 7e-8, with 16 more cells out of 370,636.

**Integration suite, Metal** (the other 52 or 53 metrics are at their safe-mode values, median ratio 1.0–1.05):

| build | result |
|---|---|
| safe (default) | 25/25 |
| relaxed, with the NaN fix | 25/25 (worst metric 2.2× its safe value) |
| relaxed, solver exact | 24/25: `gpu_backend_sources`, strings viscous, conservation 1.05e-3 (limit 1e-4), from a jump of +0.106 GeV in the step at τ = 1.36 |

**Speed of the variants** (wall, 1.27M cells, 101 steps):

| case | safe | relaxed | relaxed, solver exact |
|---|---|---|---|
| viscous, EoS 9 | 7.86 s | 6.07 s | 6.61 s |
| viscous, EoS 20 | 12.91 s | 8.57 s | 9.86 s |
| ideal, EoS 9 | 3.57 s | 2.37 s | 3.10 s |

**Decision:** the default stays `MTLMathModeSafe`, and the source is unchanged. Making relaxed math safe would need a primitive solve whose decisions are robust to rounding by design (section 4, item 1).
- **The patch of the attempt:** `doc/gpu_relaxed_math_attempt.patch`. It contains the bit-level NaN tests, the `MD_EXACT_FP` pragma on the solver, and the `gpu_relaxed_math` input parameter (Metal only, default 0).

## 8. Done: the cell kernel per RK stage, and the freeze-out copies (items 3 and 6)

Both changes keep every output byte-identical, except `global_conservation_laws.dat` of runs with freeze-out (below). Backend 1 stays bit for bit with the CPU code.

**Wall time, 7afc004 against this change** (Metal; two interleaved runs of each build; 101 steps on 1.27M cells, freeze-out 381 steps on 418k cells):

| case | 7afc004 | now | change |
|---|---|---|---|
| viscous, EoS 9 | 8.16 / 7.54 s | 7.29 / 7.28 s | −7% |
| viscous, EoS 20 | 12.95 / 12.53 s | 11.96 / 11.98 s | −6% |
| viscous, EoS 14, diffusion | 11.56 / 11.49 s | 10.91 / 10.82 s | −6% |
| ideal, EoS 9 | 3.60 / 3.45 s | 3.47 / 3.47 s | 0 (no cellVisc) |
| viscous, EoS 9, freeze-out to the end | 11.33 / 11.33 s | 10.60 / 10.46 s | −7% |

### 8.1 The stage-2 fixed-point loop (item 3)

**Measured first** (`MUSIC_METAL_PROFILE=1`, the kernel time split by stage):

| cellVisc GPU time | stage 1 | stage 2 |
|---|---|---|
| EoS 9 | 1.47 s | 1.94 s |
| EoS 20 | 2.23 s | 2.94 s |

**Passes per cell** (a counter on backend 1 in float, over all cells and steps of the EoS 9 run): 52% of the cells need no extra pass, 23% one, 21% two, 3.2% three and 0.5% more. On average a cell needs 0.90 extra passes, but its 32-lane SIMD group runs 1.92 passes per lane, so divergence doubles the cost of the loop.

**What compaction could gain:**
- With the loop capped at zero passes (wrong results, timing only), stage 2 takes 1.57 s instead of 1.94 s, so the whole loop costs 0.37 s.
  - Caps of one to three passes are slower than the full loop, because every capped cell also runs the fallback solve from its best iterate.
- With no divergence at all, about 0.2 s would remain to gain.
- **A compacted second pass gives back most of that 0.2 s:**
  - deferring the cells after two passes means recomputing their stage for 3.7% of the cells;
  - deferring them after one pass means recomputing it for 25% of the cells, or storing about 80 values per deferred cell;
  - either way the estimated net gain is under 1%.
- **Not built.**

**What the measurement showed instead:**
- With EoS 20, the stage-1 kernel was 0.2 s (9%) faster when the loop's code was absent, although stage 1 never runs it. This held in two repeats; EoS 9 showed no difference. The loop's code costs registers, or code layout, in both stages.
- **The fix:** `kernelCellVisc<Real, kStage>` (`device_visc.h`).
  - kStage = −1 reads the stage from `P.rkFlag`, as before; the host and CUDA use it.
  - Metal compiles `musicCellVisc0` and `musicCellVisc1` (`metal_kernels.metal`) and dispatches `psCellVisc[rk]`.
  - kStage gates only the fixed-point branch; every arithmetic use of `P.rkFlag` stays a run-time value. So the arithmetic is the same, and the outputs are byte-identical.
- **Result:** the stage-1 cellVisc falls from 1.47 to 1.25 s (EoS 9), and from 2.23 to 1.64 s (EoS 20). The larger fall than the probe's 0.2 s suggests that the stage-1 kernel loses more than the loop alone. Stage 2 is unchanged.
- The profile now lists `cellVisc_rk0` and `cellVisc_rk1` separately.

### 8.2 The freeze-out steps (item 6)

**Host work per freeze-out step** (76 steps of the freeze-out run, a temporary timer):

| | per step | total | |
|---|---|---|---|
| full sync of both slices | 6.9 ms | 0.52 s | the GPU waits |
| Cornelius search | 8.8 ms | 0.67 s | during the next step (22 ms of kernels) |
| copy of the two stored slices | 6.7 ms | 0.51 s | during the next step |

Only the sync was on the critical path.

**What the search reads:**
- the current slice: e at every cube corner, and the full cells of the intersecting cubes;
- the previous slice and its stored copy: only for the vorticity output.

**The change:**
- **Evolve:** on a freeze-out step that is not also an output, start, check or vorticity step, it no longer syncs before the step. After `beginStep` it calls `DeviceHydro::syncCurrentDuringStep(fpCurr)`.
- **Why it is safe:** after `beginStep`, the slice of τ_n is the backend's previous slice (`state(1)`). Both stages read it and none writes it; the next step's second stage is the first to overwrite it.
- **The backend API:** the new entry `fetchUnwritten` serves only that slice.
  - Metal returns at once (shared memory), so the copy overlaps the step.
  - CUDA does its usual blocking fetch (its stages are synchronous).
  - The CPU backend has nothing to copy.
- The `device:` log line counts these copies separately.

**Result:** in the freeze-out run, the copy now takes place during the step, and the run is 7% faster.

**The one changed output:**
- **Cause:** on a step with a full sync, the host Fields are current, so the conservation sums came from the host Fields in double. These freeze-out steps now use the device rows in float, as every other step already did.
- **Size:** `global_conservation_laws.dat` changes on those 74 lines by at most 5e-8 of E (momenta 1e-8, N_B 4e-9).
- **Against the CPU code:** the difference is unchanged at 9.8e-7 of E, and the freeze-out surface is byte-identical.

### 8.3 Tests

- **New unit test:** "Device copy of the current slice during a step equals syncHost", on backend 1, ideal and viscous, in double and float.
  - It overwrites the host Fields, copies during the step, and compares with `syncHost` after the step, ghost cells included.
  - With the copy reading slice 0 instead of slice 1, it fails in all four configurations.
- **Unit tests (bgunit, clang):** all pass:
  - unittest_bugfix: 44 cases;
  - cell, grid, minmod and reconst: all pass.
- **CUDA:** the syntax check (host and device passes) passes.
- **Integration (bmetal):** 25/25 pass.
- **Byte identity with fe606d7:** the six small cases (CPU code, backend 1 in double and float, Metal viscous, ideal, freeze-out) and the 1.27M-cell runs are identical, apart from the conservation file of freeze-out runs.
