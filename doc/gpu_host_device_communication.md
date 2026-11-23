# GPU: communication between device and host — measurements and options

Branch `gpu` at 8aa8b51. 2026-11-23. Apple M4 (10 GPU cores), Metal, float.
- **Sections 1–4:** the investigation.
- **Section 5:** step 1 of its recommendation, implemented.

## Summary

- **The GPU is idle 6–9% of the time on production grids, and 18% on a small one.** All of that idle time is host work done between the end of one step and the start of the next. The host work done during a step already hides behind the kernels.
- **The conservation check:**
  - **Its sums already run on the device.** The host reads back 12 values per (η, x) row.
  - **The cost is the round trip,** not the arithmetic. It costs 0.26 s per 101-step run (3%) because the check is a separate launch and wait before each step.
  - **Encoding the conservation kernel into the step's own command buffer** removes the round trip and keeps the output bit-identical.
  - **Computing it less often, or only at freeze-out, is not needed.** It would also drop lines from `global_conservation_laws.dat`, and the grid maxima behind the per-step "maximum energy density increased" check.
- **The η-window copy and the outflow-face gather** read only τ_n, which no stage of the next step writes. They can move into the step, as the freeze-out copy did in 8aa8b51.
- **Together these three changes would remove nearly all GPU idle time.** Estimate: about −0.6 s per run (−5% to −8% of the wall time, −15% on small grids), with every output bit-identical. This is the recommended next step.
- **Step 1 is done (section 5).** Wall time −3% to −5% on 1.27M cells and −23% on the small grid, with every output byte-identical to 8aa8b51. The remaining idle time, about 3 ms per step, is the launch floor: encoding stage 1, and the wake-up after the wait.
- **Freeze-out detection on the device:**
  - **Here, it is not on the critical path.** The intersecting cubes touch 2.25% of the cells, and their detection on the device would remove 14 of the 22 ms of host work per freeze-out step. On the M4 that host work already overlaps the next step (24.6 ms), so the wall-time gain here is small (estimated 1–2%).
  - **It becomes worthwhile when the host is the bottleneck:**
    - after the change above (freeze-out steps would be host-bound by 1–2 ms);
    - on CUDA, where the copy of the current slice blocks and goes over PCIe;
    - on GPUs several times faster than the M4;
    - with several freeze-out energy densities;
    - to save the host memory of the two stored full-grid Fields.

  It is the second step. It is larger (new kernels on three backends, device-resident stored slices, care with beastMode regridding), and it can stay bit-identical.

## 1. What is exchanged, and when

On the device path (`gpu_backend` > 0), each step of Evolve exchanges the following with the device, in beastMode 0:

| exchange | when | reads | placement |
|---|---|---|---|
| η window \|η\| < 0.5 of fpCurr (`output_momentum_anisotropy_vs_tau`) | every step | τ_n | before `beginStep` |
| conservation rows: kernel, wait, readback of 12 values per (η, x) row | every step | τ_n with τ_{n−1} | before `beginStep` |
| outflow faces of fpCurr (gathered from the device state) | every step | τ_n | before `beginStep` |
| full sync of both slices | output, freeze-start and check steps | τ_n, τ_{n−1} | before `beginStep` |
| current slice (8aa8b51) | freeze-out steps | τ_n | during the step |

In beastMode 1 to 3 the per-step diagnostics are off, and only the full syncs and the beastMode 3 edge band remain.

## 2. Measurements

**Two temporary probes** (removed):
- a Metal command-buffer timeline: the GPU start and end time of every command buffer, which gives the GPU's idle gaps and what each gap ended with;
- host timers around each exchange in Evolve.

**GPU timeline** (one run each):

| case | GPU span | busy: stage 1 / stage 2 / conservation | idle before the conservation kernel | idle before stage 1 | idle share |
|---|---|---|---|---|---|
| viscous EoS 9, 161×161×49, 101 steps | 7.02 s | 2.74 / 3.43 / 0.18 s | 0.35 s | 0.32 s | 9% |
| viscous EoS 20, same grid | 10.69 s | 4.27 / 5.55 / 0.19 s | 0.35 s | 0.34 s | 6% |
| freeze-out run, 101×101×41, 381 steps | 10.27 s | 3.77 / 5.61 / 0.18 s | 0.40 s | 0.32 s | 7% |
| viscous EoS 9, 41×41×17, 101 steps | 0.40 s | 0.13 / 0.19 / 0.01 s | 0.04 s | 0.03 s | 18% |

Stage 2 starts right after stage 1 in every run: the GPU never waits between the two stages.

**Host time per run, viscous EoS 9 on 1.27M cells:**

| | total | per step |
|---|---|---|
| *before `beginStep` (the GPU waits)* | **0.68 s** | 6.7 ms |
| conservation rows (launch, wait, readback) | 0.26 s | 2.7 ms |
| η-window copy | 0.16 s | 1.7 ms |
| outflow-face gather | 0.14 s | 1.5 ms |
| 3 full syncs (output steps) | 0.09 s | |
| *after `beginStep` (overlaps the step)* | 1.68 s | 17 ms |
| of which `output_momentum_anisotropy_vs_tau` | 1.03 s | 10 ms |
| of which `check_conservation_law` | 0.27 s | 2.7 ms |
| waiting for the step | 4.73 s | 47 ms |

**Freeze-out steps** (76 of the 381 steps; freeze-out every 5 steps, `fac_x = fac_y = 1`):

| host work per freeze-out step | ms |
|---|---|
| copy of the current slice (during the step) | 4.6 |
| intersection scan over all cubes (in the Cornelius search) | 2.4 |
| surface work of the 4,877 intersecting cubes: Cornelius, interpolation, EoS, writing | 8.3 |
| the two full copies of the stored slices | 7.0 |
| other host work of the step | ≈ 0.5 |
| **total, against a GPU step of 24.6 ms** | **≈ 22** |

On freeze-out steps `AdvanceRK` waits only 1.4 ms for the GPU: these steps are almost host-bound.
- The intersecting cubes touch 9,400 distinct cells per step, 2.25% of the grid.
- The scan was timed by skipping the surface work (`continue` after the test).

## 3. Options

### 3.1 The conservation check (recommended, bit-identical)

The sums already run on the device (`kernelConservationRow`, phase 5). What costs is the round trip:
- **The order of the work:** after the step finishes, the host launches the kernel, waits for it, and reads back 7889 × 12 values. Only then does it start the next step.
- **The idle time:** the GPU idles before the kernel, during the η-window copy, and again before the next stage 1, during the readback and the outflow gather.

**Proposal:** encode `kernelConservationRow` of the new slice with the old one at the end of the step's second command buffer (`stage(1)`). The rows are then ready when `finish` returns.
- **Fallback:** the first step, and steps after a host upload, keep the separate call.
- **beastMode > 0:** the kernel is not encoded.
- **Output:** the same kernel and the same inputs, so the rows are bit-identical.
- **Cost that remains:** the kernel's 0.18 s of GPU time (2.5%). It is bound by the serial sum over y in each row. A two-level reduction would bring it to about 0.02 s, but it would change the float summation order: `global_conservation_laws.dat` would change at 1e-7. That is optional.

**Against computing it less often, or only at freeze-out:**
- **The data is used every step:** the rows give each line of `global_conservation_laws.dat` and the grid maxima. After the freeze-out starts, the maxima feed the per-step "maximum energy density increased" check, which exits on a runaway.
- **No speed reason:** once the round trip is gone, the check costs nothing on the critical path.
- **If fewer lines are wanted anyway,** an input parameter (every N steps) is easy. It is an output choice, not a performance need.

### 3.2 The η-window copy and the outflow gather (recommended, bit-identical)

Both read only τ_n. While the next step runs, τ_n is the backend's `state(1)`, read by both stages and written by neither. This is the property that 8aa8b51 uses for the freeze-out copy. The outflow faces of τ_{n−1} come from the previous step's gather and need no read.
- **The move:** `syncEtaWindow(fpCurr)` and `outflowFaces(fpCurr, 0)` go after `beginStep`, reading `state(1)` through `fetchUnwritten`.
- **What moves off the critical path:** 3.2 ms per step.

### 3.3 Expected effect of 3.1 and 3.2

**Estimate:** the GPU idle time falls from about 0.66–0.71 s per run to the launch latency:

| case | wall now | estimate |
|---|---|---|
| viscous EoS 9, 1.27M cells | 7.3 s | −0.6 s (−8%) |
| viscous EoS 20, 1.27M cells | 12.0 s | −0.6 s (−5%) |
| freeze-out run | 10.5 s | −0.6 s (−6%), less if its freeze-out steps become host-bound (3.4) |
| 41×41×17 | 0.47 s | −0.07 s (−15%) |

The work moved into the step adds about 3 ms to the host work during each step. Normal steps have plenty of margin (17 → 20 ms against 62 ms of GPU). Freeze-out steps do not (22 → 25 ms against 24.6 ms of GPU).

### 3.4 Freeze-out detection on the device (worthwhile later)

**What the host reads now on a freeze-out step:**
- the e of every cube corner of the current slice and of the stored slice (the intersection test);
- the full cells at the corners of the intersecting cubes (2.25% of the cells);
- the previous slices only for the vorticity output.

The stored slices are full host copies of the Fields (7 ms per freeze-out step, and two full-grid Fields of host memory).

**Design:**
- **Device-resident stored slices:** at a freeze-out step, the device keeps a copy of τ_n, a GPU blit of about 0.5 ms (and of τ_{n−1} for the vorticity).
- **Marking kernel:** a kernel tests the cubes against the stored and current e. It writes the intersecting cubes, and the union of their corner cells, to compact lists.
- **The copy:** the host copies only those cells of both slices.
- **Unchanged output:** the host runs Cornelius over the cube list, sorted in the host's (η, x, y) order, so the surface files stay byte-identical.

**What it saves:** about 14 of the 22 ms of host work per freeze-out step (the copy, the scan and the stored copies). The surface work (8.3 ms) stays on the host.

**When it matters:**
- **On the M4, as measured:** the freeze-out host work already overlaps the next step, so the gain is the 1.4 ms margin, about 0.1–0.2 s per run. After 3.1 and 3.2 the freeze-out steps would be host-bound by about 1–2 ms each, and this change would recover that.
- **On CUDA** (estimates, not measured; there is no NVIDIA hardware here):
  - **The copy blocks:** `fetchUnwritten` waits for the step, as CUDA's fetch must, so the copy of the current slice is not overlapped. It also crosses PCIe: 52 MB per freeze-out step for 418k cells and 157 MB for 1.27M cells in float, about 2–6 ms at 25 GB/s, plus the float→double conversion (4.6 ms here for 418k cells).
  - **The GPU step is shorter:** a data-centre GPU would run the step several times faster than the M4, and the 22 ms of host work would then dominate the freeze-out steps.
  - **Device-side detection:** it would cut the copy by about 98%. An asynchronous copy on a second stream would hide the copy, but not the host work.
- **Several freeze-out energy densities:** each surface repeats the scan.

**Cost:**
- a marking kernel and a gather kernel on Metal, CUDA and the CPU backend;
- two extra device slices (157 MB each in float for 1.27M cells);
- the beastMode paths that regrid or extrapolate the stored slices, which would have to sync them first;
- the boost-invariant and low-temperature freeze-out variants.

### 3.5 Not worth doing now

- **`output_momentum_anisotropy_vs_tau` on the device:** it is the largest host item (10 ms per step), but it overlaps the step with a large margin. Move it to the device, or make it optional, only if a faster GPU makes the host the bottleneck.
- **The full syncs of output steps:** three per run here. The τ_n part could also move into the step, but τ_{n−1} is overwritten by stage 2 and would need a wait for stage 1.

## 4. Recommendation

1. **Do 3.1 and 3.2 together:** conservation rows in the step's command buffer, and the τ_n reads during the step. Expected −5% to −8% on production grids and −15% on small ones, bit-identical.
2. **Then 3.4:** freeze-out detection on the device, with device-resident stored slices. It is mainly for CUDA and faster GPUs, and it keeps the freeze-out steps from becoming host-bound after step 1.
3. **Keep the conservation check every step.** With 3.1 it costs nothing on the critical path.

## 5. Done: step 1 (the conservation rows after the step, the τ_n reads during it)

### 5.1 The changes

**The conservation rows are computed right after the step that made the slice.**
- `DeviceHydro::conservationAhead(true)`, which Evolve sets in beastMode 0, makes `beginStep` queue `kernelConservationRow` of the new slice with τ_n at τ_n + `delta_tau`, behind the second stage.
- **Backend API:** the new pair `conservationBegin` and `conservationEnd`.
  - Metal commits a command buffer without waiting; `waitPending`, and thus `finish`, also waits for it.
  - The CPU backend and CUDA, whose stages are synchronous, compute the rows at once (the `Backend` default; CUDA's entries are NULL).
- **Use:** `conservationRows(τ)` takes the queued rows when τ is exactly the τ they were queued for. Otherwise, after an upload, on the first step, or when Evolve's τ differs, it computes them as before. The kernel and its parameters are the same, so the rows are bit-identical.

**The τ_n reads happen while the next step runs.**
- `syncEtaWindow` and `outflowFaces(curr, 0)` may now be called during a step: they read the backend's slice 1 through `fetchUnwritten`.
- The outflow faces of τ_{n−1} are still taken before `beginStep`: they come from the cache, which the previous step's gather filled. During a step, `outflowFaces(prev, 1)` fails unless they are cached.
- On freeze-out steps the window copy is skipped: `syncCurrentDuringStep` has copied all of fpCurr.

**The log line** reports how many conservation calls used the queued rows.

### 5.2 Results

**Wall time, 8aa8b51 against step 1** (Metal; two interleaved runs of each build):

| case | 8aa8b51 | step 1 | change |
|---|---|---|---|
| viscous, EoS 9, 161×161×49 | 8.01 / 7.34 s | 6.96 / 6.96 s | −5% |
| viscous, EoS 20, 161×161×49 | 12.28 / 11.94 s | 11.56 / 11.52 s | −4% |
| viscous, EoS 14, diffusion, 161×161×49 | 10.86 / 10.71 s | 10.63 / 10.66 s | −1% |
| ideal, EoS 9, 161×161×49 | 3.62 / 3.47 s | 3.38 / 3.39 s | −3% |
| freeze-out run, 101×101×41 | 10.45 / 10.51 s | 10.08 / 10.10 s | −4% |
| viscous, EoS 9, 81×81×25 | 1.16 / 1.17 s | 1.09 / 1.10 s | −6% |
| viscous, EoS 9, 41×41×17 | 0.49 / 0.48 s | 0.38 / 0.37 s | −23% |

**Every output file of these seven cases is byte-identical** to 8aa8b51. The six small regression cases are byte-identical to fe606d7, apart from the freeze-out conservation file that 8aa8b51 already changed.

**GPU timeline** (the same temporary probe as in section 2):

| | 8aa8b51 | step 1 |
|---|---|---|
| EoS 9: idle before the conservation kernel / before stage 1 | 0.37 / 0.31 s | 0 / 0.32 s |
| diffusion: idle before the conservation kernel / before stage 1 | 0.35 / 0.33 s | 0 / 0.27 s |

The kernel times are unchanged, so the copies during the step do not slow the GPU.

**The idle time before stage 1 is still about 3 ms per step.** A host timer between `advanceStep` and the commit of stage 1 shows where it goes:
- the three full syncs of output steps (0.09 s per run);
- `beginStep` (1.4 ms per step: the parameters and the encoding of both stages; stage 1 is committed halfway);
- the conservation readback and the outflow cache (0.1 ms);
- the wake-up of the wait and the GPU's start-up latency.

Removing it would need the next step to be encoded before the host has finished with the current one. For about 0.2 s per run that is not done.

**Why the diffusion run gains least:**
- Its host wait fell by 0.8 s, and its GPU idle before the conservation kernel by 0.35 s, as for EoS 9.
- In 8aa8b51 its wall time was already shorter than the sum of its parts would suggest, and the repeated runs differ by more than for the other cases.
- The cause of the smaller wall-time gain is not established.

### 5.3 Tests

- **New unit test:** "Device reads of τ_n during a step equal those after it", on backend 1, ideal and viscous, in double and float. It checks that:
  - the queued rows equal a direct `conservationRows` call, and that they were used;
  - the η window copied during the step equals `syncHost` of the same slice after it;
  - the outflow faces taken during the step give the outflow rate of the Fields of that slice, and carry over as the previous slice's faces.

  **Mutations:** reading slice 0 during the step fails 24 of its 76 assertions, and queuing the rows at τ_n fails 2.
- **Unit tests (bgunit, clang):** all pass:
  - unittest_bugfix: 45 cases;
  - cell, grid, minmod and reconst: all pass.
- **CUDA:** the syntax check (host and device passes) passes.
- **Integration (bmetal, with the bgunit memory driver):** 25/25 pass.
