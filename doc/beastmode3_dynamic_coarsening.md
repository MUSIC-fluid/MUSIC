# beastMode 3: deciding dynamically when to coarsen the grid

Investigation on branch `gpu` (after 82cc757). The prototype is in the working tree and is not committed.

## Summary

- **Current trigger.** Every `NtauBlock` steps (at least 192), beastMode 3 doubles Δτ and coarsens dx and dy while Δτ < 0.05 fm/c. The trigger counts steps and ignores the fluid, and the number of coarsenings has no limit.
- **Proposed trigger.** Coarsen when an estimate of the error the next coarsening would make falls below a tolerance. Δτ, dx and dy each grow by at most a factor of 4 over the run.
- **Indicator.** For each cell above the freeze-out energy density, the relative second difference of e at twice the current spacing, |e(i+2) − 2e(i) + e(i−2)| / e(i), maximised over x and y. The trigger uses its energy-weighted mean over those cells. This is the curvature that a grid with twice the spacing would have to represent. The scheme's truncation error and the content the 1/4-1/2-1/4 coarsening filter removes are both proportional to it.
- **Results** on three lumpy (2+1)D viscous events, compared with a run that is never coarsened:
  - **Tolerance 0.01:** surface E and S within 3 × 10⁻⁴, ⟨u_⊥⟩ within 8 × 10⁻⁴ and the v2 proxy within 0.12 %, for 1.3–4.3× fewer cell updates.
  - **Current schedule:** 3.0–4.4× fewer cell updates, but up to 1.1 % error in the v2 proxy and 1.4 × 10⁻³ in S on the central lumpy event, because it coarsens at τ = 2.3 fm/c regardless of the profile.
  - **The indicator follows the profile.** It coarsens the smooth event early (τ = 1.8 fm/c at tolerance 0.01) and the lumpy events late (τ = 5.9–6.4 fm/c).
- **(3+1)D string event** (see the last section): the results agree with the 2D ones.
  - Tolerance 0.01: every observable within 7 × 10⁻⁵, and dS/dη within 7 × 10⁻⁴, for 1.74× fewer cell updates.
  - Current schedule (with the cap): 2.64× fewer, but 0.3 % in the v2 proxy and 0.4 % in dS/dη.
  - Without the cap, the current schedule doubles Δτ a third time, to 8×. That costs 0.5 % in E and S and 1 % in dS/dη and saves only 4 % more, which supports the factor-4 cap.
- **Recommendation.** Use the energy-weighted mean of the relative e curvature, with tolerance 0.01 as the default and 0.02 as an option when speed matters more (v2 proxy within 0.4 % in 2D, 0.05 % in 3D). Before making it the default, give the device a reduction kernel for the indicator, so the check doesn't need a full sync every 24 steps.

## Prototype (working tree)

| parameter | default | meaning |
|---|---|---|
| `beastModeCoarsenTolerance` | 0.01 | > 0: beastMode 3 coarsens when the indicator's energy-weighted mean drops below it. 0 keeps the old step-count schedule. |
| `beastModeMaxCoarsenFactor` | 4 | Δτ, dx and dy each grow at most by this factor over the run. Applies to beastMode 2 and 3, with either trigger. |

- **Where** (`src/evolve.cpp`): `Evolve::coarseningError` computes the indicator: energy-weighted mean and maximum, for e and for u^x and u^y. The main loop checks it every 24 steps, after the string sources are done, and logs it.
- **When the grid coarsens, Δτ doubles with it** if its cap allows. Δτ still doubles alone on the old CFL margin, dx / (2 Δτ) > 5, every `NtauBlock` steps.
- **Behaviour change:** the factor-4 cap also limits the existing schedule. Before, it had no limit: for example, beastMode 2 went from Δτ 0.005 to 0.04, a factor of 8, in the 2026-11-22 test.
- All 46 unit tests pass.

## Test setup

- **Initial states** (Initial_profile 9, written with the helpers of `tests/integration`): (2+1)D, τ0 = 0.4 fm/c, a 241² grid with dx = 0.1 fm (24 fm wide), Δτ = 0.01 fm/c.
  - **cent:** 300 hot spots of width 0.4 fm, at participant-like positions of two Woods–Saxon discs with b = 4 fm. Gamma-distributed weights, e_max = 67 GeV/fm³.
  - **peri:** the same with b = 9 fm and 120 spots.
  - **smooth:** the same as cent but with spots 0.8 fm wide, which gives a smooth profile.
- **Physics:** EoS 9, η/s = 0.12, second-order terms, freeze-out at e = 0.18 GeV/fm³, gridPadding 2 fm. The fluid never reached the grid edge, so no run enlarged the grid.
- **Reference:** beastMode 1, the same face scheme as beastMode 3 but with no grid or Δτ changes.
- **Observables:** the freeze-out surface's ideal-fluid energy and entropy fluxes, ⟨u_⊥⟩, and a v2 proxy Σ s dΣ·u (u_x² − u_y²) / Σ s dΣ·u (u_x² + u_y²).
- **Cost:** counted in cell updates, Σ_steps nx·ny. The 18 runs ran 8 at a time on 4 performance and 6 efficiency cores, so the wall times are noisy, though they agree with the cell counts.

## How the indicator evolves (cent, never coarsened)

| τ [fm/c] | 0.64 | 1.6 | 2.56 | 3.52 | 4.48 | 5.44 | 6.4 | 7.36 | 8.32 | 9.28 | 10.24 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| e mean | 0.124 | 0.054 | 0.045 | 0.031 | 0.021 | 0.014 | 0.0091 | 0.0059 | 0.0041 | 0.0032 | 0.0028 |
| e max | 2.3 | 0.74 | 0.67 | 0.38 | 0.28 | 0.18 | 0.14 | 0.034 | 0.016 | 0.0069 | 0.0040 |
| u mean | 0.025 | 0.027 | 0.019 | 0.014 | 0.0097 | 0.0062 | 0.0039 | 0.0023 | 0.0014 | 0.00085 | 0.00030 |

- **The mean falls smoothly and monotonically**, so it can be used as a trigger.
- **The maximum isn't usable.** It is dominated by a few cells at the edges of hot spots and of the fireball, and it jumps by factors of 4 between checks.
- **The velocity curvature adds little:** it tracks the e curvature at about 0.3–0.5 of its value.
- **The current schedule coarsens at τ = 2.3 fm/c**, where the mean is still 0.047, and again at τ = 6.1.
- **After a coarsening the indicator is about 4 times larger**, because it is now measured at four times the original spacing. A second coarsening therefore needs the profile to smooth out by another factor of 4, so tight tolerances rarely coarsen twice.

## Results

Relative differences to the reference. "When" lists the τ (fm/c) of each coarsening; the cap allows at most two.

| event | trigger | when | cell updates saved | E | S | ⟨u_⊥⟩ | v2 proxy |
|---|---|---|---|---|---|---|---|
| cent | current schedule | 2.31, 6.14 | 4.31× | +1.2e-3 | +1.4e-3 | −1.3e-4 | **+1.1e-2** |
| cent | tol 0.04 | 3.04, 5.92 | 3.47× | +4.1e-4 | +6.3e-4 | −2.3e-4 | +8.7e-3 |
| cent | tol 0.02 | 4.72, 7.60 | 2.26× | −7.7e-4 | −5.5e-4 | −6.7e-4 | +3.9e-3 |
| cent | tol 0.01 | 6.40 | 1.62× | +6.5e-5 | +1.4e-4 | −1.6e-4 | +5.2e-4 |
| cent | tol 0.005 | 7.84 | 1.36× | −2.6e-4 | −2.3e-4 | −1.2e-4 | +5.0e-4 |
| peri | current schedule | 2.31, 6.14 | 3.01× | −7.1e-4 | −3.2e-4 | −9.1e-4 | −2.1e-3 |
| peri | tol 0.04 | 3.04, 5.92 | 2.41× | −4.6e-4 | +2.9e-4 | −1.7e-3 | −1.6e-3 |
| peri | tol 0.02 | 4.72, 7.60 | 1.56× | −1.7e-3 | −1.9e-3 | +3.8e-4 | −3.9e-4 |
| peri | tol 0.01 | 5.92 | 1.27× | −2.5e-4 | −1.1e-4 | −3.8e-4 | −1.2e-3 |
| peri | tol 0.005 | 7.60 | 1.01× | −3.4e-4 | −4.1e-4 | +1.9e-4 | −1.9e-4 |
| smooth | current schedule | 2.31, 6.14 | 4.37× | +3.1e-4 | +3.4e-4 | +1.4e-4 | +4.7e-4 |
| smooth | tol 0.04 | 0.64, 2.08 | 19.1× | +3.7e-3 | +3.4e-3 | +2.2e-3 | −4.7e-3 |
| smooth | tol 0.02 | 0.88, 6.64 | 8.48× | +6.1e-4 | +5.7e-4 | +3.4e-4 | +1.7e-3 |
| smooth | tol 0.01 | 1.84, 10.48 | 4.27× | +1.3e-4 | −7.1e-5 | +8.0e-4 | −1.1e-3 |
| smooth | tol 0.005 | 6.64 | 1.59× | −1.0e-4 | −7.3e-5 | −8.1e-5 | +8.0e-4 |

Reference values (beastMode 1):

| event | E | S | ⟨u_⊥⟩ | v2 proxy | freeze-out τ |
|---|---|---|---|---|---|
| cent | 582.2 GeV | 3762 | 0.524 | 0.087 | 11.1 fm/c |
| peri | 231.8 GeV | 1418 | 0.622 | 0.229 | 7.7 fm/c |
| smooth | 610.0 GeV | 4012 | 0.498 | 0.132 | 11.3 fm/c |

What the results show:
- **The error follows the indicator, not the schedule.** On the central lumpy event the v2-proxy error falls steadily, from 1.1 % (current schedule) to 0.87 %, 0.39 % and 0.05 % (tol 0.04, 0.02, 0.01). The smooth event's v2 changes by only 0.05 % under the current schedule, even though that schedule coarsens it at the same times as the lumpy one.
- **At the same speed-up, the dynamic trigger is more accurate.** The current schedule gives the smooth event 4.4× fewer cell updates with 3–5 × 10⁻⁴ errors. Tolerance 0.01 gives 4.3× fewer, also with errors of a few 10⁻⁴, and it does that without coarsening the lumpy events early.
- **Tolerance 0.04 is too loose:** it coarsens the smooth event at τ = 0.64 fm/c, which costs 0.4 % in E and S.
- **Peripheral events gain little with tight tolerances.** They freeze out at τ ≈ 7.7, before the profile is smooth enough: tolerance 0.005 coarsens at τ = 7.6 and saves only 1 %.
- **Only E and S are near a noise floor.** Their differences of a few 10⁻⁴ at tolerances 0.005–0.01 don't fall with the tolerance: they are about the size of the jumps that the Δτ switch and the remap make in conservation (2026-11-22 report). The v2 proxy and ⟨u_⊥⟩ are the sensitive observables.

## Open points before making it the default

1. **(3+1)D.** Tested on one string event (next section). Only one event was run, and it freezes out at τ = 6 fm/c, so a longer-lived (3+1)D event, a central one at higher energy, would show more of the grid coarsening.
2. **Device cost.** Each check calls `syncAll` and computes the indicator on the host, every 24 steps. That is acceptable on the CPU (a single pass over the grid). On the GPU it would need a reduction kernel, or a check every `NtauBlock / 2` steps.
3. **The tolerance is relative.** It depends on how lumpy the initial state is and on the starting dx. A run that starts at dx = 0.2 measures its first check at 0.4 fm, so it coarsens later at the same tolerance. That is the intent: the tolerance limits the error of the coarse grid, not the coarsening factor.
4. **The cap changes the existing behaviour** for beastMode 2 and 3 runs that would have doubled more than twice. Either keep it, or default it to a large value when `beastModeCoarsenTolerance = 0`.
5. **Bulk viscosity and baryon diffusion** were not part of this test.

## Files

- Prototype: `src/evolve.cpp`, `src/evolve.h`, `src/data.h` and `src/read_in_parameters.cpp` in the working tree.
- Driver and logs: `common.py`, `probe.py`, `one.py`, `results.jsonl`, `cost.json` and one directory per run, in the session scratchpad `/private/tmp/claude-503/-Users-gf8206-Desktop-codes-MUSIC/8470541e-80e8-49d9-af7d-f4144ff84680/scratchpad/dyn/`. The scratchpad is temporary, so copy what you want to keep.

## (3+1)D test with the strings event

### Setup

- **Input:** `example_inputfiles/3D_dynamical/music_input_mode_2` with `strings_event_0.dat`, Initial_profile 13. The beam rapidity from the header is 5.36, and the string sources end at τ = 0.5 fm/c, so all of them are deposited before the first check.
- **Grid:** the dynamic grid is 120 × 120 × 107, with dx = 0.111 fm, dy = 0.129 fm and Δη = 0.119. It starts at τ0 = 0.48 fm/c with Δτ = 0.0057 fm/c, and all runs freeze out at τ = 6.0 fm/c.
- **Runs:** a g++-15 OpenMP build on 10 threads, one run at a time, so the wall times can be compared.
- **Reference:** beastMode 1.
- **Observables:** as in 2D, plus the net baryon number through the surface, N_B, and dS/dη in 24 bins over the surface. The dS/dη column gives the largest relative difference over the bins that carry more than 20 units of entropy.

### How the indicator evolves (never coarsened)

| τ [fm/c] | 0.62 | 0.89 | 1.17 | 1.58 | 2.13 | 2.68 | 3.23 | 3.78 | 4.33 | 4.88 | 5.43 | 5.70 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| e mean | 0.128 | 0.105 | 0.074 | 0.050 | 0.050 | 0.047 | 0.040 | 0.032 | 0.023 | 0.017 | 0.013 | 0.0099 |
| cells above freeze-out | 235k | 221k | 215k | 210k | 199k | 188k | 177k | 161k | 126k | 86k | 36k | 5k |

- **The curve has the same shape as in 2D:** the mean falls fast at first, then plateaus at about 0.05 for τ = 1.5–2.7 fm/c.
- **It reaches 0.01 only at τ = 5.7 fm/c,** when the fireball is almost frozen out.

### Results

"When" lists the τ (fm/c) of each change: t for Δτ doubling alone, g for a grid coarsening (with Δτ doubling, if the cap allows).

| trigger | when | cell updates saved | wall time saved | E | S | N_B | ⟨u_⊥⟩ | v2 proxy | dS/dη max |
|---|---|---|---|---|---|---|---|---|---|
| current schedule, cap 4 | 1.57t, 3.21g | 2.64× | 2.49× | −6.4e-4 | −3.0e-4 | −1.8e-3 | −2.3e-4 | +2.9e-3 | 3.9e-3 |
| current schedule, no cap | 1.57t, 3.21g, 4.32g | 2.76× | 2.65× | **−5.4e-3** | **−5.0e-3** | **−8.5e-3** | −2.5e-3 | +4.8e-3 | **1.0e-2** |
| tol 0.04 | 1.57t, 3.49g, 5.68g (dx only) | 2.51× | 2.29× | −9.7e-4 | −2.0e-3 | −1.9e-3 | −3.5e-5 | +1.1e-3 | 5.2e-3 |
| tol 0.02 | 1.57t, 4.59g | 2.05× | 2.10× | −3.5e-4 | −3.2e-4 | −8.3e-4 | −1.3e-4 | −4.9e-4 | 1.3e-3 |
| tol 0.01 | 1.57t, 5.69g | 1.74× | 1.82× | −2.2e-5 | −6.7e-5 | −1.6e-5 | +6.3e-5 | +5.0e-5 | 7.3e-4 |
| tol 0.005 | 1.57t | 1.67× | 1.72× | −1.0e-5 | −5.3e-6 | −2.8e-6 | +7.1e-5 | +6.7e-5 | 8.7e-4 |

Reference values (beastMode 1): E = 4041 GeV, S = 3944, N_B = 39.5, ⟨u_⊥⟩ = 0.491, v2 proxy = 0.273, 814 s wall time.

What the results show:
- **Doubling Δτ alone is almost free and gives most of the gain.** At τ = 1.57 fm/c the old CFL margin, dx / (2 Δτ) > 5, triggers in every run. That doubling alone gives 1.67× fewer cell updates, and the tolerance 0.005 run, which never coarsens the grid, changes every observable by less than 10⁻⁴.
- **Coarsening the grid is where the error comes from,** and it follows the indicator, as in 2D:
  - Coarsening at τ = 3.2 fm/c, where the indicator is 0.040 (current schedule): 0.3 % in the v2 proxy and 0.4 % in dS/dη.
  - At τ = 4.6 fm/c (tolerance 0.02): 0.05 % in the v2 proxy and 0.13 % in dS/dη.
  - At τ = 5.7 fm/c (tolerance 0.01): no change above the Δτ-switch floor.
- **N_B is the most sensitive conserved quantity:** −0.19 % at tolerance 0.04 and −0.85 % without the cap. Coarsening smooths the steep N_B gradients at forward rapidity.
- **The cap matters.** Without it the current schedule doubles Δτ a third time at τ = 4.3 fm/c, to 0.046 fm/c, or 8 Δτ0. With dx = 0.44 fm and dy = 0.52 fm that is still stable, but E and S change by 0.5 % and dS/dη by 1 %, ten times the error of the capped schedule, for 4 % fewer cell updates.
- **The 0.04 run's second coarsening is dx only:** Δτ had already reached its cap of 4 Δτ0. It came at τ = 5.68 fm/c, when only 5 k cells were above freeze-out, and it doubled the S error, from about 3e-4 to 2e-3. Coarsening this close to freeze-out is cheap in cell updates but costs accuracy at the surface. The tighter tolerances avoid it.
- **Wall times agree with the cell counts within 10 %.** The extra syncs for the check every 24 steps don't show at this size on the CPU.
- **Memory:** the current-RSS report falls from 2.1 GB to 0.86 GB after the two coarsenings in the 0.04 run.

### Files

- Driver: `run3d.py`, `queue.sh`, `ana.py`. Results: `results.jsonl`, `probe/log.txt` (indicator series) and one `s3d_*` directory per run. All of them are in the session scratchpad `/private/tmp/claude-503/-Users-gf8206-Desktop-codes-MUSIC/8470541e-80e8-49d9-af7d-f4144ff84680/scratchpad/dyn3d/`.
- The OpenMP build needs the uncommitted fix in `src/pretty_ostream.cpp`: 6031ba6 includes `<mach/mach.h>`, which GCC cannot parse on macOS, so the fix declares `proc_pid_rusage` directly.
