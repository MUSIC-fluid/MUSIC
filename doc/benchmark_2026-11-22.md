# MUSIC overnight benchmark, 2026-11-21/22

Branch `optimization` at 201ccad **plus the uncommitted items 1–4** (fused W flux/source pass, EoS lookup reuse, cheaper Newton start, inlined index helpers). Apple M4 (4 performance + 6 efficiency cores), g++-15 -O3 OpenMP, every run under `caffeinate`. 226 runs from 22:58 to 04:09 EST, plus six diagnostic runs afterwards.

## Verdict

- **Items 1–4 are safe to commit.**
  - Unit tests pass and the integration suite gives 18 passed, 1 skipped.
  - With item 3 (the Newton start) off, all 760 outputs are bit-identical to 201ccad.
  - Item 3 changes results only at the level of the 1e-14 Newton tolerance. Conservation and freeze-out surfaces match 201ccad to every printed digit; see the baseline tables.
  - Speed-up against 201ccad: 3–11% at 10 threads in back-to-back tests (1–4% in this matrix). At 1 thread it is 1–4% in clang back-to-back tests and −2% to +3% (noise) in this matrix.
- **Stability.**
  - 209 of the 210 runs of the current code froze out cleanly, with no NaN, no error and N_B ≥ 0. The exception is the dense beastMode 3 run below. The other 4 failures are public_stable runs.
  - The only warning, in 97 runs, is the input check "Delta_Tau maybe too large", which resets Δτ.
- **Pre-existing bug: beastMode 3 is unstable after a mid-run grid coarsening in long, dense (3+1)D evolutions** (see the beastMode 3 section).
  - It was not introduced by the optimisation branch: the base 227d813 fails the same way, and so do 201ccad and face_reconstruction 0.
  - The coarse grid is stable when used from the start, so the coarsening/Δτ switch itself is the suspect.
- **Conservation in (3+1)D at beastMode 0:** max |E/E0 − 1| = 2–10 × 10⁻⁶, max |N_B/N_B0 − 1| = 1.6–3.5 × 10⁻⁵, and |ΔP_z|/E0 ≤ 1 × 10⁻⁶. beastMode 1 is just as good. beastMode 3 is 10–100× worse, because of the remap at the coarsening.
- **Threads:** 1 and 10 threads give bit-identical output files in all 84 configurations. Surface flux sums differ by ≤ 6 × 10⁻¹⁴ (relative), because the surface cells are summed in a different order.
- **Speed (3+1)D, 10 threads:**
  - 8 s for EoS 0 ideal up to 52 s for EoS 20 with shear+bulk+diffusion, on 61×61×31.
  - 10-thread scaling is 6.1–6.6×; 4 P + 6 E cores give about 6× at most.
  - beastMode 1 is 1.6–3.4× faster than beastMode 0.
  - Against public_stable: 3.2× (EoS 0) and 2.6× (EoS 9 shear+bulk). public_stable cannot run the other two configurations.

## Setup

| | (3+1)D | (2+1)D, boost invariant |
|---|---|---|
| grid | 61×61×31, 12 fm × 12 fm × 8 units (dx = 0.2 fm) | 201×201, 15 fm (dx = 0.075 fm) |
| initial state | lumpy TATB (integration-test generator), τ0 = 1 fm/c | same |
| freeze-out | e = 0.3 GeV/fm³ | e = 0.15 GeV/fm³ |
| viscosity | η/s = 0.12, bulk on by parametrisation, κ coefficient 0.2 | same |

- **Matrix:** EoS 0, 9, 14, 20 × {ideal, shear, shear+bulk, shear+bulk+diffusion} × beastMode {0, 1, 3} × {1, 10 threads} × {2+1, 3+1}. Diffusion is run only for the baryon EoS 14 and 20. That makes 168 runs.
- **beastMode settings:** beastMode > 0 uses face_reconstruction 1 by default. beastMode 3 also uses gridPadding = 2 fm.
- **beastMode conservation:** beastMode 1 and 3 skip every diagnostic. Their conservation numbers come from a measurement-only build that keeps `check_conservation_law` on.
- **Baselines:** public_stable (1fd9f40, 445 commits older), 201ccad, and the current tree.
- **Long tests:**
  - (3+1)D at 121×121×61;
  - a dense, hot state (197 nucleons per nucleus, √s = 200 GeV, e_max ≈ 2.9 GeV/fm³) with freeze-out at 0.1 GeV/fm³, evolved to τ ≈ 20 fm/c in (3+1)D.

## Run time, current code [s] (1 thread / 10 threads, speed-up)

### (3+1)D, 61×61×31

| EoS | viscosity | beastMode 0: 1 / 10 threads [s] | beastMode 1 | beastMode 3 |
|---|---|---|---|---|
| 0 | ideal | 51.3 / 8.2 (6.3×) | 14.7 / 3.2 (4.6×) | 24.0 / 5.4 (4.4×) |
| 0 | shear | 66.8 / 10.9 (6.1×) | 33.3 / 6.3 (5.3×) | 54.4 / 10.2 (5.3×) |
| 0 | shear+bulk | 72.9 / 12.0 (6.1×) | 38.5 / 7.2 (5.3×) | 63.2 / 11.6 (5.5×) |
| 9 | ideal | 90.7 / 14.2 (6.4×) | 24.0 / 5.0 (4.8×) | 28.8 / 6.3 (4.6×) |
| 9 | shear | 113.1 / 17.9 (6.3×) | 50.6 / 9.5 (5.3×) | 68.6 / 12.9 (5.3×) |
| 9 | shear+bulk | 137.3 / 21.9 (6.3×) | 66.1 / 12.3 (5.4×) | 77.4 / 14.5 (5.3×) |
| 14 | ideal | 123.8 / 18.7 (6.6×) | 30.2 / 6.1 (5.0×) | 37.1 / 7.7 (4.8×) |
| 14 | shear | 146.4 / 22.5 (6.5×) | 57.9 / 11.0 (5.3×) | 79.4 / 15.0 (5.3×) |
| 14 | shear+bulk | 176.8 / 27.7 (6.4×) | 75.9 / 14.3 (5.3×) | 89.5 / 16.8 (5.3×) |
| 14 | shear+bulk+diff | 219.5 / 35.0 (6.3×) | 117.9 / 21.5 (5.5×) | 135.6 / 24.9 (5.4×) |
| 20 | ideal | 207.6 / 31.6 (6.6×) | 51.6 / 9.4 (5.5×) | 64.4 / 12.4 (5.2×) |
| 20 | shear | 226.8 / 35.0 (6.5×) | 86.9 / 15.4 (5.6×) | 118.1 / 21.4 (5.5×) |
| 20 | shear+bulk | 272.4 / 42.1 (6.5×) | 110.1 / 19.4 (5.7×) | 128.2 / 23.1 (5.5×) |
| 20 | shear+bulk+diff | 332.4 / 52.2 (6.4×) | 168.1 / 29.5 (5.7×) | 193.5 / 34.5 (5.6×) |

### (2+1)D, 201×201

| EoS | viscosity | beastMode 0: 1 / 10 threads [s] | beastMode 1 | beastMode 3 |
|---|---|---|---|---|
| 0 | ideal | 17.3 / 3.0 (5.8×) | 4.9 / 1.3 (3.7×) | 4.7 / 1.3 (3.5×) |
| 0 | shear | 22.6 / 3.8 (5.9×) | 11.2 / 2.3 (4.9×) | 11.0 / 2.3 (4.8×) |
| 0 | shear+bulk | 23.6 / 4.1 (5.8×) | 12.7 / 2.5 (5.1×) | 12.5 / 2.5 (5.0×) |
| 9 | ideal | 26.0 / 4.5 (5.8×) | 7.6 / 1.9 (4.0×) | 5.5 / 1.5 (3.7×) |
| 9 | shear | 33.4 / 5.7 (5.9×) | 16.3 / 3.2 (5.1×) | 12.6 / 2.6 (4.8×) |
| 9 | shear+bulk | 34.7 / 5.9 (5.9×) | 18.4 / 3.5 (5.2×) | 14.1 / 2.8 (5.0×) |
| 14 | ideal | 54.4 / 8.9 (6.1×) | 13.9 / 3.3 (4.2×) | 7.5 / 2.2 (3.5×) |
| 14 | shear | 43.6 / 7.5 (5.8×) | 19.0 / 3.9 (4.8×) | 14.8 / 3.3 (4.4×) |
| 14 | shear+bulk | 45.3 / 7.8 (5.8×) | 21.5 / 4.3 (5.0×) | 16.1 / 3.5 (4.6×) |
| 14 | shear+bulk+diff | 56.6 / 9.7 (5.8×) | 32.2 / 6.3 (5.1×) | 24.3 / 5.0 (4.9×) |
| 20 | ideal | 61.4 / 10.7 (5.8×) | 25.4 / 5.4 (4.7×) | 14.0 / 3.5 (4.0×) |
| 20 | shear | 67.9 / 11.7 (5.8×) | 29.9 / 5.9 (5.0×) | 22.9 / 4.9 (4.7×) |
| 20 | shear+bulk | 70.8 / 12.2 (5.8×) | 32.9 / 6.4 (5.1×) | 24.4 / 5.1 (4.7×) |
| 20 | shear+bulk+diff | 87.2 / 15.0 (5.8×) | 49.0 / 9.1 (5.4×) | 36.4 / 7.2 (5.1×) |

- **beastMode 3 against beastMode 1 in short (3+1)D runs.** beastMode 3 is *slower* here, because it enlarges the grid from 61² to 81² at the start (gridPadding). The Δτ doubling and coarsening only trigger at step 192, close to freeze-out.
- **Long runs.** At 121×121×61 the trade-off reverses: EoS 20 takes 779 s at beastMode 0, 428 s at beastMode 1 and 301 s at beastMode 3. beastMode 3 also uses about 3× the memory of beastMode 1 in the short runs (549 MB against 173 MB).
- **(2+1)D EoS 14 ideal** runs longer than shear (54 s against 44 s at 1 thread) only because the ideal fluid lives longer (τ_end 4.2 against 3.3 fm/c).

## Baselines (3+1)D, beastMode 0 [s]

| EoS | viscosity | 1 thread: public_stable | 201ccad | now | speed-up | 10 threads: public_stable | 201ccad | now | speed-up |
|---|---|---|---|---|---|---|---|---|---|
| 0 | ideal | 176.8 | 51.1 | 51.3 | 3.45× | 25.9 | 8.4 | 8.2 | 3.16× |
| 9 | shear+bulk | 350.2 | 138.9 | 137.3 | 2.55× | 56.1 | 22.9 | 21.9 | 2.56× |
| 14 | shear+bulk+diff | – | 216.2 | 219.5 | – | – | 35.3 | 35.0 | – |
| 20 | shear+bulk+diff | – | 341.8 | 332.4 | – | – | 54.3 | 52.2 | – |

| EoS | viscosity | build | exit | max|E/E0-1| | max|N_B/N_B0-1| | surface E [GeV] | surface S |
|---|---|---|---|---|---|---|---|
| 0 | ideal | orig | 0 | 1.6e-01 | – | – | – |
| 0 | ideal | 201 | 0 | 4.9e-06 | – | 545.93 | 2377.4 |
| 0 | ideal | new | 0 | 4.9e-06 | – | 545.93 | 2377.4 |
| 9 | shear+bulk | orig | 0 | 1.1e-01 | – | – | – |
| 9 | shear+bulk | 201 | 0 | 7.4e-06 | – | 578.71 | 2124.3 |
| 9 | shear+bulk | new | 0 | 7.4e-06 | – | 578.71 | 2124.3 |
| 14 | shear+bulk+diff | orig | 1 | 1.9e-02 | 2.0e+08 | – | – |
| 14 | shear+bulk+diff | 201 | 0 | 9.6e-06 | 1.5e-05 | 572.43 | 2176.4 |
| 14 | shear+bulk+diff | new | 0 | 9.6e-06 | 1.6e-05 | 572.43 | 2176.4 |
| 20 | shear+bulk+diff | orig | 1 | – | – | – | – |
| 20 | shear+bulk+diff | 201 | 0 | 9.0e-06 | 3.0e-05 | 571.20 | 2188.1 |
| 20 | shear+bulk+diff | new | 0 | 9.0e-06 | 3.0e-05 | 571.20 | 2188.1 |

- **public_stable failures.** It has no EoS 20 ("EOS_to_use unspecified or invalid option: 20"). EoS 14 with diffusion blows up: N_B grows by 2 × 10⁸ and the run exits with NaN in MaxSpeed. This is the lumpy-diffusion instability fixed later on this branch.
- **public_stable conservation.** Its 11–16% energy "loss" is not comparable: its conservation diagnostic did not account for the flux leaving the grid.
- **201ccad against now:** identical physics to every printed digit, with the 1-thread times within noise (−2% to +3%). At 10 threads the current tree is faster for EoS 9 (22.9 → 21.9 s) and EoS 20 (54.3 → 52.2 s).

## Conservation, (3+1)D, 10 threads

Each cell gives max |E/E0 − 1| / max |ΔP_z|/E0 / max |N_B/N_B0 − 1| over the evolution. (2+1)D runs write no conservation file, because dE/dη is not conserved under longitudinal work.

| EoS | viscosity | beastMode 0 | beastMode 1 | beastMode 3 |
|---|---|---|---|---|
| 0 | ideal | 4.9e-06 / 6.8e-07 / – | 3.6e-06 / 6.7e-07 / – | 2.7e-04 / 2.5e-05 / – |
| 0 | shear | 7.1e-06 / 1.0e-06 / – | 5.6e-06 / 1.1e-06 / – | 4.0e-06 / 3.7e-07 / – |
| 0 | shear+bulk | 7.5e-06 / 1.0e-06 / – | 6.1e-06 / 1.1e-06 / – | 4.1e-06 / 3.0e-07 / – |
| 9 | ideal | 2.6e-06 / 4.9e-07 / – | 9.5e-07 / 4.1e-07 / – | 1.7e-04 / 2.2e-05 / – |
| 9 | shear | 5.3e-06 / 3.8e-07 / – | 3.5e-06 / 4.2e-07 / – | 1.5e-04 / 2.6e-05 / – |
| 9 | shear+bulk | 7.4e-06 / 5.4e-07 / – | 5.4e-06 / 4.1e-07 / – | 9.0e-05 / 2.1e-05 / – |
| 14 | ideal | 2.3e-06 / 7.2e-07 / 2.3e-05 | 7.7e-07 / 3.5e-07 / 2.1e-05 | 1.7e-04 / 2.1e-05 / 6.1e-05 |
| 14 | shear | 5.1e-06 / 3.8e-07 / 2.4e-05 | 3.4e-06 / 3.9e-07 / 1.5e-05 | 1.5e-04 / 2.6e-05 / 5.8e-05 |
| 14 | shear+bulk | 9.4e-06 / 2.4e-07 / 2.0e-05 | 5.4e-06 / 4.1e-07 / 1.3e-05 | 8.9e-05 / 2.2e-05 / 3.3e-05 |
| 14 | shear+bulk+diff | 9.6e-06 / 2.5e-07 / 1.6e-05 | 5.4e-06 / 4.1e-07 / 1.3e-05 | 8.8e-05 / 2.2e-05 / 3.6e-05 |
| 20 | ideal | 2.5e-06 / 4.1e-07 / 2.5e-05 | 8.0e-07 / 2.9e-07 / 2.3e-05 | 1.7e-04 / 2.1e-05 / 6.9e-05 |
| 20 | shear | 5.1e-06 / 3.5e-07 / 3.0e-05 | 3.3e-06 / 4.2e-07 / 2.5e-05 | 1.5e-04 / 2.5e-05 / 5.7e-05 |
| 20 | shear+bulk | 8.8e-06 / 3.5e-07 / 3.5e-05 | 5.5e-06 / 4.0e-07 / 2.7e-05 | 8.8e-05 / 2.1e-05 / 3.3e-05 |
| 20 | shear+bulk+diff | 9.0e-06 / 3.8e-07 / 3.0e-05 | 5.6e-06 / 4.0e-07 / 2.3e-05 | 8.8e-05 / 2.2e-05 / 3.9e-05 |

- beastMode 0 and 1 conserve E to about 10⁻⁵ and N_B to about 3 × 10⁻⁵, consistent with the KT outflow accounting.
- **beastMode 3 is 10–100× worse:** a jump of about 1e-4 (up to 1e-3 in the dense run) at the coarsening remap. The remap filter conserves the integrals in the interior, so the loss most likely comes from the grid edge or from the time-derivative bookkeeping (see below).

## beastMode accuracy

Relative change of the freeze-out surface fluxes against beastMode 0: E / S / ⟨u_⊥⟩ / N_B (N_B only for EoS 14 and 20).

### (3+1)D

| EoS | viscosity | beastMode 1 | beastMode 3 |
|---|---|---|---|
| 0 | ideal | +4.0e-03 / +9.9e-04 / +1.1e-02 | +4.0e-03 / +9.9e-04 / +1.1e-02 |
| 0 | shear | +4.6e-03 / +1.7e-03 / +8.5e-03 | +4.6e-03 / +1.7e-03 / +8.5e-03 |
| 0 | shear+bulk | +4.7e-03 / +1.8e-03 / +8.7e-03 | +4.7e-03 / +1.8e-03 / +8.7e-03 |
| 9 | ideal | +3.2e-03 / +1.7e-03 / +9.6e-03 | -1.4e-03 / -3.7e-03 / +5.5e-03 |
| 9 | shear | +3.4e-03 / +2.1e-03 / +7.3e-03 | +1.4e-03 / -8.1e-04 / +6.9e-03 |
| 9 | shear+bulk | +5.3e-03 / +4.0e-03 / +1.2e-02 | +7.7e-04 / -9.7e-04 / +9.6e-03 |
| 14 | ideal | +2.1e-03 / +9.3e-04 / +7.1e-03 / +5.9e-04 | -4.6e-03 / -8.0e-03 / +2.6e-03 / -1.3e-02 |
| 14 | shear | +2.6e-03 / +1.8e-03 / +5.4e-03 / +1.9e-03 | +3.3e-04 / -1.7e-03 / +5.0e-03 / -3.2e-03 |
| 14 | shear+bulk | +3.6e-03 / +2.6e-03 / +6.7e-03 / +2.1e-03 | -7.8e-04 / -3.1e-03 / +5.9e-03 / -5.7e-03 |
| 14 | shear+bulk+diff | +3.5e-03 / +2.6e-03 / +6.2e-03 / +2.4e-03 | -5.4e-04 / -2.6e-03 / +5.4e-03 / -3.6e-03 |
| 20 | ideal | +1.9e-03 / +7.8e-04 / +6.5e-03 / +5.8e-04 | -5.3e-03 / -8.6e-03 / +1.9e-03 / -1.3e-02 |
| 20 | shear | +2.5e-03 / +1.7e-03 / +5.3e-03 / +1.7e-03 | +5.2e-05 / -2.0e-03 / +4.9e-03 / -3.2e-03 |
| 20 | shear+bulk | +3.3e-03 / +2.4e-03 / +6.2e-03 / +1.9e-03 | -1.2e-03 / -3.6e-03 / +5.2e-03 / -6.0e-03 |
| 20 | shear+bulk+diff | +3.2e-03 / +2.3e-03 / +5.9e-03 / +1.8e-03 | -9.9e-04 / -3.1e-03 / +4.8e-03 / -4.9e-03 |

### (2+1)D

| EoS | viscosity | beastMode 1 | beastMode 3 |
|---|---|---|---|
| 0 | ideal | +3.6e-04 / +3.0e-04 / +1.1e-03 | -4.3e-03 / -4.6e-03 / +8.3e-04 |
| 0 | shear | +4.0e-04 / +3.2e-04 / +9.9e-04 | +4.0e-04 / +3.2e-04 / +9.9e-04 |
| 0 | shear+bulk | +4.0e-04 / +3.2e-04 / +9.8e-04 | +4.0e-04 / +3.2e-04 / +9.8e-04 |
| 9 | ideal | +3.3e-04 / +3.0e-04 / +1.0e-03 | -2.3e-03 / -2.3e-03 / -2.4e-03 |
| 9 | shear | +3.6e-04 / +3.2e-04 / +8.8e-04 | -1.9e-03 / -1.9e-03 / -1.2e-03 |
| 9 | shear+bulk | +3.4e-04 / +2.9e-04 / +1.1e-03 | -2.3e-03 / -2.2e-03 / -2.2e-03 |
| 14 | ideal | +2.7e-04 / +2.3e-04 / +1.1e-03 / -3.2e-04 | -4.4e-03 / -5.0e-03 / -5.5e-03 / -1.0e-02 |
| 14 | shear | +3.9e-04 / +3.4e-04 / +1.1e-03 / +2.9e-04 | -2.0e-03 / -2.0e-03 / -1.0e-03 / -1.6e-03 |
| 14 | shear+bulk | +3.4e-04 / +2.9e-04 / +1.1e-03 / +2.5e-04 | -2.3e-03 / -2.2e-03 / -2.2e-03 / -1.7e-03 |
| 14 | shear+bulk+diff | +3.2e-04 / +2.7e-04 / +1.0e-03 / +2.6e-04 | -2.3e-03 / -2.2e-03 / -2.4e-03 / -1.8e-03 |
| 20 | ideal | +2.9e-04 / +2.7e-04 / +1.2e-03 / -1.8e-04 | -2.7e-03 / -2.7e-03 / -2.7e-03 / -2.5e-03 |
| 20 | shear | +4.0e-04 / +3.4e-04 / +1.1e-03 / +2.9e-04 | -2.3e-03 / -2.3e-03 / -1.1e-03 / -1.9e-03 |
| 20 | shear+bulk | +3.5e-04 / +3.0e-04 / +1.1e-03 / +2.5e-04 | -2.4e-03 / -2.3e-03 / -2.2e-03 / -1.8e-03 |
| 20 | shear+bulk+diff | +3.2e-04 / +2.7e-04 / +1.1e-03 / +2.9e-04 | -2.5e-03 / -2.4e-03 / -2.4e-03 / -2.2e-03 |

- beastMode 1 changes the surface integrals by 0.1–0.5% (E, S) and 0.5–1.2% (⟨u_⊥⟩) in (3+1)D, and by 0.03–0.1% in (2+1)D. This comes mostly from face_reconstruction 1 and the 1e-6 Newton tolerance.
- beastMode 3 changes them by up to 0.9% (S) and 1.3% (N_B).
- In the dense (2+1)D long run (coarsened twice, to dx = 0.4 fm), beastMode 3 is 12× faster than beastMode 0 (5 s against 65 s), but surface E is −3.9%, S −4.2% and N_B −3.2%.

## Long tests

| run | exit | wall [s] | tau_end | max|E/E0-1| | max|dPz|/E0 | max|N_B/N_B0-1| | frozen out | warnings |
|---|---|---|---|---|---|---|---|---|
| long_hires_d3_eos20_shear_bulk_diff | 0 | 778.5 | 7.20 | 1.8e-06 | 8.8e-08 | 1.2e-05 | True | 1 |
| long_hires_d3_eos14_shear_bulk_diff | 0 | 517.1 | 7.20 | 1.7e-06 | 1.1e-07 | 6.3e-06 | True | 1 |
| long_hires_d3_eos9_shear_bulk | 0 | 326.7 | 7.10 | 1.6e-06 | 9.4e-08 | – | True | 1 |
| long_cool_d2_eos20_shear_bulk_diff | 0 | 7.6 | – | – | – | – | True | 1 |
| long_cool_d2_eos9_shear_bulk | 0 | 3.1 | – | – | – | – | True | 1 |
| long_hires_d3_eos20_shear_bulk_diff_bm1 | 0 | 428.2 | – | – | – | – | True | 1 |
| long_hires_d3_eos20_shear_bulk_diff_bm3 | 0 | 300.9 | – | – | – | – | True | 1 |
| long_hot_d2_eos9_shear_bulk | 0 | 29.2 | – | – | – | – | True | 1 |
| long_hot_d2_eos14_shear_bulk_diff | 0 | 44.1 | – | – | – | – | True | 1 |
| long_hot_d2_eos20_shear_bulk_diff | 0 | 65.1 | – | – | – | – | True | 1 |
| long_hot_d2_eos20_shear_bulk_diff_bm3 | 0 | 5.2 | – | – | – | – | True | 1 |
| long_hot_d3_eos20_shear_bulk_diff | 0 | 5878.6 | 19.80 | 3.6e-06 | 5.6e-07 | 3.1e-06 | True | 0 |
| long_hot_d3_eos9_shear_bulk | 0 | 1546.1 | 19.70 | 3.4e-06 | 5.4e-07 | – | True | 0 |
| long_hot_d3_eos20_shear_bulk_diff_bm3 | 1 | 111.1 | 7.04 | 1.3e+01 | 1.3e+01 | 6.9e-01 | False | 0 |

- **Results:**
  - The dense (3+1)D runs at beastMode 0 lasted 98 min (EoS 20) and 26 min (EoS 9), to τ ≈ 20 fm/c, with |ΔE/E0| ≤ 3.6 × 10⁻⁶ and |ΔN_B/N_B0| = 3.1 × 10⁻⁶. They produced no warnings.
  - The "cool" (2+1)D runs froze out within seconds, so the dense (2+1)D runs (long_hot_d2) were added.
- **Failure: long_hot_d3_eos20_shear_bulk_diff_bm3, a pre-existing beastMode 3 bug.**
  - Before the coarsening, the energy is conserved to the printed digit through 3 grid enlargements (61² → 121² at dx = 0.2 fm).
  - At step 192 (τ = 4.9 fm/c) Δτ doubles to 0.04 and the grid coarsens to dx = 0.4 fm. E jumps by 1 × 10⁻³ and keeps drifting: 38227 → 38769 GeV by τ = 5.4.
  - After that it grows exponentially (4 × 10⁵ GeV at τ = 6.8) until N_B < 0 at τ = 7.04 and the run aborts. The grid keeps being enlarged, up to 44 fm, because the growing energy reaches the edges.
- **Narrowing runs (same input, failure time in fm/c):**

| variant | result |
|---|---|
| current tree (OpenMP, no diagnostics) | blows up, dies at τ = 8.8 |
| 201ccad | blows up, dies at τ = 8.9 |
| **227d813 (before any optimisation commit)** | **blows up, dies at τ = 10.3** |
| current, face_reconstruction 0 | blows up, dies at τ = 10.6 |
| current, EoS 9 (no baryons) | blows up after the coarsening |
| current, gridPadding 0 (coarsening only, no enlargement) | blows up after the coarsening, τ ≈ 7.7 |
| **current, beastMode 0, dx = 0.4 fm and Δτ = 0.04 from τ0** | **stable**, E conserved to 1e-5 past τ = 16.8 |
| current, beastMode 0, dx = 0.4 fm, Δτ = 0.02 | stable to τ = 8.4 (stopped) |

  - The bug does not depend on the EoS, on baryons, on the face reconstruction or on the grid enlargement.
  - It exists before the optimisation branch.
  - The coarse resolution is stable by itself, so the culprit is the mid-run switch: `regridTransverse(..., coarsen=true)` together with `extrapolatePrevForNewDtau` in `evolve.cpp` (from 227d813 / 1785651).
- **Suspects to check:**
  - The previous slices used in the time derivatives of the viscous sources. `fpPrev` is coarsened and then extrapolated linearly to τ − 2Δτ, which on a freshly filtered grid can make ∂_τ u inconsistent.
  - The dissipative fields after the projection in the remap.
  - The freeze-out prev/curr slices.
- The short and 121×121×61 beastMode 3 runs survived because they froze out soon after the coarsening.

## Thread determinism

- 84 configurations were run with 1 and 10 threads: every output .dat file is bit-identical.
- The surface flux sums differ by at most 6 × 10⁻¹⁴ relative, which is summation order.

## Caveats

- Timings are wall-clock for one run each. The machine was otherwise idle, and earlier back-to-back repeats varied by about 1–2%.
- Surface N_B cannot be compared with the initial N_B: cells below e_fo at τ0 never cross the surface. Surface fluxes are used only to compare runs.
- beastMode 3 conservation and the narrowing runs used a measurement build (current tree plus the conservation check at beastMode > 0). The failing case was reproduced with the regular build as well.

## Files

In `~/Desktop/codes/MUSIC_benchmark_2026-11-22/` (the data stays outside the repository):
- `results.jsonl`: one line per run.
- `logs/`: the log, music_input and conservation file of every run.
- `nightly.py` and `report.py`: the driver and the table generator.
- `chain.log`: the phase timings.

## Follow-up 2026-11-22: beastMode 3 coarsening fix and beastMode 2 test

Committed as d966c76 on `optimization`. Raw data is in `~/Desktop/codes/MUSIC_benchmark_2026-11-22/bm3fix/` (`results.jsonl`, `fixbench.py`, `nopad.py`, logs).

- **Cause.** The conservative coarsening remap at the dilute fluid edge reconstructed a fast u from the averaged T^{τμ}. Projecting the averaged π^{μν} onto that u then created energy (one coarse cell went from e 0.08 to 1.2 with π ~ 270). That seeded an exponential instability.
- **Fix** (`coarseGrainAndEnlargeGrid`). If the reconstructed cell is unphysical, its projected dissipative part exceeds 2× the parents' magnitude, or its u^τ exceeds both the parents' maximum and that of the ideal-part average, the cell keeps only the conserving average of the ideal part. Ideal runs never use this fallback. Viscous runs use it for about 100 cells per field, and the dense run for 300–900.
- **Diagnostic fix** (`check_conservation_law`). On the step where Δτ changes, the boundary outflow is now integrated over the true interval. beastMode 0 output is bit-identical.
- **Regression test.** "Evolve: coarsening a dilute, viscous edge gives physical cells" in `bugfix_unittest.cpp` fails on the old code and passes now.

**Results** (90 runs plus 1 extra; all exit 0, all freeze out, no NaN or error, N_B ≥ 0):

| run | before | after |
|---|---|---|
| long_hot_d3 EoS 20 shear+bulk+diff, beastMode 3 | crashed at τ = 7.04 | freezes out at τ = 19.6, max \|dE\| 3.1e-4; surface E −1.0 %, S −0.5 % vs beastMode 0 |
| long_hot_d3 EoS 9 shear+bulk, beastMode 3 | – | freezes out, max \|dE\| 3.3e-4 |
| same EoS 20 run, fixed 40 fm grid (no enlargement) | – | max \|dE\| 4.9e-5, max \|dN_B\| 4.1e-5 |
| (3+1)D matrix, beastMode 3, conservation | 1e-4 to 1e-3 | 4e-6 to 3e-5 |
| (3+1)D and (2+1)D matrix, beastMode 3, surface fluxes | – | identical to the overnight run (≤ 1e-10) |

- **The remaining 3e-4 in the dense beastMode 3 runs** comes from grid enlargement, which predates this fix. It shows up as single-step jumps at the enlargements at τ = 11.1 and 16.1. Padding cells get Gaussian-damped copies of the edge fluid, which is above e_FO there. On that step the outflow is also counted at the new edge. Without enlargement only the ~3e-5 jumps at the Δτ switches remain; they come from the lagged viscous term of the scheme.
- **beastMode 2 is not affected.** It doubles Δτ but never remaps the grid. It was tested with Δτ 0.005 in (3+1)D and 0.0025 in (2+1)D, so the doubling actually triggers. All runs are stable. Conservation in (3+1)D is 3e-6 (ideal) to 6e-5. Surface fluxes are within +0.3 to +0.6 % for E and S of beastMode 0 in (3+1)D, and within 0.05 % in (2+1)D. The dense long run (101×101×31, Δτ 0.005 → 0.02) gives max \|dE\| 4.5e-5.

## Follow-up 2026-11-22 (2): conservative grid enlargement

Committed as 5f82266 on `optimization`. Raw data is in `~/Desktop/codes/MUSIC_benchmark_2026-11-22/padfix/`.

- **Change.** Padding cells used to be full, damped copies of the nearest edge cell, added on top of it. Now the edge cell and its padding copies *share* the edge cell's T^{τμ} and J^τ, with weights damp / Σ damp. They go through the same remap as coarse cells: u is kept, π^{μν} and q^μ are projected, and the ideal-part fallback still applies. The conservation diagnostic now integrates the outflow at a remap step with the rate through the old grid faces, recorded just before the remap (`Cell_info::noteRegrid`).
- **The old padding added real content when the grid was enlarged at the start.** The long beastMode 3 tests start from 10–12 fm grids that cut off the initial condition. The old code's first enlargement (before the first time step and the first conservation record) added 2.2 % of the energy (E(τ0) 38227 against 37405 GeV on the grid). That is why it looked closer to beastMode 0. The full initial condition on a 40 fm grid is 39041 GeV.
- **Mid-run enlargement, from a 20 fm grid that holds the initial condition** ((3+1)D EoS 20 shear+bulk+diff, 5 enlargements): the surface fluxes of the new and old code agree to 5e-5 to 1e-4. Max |dE| is 4.6e-5. The only jumps left are at the Δτ switches (3e-5 and 1.4e-5). In (2+1)D at 20 fm no enlargement triggered.
- **Conservation** in the beastMode 3 matrix: ideal cases 6–9 × 10⁻⁷ (was about 8e-6); viscous cases 2–3 × 10⁻⁵. Dense long runs: 5e-5 (was 3e-4).
- **Matrix surface fluxes** are unchanged from before (≤ 2e-6 relative). Every run froze out cleanly.
- **Integration tests:** 18 passed. Only the 12 output files of the beastMode 3 test changed, by about 0.1 %, because its first enlargement no longer adds content.
