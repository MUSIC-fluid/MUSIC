# Performance optimizations on `optimization` compared with `chun_dev`

2026-11-22. `chun_dev` is `origin/chun_dev` = 6e19c80 (2026-11-20). `optimization` (d2964a1) contains all of it and adds 107 commits:

- **The `main_4D` merge** (6f22cfd): the 4D EoS and the evolution of n_Q and n_S. It also brings older speed-ups of the velocity solve from the `greg_dev` history.
- **Physics and conservation work** (e624038, 2a27154, 4bd855a, the beastMode fixes). It costs some time; see [Costs](#costs-added-on-the-branch).
- **14 performance commits**, 40c37ff … 0772bf4, on top of the 227d813 merge. They include the optional `face_reconstruction = 1`.

227d813 is therefore the natural baseline for the performance commits. `chun_dev` is the baseline for the branch as a whole.

## Summary

Wall time to freeze-out of one lumpy (3+1)D TATB event:
- 61 × 61 × 31 cells on 12 fm, Δτ = 0.02, η/s = 0.12, freeze-out at 0.3 GeV/fm³;
- g++-15 -O3 with OpenMP on the same machine, run back to back.
- beastMode 0 is the default (`face_reconstruction = 0`); beastMode 1 uses `face_reconstruction = 1`.

| case | threads | `chun_dev` | 227d813 | branch, beastMode 0 | speed-up vs `chun_dev` (vs 227d813) | branch, beastMode 1 | speed-up vs `chun_dev` |
|---|---|---|---|---|---|---|---|
| EoS 9, shear + bulk | 1 | 221.2 s | 277.9 s | **135.6 s** | 1.63× (2.05×) | 66.2 s | 3.34× |
| EoS 14, shear + bulk + diffusion | 1 | 316.4 s | 405.3 s | **213.2 s** | 1.48× (1.90×) | 118.6 s | 2.67× |
| EoS 0, ideal | 10 | 13.9 s | 14.7 s | **8.2 s** | 1.70× (1.79×) | – | – |
| EoS 9, shear + bulk | 10 | 34.7 s | 45.8 s | **21.8 s** | 1.59× (2.10×) | 12.1 s | 2.87× |
| EoS 14, shear + bulk + diffusion | 10 | 53.5 s | 70.2 s | **35.3 s** | 1.52× (1.99×) | 21.6 s | 2.48× |
| EoS 20, shear + bulk + diffusion | 10 | not available | 124.3 s | **52.4 s** | – (2.37×) | 29.7 s | 4.19× |

- **227d813 is 1.06–1.32× slower than `chun_dev`.** It carries the stage-2 viscous fixed point (about 1.2×) and the four-density EoS interface. The performance commits more than recover that.
- **At the same settings the branch is 1.5–1.7× faster than `chun_dev`** and 1.8–2.4× faster than 227d813.
- **beastMode 1 is another 1.6–2.0× faster.** Its `face_reconstruction = 1` changes the face states, so its results are not bit-identical.
- **The optimizations keep the results.** At beastMode 0 the branch reproduces the surface energy and entropy of 227d813 to every printed digit for EoS 0, 9 and 14. For EoS 20 they differ by 1.8e-4 in E and 1e-3 in S, from the physics correction to dP/de of the 4D EoS (4bd855a). The energy and N_B conservation are unchanged. `chun_dev` differs physically: it lacks the stage-2 fixed point and the other physics changes, and its energy drift is about 2× larger (1.3–1.5e-5 against 7–10e-6).

## Ground rules the optimizations followed

- **Bit-identical, or stated.** Each commit was checked against its parent, almost always for bit-identical output (clang; g++ at round-off where FMA contraction changed). The exceptions are called out below. They change results only at a solver tolerance (≤ 1e-8), or by an opt-in mode (`face_reconstruction = 1`).
- **Conservative by construction.** Fluxes are computed once per face and differenced, never once per cell side (see [the conservation report](conservation_treatments_vs_origin_main.md), section 1.1).
- **Thread-count independent.** Every reduction that writes output sums row partial sums in a fixed order. The integration suite compares 1 and 4 threads, and the 2026-11-22 benchmark found 1- and 10-thread runs identical in every output file.
- **Tested.** The unit tests and the 18-test integration suite (760 output files compared bit by bit against the previous commit) ran after each step.

## The optimizations

Each speed-up was measured when the commit was made, on its own, against its parent: "3D shear+bulk" is the 61 × 61 × 31 benchmark above, and "EoS 20" the same with the 4D EoS.

### 1. Do each piece of work once

| commit | what was repeated | now | speed-up |
|---|---|---|---|
| 40c37ff | the spatial part of the W update (`Make_uWRHS`) recomputed in every stage-2 fixed-point iteration | once per cell and stage, passed to the iterations | 7 % |
| a7e96e5 | 13 EoS calls per cell (the four neighbours in each direction, plus the cell) | one pressure array per stage, read by the stencils | 2.4 % (EoS 20) |
| c4065ec | each KT interface reconstructed by both of its cells (12 velocity solves per cell and stage) | a `FaceFluxes` pass over the faces; cells take the difference | 21–32 % |
| d1b3c5e | the u^i W face flux built for both faces of every cell | `computeWFluxDivergence`, once per face | part of −4 % / −11 % (1 / 10 threads) |
| 201ccad | the face averages of W in the W source, per cell | once per face (`computeWSourceDivergence`, bit-identical to `MakeWSource`) | part of 201ccad |
| 0772bf4 | separate passes for the W flux and W source divergences | one fused per-line pass (x lines, 32-cell blocks along y and η) | part of 0772bf4 |
| 7015f66 | the EoS at the cell and its six neighbours for ∇P in `MakedU`; the thermal variables in every stage-2 iteration; T and μ_B in `MakeDTau` without diffusion; ghost updates of unchanged slices | pressure array reused; thermal variables looked up once per cell; lookups skipped when unused; ghosts updated only after remaps | 5–8 % (EoS 9), 7–14 % (EoS 20) |
| 0772bf4 | the stage-1 pressure recomputed for stage 2 | reused | part of 0772bf4 |
| 201ccad | separate sweeps for conservation, grid maxima and cosh/sinh per cell | one sweep; cosh/sinh tabulated per η slice; first-stage pressure reused | part of 201ccad |
| 201ccad | copies of the current cell in `FirstRKStepW` and every stage-2 correction | one shared copy; the next-slice cell skips the W and Π loads | part of 201ccad |

### 2. Memory layout

| commit | change | speed-up |
|---|---|---|
| 361011e | `Fields`: its 24 components move from separate `std::vector`s (and vectors of vectors) to one 64-byte-aligned structure-of-arrays buffer, with two ghost layers per direction instead of clamped neighbour indices. A direction with one cell gets stride 0. `-DMUSIC_CHECK_GHOSTS` checks every halo read against its edge cell. | 9–21 %; memory +19 % on 101 × 101 × 41, +8 % on 200 × 200 × 100 |
| b1bd88b | the thermodynamic variables in a fixed-size `ThermalVec` (`std::array<double, 15>`) instead of a `std::vector` built for every call; `Eigen::Matrix4d` instead of `MatrixXd` in the causality check | 4–8 % |
| 7609fb7 | EoS tables: one contiguous buffer per quantity (`EosTable`, with offsets for a later device copy) instead of `double***` with one heap block per n_B row; the 4D corners read without bounds-checked `.at()` | 6 % (EoS 20) |
| 0772bf4 | `map_2d_idx_to_1d`, `map_1d_idx_to_2d` inline in `util.h` | part of 0772bf4 |

### 3. Cheaper EoS evaluation

| commit | change | speed-up |
|---|---|---|
| 58be62b | `EOS_4D::FourDLInterp`: T̃⁵ differences of each table cell precomputed at load; T̃⁴ = (T̃²)²; the 16 corners read as a base index plus fixed strides | 6 % (EoS 20, 1 thread), bit-identical |
| 58be62b | `ReconstIt_shell` returns the EoS at its solution (`ReconstEos`): the face reconstruction uses `get_pressure_with_gradients_and_cs2`, and `MaxSpeed` reuses it instead of a second EoS call | round-off change |
| 0772bf4 | `EOS_4D::getThermalVariables` sets up the 4D interpolation cell once for p, T and the chemical potentials | part of 0772bf4 |
| 2a27154 | c_s² from dP/de and dP/dn; the `neos4d_cs` tables are no longer read. Mainly a physics fix (the tables gave 0.01 in cool baryon-rich cells) | – |
| 13e6c58, 53a26fe, 4661ece (`main_4D`) | c_s² not computed inside the velocity solve; df/dv and df/du⁰ without extra pressure calls; `solve_u0_hybrid` skipped for e < 1e-6 | from the 2024 `greg_dev` work |

### 4. Fewer and cheaper solver iterations

| commit | change | effect |
|---|---|---|
| 7653ee5 | the velocity solve starts with a secant iteration from the previous state; the bracketing hybrid solver only runs if that fails (0.3 % of solves with EoS 20) | pressure lookups per solve 7.6 → 5.4 (EoS 9), 8.4 → 6.2 (EoS 20); 3–9 %; results change ≤ 5e-10 |
| d1b3c5e | convergence checked before evaluating f(v_next), which only the next step needs | one lookup per solve; bit-identical |
| 0772bf4 | Newton start v = M / (T⁰⁰ + P_guess), with the cell or face pressure as the guess, and an analytic first slope 1 − 0.2 v² instead of a finite difference | differences at the 1e-14 solve tolerance |
| 201ccad | stage-2 viscous fixed point: early exit when already consistent, and a depth-1 Anderson (secant) step after the first correction | corrections per call 1.53 → 1.22; results change at the 1e-8 tolerance |

### 5. Parallelism

| commit | change | effect |
|---|---|---|
| d1b3c5e | `updateGhosts` in two parallel phases; `output_momentum_anisotropy_vs_etas` and `_vs_tau` parallel with fixed-order sums; `get_maximum_energy_density` puts every maximum in the OpenMP reduction (`rhoq_max` and `rhos_max` were a data race) | −11 % at 10 threads (EoS 9) |
| 58be62b | cell, face and W-flux loops use `schedule(dynamic, 64)` instead of `guided` | −7 % at 10 threads, −3 % at 4 threads; bit-identical |
| 201ccad | `check_conservation_law` sums rows in a fixed order | thread-count independent output |

Thread scaling on this grid with 10 threads: 6.0–6.5× at beastMode 0 and 5.5× at beastMode 1. The 2026-11-22 benchmark has the full 1/10-thread table.

### 6. Optional: primitive-variable face reconstruction (816e6d5)

`face_reconstruction = 1` limits e, n_B, n_Q, n_S and u^i at the KT faces and sets u^τ from the normalisation. That needs one EoS call per face side instead of a velocity root solve. The flux is still one per face, so it is exactly conservative.

- **Speed:** 34–42 % faster on one thread and 27 % at 10 threads when introduced. It is the default whenever beastMode > 0. `face_reconstruction = 0` stays the default otherwise and is bit-identical to before.
- **Accuracy:**
  - Against the Gubser solutions its L1 errors are 2–5× smaller (ideal with n_B) and about 2× smaller (viscous).
  - Both modes converge at the same order, and their difference falls 4–6× per refinement.
  - The ideal-fireball freeze-out balance is 2.7e-3 at dx = 0.3 fm (mode 0: 1e-5), falling to 7e-4 at dx = 0.15 fm.
- **Effect on surface results:** in the table above beastMode 1 shifts the surface E and S by 0.3–0.5 %.

## Costs added on the branch

These are not optimizations, but they explain why 227d813 is slower than `chun_dev`:

- **e624038, stage-2 viscous fixed point:** about 1.2× run time at introduction, for exact telescoping of the viscous part (energy drift of smooth events 3–9e-6 → 1–2e-7). 201ccad cut its cost by 20 %.
- **The `main_4D` merge:** every EoS call takes four densities (e, n_B, n_Q, n_S), and n_Q and n_S are evolved when enabled.
- **Ghost cells (361011e):** +8–19 % memory.
- **Face fluxes (c4065ec):** 21 doubles per cell.

## Where the time goes now

From a sampled profile at 201ccad, just before 0772bf4. 3D EoS 9 shear+bulk profiling setup (12 fm grid, 41 η cells, evolved for 0.6 fm/c), wall seconds:

| part | 1 thread, mode 0 | 1 thread, mode 1 | 10 threads, mode 0 | 10 threads, mode 1 |
|---|---|---|---|---|
| KT face fluxes (velocity solves) | 38.5 | 6.2 | 4.8 | 1.1 |
| per-cell loop (gradients, W relaxation, stage-2 fixed point) | 28.0 | 27.9 | 4.4 | 4.5 |
| W flux + source divergences | 7.2 | 7.2 | 2.0 | 2.0 |
| diagnostics, ghosts, pressure | 1.9 | 1.9 | 0.6 | 0.6 |
| total | 76.0 | 43.3 | 11.8 | 8.3 |

With EoS 20 in mode 1 the cell loop is 50 of 80 s, mostly EoS lookups in the W relaxation and the fixed point. In mode 0 the face velocity solves still dominate. In mode 1 the per-cell loop does, and that is where further work would go.

## Reproducing

- **Builds:** `git archive <rev>` built with `cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-15` (OpenMP), for 6e19c80, 227d813 and the branch.
- **Driver:** `vscd.py` (on the 2026-11-22 benchmark's `nightly.py`), with results in `~/Desktop/codes/MUSIC_benchmark_2026-11-22/vscd/`.
- **More data:** the 2026-11-22 benchmark ([benchmark_2026-11-22.md](benchmark_2026-11-22.md)) has the full EoS × viscosity × beastMode timing matrix in (3+1)D and (2+1)D, at 1 and 10 threads, with the `public_stable` baseline.
