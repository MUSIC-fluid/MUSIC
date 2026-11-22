# Conservation-law treatments on `optimization` compared with `origin/main`

2026-11-22. `origin/main` is febb2ab (2026-07-24), the merge base. The `optimization` branch (2e59337) is 145 commits ahead and 0 behind. This report covers the commits that change what the scheme conserves, or how it is measured. The pure performance commits are left out, but several of them had to keep a property described here (for example one flux per face), and those are noted.

## Summary

MUSIC evolves τT^{τμ} and τJ^τ with a Kurganov–Tadmor (KT) finite-volume update and RK2 (Heun) time stepping. Exact conservation then needs three things:

1. the two cells that share a face use the same flux;
2. the time integration of the viscous part telescopes;
3. every other change of the fields leaves the integrals unchanged (reconstruction failures, regulation, remaps of the grid).

On `origin/main` all three were violated somewhere, and its diagnostic could not tell. The branch fixes all three and rebuilds the diagnostic, so that it measures the quantity the scheme conserves.

| area | `origin/main` | `optimization` | commits |
|---|---|---|---|
| KT face flux | each cell reconstructed both sides of its faces with *its own* state as guess and fallback; a rejected solve gave the two neighbours different fluxes | each half-way state uses the cell on its own side as guess and fallback; since c4065ec the flux is computed once per face and differenced | 149d049, c4065ec |
| superluminal T^{τμ} | the cell was reverted to its previous state (content lost or created) | T^{ττ} and J^τ kept; only T^{τx}, T^{τy} scaled to \|M\| = 0.999 T^{ττ}, so E and P_z are exact | 149d049, cb5de43 |
| regulated dilute cells | n_B set to 0 | n_B, n_Q, n_S kept (= J^τ) | 6d6a2cc, 6f22cfd |
| viscous part in time | the two RK stages left ½ τ_{n+2}(V_{n+1} − V\*) per step, an O(Δτ) drift | local fixed point in stage 2 makes the average telescope (tolerance 1e-8) | d4801e0 → e624038, 201ccad |
| RK stage 2 | derivatives and transport coefficients at τ_n | at its own τ_{n+1} and state | 9bd9abd |
| bulk Π_chem | entered T^{μν} only with `turn_on_bulk = 1` | always with bulk or chemical bulk | 85cdc2c |
| ill-posed sources | κ_B < 0 (anti-diffusion), noisy D(μ_B/T), unreliable dilute 4D EoS | κ_B ≥ 0, D(μ_B/T) from the ideal equations of motion, regularised EOS_4D, Du from the ideal EoM in dilute cells | a7e141e, 7b6f772, 2a27154, 0b4ed3a |
| conserved charges | N_B only | N_B, N_Q, N_S, with sources, regulation and diagnostics | 4c194c0, 6f22cfd |
| beastMode Δτ doubling | backward time differences spanned the old Δτ | the previous slice is extrapolated to τ − Δτ_new | 1785651 |
| beastMode coarsening | e and u averaged (T^{ττ} 14/3 → 2 in the unit test) | filtered T^{τμ}, J^τ remapped and reconstructed; ideal-part fallback at the dilute edge | e5078e9, a24c3d4, d966c76 |
| beastMode enlargement | padding cells were added copies of the edge cells (+2 % of E at the start, ~3e-4 per enlargement) | the edge cell and its copies share its content | 5f82266 |
| diagnostic | outflow = change of the edge-cell content; viscous part at one lagged level only | time-integrated flux through the grid faces; the two-level invariant of the scheme; correct intervals at Δτ changes and remaps | 2cd2385, d4801e0, d966c76, 5f82266 |

Measured end to end on the same events (lumpy (3+1)D TATB, 101 × 101 × 31 on 24 fm, EoS 9 and 14, beastMode 0; see [Measurements](#measurements)):

| event | max \|E/E0 − 1\| | max \|ΔP_z\|/E0 | max \|N_B/N_B0 − 1\| |
|---|---|---|---|
| EoS 9 ideal | 8.1e-2 → **2.9e-7** | 2.6e-3 → 1.3e-7 | – |
| EoS 14 ideal | 7.8e-2 → **4.4e-7** | 2.1e-3 → 2.0e-7 | 3.0e-2 → **1.0e-6** |
| EoS 9 shear+bulk | 2.7e-2\* → **5.1e-6** | 2.2e-3 → 1.7e-7 | – |
| EoS 14 shear+bulk+diff | `origin/main` aborts at τ = 3.3 (N_B < 0) → **7.5e-6** | 1.4e-3 → 8.5e-7 | 1.4e-1 (then N_B < 0) → **5.3e-6** |

Each cell gives `origin/main` → branch. \* For viscous runs `origin/main`'s diagnostic mixes time levels (section 6), so this is an upper bound on the real drift. `origin/main` *loses* energy (E at freeze-out −8.1 % ideal, −2.7 % viscous).

## 1. Spatial fluxes (KT)

### 1.1 One flux per face (149d049, c4065ec)

The KT update of cell i is the difference of the fluxes through its faces i ± ½. The integrals telescope only if cell i and cell i + 1 use the same flux at i + ½.

On `origin/main`, `MakeDeltaQI` reconstructed all four half-way states of cell i (the left and right side of both faces) with cell i as the guess for the velocity solve. When a solve was rejected (it reverted to the guess, or u^τ moved by more than 100× from it), the result depended on the guess, so cell i + 1 computed a different flux for the same face. That broke conservation in fast, dilute regions.

- **149d049.** Each half-way state uses the cell on its own side of the face as the guess and fallback (q_{i+½,R} → cell i + 1, q_{i−½,L} → cell i − 1). Both neighbours then compute an identical flux. In a 3D ideal TATB run at 19.6 GeV the energy drift fell from −3.2 % to −6e-4, and the net-baryon drift from +13 % to the size of the initial-state artifact. At 200 GeV: −10 % / +52 % → −1.3 % / −0.4 %.
- **c4065ec.** The face fluxes are now computed once per face, in a pass over the faces, and each cell takes their difference. Equal fluxes on both sides are then guaranteed by construction, not by matching arithmetic.
- **d1b3c5e, 0772bf4.** The same holds for the W fluxes (`computeWFluxDivergence`, then `computeWDivergences`).

### 1.2 Primitive-variable limiting at the faces (816e6d5, optional)

`face_reconstruction = 1` limits e, n_B, n_Q, n_S and u^i at the faces instead of τT^{τμ} and τJ^τ. The flux is still one per face, so the update stays exactly conservative. Only the face states change.

- It is the default with beastMode > 0. Mode 0 is unchanged and bit-identical.
- Against the Gubser solutions its L1 errors are 2–5× smaller (ideal, with n_B) and about 2× smaller (viscous).
- The freeze-out balance of the ideal 3D fireball (energy through the surface against the initial energy) is 2.7e-3 at dx = 0.3 fm and 7e-4 at dx = 0.15 fm, against 1e-5 for mode 0. This is an accuracy difference, not a conservation one.

### 1.3 Boundary (361011e)

The zero-gradient boundary (clamped neighbour indices) is now stored as two ghost layers that copy the edge cells. The output is bit-identical. The diagnostic's outflow (section 6) uses the same convention: the flux through a boundary face is the flux of the edge cell.

## 2. Reconstruction failures and regulation

### 2.1 Superluminal T^{τμ} (149d049, cb5de43)

When |M| ≥ T^{ττ} (M = (T^{τx}, T^{τy}, T^{τη})), no (e, u) solves the reconstruction.

- **`origin/main`:** reverted the cell to its previous state, which does not conserve anything.
- **149d049:** keeps T^{ττ} and J^τ and caps |M| at 0.999 T^{ττ}.
- **cb5de43:** caps through T^{τx} and T^{τy} only. The Cartesian E = cosh η T^{ττ} + sinh η τT^{τη} and P_z then stay exact. All of M is scaled only if T^{τη} alone is superluminal (under-resolved η gradients).

Scaling T^{τη} as well had removed energy from every capped cell flowing outward in η. Dilute edge cells of lumpy events are capped every step. Lumpy ideal run: E drift 3.2e-4 → 1.6e-5, P_z 3.9e-5 → 9.9e-6.

### 2.2 Regulated cells keep their charges (6d6a2cc, 6f22cfd)

Where T^{ττ} drops below `abs_err`, `regulate_grid` replaces the cell.

- **`origin/main`:** set n_B = 0, which removed up to 2.3e-3 of N_B once dilute, baryon-rich cells became viscous.
- **Now:** a regulated cell keeps n_B = J^τ, and after the 4D merge also n_Q and n_S. The N_B drift of lumpy events falls from 8.9e-3 to 1.7e-4 (neos_2) and from 5.5e-3 to 4.3e-4 (neos_bqs).

### 2.3 Charges in nearly empty cells

- The EoS 0 ideal gas is now a consistent free QGP (6d6a2cc): s > 0 at any n_B, and the state solve takes n_B into account.
- EOS_4D (EoS 20) is regularised below T̃ = 40 MeV, with consistent gradients (2a27154). This alone improved energy conservation of smooth TATB runs from 2.4e-5 to ≤ 1.5e-7, and made viscous EoS 20 string runs survive.

## 3. Time integration of the viscous part

### 3.1 The telescoping condition

The viscous part V^α = W^{τα} + Π(g^{τα} + u^τ u^α) (and q^τ for J^τ) enters each RK stage through the backward difference (τ + Δτ)V_stage − τ V_prev in `Diss::MakeWSource`:

- the first stage subtracts τ_{n+1}V_n − τ_n V_{n−1};
- the second subtracts τ_{n+2}V\* − τ_{n+1}V_n, with V\* the predictor.

The RK average telescopes only if V\* = V_{n+1}. The remainder ½ τ_{n+2}(V_{n+1} − V\*) is O(Δτ²) per step, so O(Δτ) over a run. It is largest where V changes quickly, for example Π right after string deposition, where it reached 0.5 % of E.

### 3.2 Deferred correction (d4801e0), then a local fixed point (e624038, 201ccad)

- **d4801e0** first removed the remainder in the next step's first stage. That correction had to be skipped where V is not smooth in time, because it made IPGlasma 2D unstable. 3D viscous TATB: E drift 2.8e-4 → 2.0e-6 at Δτ = 0.02, now about Δτ².
- **e624038** replaced it with `Advance::stage2ViscousFixedPoint`. The remainder involves only the cell itself. So the stage-2 result is shifted by s, reconstructed, and W updated again (fluxes and gradients don't depend on s), until s = ½ τ_{n+2}(V_{n+1}(s) − V\*).
  - Full steps converge in one or two updates in almost all cells.
  - The step is halved where the miss grows (stiff diffusion, Π ≈ −(e + P)).
  - A cell that does not converge in 16 updates keeps the visited shift with the smallest miss (at worst s = 0, the predictor).
  - Smooth-event E drift fell from 3–9e-6 to 1–2e-7; N_B with diffusion by about 100×.
- **201ccad** adds an early exit when stage 2 is already consistent, and a depth-1 Anderson (secant) step. Corrections per call fell from 1.53 to 1.22; results change at the 1e-8 tolerance.

The scheme then conserves exactly the two-level quantity

  τ_n T_ideal^{τμ}(n) + ½ [τ_{n+1} V(n) + τ_n V(n−1)],

which is what the diagnostic now reports (section 6).

### 3.3 Second stage at its own time (9bd9abd)

On `origin/main` the second RK stage evaluated θ, Du, σ, the Christoffel terms, τΔη and the relaxation transport coefficients at τ_n, although its state is at τ_{n+1}. That made the time integration first order. In uniform viscous Bjorken flow T changed by 6e-4 between Δτ = 0.02 and 0.01; now it is the same to 7 digits from 0.02 to 0.0025. In lumpy runs this halved the E drift.

### 3.4 Π_chem in T^{μν} (85cdc2c)

`MakeWSource` added Π only when `turn_on_bulk = 1`, so `turn_on_bulk_chem = 1` with `turn_on_bulk = 0` (including the profile 9 initial bulk pressure) never fed back into T^{μν}. Now the total Π_bulk + Π_chem enters whenever either is on. The fixed point, the remaps and the diagnostic use the same sum.

### 3.5 Inputs to the viscous terms that drove drifts or blow-ups

These are not conservation fixes in themselves. They made V smooth or well-posed, so that 3.2 can be applied and the drift converges in Δτ.

- **Du from the ideal equations of motion in dilute cells (0b4ed3a).** Below `Du_from_ideal_EoM_e_dilute` (0.003 GeV/fm³, far below freeze-out), ∂_τ u^i comes from the ideal EoM, blended in log e up to 5× that. Lumpy viscous E drift: 2.8e-4 → 5e-6; lumpy diffusion N_B: 7.4e-4 → 8e-5. Flow observables shift by ≤ 1 %. Gubser at η/s = 0.2: u^x error 0.58 % against 0.35 % with the backward difference everywhere. Set `Du_from_ideal_EoM = 0` for the old behaviour.
- **κ_B ≥ 0 (a7e141e).** The kinetic-theory κ_B turns negative in cold, baryon-rich cells. Anti-diffusion is ill-posed: N_B swung between 2 and 609 (initial 76) and the run aborted. Smooth events are untouched.
- **D(μ_B/T) from the ideal EoM (7b6f772).** The backward-difference form carried time-step-scale noise where μ_B/T is large, and the noise grew as Δτ → 0. Lumpy-event N_B drift: 4.2e-3 → 8.3e-4 at Δτ = 0.02, now falling with Δτ.

## 4. Additional conserved charges (4c194c0 … 6f22cfd)

n_Q and n_S are evolved alongside n_B, with the 4D EoS (EoS 20). They are conserved by the same KT update. The regulation (2.2), Du, the beastMode remaps and the diagnostic carry all three. N_Q and N_S, with their boundary outflow, are columns 7 and 8 of `global_conservation_laws.dat`. There are no diffusion currents for Q and S.

## 5. beastMode grid changes

### 5.1 Δτ doubling (1785651)

When beastMode doubles Δτ, the backward time differences (curr − prev)/Δτ would mix the old and new step. The previous slice is linearly extrapolated to τ − Δτ_new, ghost cells included. The diagnostic then sees a jump of about 3e-5 at each switch. It comes from the lagged viscous invariant (V(n−1) is extrapolated, not computed), and is accepted.

### 5.2 Conservative coarsening (e5078e9, a24c3d4, d966c76)

Coarsening doubles dx and dy.

- **`origin/main`:** averaged e, n_B and u. In the unit test, two cells e = 1, u^x = ±1 became a rest cell with T^{ττ} = 2 instead of 14/3.
- **e5078e9:** filters the *total* T^{τμ} and J^τ with a 1/4–1/2–1/4 stencil, reconstructs e, n and u from them, and projects π^{μν}, q^μ onto the new u (transverse, traceless). It iterates, since the projected W depends on u.
- **a24c3d4:** drops out-of-grid neighbours instead of clamping them, so odd grids conserve their content; single-cell directions are not coarsened.
- **d966c76:** at the dilute fluid edge, projecting the averaged π onto a fast reconstructed u created energy. One coarse cell went from e 0.08 to 1.2 with |π| ~ 270, which seeded an instability that crashed dense runs. A coarse cell that is unphysical keeps the conserving average of the ideal part only. Unphysical means any of:
  - its projected dissipative part exceeds twice its parents';
  - its u^τ exceeds both its parents' and that of the ideal average.

  Ideal runs never use this fallback; viscous runs use it for about 100 cells per field.

### 5.3 Conservative enlargement (5f82266)

- **`origin/main` and the branch before 5f82266:** padding cells were Gaussian-damped copies of the nearest edge cell, added on top of it.
  - A grid that cut off the initial condition gained 2.2 % of its energy at the first enlargement, which happens before the first conservation record.
  - Each mid-run enlargement in a dense run added about 3e-4.
- **Now:** the edge cell and its copies share its T^{τμ} and J^τ, with weights damp / Σ damp, and go through the remap of 5.2.
  - From a grid that holds the initial condition, the surface fluxes agree with the old padding to 5e-5–1e-4.
  - A run whose initial grid cuts off the initial state now keeps that deficit instead of inventing content. Start beastMode 3 on a grid that holds the initial condition.

## 6. The conservation diagnostic

`Cell_info::check_conservation_law` writes `global_conservation_laws.dat`: E, P_x, P_y, P_z, N_B (and N_Q, N_S).

| | `origin/main` | `optimization` |
|---|---|---|
| content | ideal T at n, plus W and Π at the lagged level only | the two-level invariant of 3.2: τ_n T_ideal(n) + ½[τ_{n+1}V(n) + τ_n V(n−1)], with Π_bulk + Π_chem |
| outflow | the change of the content of the edge cells, counted as outflow, so deposits, inflow and p dV work there were double-counted (profile 17 looked 41 % too energetic) | trapezoidal time integral of the flux of T^{μν} and J^μ through the x, y and η grid faces, viscous and diffusion parts included (2cd2385). With fluid leaving an 8 fm box E is conserved to 8e-6 instead of 19 % |
| Δτ change | – | the interval is the true elapsed τ, using the rate stored at the previous record (d966c76) |
| remap | – | the rate through the old grid faces, recorded just before the remap (`noteRegrid`, 5f82266) |
| N_B check | exit unless 0 ≤ N_B < 500 | the same, now applied to the corrected total, and N_Q, N_S reported |
| summation | OpenMP reduction (thread-count dependent at round-off) | row sums added in a fixed order (201ccad) |

In the repository build the diagnostic runs only with beastMode 0. The beastMode measurements below used a diagnostic build that also runs it at beastMode > 0.

## Measurements

### This report: `origin/main` against the branch on the same events

The setup:
- Lumpy (3+1)D TATB events, 101 × 101 × 31 cells on 24 fm, η grid size 8, Δτ 0.02, beastMode 0, 10 threads, evolved to freeze-out at 0.3 GeV/fm³.
- On this grid the freeze-out surface stays within r = 5.4 fm, so nothing reaches the edge. Both diagnostics then measure the content of the grid, and the ideal rows are directly comparable.
- For viscous runs `origin/main`'s diagnostic adds W and Π at the lagged level only. Part of its number is that mismatch rather than real non-conservation.
- `origin/main` deletes `*.dat` in its run directory at start-up, so its TA/TB inputs were renamed to `.in`. The initial states are identical.

| run | `origin/main` | branch |
|---|---|---|
| E(τ0) [GeV] | 742.736 | 742.736 |
| EoS 9 ideal: E at freeze-out (τ = 6.65) | −8.1 % | −2.2e-7 |
| same on a 36 fm grid (151², `origin/main`) | −8.1 % (8.0885e-2 against 8.0886e-2) | – |
| EoS 14 ideal: E / N_B at freeze-out | −7.8 % / −3.0 % | −3.4e-7 / +1.0e-6 |
| EoS 9 shear+bulk: E at freeze-out (τ = 7.15) | −2.7 %\* | +4.8e-6 |
| EoS 14 shear+bulk+diff | N_B grows from 80.3 to 90.8 by τ = 3.3, then −5.8 (ρ_B max 358 fm⁻³): error exit, status 1 | froze out at τ = 7.25; E +7.3e-6, N_B −5.3e-6 |

- **The `origin/main` loss is in the scheme, not at the boundary.** Its loss is the same on a 36 fm grid as on 24 fm. Most of it comes from section 1.1 (different fluxes at a face after a rejected solve) and section 2.1 (reverted superluminal cells). In an ideal run no viscous term is involved.
- **The `origin/main` diffusion failure** is the κ_B < 0 anti-diffusion of section 3.5 (a7e141e).
- **The branch's surface fluxes:** E 543.5 / 541.3 / 575.8 / 569.6 GeV and S 1840 / 1934 / 2119 / 2171 for the four events. `origin/main` writes an older surface format, which the analysis script does not read.
- Scripts and results: `vsmain.py`, `vsmain2.py` and `vsmain3.py` under `~/Desktop/codes/MUSIC_benchmark_2026-11-22/vsmain/`.

### From the 2026-11-22 benchmark ([benchmark_2026-11-22.md](benchmark_2026-11-22.md))

These are current-code numbers on 61 × 61 × 31 grids, where fluid does leave the grid:

| | max \|E/E0 − 1\| | max \|N_B/N_B0 − 1\| |
|---|---|---|
| beastMode 0 | 2–10e-6 | 1.6–3.5e-5 |
| beastMode 1 | similar to beastMode 0 | similar to beastMode 0 |
| beastMode 3, after 5f82266 | 6e-7–3e-5 | ≤ 7e-5 |
| beastMode 3, dense long runs | 5e-5 | – |

|ΔP_z|/E0 is ≤ 1e-6 at beastMode 0.

## Remaining known errors

- **Δτ switches** (beastMode 2 and 3): about 3e-5 of E per switch, from the extrapolated lagged viscous level (5.1).
- **The stage-2 fixed point** stops at a relative miss of 1e-8. A cell that does not converge in 16 updates keeps its best shift, so a residual of that cell's remainder stays.
- **The ideal-part fallback of the remap** (5.2, 5.3) drops the viscous part of the content of the cells where it is used.
- **Boundary:** the zero-gradient boundary lets fluid leave. The diagnostic accounts for it, but the content on the grid is not conserved by itself.
- **Small cell-local changes:** density floors (`small_eps` in padding cells) and the |M| cap at 0.999 T^{ττ} change T^{τx}, T^{τy} by tiny amounts in dilute cells.
