# Skipping near-vacuum cells (`vacuumSkipEnergy`)

2026-11-25, branch `gpu` at 3be45fe plus the uncommitted changes described
here. Apple M4, 10 OpenMP threads, CPU hydro (gpu_backend 0) unless
noted.

## What it does

`vacuumSkipEnergy` (GeV/fm^3; 0 turns it off; the default depends on
`beastMode`, see Defaults) is the threshold below which a cell counts as vacuum. In each RK stage
`Advance::markVacuum` marks the vacuum cells of the current slice:

* **Faces.** An interior face whose four KT stencil cells (ILL, IL, IR,
  IRR) are all vacuum gets zero flux. `computeFaceFluxes` skips its
  reconstruction and the Riemann solve. The faces on the grid boundary are
  always computed, because `check_conservation_law` computes the boundary
  outflow from them.
* **Cells.** A cell is *frozen* if it is not at the grid edge and it and
  its neighbours at +-1 and +-2 along each axis are all vacuum. A frozen
  cell still gets `FirstRKStepT` (its fluxes are zero, so this is the
  conservative tau update of T^{tau mu} and J^tau). It skips the velocity
  derivatives, the dissipative sources and the W divergences. Its W^{mu nu}
  and Pi are only diluted by tau, `tau / (tau + dtau)` in the first stage
  (`Advance::frozenViscousUpdate`). In the second stage they are copied,
  so the backward differences of tau W telescope as they do when
  V* = V_{n+1}.
* **mu_B / T cache.** With diffusion on, mu_B / T is only filled in the
  cells that are read: those that are not frozen, and their neighbours.
* **Sources.** Skipping stays on while sources deposit, but the cells
  they heat are not vacuum (see "With sources" below).
* **Log.** At the end the run logs the fraction of cell updates and faces
  that were skipped.

On the CPU the logic is in `Advance::AdvanceIt`. The device backends do
the same (see "Device backends" below).

## Test problem

The setup is the full event of `doc/openmp_cpu_scaling.md`: TATB,
96 x 96 x 64 cells, EoS 14 (NEOS-BQS), shear + bulk + second order + baryon
diffusion, 295 steps, tau = 1.35 to 6 fm, and freeze-out at
e_FO = 0.3 GeV/fm^3. The conservation numbers are the maximum over the run
of `global_conservation_laws.dat`, boundary outflow included.

## Viscous hydro: speed against conservation

| e_vac (GeV/fm^3) | wall (s) | cell updates skipped | faces skipped | max dE/E0 | max dPz/E0 | max dN_B/N_B0 |
|---|---|---|---|---|---|---|
| 0 (off) | 119.1 | | | 1.6e-7 | 2.2e-7 | 2.2e-7 |
| 1e-10 | 105.8 | 21 % | 27 % | 1.7e-7 | 2.6e-7 | 2.2e-7 |
| 1e-8 | 100.0 | 31 % | 39 % | 1.2e-7 | 2.0e-7 | 2.3e-7 |
| **1e-6** | **84.1** | **46 %** | **54 %** | **2.2e-6** | **7.4e-7** | **2.8e-7** |
| 1e-4 | 64.9 | 64 % | 71 % | 8.5e-5 | 6.0e-5 | 7.8e-5 |
| 1e-3 | 55.7 | 72 % | 79 % | 1.6e-3 | 1.1e-2 | 3.9e-4 |

* **Up to 1e-8, conservation is unchanged.** The deviations stay at the
  1e-7 level of the run without skipping.
* **At 1e-6, the run is 29 % faster** and dE/E0 is 2e-6. That is ten times
  the baseline, but still at the O(dtau^2) level that the integration
  tests allow (2e-5).
* **At 1e-4 and above, the deviations grow** smoothly at late times
  (tau > 4.5 fm), when the fireball has spread into the dilute region.

### Where the leak comes from

* **Ideal hydro conserves exactly at any threshold.** The ideal run of the
  same event gives dE/E0 = 5.9e-9 at both 1e-4 and 1e-3 (1.6e-8 without
  skipping), and it is 2 - 2.3x faster (36.8 s against 19.0 s and 16.2 s).
  The zero-flux faces and the frozen cells conserve exactly by
  construction.
* **The viscous leak is not in the frozen W update.** Skipping only the
  faces, while every cell keeps its full viscous update, gives the same
  dE/E0 at 1e-4 (1.2e-4). Copying W in the second stage instead of
  diluting it changes nothing either.
* **It comes from zero ideal fluxes in contact with viscous dynamics.** The
  cells just outside the frozen region see a stalled neighbour. The
  velocity gradient across that boundary is large, so W there is large
  compared with e. `QuestRevert` then limits W, and it does not conserve
  T^{mu nu}. The leak therefore scales with how much energy sits in the
  cells near the threshold, and that grows quickly above 1e-6.
* **Two errors fixed along the way:**
  * Copying the frozen cell state, instead of the zero-flux update, made
    tau T^{tau tau} grow by dtau / tau every step (2.4e-3 at 1e-4).
  * Zeroing the boundary faces broke the outflow term of the
    conservation check.

## Observables

These come from the freeze-out surface at e_FO = 0.3 GeV/fm^3, all
per-thread files summed. The table shows the relative change against the
run without skipping. The weights are u^mu dsigma_mu; "v2 proxy" is
<u_x^2 - u_y^2> / <u_x^2 + u_y^2>.

| e_vac | surface cells | V = sum u.dsigma | <u_perp> | v2 proxy | N_B on the surface | <T> |
|---|---|---|---|---|---|---|
| (off: values) | 165895 | 114.49 | 0.2023 | 0.0813 | 7.262 | 0.747 |
| 1e-10 | same | -2e-12 | 8e-11 | 2e-9 | 2e-12 | -7e-13 |
| 1e-8 | same | 1e-11 | 7e-12 | -8e-11 | 1e-11 | 3e-13 |
| 1e-6 | same | 1e-12 | 3e-12 | 8e-11 | -8e-13 | 2e-13 |
| 1e-4 | same | 6e-12 | 1e-10 | -2e-9 | -1e-12 | 2e-12 |
| 1e-3 | same | 6e-9 | 5e-8 | 1e-6 | 6e-8 | -5e-9 |

The surface is the same to the float precision of the output (about 1e-7)
for every threshold up to 1e-4. Even at 1e-3 the changes are around 1e-6.
The conservation leak is in the dilute matter outside the surface, which
is energy that particlization never sees. A threshold more than 2 - 3
orders of magnitude below e_FO cannot reach the surface within the event.

## Also fixed: the ideal-gas mu_B / T solve in dilute cells

The integration test (EoS 0, the ideal gas with baryons) was first
**slower** with skipping: 14.2 s against 5.1 s. That had two causes.

1. **`EOS_idealgas::get_muB_over_T` underflow.** When n_B^4 / e^3
   underflows (n_B ~ 1e-80 in the dilute tails), the table guess is 0. The
   code then called the robust solver `x_of_rho`, which started from
   x = 1e-300 and took its full 60 log/exp Newton steps without reaching
   the answer. With skipping, the frozen cells keep their tiny n_B instead
   of being diluted further, so 15 M calls took this path, against 0.27 M
   without skipping. In this limit x << 1 and the exact small-x solution
   is x = 27 n (3a / e)^{3/4} / b. The code now returns it directly. The
   outputs are bit-identical (conservation file, with and without
   skipping), because the old path converged to the same value or to
   something negligible in these cells. This also speeds up every EoS 0
   run with baryons in its tails.
2. **Load imbalance in `computeMuBOverT` and `computePressure`.** Both
   used a static schedule, while the EoS cost differs from cell to cell.
   They now use `schedule(dynamic, 512)`, which is bit-identical.

With both fixes the integration test takes **3.8 s with skipping against
5.1 s without**. The EoS 14 full event is not affected: its table lookups
cost the same in every cell.

## Defaults

When `vacuumSkipEnergy` is not in the input file (or is negative), the
default follows `beastMode` (`read_in_parameters.cpp`). The beastMode
used is the final one, after the reset of beastMode >= 2 to 1 for the
evolution output:

| beastMode | default vacuumSkipEnergy (GeV/fm^3) |
|---|---|
| 0 | 1e-8 |
| 1 | 1e-6 |
| >= 2 | 1e-4 |

A value given in the input overrides the default, and 0 turns skipping
off. With these defaults all 25 integration tests and the unit tests
pass. The `vacuum_skip` test sets 0 explicitly for its reference run.

## Recommendation

* **Viscous hydro: `vacuumSkipEnergy 1e-6`.** That is 5.5 orders of
  magnitude below e_FO = 0.3. It gives 1.4x on this grid, with energy
  conserved to 2e-6 and an unchanged freeze-out surface. 1e-8 is the safe
  choice when conservation at the 1e-7 level matters, and gives 1.2x.
* **Ideal hydro: any threshold well below e_FO**, for example 1e-4, which
  gives 2x. Conservation is exact.
* The gain grows with the vacuum fraction of the grid, so it is larger on
  generous grids and late in the event. It is 0 % at the start of a tight
  grid.
* Set it to 0 for the tests that check the dilute dynamics themselves.

## Tests

* **New integration test `vacuum_skip`.** It runs a 41 x 41 x 49 TATB event
  with shear, bulk and diffusion (EoS 0), with and without
  `vacuumSkipEnergy 1e-6`. It checks that:
  * more than 10 % of the cell updates are skipped (35 % here);
  * E and N_B are conserved to 2e-5 (measured: 1e-7);
  * E differs from the run without skipping by less than 1e-6 (measured:
    9e-10);
  * the freeze-out surface has the same number of cells and the same
    volume (measured: identical).
* **All 25 integration tests pass** (memory_history is skipped without its
  driver), and the unit tests pass. With `vacuumSkipEnergy 0` the outputs
  are bit-identical to before.

## Device backends

The device kernels skip the same cells. The host passes the threshold in
`StepParams::eVacuum` (1/fm^4, 0 when skipping is off). There is no mask
buffer: each thread tests the e of its stencil itself, and with sources in
the stage also the source (`isVacuumIdx`, `isVacuum`, `isFrozenCell` in
`device_ideal.h`).

* **`kernelFace`.** An interior face whose four stencil cells are vacuum
  writes zero flux and returns before the reconstruction.
* **`kernelCellVisc`.** A frozen cell still does `FirstRKStepT`, with the
  W divergences. It then replaces the gradients, `firstRKStepW` and the
  second-stage fixed point with the dilution of W, Pi and Pi_chem, as
  `Advance::frozenViscousUpdate` does.
* **Unchanged.** The ideal `kernelCell` (the CPU code does the full ideal
  update of a frozen cell too), `kernelWFace` and the mu_B / T fill. The
  mu_B / T fill costs under 1 % on Metal, and only computing the cells
  that are read would need the frozen mask of the neighbours.
* **Logging.** The device run logs "skipping the cells below ...
  GeV/fm^3" in its "hydro on the device" line. It does not count the
  skipped cells.

The frozen regions are large and contiguous, so the threads of a GPU
SIMD group mostly take the same branch.

### Results (96 x 96 x 64 full event, as above)

| backend | e_vac | wall (s) | max dE/E0 | max dPz/E0 | max dN_B/N_B0 |
|---|---|---|---|---|---|
| 1 (device kernels on the CPU, double) | 0 | 133.4 | 1.6e-7 | 2.2e-7 | 2.2e-7 |
| 1 | 1e-6 | 100.5 | 2.2e-6 | 7.4e-7 | 2.8e-7 |
| 1 | 1e-4 | 83.3 | 8.5e-5 | 6.0e-5 | 7.8e-5 |
| 2 (Metal, float) | 0 | 14.0 | 4.3e-6 | 1.1e-6 | 1.7e-7 |
| 2 | 1e-6 | 11.6 | 4.6e-6 | 1.1e-6 | 1.6e-7 |
| 2 | 1e-4 | 10.8 | 4.2e-5 | 4.3e-5 | 1.7e-5 |

* **gpu_backend 1 is bit-identical to the CPU code with skipping** at
  both thresholds (`global_conservation_laws.dat`).
* **Metal is 1.2x faster at 1e-6 and 1.3x at 1e-4.** At 1e-6 the
  conservation stays at the float level of the run without skipping.
* **The freeze-out surface is unchanged.** It has the same number of
  cells, and the volume changes by under 1e-9 on Metal.
* **Ideal hydro on Metal** (EoS 9, 215 steps): 3.5 s, against 2.2 s at
  1e-4. dE/E0 goes from 3.5e-6 to 4.2e-6 and dPz/E0 from 1.4e-6 to
  6.7e-6, both at the float level.
* **The gain is smaller than on the CPU** because the Metal step is
  dominated by kernels that still run on every cell: the face kernel of
  non-vacuum faces, `kernelWFace` and the ideal part of the cell kernel.
* **CUDA** uses the same kernel functions and needs no other change, but
  it is **not compiled or tested** (no nvcc on this machine).
* **Test.** The `vacuum_skip` integration test also runs backends 1 - 4
  with 1e-6 and checks them against the CPU run with skipping, to the
  `gpu_backend` limits:
  * backend 1: 0 difference;
  * Metal: 1.6e-6 in the conservation and 1.2e-7 in the surface volume.

## With sources

Skipping used to be off in every stage until the sources were over
(`get_source_tau_end() + dtau`), because a source can heat any cell. Now
the cells away from the sources are skipped while the sources deposit.

* **CPU code.** In a stage with sources (`flag_add_hydro_source`, before
  `get_source_tau_end() + dtau`, and `is_active_at`), `AdvanceIt` first
  evaluates the source of every cell in one parallel pass
  (`computeSourceCache`, the same `cellSource` that `FirstRKStepT` used to
  call per cell). `FirstRKStepT` reads the result, so the arithmetic is
  unchanged and runs without skipping are bit-identical. `markVacuum`
  counts a cell as vacuum only if it stays below the threshold with the
  energy of this step: e + |J^tau| dtau < e_vac (qi_source = tau J). The
  faces and frozen cells then follow from this mask as before.
* **Device.** `isVacuumIdx` applies the same test with the stage's source
  buffer (the strings kernel or the host's source terms, filled before
  the faces). The face kernel now also reads the source buffer.
* **A frozen cell still gets its source** in `FirstRKStepT`, so the
  energy it receives below the threshold is conserved. It is only not
  transported until the cell rises above the threshold, and the mask is
  rebuilt every stage.
* **The ghost cells keep the test on e.** A face or stencil that reads a
  ghost cell also reads the edge cell it copies, with the source test, so
  the CPU and device masks agree.
* **Memory.** The cache is 7 doubles per cell of the grid, only during
  the source stages, and freed when they end. That is 84 MiB for
  128 x 128 x 96.

### Results: strings event

The test case is the strings event of the integration tests (40 strings,
deposits from tau = 0.5 to 2.0 fm). It runs on 61 x 61 x 49 cells, with
shear and diffusion, EoS 0, to tau = 2.48 fm, on 10 threads. The
deviations are the largest difference from the run without skipping,
relative to the largest E of the run.

| backend | e_vac | wall (s) | cells skipped | max dE / E_max | max dN_B / N_B |
|---|---|---|---|---|---|
| CPU | 0 | 85.4 | | | |
| CPU, before (skip only after the sources) | 1e-6 | 83.3 | 37 % | 7.8e-8 | 3.6e-8 |
| CPU | 1e-8 | 71.0 | 36 % | 1.3e-6 | 3.3e-7 |
| CPU | 1e-6 | 66.2 | 45 % | 2.2e-5 | 4.1e-7 |
| CPU | 1e-4 | 59.5 | 55 % | 9.2e-5 | 2.3e-6 |
| Metal | 0 | 30.2 | | | |
| Metal, before | 1e-6 | 29.7 | | 1.7e-6 | 1.2e-7 |
| Metal | 1e-8 | 28.3 | | 5.1e-5 (see below) | 4.2e-7 |
| Metal | 1e-6 | 27.2 | | 1.6e-5 | 1.4e-7 |
| Metal | 1e-4 | 25.9 | | 9.1e-5 | 1.5e-7 |

* **The CPU code is 1.2 - 1.4x faster** than before. Before, almost
  nothing was saved, because the cells were only skipped after the
  deposits (tau > 2.0 fm, about a quarter of the 338 steps). The "cells
  skipped" column counts only the stages in which skipping is on.
* **gpu_backend 1 is bit-identical to the CPU code** at 1e-6.
* **With skipping off the results are bit-identical** to before, since
  the source cache is not used then.
* **The deviation grows smoothly with tau**, from the first deposits on,
  as in the runs without sources. It is larger than there (2.2e-5 against
  2.2e-6 at 1e-6), because the thin string deposits leave much of their
  energy in dilute cells near the threshold.
* **The Metal value at 1e-8 is a single event in float.** A jump of
  +0.16 GeV in two steps at tau = 1.36 fm causes it. It is deterministic
  and the same with one command buffer per kernel, so it is not a race.
  The thresholds 3e-9, 1.01e-8, 3e-8, 1e-7 and 3e-7 give 2.4e-6, 2.2e-6,
  2.8e-6, 3.9e-6 and 8.4e-6. So it is a float-sensitive event set off by
  a small perturbation, not a systematic leak. For comparison, Metal
  differs from the CPU code by 4.5e-6 without skipping.
* **Metal gains less** (1.07 - 1.17x), because the strings source kernel
  and the kernels that run on every cell dominate its step.
* **Test.** The `vacuum_skip` integration test now also runs a strings
  event that ends while the strings deposit (tau = 1.6 fm). With 1e-6 it
  skips 43 % of the cell updates, with E within 9.7e-6 and N_B within
  9e-7 of the run without skipping. Backend 1 differs by 0 and Metal by
  4.7e-6.

## Not done

* **Device compaction.** The skipped threads still launch. A list of
  active tiles (for example 8 x 8 x 1) would also remove the launch and
  test cost of the frozen cells, and allow skipping `kernelWFace` there.
* **Tiles on the CPU.** The mask is built per cell in every stage (two
  passes over the grid, about 1 % of the time). An active-tile list would
  also skip the loop overhead of the frozen cells. That is a further gain,
  but a small one.
