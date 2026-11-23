# GPU port, phase 5: the strings source on the device, and fewer host copies

Branch `gpu`, from 1083265 (phase 4, [gpu_porting_phase4_diffusion_chem_causality.md](gpu_porting_phase4_diffusion_chem_causality.md)). 2026-11-22. Machine: Apple M4, as before.

## Summary

- **Item 5: the strings source (`HydroSourceStrings`, Initial_profile 13/131) is now evaluated on the device.**
  - Once per step, the host packs the strings of the τ frame into a buffer.
  - Each RK stage then runs `kernelStringSource` before the update.
  - On the strings event, the host source cost goes from 2.24 s to 0.6 ms, and the Metal run from 3.7 s to **1.8 s** (CPU code: 5.8 s).
- **Item 6: the host Fields are no longer filled after every step.**
  - Evolve copies from the device only what it reads:
    - all of both slices on output and freeze-out steps;
    - the η window of the per-step diagnostics;
    - the outermost cells for the boundary outflow and the beastMode 3 edge check.
  - The conservation sums and the grid maxima are computed on the device.
  - The host pressure (`preparePressure`) is skipped.
  - On the 1.27M-cell viscous benchmark, Metal wall time drops by **15% (EoS 9) and 10% (EoS 20)**.
- **Validation:**
  - **Backend 1 (double) reproduces the CPU code bit for bit** with clang:
    - the strings source (ideal and viscous, energy and n_B);
    - the conservation sums and maxima (ideal/viscous × EoS 0/9).
  - **New integration test `gpu_backend_freezeout`**, run on backends 1 and 2:
    - freeze-out and evolution output of a 3D viscous fireball;
    - beastMode 3 with grid enlargement and coarsening.

    The surfaces are identical on backend 1, and agree to 4e-6 on Metal with the same cell count.
- **Results:**
  - With `gpu_backend = 0` (the default) nothing changes: the CPU path through Evolve and `check_conservation_law` is the same code.
  - With a device backend, `global_conservation_laws.dat` changes at float rounding on Metal (1.2e-7 of E), because the sums are now taken in float on the device. Every other output file is byte-identical to phase 4.
- **Found on the way:** phase 4 made the Metal viscous kernels **25% slower** on this benchmark (7.36 s → 9.22 s for EoS 9), even with diffusion, Π_chem and the limiter off. This phase did not cause it (HEAD of phase 4 measures the same 9.22 s). It is a target for item 4 (performance).
- **Tests:** all unit tests and all integration tests pass (section 5).

## 1. The strings source on the device

| file | what |
|---|---|
| `src/device/device_strings.h` (new) | `kernelStringSource`: J^μ (energy and remnants) and the n_B and n_Q sources of one cell, a statement-by-statement copy of `get_hydro_energy_source`, `get_hydro_rhob_source` and `get_hydro_rhoQ_source`. The u comes from the stage's input slice (the predictor in the second stage). A NaN sets `kFlagSourceNan` (fatal, as `qi_source is nan` on the host). |
| `src/device/device_hydro.cpp` | `packStrings`: the union of the four τ-frame lists (energy, remnant, baryon, electric) in the order of the full list. Each record has membership flags and the string-only prefactors, computed on the host in double exactly as in the host loops. It is sent once per step (`setStrings`). |
| `src/hydro_source_strings.h` | read-only accessors for the lists, `preEqFlowFactor_` and `parton_quench_factor` |
| `src/device/device_defs.h` | `derf`: `std::erf` on the host, `erf` in CUDA. MSL has no erf, so Metal uses the Numerical Recipes `erfcc` approximation (relative error < 1.2e-7, below float precision). |
| `metal_backend.m`, `cuda_backend.cu` | the kernel, the strings buffer, `setStrings`. `StepParams::stringSource` sits at a fixed byte offset (`MUSIC_DEVICE_STRINGS_OFFSET`, checked by `static_assert`) for the Objective-C side. |

Other sources (TATB, AMPT, Tmunu) still use the phase 3 path: the host evaluates them per stage from the device's u, and the device adds them.

**Strings event** (`gpu_backend_sources`: 40 strings, 41×41×49, 101 steps, viscous with n_B):

| | phase 3 | phase 5 |
|---|---|---|
| CPU code | 5.8 s | 5.8 s |
| Metal wall | 3.7 s | **1.8 s** |
| Metal kernels | 0.75 s | 1.02 s (includes the strings kernel) |
| host strings source / packing | 2.24 s | 0.0006 s |
| Metal vs CPU: conservation / E, N_B | 5.8e-6, 7.6e-6 | 5.7e-6, 7.7e-6 |
| backend 1 (g++) vs CPU | 8.4e-13 | 8.4e-12, N_B 0 |

## 2. The host copies only what it reads

`DeviceHydro::advanceStep` no longer downloads. Evolve calls:

| call | when | what |
|---|---|---|
| `syncHost(curr, prev)` | any output step (`output_evolution_data` 1–4 and 7, memory output, Knudsen/Reynolds); the freeze-start step and every freeze-out step; `output_momentum_anisotropy_vs_etas` steps; Gubser/1+1D check times; vorticity, L^{μν} and monitor output; before beastMode changes Δτ or the grid; after a frozen step; before every return | both slices, all cells, ghosts updated |
| `syncEtaWindow(curr, η_min, η_max)` | every step in beastMode 0: `output_momentum_anisotropy_vs_tau` (−0.5, 0.5); the phase-diagram trajectories of profiles 13/131 (−0.5, 5) | the η slices of the window, with the same η formula as `Cell_info` |
| `syncShell(curr, prev, 1, 1, 1)` | every step: the boundary outflow in `check_conservation_law` | the outermost cells in x, y, η |
| `syncShell(curr, -, bandX, bandY, 0)` | beastMode 3: `fluidNearGridEdge` | its transverse band |

- **The decision is made once per step**, before any reader, so a full sync is never followed by stale partial reads.
- **A guard checks the uploads:** `advanceStep` refuses to upload host Fields that are older than the device state. Every path that changes the Fields (`regridTransverse`, `extrapolatePrevForNewDtau`) syncs first.

**Conservation diagnostic on the device** (`src/device/device_diag.h`):
- `kernelConservationRow` is the body of the row loop of `check_conservation_law`, run once per (η, x) row. It returns the 7 sums and the 5 maxima (e, n_B, n_Q, n_S, T).
- The host adds the rows in the host's order and continues with the units, the outflow and the output as before. `check_conservation_law` takes the rows as an optional argument.
- The maxima fill `GridMaxima` (also in boost-invariant runs), so `get_maximum_energy_density` does not sweep the host Fields.
- For T of ideal runs, the thermal tables are now uploaded for ideal hydro as well.

**Viscous benchmark** (the phase 2/3 event: lumpy 3D, 161×161×49 = 1.27M cells, shear + bulk + second order, n_B, 101 steps; repeated runs agree to 1%):

| Metal | EoS 9, phase 4 | EoS 9, phase 5 | EoS 20, phase 4 | EoS 20, phase 5 |
|---|---|---|---|---|
| wall | 13.30 s | **11.37 s** | 18.49 s | **16.73 s** |
| kernels | 9.22 s | 9.23 s | 12.92 s | 12.95 s |
| copies to the host | 1.06 s (101 full) | 0.52 s (4 full, 194 partial) | 1.04 s | 0.55 s |
| conservation sums on the device | — | 0.24 s | — | 0.25 s |

- **Where the saved time comes from:** half of the copy time, plus the host `preparePressure` and the host conservation sweep, about 1.3 s.
- **What is left:** the copies are now mostly the η window of `output_momentum_anisotropy_vs_tau` (every step). The conservation kernel runs only 7889 threads, one per row, which leaves the GPU underused; a two-level reduction would make it cheaper.

**Phase 4 kernel time:** phase 3 measured 7.36 s (EoS 9) and 11.0 s (EoS 20) of kernels on this event, and phase 4's HEAD measures 9.22 s and 12.92 s. The viscous kernels now hold the diffusion, Π_chem and causality code behind runtime flags. The extra register pressure or code size probably costs even when the flags are off. Metal function constants (specialisation at pipeline creation) are the natural fix.

## 3. CUDA

- **Changes:** the strings and conservation kernels, `setStrings` and `conservation`.
- **`fetch` copies whole components,** so a partial sync still moves the full slice over PCIe. The η windows are contiguous ranges, so sub-range copies are a later improvement.
- **Status:** it type-checks with clang's CUDA front end (host and device passes); it has not been run, since there is no NVIDIA hardware here.

## 4. Tests

- **Unit tests** (`src/bugfix_unittest.cpp`, backend 1, clang):

  | test | result |
  |---|---|
  | **strings source** (new), ideal and viscous: a strings event made from `tests/test_source_terms/strings_event_0.dat` with τ_0 = 0.2 fm/c, depositing E and n_B on the grid in 4 steps | 0 differing values; float L1 7e-7 (E) and 5e-7 (n_B) |
  | **conservation sums** (new): `check_conservation_law` from the device rows, with only the outermost cells on the host, against the full host sweep; ideal/viscous × EoS 0/9 × double/float | double: identical lines and maxima; float: 3e-8 (sums) and 5e-8 (maxima) |

  The existing device tests call `Evolve::syncHost` after their steps.
- **Integration tests:**
  - `gpu_backend_sources` accepts the new log line.
  - New `gpu_backend_freezeout`, backends 0/1/2:

    | case | backend 1 | Metal |
    |---|---|---|
    | 3D viscous fireball, freeze-out every 5 steps, evolution output every 10 | surface, cell count and conservation identical; evolution file 2e-11 | surface fluxes ≤ 1.4e-7, same 98,990 cells, conservation 1.5e-7, evolution file 4e-7 |
    | beastMode 3 (grid enlarged and coarsened) | surface identical | 4e-6, same 6,592 cells |

## 5. Regression

- **Unit tests (bgunit, clang):** all pass:
  - unittest_bugfix: 42 cases, 1,519,832 assertions;
  - cell, grid, minmod and reconst: all pass.
- **CUDA:** syntax check passes.
- **Integration (bmetal, with the bgunit memory driver):** 25/25 pass. The full run had 24 tests and passed; the new `gpu_backend_freezeout` was added after it and passed on its own.

Next: [phase 6](gpu_porting_phase6_performance.md), performance.
