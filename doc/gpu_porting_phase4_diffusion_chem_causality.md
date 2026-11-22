# GPU port, phase 4: the rest of the viscous update

Branch `gpu`, from 6d85b61 (phase 3, [gpu_porting_phase3_sources_eos4d.md](gpu_porting_phase3_sources_eos4d.md)). 2026-11-22. Machine: Apple M4, as before.

## Summary

- **Phase 4 finishes the leftovers of phase 2.** `gpu_backend` now also runs:
  - viscous hydro with **NEOS** (EoS 10–14);
  - the **causality limiter** (`FlagResetCausality = 1`);
  - **baryon diffusion** (`turn_on_baryon_diffusion = 1`);
  - the **chemical bulk pressure** (`turn_on_bulk_chem = 1`).

  Only two settings still fall back to the CPU code: `viscosity_flag` other than 0 or 1, and transport coefficients from files (type 23).
- **Validation:**
  - **Backend 1 (device kernels on the CPU, double) reproduces `Advance` bit for bit** with clang (unit tests). Every new path has 0 differing values:
    - viscous NEOS: 4 configurations;
    - diffusion, Π_chem, and both together: EoS 0, 14 and 20, 2 configurations each, 18 cases.
  - **The causality limiter agrees to rounding, not bit for bit.** The host takes the eigenvalues of W^μ_ν from Eigen, the device from a closed form. They differ by at most 8.5e-12 of the largest value after 4 strongly limited steps.
  - **Metal (float)** agrees with the CPU code on the lumpy TATB event at τ = 2 fm/c, as closely as the cases without the new features (table below).
  - **Viscous Gubser flow** (tests/Gubser_flow, η/s = 0.2): on Metal, the L1 error against the semi-analytic solution equals the CPU code's to 3e-4 (relative) for T, u and every π^{μν} component, at τ = 1.2, 1.5 and 2 fm/c.
- **Results:** with `gpu_backend = 0` (the default) nothing changes, because no host code changed.
- **Tests:** all unit tests and all integration tests pass (the full run is in section 4).

## 1. Code

| file | what |
|---|---|
| `src/device/device_visc.h` | **Diffusion:**<br>• ∂^μ(μ_B/T): minmod in space, backward difference in τ (`MakeDSpatial`, `MakeDTau`);<br>• `uqSource` (= `Diss::Make_uqSource`) with `dMuBOverTIdeal` (= `D_muB_over_T_ideal`);<br>• the q RK average and q^τ from transversality;<br>• `questRevertQmu` (= `QuestRevert_qmu`).<br>**Chemical bulk:** `uPiChemSource` (= `Make_uPiChemSource`), its W-rhs terms and its RK average.<br>**Causality:** `wEigenvalues` (closed form, see below) and `causalityConstraints` (= `nCausalityConstraints`).<br>**Shared:** the QuestRevert factor moved into `questFactor` (the same expression as before). Π = Π_bulk + Π_chem enters T^{μν} whenever either one is evolved, as on the host. |
| `src/device/device_eos.h` | `eosTMuB` (T and μ_B for the μ_B/T gradients, all four EoS types); `eosCs2` (= `calculate_velocity_of_sound_sq`, for the limiter); the NEOS thermal variables (T, μ_B, μ_S, μ_Q from the 2D tables, n_B capped at the table edge in the entropy) |
| `src/device/device_ideal.h`, `device_hydro.cpp` | new `StepParams`: `causality`, `bulkChem`, `diff`, `profile1`, `kappaCoeff`, `chemRateC`; `nWComp` follows `Diss::nWComponents` (10 with diffusion). The NEOS thermal tables are uploaded for viscous runs. The fallbacks for these features are removed. |
| `src/device/device_defs.h`, `cuda_backend.cu` | `acos`, `cos`, `tanh` for MSL and CUDA |
| `src/bugfix_unittest.cpp` | bit-for-bit tests:<br>• viscous NEOS (EoS 14 added to the viscous test);<br>• the causality limiter;<br>• diffusion and Π_chem. The lumpy test state gets a transverse q^μ and Π_chem. |
| `tests/integration/music_integration_tests.py` | `gpu_backend_viscous`: EoS 14 and the cases diffusion and causality; new `gpu_backend_pi_chem` (IP-Glasma start, hotQCD) and `gpu_backend_gubser_viscous` |

### The eigenvalues of W^μ_ν

`Advance::nCausalityConstraints` takes the eigenvalues of W^μ_ν from `Eigen::EigenSolver`, which the device cannot call. W^μ_ν is transverse to u, so one eigenvalue is 0. The other three are the roots of a cubic whose coefficients come from tr A, tr A² and tr A³. The device solves the cubic by the trigonometric method; the roots are real because W is symmetric in the local rest frame.

### D(μ_B/T) in float

The host takes D(μ_B/T) = −θ[(e+P) ∂α/∂e + n ∂α/∂n], with α = μ_B/T, from two central differences with step 1e-4. For the conformal ideal gas, α depends on n/s alone, and the ideal flow conserves n/s, so the two terms cancel exactly. In float, each difference keeps only about 3 of its 7 digits, so only noise is left.

For the ideal gas in float, the device instead takes one difference along (e+P, n) with step 1e-3. This is the same derivative, without the cancellation, and it cuts the float error of q by up to 20× in the unit test. Double keeps the host's expression bit for bit. The tabulated EoSs keep it in float too: their tables are piecewise linear, and a larger step crosses the kinks (tried: it made q worse).

## 2. Unit tests (backend 1, clang)

| test | cases | differing values | float (L1 of q, Π_chem, e) |
|---|---|---|---|
| viscous update, EoS 0, 9, 14 | 12 | 0 | W, Π ≤ 4.4e-4 (NEOS, QuestRevert) |
| causality limiter, EoS 0, 9 | 4 | ≤ 8.5e-12 (rounding) | — |
| diffusion / Π_chem / both, EoS 0, 14, 20 | 18 | 0 | q 8e-5 to 2.5e-3; Π_chem ≤ 1.8e-5; e ≤ 4.2e-5 |

**Float error of q in the unit test.** In this rough test state, q relaxes in τ_ρ = 3Δτ (clamped). Its explicit update, with the q·ω and q·σ terms, amplifies rounding:
- after two steps, a few cells differ by O(1) of their q, and this happens in the host's own evolution as well;
- turning off the vorticity term brings EoS 0 down to 5e-5.

The float bounds there are therefore loose (q < 3e-2), and float is judged on the realistic event below.

## 3. Metal (float) against the CPU code

Lumpy TATB event (41×41×17, τ = 1 → 2 fm/c), shear and bulk with second-order terms. The W, Π, Π_chem and μ_B columns are L1 differences over the cells with T > 0.1 GeV.

| case | E (max) | N_B (max) | W/(e+P) | Π | Π_chem | μ_B |
|---|---|---|---|---|---|---|
| EoS 0 | 8.3e-7 | 6.2e-8 | 2.0e-6 | 3.1e-6 | — | 4.2e-7 |
| EoS 9 | 6.6e-7 | 5.7e-8 | 5.3e-6 | 2.6e-6 | — | — |
| EoS 14 (NEOS, new) | 1.0e-6 | 6.6e-8 | 3.1e-6 | 2.4e-6 | — | 7.4e-7 |
| EoS 20 | 6.4e-7 | 6.3e-8 | 5.1e-6 | 3.1e-6 | — | 1.9e-6 |
| diffusion, EoS 0 | 8.3e-7 | 1.1e-7 | 1.9e-6 | 3.2e-6 | — | 1.4e-6 |
| causality, EoS 0 | 8.3e-7 | 6.1e-8 | 2.0e-6 | 3.1e-6 | — | 4.2e-7 |

Backend 1 with g++ (the Metal build) gives 0 to 1e-16 for every case except EoS 14 (8e-9, the known g++ contraction noise of NEOS) and EoS 20 (2e-11).

**Chemical bulk pressure (new test `gpu_backend_pi_chem`).** Π_chem = 0 is a fixed point of its equation: √Y_q = 1 gives ζ_chem = 0. On the TATB event, which starts at Π_chem = 0, the case is therefore identical to the run without it. With P = e/3, Π_chem stays 0 in any case.

The test therefore starts from IP-Glasma (profile 9), which sets Π_chem = e/3 − P at τ_0, and uses hotQCD. It runs 41×41×9 cells from τ = 0.6 to 1.6 fm/c, and Π_chem reaches 8× Π_bulk. Results:

| backend | E, P, N_B (max) | Π_bulk (L1) | Π_chem (L1) |
|---|---|---|---|
| 1 | 2.9e-16 | 0 | 0 |
| 2 (Metal) | 1.1e-7 | 1.7e-6 | 4.1e-6 |

**Viscous Gubser flow, L1 error against the semi-analytic solution (y = x):**

| τ (fm/c) | T: CPU / Metal | u^x: CPU / Metal | π^xx: CPU / Metal | π^ηη: CPU / Metal |
|---|---|---|---|---|
| 1.2 | 2.6185e-4 / 2.6184e-4 | 5.6888e-4 / 5.6899e-4 | 3.3671e-3 / 3.3672e-3 | 1.7215e-3 / 1.7219e-3 |
| 1.5 | 7.6377e-4 / 7.6379e-4 | 1.8756e-3 / 1.8757e-3 | 9.9532e-3 / 9.9528e-3 | 6.2150e-3 / 6.2149e-3 |
| 2.0 | 8.3707e-4 / 8.3716e-4 | 4.9904e-3 / 4.9906e-3 | 3.1215e-2 / 3.1213e-2 | 9.8376e-3 / 9.8382e-3 |

Metal and the CPU code have the same error to 3e-4 (relative); the float rounding is far below the discretization error. Run time: CPU 5.4 s, Metal 2.0 s (201×201, 200 steps).

The host's "pixy_diff" compares π^xy with the analytic π^xx (`Cell_info::Gubser_flow_check_file`), so its 0.5–0.7 is not an error of the evolution; the test only compares Metal with the CPU code there.

## 4. Regression

- **Unit tests (bgunit, clang):** all pass:
  - unittest_bugfix: 40 cases, 1,491,722 assertions;
  - cell, grid, minmod and reconst: all pass.
- **CUDA:** syntax check (host and device passes) passes.
- **Integration (bmetal, with the bgunit memory driver):** 24/24 pass, including the new `gpu_backend_pi_chem` and `gpu_backend_gubser_viscous`. (The full run had 23 tests; `gpu_backend_pi_chem` was added after it and rerun on its own together with `gpu_backend_viscous`.)
