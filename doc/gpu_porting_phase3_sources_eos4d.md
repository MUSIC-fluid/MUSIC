# GPU port, phase 3: source terms and the 4D EoS on the device

Branch `gpu`, from cebd857 (phase 2, [gpu_porting_phase2_viscous.md](gpu_porting_phase2_viscous.md)). 2026-11-22. Machine: Apple M4 (10 CPU cores, 10-core GPU, 16 GB unified memory); no NVIDIA hardware.

## Summary

- **`gpu_backend` now runs with source terms**, ideal and viscous hydro:
  - tested with TATB (112) and strings (13);
  - AMPT (30), T^{μν} (17) and the profiles 113 and 131 take the same path (any `HydroSourceBase`) but are not tested here.

  For each RK stage the host evaluates the sources from the u of the device slice (the predictor in the second stage), exactly as `Advance` does. The device adds them to q at the same point as `FirstRKStepT`.
- **`gpu_backend` now runs with the 4D EoS (EoS 20)**, ideal and viscous:
  - the device reads the same float tables as `EOS_4D`;
  - the device repeats its float/double arithmetic operation by operation.

  n_Q and n_S now reach every pressure call on the device. Before, only the EoS gradients got them, because no device EoS used them.
- **Validation:**
  - **Backend 1 (device kernels on the CPU, double) reproduces `Advance` bit for bit** with clang (unit tests). It finds 0 differing values with sources in 4 configurations (ideal and viscous, 3D and boost-invariant). The source depends on u, as the strings and AMPT sources do. With EoS 20 and nonzero n_Q, n_S it also finds 0 differences, in 4 configurations including dilute cells.
  - **Metal (float) agrees with the CPU code:**

    | event | E, P (of the largest E) | N_B |
    |---|---|---|
    | TATB source, ideal | 7.5e-7 | 6e-8 |
    | TATB source, viscous | 7.1e-7 | 7e-8 |
    | strings, viscous | 5.8e-6 | 7.6e-6 |

    With EoS 20 (lumpy event, without sources), E agrees to 6–9e-7, and W/(e+P) and Π to 3–5e-6 (L1).
- **Speed:**
  - viscous EoS 20 benchmark: **Metal is 6.9× faster** than the CPU code (16.5 s against 114.8 s);
  - strings event: the host evaluation of the string source takes 3× the kernel time, so it now limits the speed-up (see [section 2](#2-source-terms)).
- **Results:** with `gpu_backend = 0` (the default) nothing changes. The only change to host code is read-only accessors in `EOS_4D` (and `diluteScale` reading its constant through one of them, the same float). All integration tests pass.

## 1. Code

| file | what |
|---|---|
| `src/device/device_ideal.h` | `cellAddSource`: q += qi_source Δτ before the RK average, as in `FirstRKStepT`. Both cell kernels take a source buffer. `Solver` carries J⁰_Q and J⁰_S, so the velocity and u⁰ solves pass n_Q = J⁰_Q √(1−v²) (or J⁰_Q/u⁰) and n_S to the EoS, as `Reconst` does. |
| `src/device/device_eos.h` | EoS type 3, `EOS_4D`: `eos4DTilde`, `eos4DDilute`, `eos4DSetup` (= `FourDSetup`), `eos4DEval` (= `FourDEval`), the pressure, the gradients with c_s², and the thermal variables. `eosPressure` takes n_Q and n_S. |
| `src/device/device_visc.h` | n_Q, n_S in the QuestRevert pressures and in the q of reverted cells |
| `src/device/device_hydro.*` | the sources per stage (`computeSource`, the loop of `Advance::AdvanceIt`); EoS 20 in `fillEos`; the staging copy of the EoS tables is released once the backend has its own |
| `device_api.h`, `metal_*`, `cuda_backend.cu` | one call per RK stage (`stage(rk)`) instead of one per step; `source()`, `uploadSource()` and `fetch()`; `state(2)`, the predictor |
| `src/eos_4D.h` | read-only accessors for the tables, the grid, α_Nf and the dilute reference |
| `src/evolve.cpp` | passes the source object to `DeviceHydro` |

## 2. Source terms

### Design

- **Why the host:** `Advance::FirstRKStepT` evaluates `get_hydro_energy_source` and `get_hydro_rhob_source` / `get_hydro_rhoQ_source` for each cell and stage. It uses the cell's position, τ of the stage and the u of the stage's input slice. The four sources differ a lot:
  - **TATB** reads TA, TB at one step and returns 0 at the others.
  - **Strings** loop over the strings of the current τ frame.
  - **AMPT** loops over partons.
  - **T^{μν}** reads a grid.

  Strings and AMPT use u in the n_B source. A device version of each would be a port of its own, so the device takes the source values as an input instead.
- **Per stage:**
  1. The host evaluates qi_source in an OpenMP loop over the cells, with the u of the device slice: τ_n in the first stage, the predictor in the second.
  2. It uploads the 7 components.
  3. The device adds qi_source Δτ where `FirstRKStepT` does.

  If every value is zero (TATB away from its deposition step, strings outside their τ range), nothing is uploaded and the kernel skips the add.
- **Stages:** a step is now two `stage(rk)` calls, so the host can read the predictor in between.
  - **Metal:** the first stage is committed without waiting. `fetch()` waits for it, and its u is already in shared memory.
  - **CUDA:** `fetch()` copies u (4 of the 31 components) of the predictor to its pinned host buffer.
  - **Without sources:** Metal commits the two stages back to back as before, with no host work in between.
- **The stage-1 τ:** the second stage evaluates the sources at τ_n + Δτ, as `AdvanceIt` does with `tau_rk`.

### Validation

- **Unit test** `Device kernels reproduce the update of Advance with sources`:
  - a smooth energy-momentum and n_B source that depends on u;
  - ideal and viscous, 3D and boost-invariant, 4 steps;
  - double: 0 differences; float: 1e-6 of e (L1);
  - the source changes the state by up to 70, far above the tolerance.
- **Integration test** `gpu_backend_sources` (new): backends 1 and 2 against 0, conservation series relative to the largest E of the run (the strings run starts from an empty grid).

  | event | backend 1 (g++ build) | Metal | Metal N_B |
  |---|---|---|---|
  | lumpy TATB source (112), ideal, 41×41×17 | 1.3e-16 | 7.5e-7 | 6.4e-8 |
  | the same, viscous (shear + bulk + second order) | 1.3e-16 | 7.1e-7 | 6.9e-8 |
  | strings (13), viscous, 41×41×49, τ 0.4 → 2.4 | 8.4e-13 | 5.8e-6 | 7.6e-6 |

  - **The strings event:** 40 copies of the three strings of `tests/test_source_terms/strings_event_0.dat`, moved in x, y and boosted in η, with collision times τ_0 ∈ [0, 0.4] fm/c so that deposition runs through the whole evolution. One string carries its baryon number at η_s = −0.4, so 13 units of N_B are deposited at mid-rapidity through the u-dependent n_B source.
  - **Backend 1:** 8e-13 there is g++ contracting multiply-adds differently in the two code paths, as in phases 1 and 2.

### Cost

| event | CPU code | Metal total | Metal kernels | host sources |
|---|---|---|---|---|
| TATB source, viscous | 0.8 s | 0.3 s | 0.10 s | 0.10 s |
| strings, viscous (101 steps) | 5.8 s | 3.7 s | 0.75 s | **2.24 s** |

- **The strings source:** it costs 11 ms per stage on 10 threads for 82k cells, because every cell loops over the strings of the τ frame. The CPU code pays the same cost, but on Metal it is 3× the kernel time.
- **What would remove it:** a device version of `HydroSourceStrings::get_hydro_energy_source` / `get_hydro_rhob_source`, i.e. a kernel over cells with the frame's string list as a buffer. The source does not depend on the hydro state except through u in the n_B part, so the energy part could also overlap with the kernels.

## 3. The 4D EoS (EoS 20)

### Design

- **Tables:**
  - `EOS_4D` keeps P, T, μ_B, μ_Q and μ_S as float tables on a grid in (T̃, μ̃_B, μ̃_Q, μ̃_S), 18.4M points (73.5 MB) each.
  - The device gets the same floats in both precisions, after `T5diff_` in the Real table buffer, and reads them through a float pointer.
  - Ideal runs upload P only (73.5 MB); viscous runs all five (367 MB).
- **Arithmetic:**
  - The host mixes float (weights, table values, interpolation, the μ̃ derivatives) and double (T̃ and μ̃ from e and n, the T̃⁵ slope, the chain rule to e and n).
  - The device repeats each operation in the same type, with Real for the host's double. That makes the kernels on the CPU in double exact, and on Metal everything is float.
- **The dilute region** (T̃ < 40 MeV): the same rescaling as `EOS_4D::diluteScale` and the same derivative formulas.
- **c_s²:** from the pressure derivatives, with e + P + ε as in `EOS_4D::get_pressure_with_gradients_and_cs2`.
- **The thermal variables of the viscous update** (`EOS_4D::getThermalVariables`): T, μ_B, μ_S, μ_Q from the tables with one set of weights, s = (e + P − μ_B n_B − μ_S n_S − μ_Q n_Q)/T.
- **n_Q and n_S:**
  - EoS 20 is the first device EoS that depends on them.
  - The pressure kernel, the face reconstruction, the velocity and u⁰ solves, QuestRevert and the q of reverted cells now pass them.
  - Every other device EoS ignores them.
- **Memory:** `DeviceHydro` used to keep its staging copy of the tables, which would be 367 MB more for EoS 20. It now frees that copy after `setEos`.

### Validation

- **Unit test** `Device kernels reproduce the update of Advance with the 4D EoS` (needs `EOS/neos4D`):
  - the lumpy state with n_Q = 0.4 n_B (1 + 0.3 sin) and n_S = 0.1 n_B cos, and cells at e = 1e-7 (the dilute branch);
  - ideal (3D mode 0, boost-invariant mode 1) and viscous (configurations 0 and 3 of phase 2), 4 steps;
  - double: 0 differences;
  - float, relative L1 difference of e: 3–9e-7.
- **Integration tests** `gpu_backend` and `gpu_backend_viscous` now include EoS 20 (lumpy event with n_B):

  | | ideal | viscous |
  |---|---|---|
  | Metal, conservation series / E | 9.3e-7 | 6.4e-7 |
  | Metal, W/(e+P), Π (hot cells, L1) | — | 5.1e-6, 3.1e-6 |
  | backend 1 (g++ build) | — | 1.9e-11 |
  | max \|E/E₀ − 1\| CPU / Metal | — | 4.8e-6 / 4.9e-6 |

### Speed

The phase 2 benchmark event: lumpy 3D with n_B, 161×161×49 = 1.27M cells, shear + bulk + second order, 101 steps, no freeze-out. The CPU code ran on 10 threads (g++-15 -O3). Both runs were repeated, and each agreed with its repeat to 0.2%.

| | EoS 9 (phase 2) | EoS 20 |
|---|---|---|
| CPU code | 60.9 s | 114.8 s |
| Metal wall | 10.8 s | **16.5 s** |
| Metal kernels | 7.3 s | 11.0 s |
| per cell and RK stage | 28 ns | 43 ns |
| speed-up (wall) | 5.6× | **6.9×** |

- **Cost of the 4D lookups:** each one reads 16 scattered corners of a 73.5 MB table. On the CPU they almost double the run time over EoS 9; on Metal they add 50% to the kernels. The viscous update reads five tables per cell.
- **Agreement on this event:** the conservation series of the two runs agree to 1.2e-6 of E; max \|E/E₀ − 1\| is 2.2e-7 (CPU) and 1.2e-6 (Metal).

## 4. CUDA

- **Changes:**
  - the per-stage entry point;
  - the source buffer (device and pinned host);
  - `fetch()` with a pinned copy of the predictor;
  - the 4D tables in the table buffer.

  It type-checks with clang's CUDA front end in the host and device passes.
- **Not run:** there is no NVIDIA hardware here.
- **Memory:** a viscous EoS 20 run needs about 367 MB of tables plus 31 × 8 bytes per cell and slice in double.

## 5. Tests

- **Unit tests (new):**
  - `Device kernels reproduce the update of Advance with sources`;
  - `Device kernels reproduce the update of Advance with the 4D EoS`.
- **Integration tests:**
  - `gpu_backend_sources` (new);
  - EoS 20 added to `gpu_backend` and `gpu_backend_viscous`.
- **Regression:**
  - all unit tests pass: `unittest_bugfix` 38 cases, 516278 assertions (clang), plus the cell, grid, minmod and reconst suites;
  - all 22 integration tests pass with the Metal build (g++-15). The memory-history driver comes from the unit-test build, because `integration_hydro_memory.e` is only built with `unittest=ON`.

## 6. Next steps

1. The strings source on the device (section 2): now the largest host cost in a strings run on Metal.
2. Baryon diffusion and the chemical bulk pressure; viscous hydro with NEOS (10–14), which needs μ_B, μ_Q, μ_S tables in `eosThermal`.
3. `FlagCausalityMethod = 1`.
4. Fewer host copies: download only on output and freeze-out steps.
5. CUDA validation and timing on NVIDIA hardware.
