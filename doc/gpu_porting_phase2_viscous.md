# GPU port, phase 2: viscous hydro on the device

Branch `gpu`, from b52193c (phase 1, [gpu_porting_phase1.md](gpu_porting_phase1.md)). 2026-11-22. Machine: Apple M4 (10 CPU cores, 10-core GPU, 16 GB unified memory); no NVIDIA hardware.

## Summary

- **`gpu_backend` now also runs viscous hydro** (`Viscosity_Flag_Yes_1_No_0 = 1`) on the device:
  - shear and bulk, with or without the second-order terms;
  - the vorticity terms, `FlagResumTransportCoeff`, QuestRevert;
  - `Du_from_ideal_EoM` 0 and 1;
  - all bulk relaxation types (0, 1, 2);
  - every η/s and ζ/s parametrization except the tabulated one (type 23): shear T types 0/1/2/3/11 and μ_B types 7/10, bulk types 1/2/3/7/8/9/10;
  - EoS 0 (ideal gas, with n_B) and 9/91 (hotQCD).
- **Not on the device yet** (the CPU code runs, with a warning saying why):
  - baryon diffusion;
  - the chemical bulk pressure;
  - `FlagCausalityMethod = 1`;
  - transport coefficients from files (type 23);
  - viscous hydro with the NEOS tables (10–14).
- **Validation:**
  - **Backend 1 (device kernels on the CPU, double) reproduces `Advance` bit for bit** with clang. The unit test checks eight configurations over 4 steps (EoS 0 and 9, 3D and boost-invariant, both face modes) and finds 0 differing values in e, n_B, u, W^{μν} and Π.
  - **Metal (float) agrees with the CPU code to 8e-7 of E**, and to 2–5e-6 (L1) in W/(e+P) and Π in the cells above 0.1 GeV/fm³. This is on the lumpy 3D event with n_B, shear + bulk + second-order terms, for both EoSs. Energy conservation is the same as on the CPU: 1.8e-6 (EoS 0) and 3.6e-6 (EoS 9).
- **Speed:** see [section 4](#4-speed).
- **Results:** with `gpu_backend = 0` (the default) nothing changes. The only change to host code is five accessors in `EOS_idealgas`. All integration tests pass.

## 1. Code

| file | what |
|---|---|
| `src/device/device_visc.h` (new) | `kernelWFace` (the face values of `Diss::computeWDivergences`) and `kernelCellVisc`, the viscous cell loop of `Advance::AdvanceIt`. The cell kernel runs: `FirstRKStepT` with ∂_m W^{mα}, then `U_derivative` (MakedU, θ, Du, σ, ω), `Make_uWRHS`, the thermal variables, `FirstRKStepW` (`Make_uWSource`, `Make_uPiSource`, QuestRevert, `QuestRevertResummedTransCoeff`), and the stage-2 fixed point. It also holds `TransportCoeffs`. |
| `src/device/device_eos.h` | `eosThermal`, the entries of `getThermalVariables` that the viscous update reads. For the ideal gas this is μ_B/T from its x table with the same Newton steps; for hotQCD, T from the T⁵ table. |
| `src/device/device_ideal.h` | `kernelCell` split into `cellQiFlux` and `cellRkAverage`, which the viscous kernel reuses. `StepParams` gains `nComp` and the viscous parameters, including g^{μν}, read at run time as on the host. |
| `src/device/device_defs.h` | the slice layout of viscous runs; `realMax`, `dexpm1`; more math functions |
| `device_hydro.*`, `device_api.h`, `metal_*`, `cuda_backend.cu` | slices of `nComp` components; `state(which)` for the current and the previous slice; the W face buffers; two more kernels per stage |
| `src/eos_idealgas.h` | read-only accessors for the x table and the coefficients |

### Slice layout and the previous slice

- **Layout:** a viscous slice has 31 components: the 15 of phase 1 (e, n_B, n_Q, n_S, u^μ, stored τT^{τμ}), then the 14 W^{μν} in `Fields::Wmunu_` order, Π and Π_chem.
- **The previous slice:** the viscous update takes ∂_τ u, ∂_τ W and ∂_τ Π against `arenaFieldsPrev`. That slice is τ_{n−1} in the first stage and τ_n in the second.
- **No fourth slice is needed.** `Evolve::AdvanceRK` writes τ_{n+1} into the τ_{n−1} slot, and the device's pointer swap does the same. So at the start of a step the backend's `next` slice already holds τ_{n−1}.
- **Uploads:** when the device copy is stale, `fpPrev` is uploaded with `fpCurr`. After each step W and Π are copied back with the primitive variables, so freeze-out, the output and the conservation diagnostic see the same fields as before.

### Kernels per RK stage

| kernel | threads |
|---|---|
| pressure | cells |
| KT face (×3 directions) | faces |
| W face (×3): KT fluxes of u^i W, face averages of τW^{iα} and of the Π part | faces |
| viscous cell | cells |

The cell kernel does everything that is local to one cell. That includes the stage-2 fixed point (up to 16 reconstructions and W updates), so nothing goes back to the host within a step.

## 2. Bit-for-bit reproduction

As in phase 1, every device function copies its host version **statement by statement**: clang fuses `a*b + c` only within a statement.

Two details mattered:

- **g^{μν} is read at run time.** `DATA.gmunu` is run-time memory on the host. If the device saw it as constants, the compiler could fold the zero and unit entries and change which products are fused. `StepParams::g` therefore carries it.
- **The stage-2 shift is applied the host's way.** The fixed point applies its shift to a copy of q_i and reconstructs from the stage's cell, as `applyShift` does. The stored τT^{τμ} (float runs) takes the q_i of the accepted shift.

**Unit test result.** `Device kernels reproduce the viscous update of Advance` (unittest_bugfix) runs 4 steps of a lumpy viscous state. The state is a traceless, transverse W of 3–8% of e, Π = −0.1 e, and a previous slice that differs by a few %. The configurations are:

| config | grid, face mode | options |
|---|---|---|
| 0 | 3D, mode 0 | second order, vorticity, resummation, QuestRevert, η/s type 3 with μ_B type 7, ζ/s type 10, bulk relaxation 0 |
| 1 | boost-invariant, mode 1 | first order, η/s type 1, ζ/s type 1, bulk relaxation 2, `Du_from_ideal_EoM = 0` |
| 2 | 3D, mode 0 | as 0 with η/s type 2, ζ/s type 3, bulk relaxation 1 |
| 3 | 3D, mode 1 | shear only (5 W components) |

For EoS 0, and for EoS 9 when `HYDROPROGRAMPATH` points to the tables: **maxDiff = 0 in all eight runs.**

With g++ (the integration test), backend 1 differs from the CPU code by 1.3e-16 of E and by 0 in the printed W and Π. The cause is g++'s cross-statement contraction (see phase 1).

## 3. Metal (float)

### Accuracy

The lumpy 3D event is 41×41×17 cells with n_B, shear (η/s = 0.12) + bulk + second-order terms, run for 1 fm/c (`gpu_backend_viscous` integration test):

| EoS | max \|Δ(E, P, N_B)\| / E | L1 Δ(W/(e+P)), e > 0.1 GeV/fm³ | L1 ΔΠ | max \|E/E₀ − 1\|: CPU, Metal |
|---|---|---|---|---|
| 0 | 8.3e-7 | 2.0e-6 | 3.1e-6 | 1.86e-6, 1.81e-6 |
| 9 | 6.6e-7 | 5.3e-6 | 2.6e-6 | 3.50e-6, 3.60e-6 |

The largest single-component differences of W/(e+P) are 1e-4 to 3e-4 of its maximum, in the first 0.2 fm/c, falling to ~4e-5 by τ = 2 fm/c.

### Float and QuestRevert

The unit test's state is far harsher than an event: cells with e = 1e-7 next to e ≈ 2, and u^τ up to 18. On it, float gives:

- e to 1e-6 and W to ~2e-6 in every configuration without QuestRevert;
- W to 3e-4 (L1) with QuestRevert.

Switching the options off one by one isolates QuestRevert: off, W agrees to 5e-6 in every configuration. The error comes from cells that start nearly empty and fill within the four steps. QuestRevert caps their W at |W| ~ e. There, π_{μν}π^{μν} is a sum of terms of size (u^τ)⁴ π_LRF² that cancel to π_LRF², so float loses ~(u^τ)⁴ of its digits. Computing it with error-free products and sums did not help, so the inputs, not the arithmetic, limit it. The unit test's float limit for W is 1e-3, with this reason next to it. The realistic event above is at 2–5e-6.

**One float-specific change** is in QuestRevert's strength factor. The host writes it as the difference of two logistic functions:

`1/(exp(−(e−0.1)/0.05)+1) − 1/(exp(2)+1)`

In dilute cells both terms are 0.119, and the difference cancels log10(0.06/e) digits: 4e-2 relative error at e = 1e-7 in float. In float the device uses the same difference without the cancellation:

`−e^{2} expm1(−e/0.05) / ((1+e^{−(e−0.1)/0.05})(1+e^{2}))`

This is accurate to 1e-7 for all e. `expm1` is Kahan's (u−1)·y/log u, since MSL has none. In double the host formula is kept, so double stays bit-identical.

### Porting fixes

- **Address spaces:** MSL requires them on references to the parameter block (`constant const StepParams<float> &P`).
- **Double literals:** MSL has no double, so `numeric_limits<double>::max()` (the fixed point's initial "best miss") and the 1e-300 floor of `x_of_rho` need per-target values (`realMax`, `#if __METAL_VERSION__`).

## 4. Speed

### Setup

The phase 1 benchmark event, now viscous:

- lumpy 3D event with n_B, 161×161×49 = 1.27M cells;
- shear (η/s = 0.12) + bulk + second-order terms;
- Δτ = 0.01, 101 steps, mode 1, no freeze-out.

The CPU code ran on 10 threads (g++-15 -O3). Metal repeated to 1%.

### Results

| | EoS 0 wall | EoS 0 kernels | EoS 9 wall | EoS 9 kernels |
|---|---|---|---|---|
| CPU code (`gpu_backend 0`) | 56.1 s | — | 60.9 s | — |
| device kernels on CPU (1) | 57.8 s | 53.7 s | 55.2 s | 51.4 s |
| Metal M4 (2) | **11.6 s** | **8.0 s** | **10.8 s** | **7.3 s** |

- **Speed-up:** Metal is **4.8× (EoS 0) and 5.6× (EoS 9) faster in wall time**, and 6.7× and 7.0× faster in the hydro kernels than the device kernels on the CPU.
- **Per cell and RK stage:**
  - Metal: 31 ns (EoS 0) and 28 ns (EoS 9). The ideal update took 12 ns.
  - Device kernels on the CPU: 209 ns and 200 ns.
- **Where the time goes:** the viscous cell kernel is ~2.5× the cost of the ideal one. It does the gradients, the sources, the transport coefficients and, in stage 2, the fixed point.
- **Host copy:** 1.0 s, i.e. 10 ms per step. That is 13% of Metal's kernel time, twice phase 1's because W and Π are copied too.
- **Outside the hydro kernels:** about 2.6 s, as in phase 1 (initialisation and host diagnostics).
- **Conservation on this event:** max |E/E₀ − 1| is 8.9e-7 (EoS 0) and 1.1e-6 (EoS 9) on Metal, against 3.0e-7 and 2.2e-7 on the CPU. The device kernels on the CPU match the CPU code to 2e-17. The time series agree to 8.8e-7 (EoS 0) and 1.1e-6 (EoS 9) of E.

### Next gains

- Skip the host copy on steps that neither freeze-out nor output read.
- Merge the W face values into the KT face kernel, which reads the same cells.

## 5. CUDA

- **Changes:** the backend has the two new kernels and per-slice pinned host copies. It type-checks with clang's CUDA front end, host and device passes (see `cuda_backend.cu`).
- **Not run:** there is no NVIDIA hardware here.
- **Launch parameters:** `StepParams` is passed by value and is now about 0.5 kB, well under the 4 kB limit.

## 6. Tests

- **Unit test** `Device kernels reproduce the viscous update of Advance` (new):
  - double: exact, 4 configurations × 2 EoSs;
  - float: e < 1e-5 and W, Π < 1e-3 (L1).
- **Integration test** `gpu_backend_viscous` (new): backends 1 and 2 against 0, for EoS 0 and 9.
  - conservation-law time series: limits 1e-6 (backend 1) and 1e-4 (Metal);
  - E conservation: below 2e-5;
  - W/(e+P) and Π in the hot cells at the last output slice: limits 1e-5 and 1e-4.
- **Regression:** all unit tests pass, and all 21 integration tests pass with the Metal build (g++-15).

## 7. Next steps

1. Baryon diffusion (`Make_uqSource`, QuestRevert_qmu, the μ_B/T gradients) and the chemical bulk pressure.
2. Viscous hydro with NEOS: μ_B, μ_Q, μ_S and the table densities in `eosThermal`.
3. Source terms (TATB, strings) as per-step device input.
4. Fewer host copies: download only on output and freeze-out steps (W makes the copy 2× larger).
5. CUDA validation and timing on NVIDIA hardware.
