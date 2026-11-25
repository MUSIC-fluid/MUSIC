# Temperature lookups of the tabulated EoS

2026-11-25, branch `gpu` at 4ceacd5. The profiles come from the
10-thread full event of `doc/openmp_cpu_scaling.md`: 96 x 96 x 64 cells,
EoS 14 (NEOS-BQS), shear + bulk + baryon diffusion.

## Why the tables store T^5

`EOS_hotQCD` and `EOS_neos` store T^5 in place of T and return
`pow(T5, 0.2)` (`eos_hotQCD.cpp:77`, `eos_neos.cpp:194` and `:227`).
This came in with 52fb69c and a4d1475 (2020-11-05). The commit message
says that linear interpolation of T made e/T^4 against T non-monotonic
at low energy density, and that T^5 against e "has a better behavior at
low energy density, even extrapolating to zero energy density".

There are two effects:

1. **Between the table points.** The tables are uniform in e, and T(e)
   is strongly concave (roughly e^{1/4} in the plasma, and steeper in
   the hadron gas). Linear interpolation of T underestimates T between
   the points. T^5 is much closer to linear in e.
2. **Below the first table point.** `interpolate1D` and `interpolate2D`
   scale the stored value as e / e0 below e0. With T stored, T would go
   to 0 linearly in e. With T^5 stored, T goes as e^{1/5}, which is
   smooth and stays positive.

### Leave-one-out test on the NEOS-BQS tables

I interpolated every odd e point from its even neighbours, storing T,
T^4, T^5 or T^6. The table shows the largest relative error in T:

| table | e (GeV/fm^3) | n_B column | T | T^4 | T^5 | T^6 |
|---|---|---|---|---|---|---|
| 1 | 3e-4 - 3.3e-3 | 0 | 2.1e-2 | 5.3e-3 | **4.7e-4** | 4.1e-3 |
| 1 | 0 - 3.6e-3 | 125 | 4.2e-1 | 4.4e-2 | **2.2e-2** | 2.2e-2 |
| 1 | 0 - 3.6e-3 | 250 | 4.3e-1 | 6.2e-2 | 2.9e-2 | 5.8e-3 |
| 2 | 3.6e-3 - 1.4e-2 | 0 | 1.7e-3 | 5.8e-4 | 2.2e-4 | 1.4e-4 |
| 3 | 0.015 - 0.045 | 0 | 3.1e-4 | 1.4e-4 | 8.3e-5 | 2.7e-5 |
| 4 | 0.045 - 0.445 | 0 | 2.4e-3 | 1.2e-3 | 8.5e-4 | 4.8e-4 |
| 5 | 0.455 - 20.3 | 0 | 1.2e-3 | 3.0e-4 | **1.0e-4** | 2.9e-4 |
| 6 | 20.4 - 218 | 0 | 1.8e-4 | 8.6e-6 | 4.9e-5 | 1.1e-4 |

At finite n_B in tables 2 - 5 the error is about the same for every power
(it comes from the n_B structure). T^5 is the best overall choice, and it
matters most in the dilute hadron gas, where T stored directly is off by
up to 40 % between the points. **Keep T^5.** There is no exact cheaper
replacement for the 5th root: a fast root or `exp(0.2 log x)` changes the
results at the rounding level for little gain.

## Where the time goes

Temperature and mu_B lookups, including `pow`, take **14.6 %** of all
thread samples. By caller:

| caller | share | what it does |
|---|---|---|
| `U_derivative::MakeDSpatial` | 8.3 % | mu_B / T at the cell and both neighbours, in each of the 3 directions (diffusion) |
| `Diss::D_muB_over_T_ideal` | 3.3 % | central differences of mu_B / T in e and n_B: 4 (T, mu_B) pairs |
| `U_derivative::MakeDTau` | 1.2 % | T and mu_B of the cell in the current and previous slices |
| rest (thermal variables, grid maxima, ...) | 1.8 % | |

The cost is in the **repetition**, not in the 5th root. In each RK stage
`MakeDSpatial` evaluates `get_muB(...) / max(get_temperature(...), eps)`
9 times per cell. The centre value is computed 3 times, and each
neighbour value is recomputed by the 6 cells around it. Each evaluation
finds the table index and the interpolation weights twice (once for
mu_B, once for T) and calls `pow`. With `output_vorticity 1`, the same
loop also computes T at three points for each of the 4 components of u,
which is 36 more T lookups per cell.

## Done: mu_B / T cached once per slice (option 1)

With viscous hydro and baryon diffusion on (`turn_on_diff 1`),
`Advance::AdvanceIt` fills `muBOverTCurr_` and `muBOverTPrev_`. These
hold `get_muB / max(get_temperature, small_eps)` of every stored cell of
the current and previous slices, ghost cells included. They are filled
in one parallel pass each (`Advance::computeMuBOverT`). In the second RK
stage the previous slice is the first stage's current slice, so its
array is swapped in instead of being recomputed (as `pressurePrev_` is).
`U_derivative::MakedU` passes the arrays on: `MakeDSpatial` reads the
centre and neighbour values from them, and `MakeDTau` reads the two
slice values and skips its T lookups unless the vorticity output needs
T. Without the arrays (other callers, unit tests) both still call the
EoS. The device kernels do the same (see below).

Each step now fills 3 arrays per cell (2 in the first stage, 1 in the
second), against 22 (mu_B, T) evaluations per cell before (9 in
`MakeDSpatial` and 2 in `MakeDTau`, in each of the 2 stages).

| run (96 x 96 x 64, visc + diffusion) | before (s) | after (s) |
|---|---|---|
| 38 steps, 10 threads | 15.6, 15.8 | 14.2, 14.3 |
| 38 steps, 1 thread | 93.6 | 83.0 |
| full event, 295 steps, 10 threads | 130.6 | 120.6 |

That is 8 - 11 % faster. The conservation laws, eccentricities, mean p_T,
momentum anisotropy, inverse Reynolds numbers and the freeze-out surface
(4.8 MB, 1 thread) are bit-identical. Unit tests and the 24 integration
tests pass.

### Memory

There are two doubles per stored cell, ghost cells included, and only
with diffusion on. For 128 x 128 x 96: (128 + 4) x (128 + 4) x (96 + 4)
= 1,742,400 stored cells, so 13.3 MiB per array and **26.6 MiB in all**.
Measured peak RSS for the 17-step 128 x 128 x 96 event of
`doc/openmp_cpu_scaling.md` (EoS 14, shear + bulk + diffusion, CPU):
2149 MB before, 2176 MB after (+27 MB, +1.3 %).

## Done: the same on the device backends

`kernelCellVisc` (`device_visc.h`) used to evaluate mu_B / T 11 times per
cell and stage: the centre and both neighbours in 3 directions, and the
cell in the current and previous slices. It now reads them from two
arrays, `aCurr` and `aPrev`, and passes them to `ViscSolver::gradients`.
A new kernel, `kernelMuBOverT`, fills them with the same expression
(`eosTMuB`, then `muB / max(T, 1e-16)`).

Each backend keeps two buffers, `aCur` and `aOld`, of n = nx ny neta
reals. The device slices have no ghost cells. The buffers are allocated
on the first stage with diffusion (`StepParams::diff`) and freed on a
resize.

* First stage: fill `aCur` from the input slice and `aOld` from
  tau_{n-1}.
* Second stage: fill `aOld` from the input slice. The previous slice is
  the first stage's input, so its values are already in `aCur`.

Without diffusion no fill kernel runs, and the cell kernel gets the
pressure buffer as a placeholder that it never reads.

* **CPU backend (gpu_backend 1):** two `omp for` loops after the
  pressure (`device_hydro.cpp`).
* **Metal:** `musicMuBOverT`, with buffers 17 and 18 of the cell
  kernels. The stage reads `diff` from the raw parameters at
  `MUSIC_DEVICE_DIFF_OFFSET_FLOAT` (`device_api.h`, checked by a
  `static_assert` in `device_hydro.cpp`, like the rkFlag and strings
  offsets).
* **CUDA:** `cudaMuBOverT`, timed as `muBOverT` in `MUSIC_CUDA_PROFILE`.
  **Not compiled or tested: there is no nvcc on this Mac.**

### Results (96 x 96 x 64, visc + diffusion, EoS 14)

* **Bit-identical outputs**, before against after, on gpu_backend 1
  (double) and on Metal (float, also over the full event): conservation
  laws, momentum anisotropy, eccentricities, mean p_T and inverse
  Reynolds numbers. All 24 integration tests pass, and the
  device-against-CPU differences of the diffusion test are unchanged.
* **Metal, full event (295 steps):** wall time 13.6 s against 12.5 s
  (two runs each). GPU time per kernel (`MUSIC_METAL_PROFILE`, 38 steps):

| kernel | before (s) | after (s) |
|---|---|---|
| cellVisc_rk0 | 2.932 | 2.089 |
| cellVisc_rk1 | 4.727 | 4.294 |
| muBOverT | | 0.075 |
| face, wface, pressure | 5.069 | 5.035 |
| total | 12.73 | 11.49 |

* **gpu_backend 1 (10 threads, 38 steps):** 17.2 s against 16.3 s.

### Device memory

The cost is two arrays of n reals, only with diffusion. For
128 x 128 x 96 (n = 1,572,864):

| backend | extra memory |
|---|---|
| Metal (float) | 12.0 MiB |
| CUDA float | 12.0 MiB |
| CUDA double, gpu_backend 1 (double) | 24.0 MiB |

GPU runs do not use the host arrays of `Advance`, so this is the whole
increase for them. The measured peak RSS of the 128 x 128 x 96 Metal
event is 2281 MB both before and after (two runs each). The Metal
buffers are in unified memory, but they are allocated after the peak,
which is set during the initialisation. Host CPU runs (gpu_backend 0)
add the 26.6 MiB above instead.

## Options

1. **Cache alpha = mu_B / T once per slice** (done, see above). Fill an array
   of alpha for the current slice, ghost cells included, in one parallel
   pass, the same way `Advance::computePressure` fills the pressure.
   `MakeDSpatial` and `MakeDTau` then read it. The previous slice's alpha
   is the current one of the step before (keep it with the same owner and
   swap logic as `pressureCurr_` / `pressurePrev_`). The arithmetic is the
   same function on the same inputs, so the results are bit-identical.
   This removes about 9.5 % of the thread samples with diffusion on, at
   the cost of one extra double per cell per slice. A T array would do
   the same for the vorticity output.
2. **Combined T and mu_B lookup.** A `get_T_muB` that finds the table
   index and weights once and interpolates both tables would roughly
   halve the remaining lookups (`D_muB_over_T_ideal`, the thermal
   variables). It is bit-identical if it reuses the same arithmetic.
3. **Analytic d(mu_B/T)/de and d(mu_B/T)/dn in D_muB_over_T_ideal**, from
   the slopes of the table cell instead of 4 finite-difference lookups.
   This changes the results, so it needs validation.

Options 2 and 3 also apply to the device kernels: `ViscSolver::dMuBOverTIdeal`
(`device_visc.h`) still makes its 4 finite-difference `eosTMuB` calls.
