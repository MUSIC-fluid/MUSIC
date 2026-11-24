# Denominator safety pass (branch `gpu`, on top of f581a30)

Goal: no division by 0 in any case. Many denominators of the form
`x + abs_err` (or `+ small_eps`, `+ 1e-16`) were meant to be safe, but are
not for every input: `x` can be negative (then `x + eps` can be 0), or the
guard was missing. Every division in `src/` (host and device) was reviewed.

## Conventions

1. **State-dependent denominators** (e + P, T, Ttilde^2, n_B^2, a bracket
   width, a secant step, ...) get a floor or a guarded branch:
   `Util::tiny_normal` (2.2250738585072014e-308) on the host,
   `tinyNormal<Real>()` (FLT_MIN on Metal) on the device
   (`device_defs.h`).
2. **Rule:** a guard changes no value that the old code computed without
   dividing by 0. It acts only where the old denominator was 0 (or
   subnormal, or the old result was inf/NaN). `x + eps` with x ≥ 0 by
   physics becomes `max(x, 0) + eps`, the same value wherever x ≥ 0.
3. **Input parameters** that the code divides by are checked at startup,
   and a bad value stops the run with `exit(1)` and a message.
4. `u0^2 - u3^2 = 1 + u1^2 + u2^2 ≥ 1`: floored at 0.5, not 1. Float
   rounding takes it slightly below 1, and a floor of 1 changed Metal
   results (one freeze-out cell at the e = 0.1 cut). It reaches 0 only for
   u0 ~ 1/sqrt(eps).

## Changes by file

| File | Change |
|---|---|
| reconst.h/.cpp, device_ideal.h | `secant_denominator(d)`: `d + abs_err`, or `-abs_err` (= d, the plain secant step) where `d + abs_err` is 0 or subnormal; used at both v-solver secant sites. A zero slope in the v and u0 hybrid solves forces bisection; a zero slope in the Newton solves fails them (the fallback solver runs). u0 hybrid: a step that does not move the root stops the solve and fails it unless converged, as the old `0/0 = NaN` did. `reconst_u0_f`/`u0F`: `sqrt(max(temp1, abs_err^2))`. `rel_error_v` uses `max(v, 0)`. |
| advance.cpp, device_visc.h | `scaleT`, `eq_size`, `factor`, `rho_q` floored (an empty cell then gives rho = inf, i.e. full removal, instead of 0/0); enthalpy floored in the causality constraints; the `alpha` of n6 only where its denominator is < 0; R_shear/R_bulk guard uses `max(e+P, 0) + 1e-16`; the cubic root `r/(q sq)` guarded; `u0^2 - u3^2` floor 0.5 |
| u_derivative.cpp | `u0^2 - u3^2` floor 0.5; `muB/T` with T ≥ small_eps; enthalpy floored; the d/dtau causal solve returns (keeps the backward difference) if `cs^2 v^2 ≥ 1` |
| dissipative.cpp | R numbers; kappa/eta/zeta output tables: T and T·s floored |
| eos_base.cpp/.h, eos_4D.cpp, eos_idealgas.cpp, device_eos.h | entropy, cs^2, enthalpy, `/T` guards with `max(.,0) + small_eps`; Ttilde^2 and `scale^2` floored (e = 0); `1/eLocal` floored; bisection stop tests `(hi - lo) > rel·mid` instead of `/mid`; ideal-gas muB/T Newton breaks to the robust solve on a slope ≤ 0 |
| evolve.cpp | coarse-grain weights `dampSum`, `wSum` guarded; freeze-out `/T` floored |
| grid_info.cpp | W/(e+P), inverse Reynolds, meanpT estimators, Gubser autocheck guarded (output only) |
| cornelius.cpp | centroid of a degenerate (zero-weight) element = the mean of its points |
| HydroinfoMUSIC.cpp | eta = 0 outside the light cone (t ≤ \|z\|) |
| init.cpp, hydro_source_Tmunu.cpp | `utau^2 - ueta^2` floor 0.5; renormalisation `max(e+P,0) + 1e-16` |
| hydro_source_ampt.cpp | rapidity with `E - pz` floored; `pT/m` with m floored |
| hydro_source_strings.cpp | beam rapidity from `acosh(E/(N_B m_N))` only for N_B > 0 and E ≥ N_B m_N; otherwise a warning and the ecm value |
| read_in_parameters.cpp | new `check_parameter_denominators` (below) |
| transport_coeffs.cpp/.h | new `check_denominators`: relaxation-time factors > 0, increasing T grid of the viscosity file, muB scale and bulk widths > 0, \|lambda_asymm\| < 1 |

Startup checks in `check_parameter_denominators`: nx, ny ≥ 2; x_size,
y_size > 0; eta_size > 0 unless boost invariant; delta_tau > 0;
dtaudxRatio > 0; tau0 > 0 (except profiles 13, 131, 30); profiles
11/111/112/113: beam rapidity > 0 (ecm > 2 m_N) and both rhob widths > 0
with turn_on_rhob; eta_fall_off > 0 for profiles 9, 91–95; 112/113:
nucleonWidth > 0 and eta_m > 0 if tilted_fraction > 0; 13/131: string
sigma_x, sigma_eta > 0.

Left as they are (safe by construction): TA/TB, the device strings and
device_hydro sites, minmod, the cornelius edge cuts, lin_int, the EoS table
constants, and the kappa_hat table values. **Not covered:** the legacy
spectra code `freeze_pseudo.cpp` and `reso_decay.cpp`.

## Guards revised during testing

Three first versions broke the "change nothing that was safe" rule. Tests
caught each one, and each was narrowed:

1. **Secant denominator.** The first version returned `-abs_err` wherever
   `|d + abs_err| < abs_err`, i.e. for d in (-2·abs_err, 0), and float
   iterates near convergence often land there. It changed EoS 14 float
   results (unittest_bugfix relE 3.2e-5 > 1e-5), Metal visc20/diff14
   observables in the last digits, the CPU visc20 bench from the first
   steps, and beastmode eccentricities at 1e-17. Now only
   `|d + abs_err| < tiny_normal` is replaced, and all of these are
   bit-identical again.
2. **u0 hybrid zero step.** The first version accepted a stalled root; in
   float that changed the Metal TATB/strings source outputs. Now the solve
   fails unless converged, as before.
3. **`u0^2 - u3^2` floor of 1.** Changed by float rounding; now 0.5.

## Impact on the tests (final source)

- Unit tests: all pass (bugfix 1520088/1520088, cell 11, grid 32,
  minmod 6, reconst 112).
- Integration suite: 25/25 pass; all 2251 output files are **bit-identical**
  to the run before the pass (after the section 9 fast division).
- Bench (AuAu 200 GeV, CPU `_0` and Metal `_2a`; ideal9, visc9 at
  41/81/121/161, visc20, diff14, fullFO9): all 16 runs are
  **bit-identical** to the bench18 reference, every output file (28–38
  per case).
- Metal wall times are unchanged (diff14 10.0 s, visc20 10.0 s,
  fullFO9 9.4 s, ideal9 3.2 s, visc9_161 6.4 s). The CPU times in this run
  were taken on a loaded machine and are not comparable.

The commit of this pass was also checked without the section 9 fast
division (doc/gpu_benchmark_6880b71_hotspots.md), that is with plain
division in device_ideal.h: unit tests and integration 25/25 pass, the
2251 integration files and the 8 Metal bench runs are bit-identical to the
same references.

So none of the tested results change. The guards act only at a zero (or
subnormal) denominator, which none of these runs reaches. Inputs that
divided by 0 before now either stop at startup with a message (bad
parameters) or give a finite, defined value (empty cells, degenerate
geometry).
