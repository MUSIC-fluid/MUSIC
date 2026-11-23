# GPU benchmark at 6880b71 and runtime hotspots (Apple M4, Metal)

Branch `gpu`, commit 6880b71. Metal backend in float (MTLMathModeSafe),
CPU reference = backend 0. All probes used for this study were temporary
and have been removed; the tree is unchanged.

## 1. Benchmark (wall time, Metal best of 2)

| case        | CPU     | Metal   | speed-up |
|-------------|---------|---------|----------|
| visc9_41    | 1.5 s   | 0.32 s  | 4.8x     |
| visc9_81    | 7.8 s   | 1.09 s  | 7.1x     |
| visc9_121   | 24.9 s  | 3.09 s  | 8.1x     |
| visc9_161   | 59.8 s  | 6.79 s  | 8.8x     |
| ideal9_161  | 39.8 s  | 3.30 s  | 12.0x    |
| visc20_161  | 116.0 s | 11.36 s | 10.2x    |
| diff14_161  | 100.7 s | 10.27 s | 9.8x     |
| fullFO9_101 | 99.7 s  | 9.97 s  | 10.0x    |

Agreement with the CPU: dE/E 0.8-2.6e-6, dN_B/N_B 5e-8 to 3.5e-7.
Freeze-out surface (fullFO9_101): the same 370,636 cells, with relative
differences of 8.3e-8 in the energy and 6.7e-8 in the entropy.

## 2. Where the GPU time goes (MUSIC_METAL_PROFILE=1, GPU seconds)

| case        | face          | cellVisc_rk1  | cellVisc_rk0  | wface         | other |
|-------------|---------------|---------------|---------------|---------------|-------|
| visc9_161   | 1.92 (30%)    | 1.96 (31%)    | 1.25 (20%)    | 1.03 (16%)    | 0.22  |
| visc20_161  | 4.15 (41%)    | 2.95 (29%)    | 1.65 (16%)    | 1.03 (10%)    | 0.3   |
| diff14_161  | 2.71          | 3.48 (36%)    | 2.01          | 1.35          |       |
| fullFO9_101 | 3.21          | 3.35          | 1.52          | 1.23          |       |
| ideal9_161  | 2.01 (71%)    | cell 0.70 (25%) |             |               |       |

The "other" column is mostly conservation (0.18 s) and pressure (0.04 s).

For visc9_161 the kernels add up to 6.4 s of the 6.8 s wall time. The
host is not the bottleneck: per-step host work overlaps the GPU, and the
GPU sits idle about 3 ms per step. The exception is freeze-out steps,
which are close to host-bound (22 of 24.6 ms; see
gpu_host_device_communication.md).

### Xcode Instruments

Metal System Trace with GPU counters (visc9_81) gave no attribution
inside the kernels:
- the shader-profiler sample table is empty (no PC samples);
- the counter set exposed only "RT Unit Active".

The attribution below therefore comes from operation counts and
ablation runs.

## 3. Operation counts per kernel call (backend 1, float, 1.27M cells)

Counters were compiled into the host build only. Each entry is calls
per face or per cell:

| kernel       | EoS 9: pressure / grad / thermal | EoS 20: pressure / grad | primitive solves | velocityF iter. |
|--------------|----------------------------------|-------------------------|------------------|-----------------|
| face         | 3.68 / 2.00 / -                  | 4.52 / 2.00             | 2.00             | 3.1 (EoS 20: 3.9) |
| cellVisc_rk0 | 3.69 / 1.85 / 1.00               | 3.71 / 0.85 / 1.00      | 1.00             | 1.4             |
| cellVisc_rk1 | 7.86 / 1.85 / 1.00               | 7.77 / 0.85 / 1.00      | 1.77             | 2.6             |

Each RK stage has 3 face calls per cell, one per direction. So for
every cell update the face kernel makes about 11 pressure lookups and 6
gradient lookups, against the cell kernel's own 3.7 (stage 1) to 7.9
(stage 2).

The root solve converges in about 1.5-2 iterations, so the initial guess
is already good. More iterations are not the problem; the cost is the
EoS lookups inside each iteration.

## 4. Ablations on Metal (GPU seconds; the outputs stay bit-identical)

### Doubling probes

Each probe computes a piece twice and adds the second result times 0,
so the output is unchanged and the extra time is roughly the cost of
that piece:

| added cost of ...       | face visc9 | face EoS 20 | rk0 visc9 | rk0 EoS 20 | rk1 visc9 | rk1 EoS 20 | face ideal9 |
|-------------------------|-----------:|------------:|----------:|-----------:|----------:|-----------:|------------:|
| baseline                | 1.91       | 4.14        | 1.26      | 1.64       | 1.96      | 2.93       | 2.03        |
| every EoS call          | +0.44      | +2.91       | -0.01     | +1.21      | +0.14     | +1.18      | +0.45       |
| primitive solve         | +1.42      | +3.02       | +0.05     | +0.12      | +0.38     | +0.89      | +1.45       |
| gradients               | -          | -           | +0.67     | +0.75      | +0.53     | +0.30      | -           |
| firstRKStepW (W update) | -          | -           | +0.01     | +0.34      | +0.54     | +0.32      | -           |

### Removal probes for cellVisc_rk0

These runs are for timing only; their outputs are wrong:

| removed                         | visc9          | EoS 20         |
|---------------------------------|----------------|----------------|
| nothing                         | 1.24           | 1.64           |
| gradients                       | 0.82 (-34%)    | 0.97 (-41%)    |
| W divergences (reads wface)     | 0.82 (-34%)    | 1.17 (-29%)    |
| firstRKStepW                    | 1.17 (-6%)     | 1.29 (-21%)    |
| solve                           | 1.26 (0%)      | 1.64 (0%)      |
| all four                        | 0.49           | 0.49           |

### face_reconstruction 1 (existing option: limit e, n, u, so no root solve at the faces)

| case       | face, mode 0 | face, mode 1 | total GPU time, mode 0 → mode 1 |
|------------|-------------:|-------------:|---------------------------------|
| visc9_161  | 1.91         | 0.62         | 6.4 → 5.1 s (-20%)              |
| visc20_161 | 4.14         | 1.32         | 10.1 → 7.3 s (-28%)             |

## 5. Hotspots, in order

1. **The face primitive solve: 74% of face time, and 23-30% of all GPU
   time.**
   - Two solves per face, each making about 2 pressure lookups.
   - With EoS 20 the EoS lookups alone are 70% of the face kernel: the
     4D table interpolation in float, with scattered global loads.
2. **Gradients in both cell stages: 34-41% of rk0, and about a quarter of
   rk1.** They read 3x2 neighbours per direction for u, e, n_B, and the
   neighbour pressures.
3. **The W divergence in cellVisc_rk0, plus the wface kernel.**
   - The divergence is 30-34% of rk0.
   - wface takes 1.0 s regardless of EoS, so it is memory-bound: 10-16%
     of GPU time.
   - Together they move the W face fluxes through memory twice: wface
     writes them and the cell kernel reads them back.
4. **EoS 20 inside cellVisc:** about 1.2 s across rk0 and rk1.
5. **The stage-2 fixed-point loop:** 1.77 solves per cell instead of 1;
   0.37 s for EoS 9 and 0.8 s for EoS 20.
6. **The host and the freeze-out steps.** They are small in the
   benchmark, but freeze-out steps stay host-bound (step 2 of
   gpu_host_device_communication.md).

## 6. Recommended next optimisations

| # | change | expected gain | changes results? |
|---|--------|---------------|------------------|
| A | Reuse the 4D-EoS cell setup (`eos4DSetup`/`eos4DTilde`) between the last pressure call of a solve and the `eosPressureGrad` at the same (e, n): one combined evaluation of p and its gradients | about 1/4-1/3 of the EoS 20 face cost, ~0.7-1 s on visc20_161 | no (same arithmetic) |
| B | Fuse wface into the face kernel. Both read the same L/LL/R/RR cells and are launched over the same faces. | most of wface's 1.0 s of memory traffic: ~0.4-0.7 s | no |
| C | Store the face-centred derivatives (or the stage-1 gradients) needed for the stage-2 gradients. Alternatively, compute the gradients in a tiled kernel with threadgroup memory, so neighbours are loaded once per tile. | 0.3-0.6 s | no |
| D | Use `face_reconstruction 1` for production runs where it is acceptable | -20% (EoS 9) to -28% (EoS 20) of GPU time | yes: a different limiter (physics choice) |
| E | Freeze-out detection on the device (step 2 of the communication plan) | freeze-out steps, about 10% of the fullFO wall time | no |

My first recommendation was A and B. Section 7 shows that neither one
pays off. D is the largest single gain, but it changes the numerics, so
it is a physics decision. Relaxed math stays rejected (section 7 of
gpu_benchmark_and_next_steps.md).

## 7. Tried: A and B (neither adopted; the source is unchanged)

### A: sharing the 4D-EoS setup

**The premise was wrong.** The `eosPressureGrad` after a solve is not
evaluated at a point the solve already evaluated:
- the final v (`vNext`) is never passed to `velocityF`;
- n_B = J^0_B / u^0, with u^0 = 1/(sqrt(1-v^2) + v absErr), differs from
  the J^0_B sqrt(1-v^2) of `velocityF`.

Sharing the setup would therefore change results. It is not the
redundant work the report assumed.

The probe iterations have their minimum already: about 2.25 pressure
lookups and 1 gradient lookup per solve. No point is evaluated twice.

**What the lookup costs** (EoS 20, face kernel, timing-only doubling):

| doubled                                          | face   |
|--------------------------------------------------|--------|
| nothing                                          | 4.16 s |
| the 16 table loads and the interpolation         | 4.65 s (+0.49) |
| `eos4DTilde`, `eos4DDilute` and the cell index   | 5.80 s (+1.64) |

The 73 MB table is not the bottleneck: the loads mostly hit the cache,
because the fluid samples a narrow region of the table. The cost is the
correctly rounded arithmetic of the setup:
- two sqrt and three divisions by T~^2 in `eos4DTilde`;
- four divisions by the grid spacings in `eos4DSetup`.

Changing their rounding is not an option either. (Correction: an earlier
version of this section said that relaxing only the EoS broke
conservation at 3.9e-4. In section 7 of gpu_benchmark_and_next_steps.md,
3.9e-4 is the run with the EoS kept exact and everything else relaxed. A
run with only the EoS relaxed was never made.)

Section 8 shows that the divisions can be made cheaper exactly.

### B: the W face values

All three variants are bit-identical to the reference (Metal outputs
equal to bench18). None is faster:

GPU seconds of face + wface + cellVisc rk0 + rk1:

| variant | visc9_161 | visc20_161 | diff14_161 |
|---|---|---|---|
| before | 6.16 | 9.77 | 9.54 |
| W face values computed in the cell kernel (no wface kernel, no face buffers) | 6.25 | 9.86 | 10.11 |
| wface fused into the face kernel (one thread does both) | 6.15 | 9.88 | 9.63 |

Wall time with face and wface dispatched concurrently
(`MTLDispatchTypeConcurrent`, barriers only between the groups
{strings, pressure}, {faces}, {cell}), against the bench18 baseline:

| case | concurrent | baseline |
|---|---|---|
| visc9_161 | 6.89 s | 6.79 s |
| visc20_161 | 11.43 s | 11.36 s |
| diff14_161 | 10.54 s | 10.27 s |

Why none of them helps:
- **Concurrent dispatch:** each dispatch already fills the GPU, so
  nothing is gained by overlapping them.
- **Computing in the cell kernel:** it saves the 76 B per face written
  and read back. But the cell then reads the η neighbours (±1 and ±2
  planes, 2.4 MB per plane) itself, and those are not in cache either.
  The added registers also cost occupancy in an already large kernel.
- **Fusion:** it saves only the u^0 and u^d loads that the two kernels
  share.

wface moves about 168 B per face. At the M4's bandwidth that comes to
about 1.1 s, which is the 1.0 s it takes. It is at the hardware limit.

The in-cell variant is not kept in the repository.

### What is left

- **D (`face_reconstruction 1`):** a physics choice; -20 to -28% of GPU
  time.
- **C (gradients):** the next exact candidate. It has the same
  caching caveat as B, so a tiled kernel with threadgroup memory is the
  form that could pay off.
- **E (freeze-out detection on the device):** helps the freeze-out
  steps only.

## 8. Done: an exact fast division for the 4D EoS (Metal)

### What the operations cost

Tested on the M4 in `MTLMathModeSafe`, with a dependent chain of 256
operations over 4M threads:

| operation | ns per operation | correctly rounded? |
|---|---|---|
| `a * b + c` | 2.8 | yes |
| `a / b` (safe) | 13.5 | yes: no error in 4.3e9 random pairs |
| `sqrt` (safe) | 3.0 | **no**: faithful only, 28% of results differ from the correctly rounded value |
| reciprocal from `fast::divide(1, b)`, one Newton step, Markstein's correction `fma(fma(-b, q, a), y, q)` | 3.1 | equal to safe `/` |

- **The square roots cost nothing to speak of.** They are as cheap as a
  multiply and are not worth rewriting. They also round differently from
  the CPU, which is one more reason Metal float cannot match backend 1
  float bit for bit.
- **The cost is in the divisions:** 4.4 times the short sequence.

### Exactness of the short sequence (equal bits to safe `/`)

- **Exhaustive:** all 2^46 pairs of float mantissas, in 1232 s on the GPU,
  with no difference.
- **Map over all 256 x 256 exponent fields,** with 4096 random mantissas
  and signs each. It found three cases that need handling:
  - a numerator with exponent field 1-25: the remainder underflows. Fixed
    by scaling it by 2^64.
  - a = -0: the correction gives +0. Fixed by returning the product.
  - inf or NaN operands and b = 0. Handled with selects in `mdDiv` only.
- **Not covered:** |b| ≥ 2^126 finite, and quotients that overflow. These
  are the preconditions.

After the fixes, neither version differs from safe `/` anywhere inside
its precondition:
- **`mdDivF`:** a not ±inf, b normal, |b| < 2^126.
- **`mdDiv`:** additionally any a, and b = 0 or ±inf.

On the CPU and CUDA backends both functions are plain `/`, so backend 1
is unchanged.

### The change

| file | what changes |
|---|---|
| `device_defs.h` | `mdDivF` and `mdDiv` |
| `device_eos.h` | every division of the 4D lookup (tilde variables, dilute scaling, cell setup, derivatives, dilute gradient branch, cs^2) |
| `device_hydro.cpp` | a check at upload: grid spacings, T_ref, T5diff and the range of the clamped T normal and below 2^126, and every table entry finite. It fails with an error otherwise. |

How each call site chooses:
- **`mdDivF`** where the operands are provably in range:
  - table constants;
  - powers of the clamped T;
  - differences of finite table values;
  - the tilde and setup divisions, behind one finiteness test per call.
- **`mdDiv`** at `e / 3`, `/ eLocal`, `/ enthalpy`, and in the dilute
  gradient branch.

**The guards cost time:**

| variant | face, visc20_161 |
|---|---|
| before | 4.15 s |
| a range test that falls back to safe `/` | 4.89 s (the compiler evaluates the slow division too) |
| `mdDiv` everywhere | 3.45 s |
| per site (adopted) | 3.27 s |
| no guards at all (not exact) | 2.75 s |

### Results

**GPU seconds, visc20_161:**

| kernel | before | after |
|---|---|---|
| face | 4.15 | 3.27 |
| cellVisc_rk1 | 2.95 | 2.44 |
| cellVisc_rk0 | 1.65 | 1.68 |

**Wall time:** visc20_161 goes from 11.36 to 10.18 s (best of 2, -10%).
The EoS 9 and 14 cases are unchanged (visc9_161, diff14_161 profiled),
as expected.

**Outputs:** `global_conservation_laws.dat` is bit-identical to bench18
for visc20_161, visc9_161 and diff14_161.

**Tests:**
- unit tests (clang, `-Dunittest=ON`): bugfix 45/45 (1,520,088
  assertions), reconst, cell, grid and minmod all pass;
- integration suite on the Metal build: 25/25. The EoS 20 metrics are at
  their safe-mode values, e.g. μ_B 1.87e-6 against 1.9e-6 before.

### Possible next step

The same functions could replace the solver's divisions (48 lines in
device_ideal.h, 176 in device_visc.h), which run for every EoS.

Each site would need the same care. On overflow, `mdDiv` gives NaN where
`/` gives inf, and some solver tests (`fabs(u0 - u0_guess) / u0_guess >
100`) treat the two differently. So an overflow rule, or a per-site
argument that overflow cannot happen, comes first.

