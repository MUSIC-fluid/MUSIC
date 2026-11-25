# OpenMP scaling of the CPU hydro (gpu_backend 0)

2026-11-25, branch `gpu` at 96167c9, Apple M4 (4 performance + 6
efficiency cores), Apple clang with Homebrew libomp unless noted.

## Test problem

TATB initial state (`Initial_profile 111`, Gaussian thickness functions,
sigma = 1.2 fm), 96 x 96 x 64 cells in a 15 fm x 15 fm x 12 eta box,
freeze-out at 0.3 GeV/fm^3 every 5 steps.

* `visc`: EOS 14 (NEOS-BQS), shear + bulk + second order + baryon
  diffusion.
* `ideal`: EOS 9 (hotQCD), no net baryon density.

Short runs are 38 steps (tau = 1.35 - 1.8 fm); the full event runs 295
steps until every cell has frozen out (tau = 6 fm).

## Scaling

| threads | visc wall (s) | speed-up | ideal wall (s) | speed-up |
|---|---|---|---|---|
| 1 | 91.9 | 1.00 | 36.3 | 1.00 |
| 2 | 49.6 | 1.85 | | |
| 4 | 26.2 | 3.51 | 10.1 | 3.59 |
| 6 | 21.1 | 4.36 | | |
| 8 | 18.5 | 4.97 | | |
| 10 | 16.3 | 5.64 | 6.1 | 5.95 |

### What limits it: the hardware, not OpenMP

* **Throughput ceiling of the machine.** Ten independent single-thread
  `ideal` runs at the same time take 67-68 s each against 36.3 s alone:
  the 10 cores deliver 5.3 single-thread equivalents (the efficiency cores
  run at about 40 % of a performance core). The 10-thread OpenMP run
  (5.95x) is already above this ceiling (one process shares the EoS tables
  and the caches).
* **Four concurrent single-thread runs** take 38.2 s each against 36.3 s
  alone: 5 % is lost to the shared caches and lower clocks with four
  performance cores busy.
* **Profile of the 10-thread runs** (`sample`, short and full events): the
  main thread is inside parallel regions 98.4 % of the time, and every
  thread waits at a barrier or sleeps 1.3 - 2.3 % of the time, evenly
  spread. So there is no load imbalance between the performance and
  efficiency cores: the `schedule(dynamic, 64)` loops balance them. The
  freeze-out finder (0.2 %), the grid maxima (1.5 %) and the anisotropy
  output (0.3 %) are parallel and small.
* **Start-up** (EoS tables, initial state) takes 1.0 s with 10 threads.
* **Runtime settings make no difference:** `KMP_BLOCKTIME=infinite`
  16.9 s, `OMP_PROC_BIND=spread` 16.7 s, default 16.3 s.
* **GCC 15 + libgomp** is slower on one thread (110.5 s) and on ten
  (17.3 s) than Apple clang + libomp (91.9 s, 16.3 s).

## Done: loop order in Advance::AdvanceIt

The cell loop of `Advance::AdvanceIt` ran `ieta, ix, iy` with `iy`
innermost, but the fields are stored with x fastest (`strideX_ = 1`,
`strideY_ = Sx_`). Each thread's chunk of 64 cells then walked along y
with a stride of `Sx_` doubles in each of the field components, and in
the pressure, face-flux and W-divergence arrays. The face-flux and
W-divergence loops already run x innermost. The loop is now
`ieta, iy, ix`. Each cell is independent, so the results are
bit-identical. The grid-maxima loop of `Cell_info::check_conservation_law`
(`grid_info.cpp`, a max reduction, so also bit-identical) is reordered
the same way.

| run (10 threads) | before (s) | after (s) |
|---|---|---|
| visc, 38 steps (two interleaved pairs) | 17.0, 17.2 | 15.6, 15.8 |
| visc, full event, 295 steps | 145.5 | 130.6 |
| ideal, 38 steps | 6.1 | 5.8 |
| visc, 4 threads | 26.2 | 24.8 |
| visc, 1 thread | 91.9 | 93.6 |

8 - 10 % with many threads and nothing on one thread: the strided access
costs memory bandwidth and cache that the threads share.
`global_conservation_laws.dat`, the momentum anisotropy, the
eccentricities and the logged grid maxima are bit-identical to before.

The other loops with y innermost are in the initial conditions (run once),
in optional outputs (`compute_Lmunu`, `compute_angular_momentum`) and in
the freeze-out finder, where reordering would change the order of the
surface cells in the output file. They cost under 1 % and are left as
they are.

## Remaining opportunities

### 1. Dilute cells cost almost as much as dense ones (largest)

The same event on a 144 x 144 x 64 grid with the same spacing (the extra
cells are all near vacuum) takes 33.0 s against 16.3 s: 2.25x the cells
for 2.02x the time. A near-vacuum cell costs about 80 % of a cell in the
fireball: it still goes through the Newton solve of the primitive
variables, the EoS lookups, the velocity derivatives and the dissipative
update. With this setup most of the grid is dilute (the fireball covers
roughly a fifth of the transverse area at the start and less than 60 % of
the eta range). Production grids are chosen with margins, so typically
half of the cells or more are dilute, and the fraction changes during the
event.

Possible approach: keep a list of active tiles (for example 8 x 8 x 1
cells) whose cells or stencil neighbours (two cells) have e above a small
threshold, and update only those. Outside the active tiles the fields
stay at the vacuum floor. This changes the near-vacuum values (they no
longer follow the Bjorken dilution exactly) and has to keep conservation
and the outflow boundary exact, so it needs its own validation. The
`fluidNearGridEdge` band check of beastMode 3 is a similar idea at the
grid edge. Expected gain: up to the dilute fraction of the cell updates,
1.5 - 3x on generous grids. It helps every thread count and the GPU
backends too.

### 2. EoS lookups and the 5th root

In the 10-thread profile, `EOS_base::interpolate2D` and
`get_table_idx` take 15 % of the samples, `pow` 6 % and
`Reconst::solve_v_from_guess` 13 %. Most `pow` calls come from
`EOS_neos::get_temperature` (`T = pow(T5, 0.2)`, `eos_neos.cpp:227`):
the NEOS tables store T^5. A T table (or a fast 5th root) and computing
the table index and weights once for P, T and mu_B of the same cell
would save part of this. Changes the results at the rounding level.
These are single-thread costs that scale with the threads.

### 3. NUMA placement on multi-socket Linux nodes (not measurable here)

`Fields::resizeFields` zero-fills the arrays with one `std::fill`, so the
first touch puts all field pages on the socket of the main thread. On a
two-socket node the threads of the other socket then read remote memory
in every step. Two options:

* a quick test with no code change: `numactl --interleave=all MUSIChydro`
  against a plain run with `OMP_PROC_BIND=close OMP_PLACES=cores`;
* first-touch the arrays in parallel (`#pragma omp parallel for
  schedule(static)` in `resizeFields`) and use a static schedule over
  eta slices in the hot loops, so that each thread mostly reads memory
  that it touched first. On homogeneous x86 cores the dynamic schedule is
  not needed for load balance; `schedule(runtime)` would let `OMP_SCHEDULE`
  choose per machine.

Worth measuring on the cluster before changing code: on a single-socket
machine (or this M4) it makes no difference.

## Also fixed: crash of gpu_backend 1 with EoS 20 in ideal hydro

The `gpu_backend` integration test crashed (SIGSEGV) in its EoS 20 case on
backend 1 (the device kernels on the CPU). This predates this session's
commits: fd7e7fd already crashes. For ideal hydro `DeviceHydro` uploads
only the P and T tables of the 4D EoS (`nTab = 2`), but `eosTMuB`, called
by the grid maxima of the conservation rows, also read the mu_B table,
past the end of the upload. On Metal the read past the buffer returned
garbage that was never used (the rows use T only); on the CPU it
faulted. `eosTMuB` now reads mu_B only when that table is uploaded
(`d.f4NTables > kEos4DMuB`) and returns 0 otherwise. All 24 integration
tests pass (memory_history skipped without its driver) and the unit tests
pass.
