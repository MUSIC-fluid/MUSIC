# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

MUSIC is a 3+1D relativistic second-order viscous hydrodynamics code for heavy-ion collisions (C++11, OpenMP), with optional GPU backends (Metal, CUDA).

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release    # options: -Dmetal=ON (macOS), -Dcuda=ON (needs nvcc), -Dopenmp=OFF, -Dunittest=ON
make -j10 && make install              # installs MUSIChydro (and libmusic) into the repo root
```

* With Apple clang, OpenMP needs Homebrew `libomp` (or `LIBOMP_ROOT`).
* `-Dunittest=ON` builds the unit-test executables and `integration_hydro_memory.e`, **instead of** `MUSIChydro`. Use a separate build directory for it.
* Formatting: `bash formatCode.sh`, which runs clang-format with the repo's `.clang-format`, 80 columns.

## Run

```bash
OMP_NUM_THREADS=4 ./MUSIChydro music_input   # input: "name value" lines ending with EndOfData
```

* MUSIC reads the EoS tables from `./EOS` relative to the working directory.
* At start-up it deletes old `*.dat` outputs in that directory.
* Every run writes `music_parameters.dat`. Its `[input]` part, with the defaults filled in, is an input file for the same run. `[unused]` lists the input entries that MUSIC does not read, and `[run]` lists every `InitData` field. `ReadInParameters::writeParameterFile` spells the fields out, and `unittest_bugfix` checks that the list is complete.
* The parameters and their defaults are in `src/read_in_parameters.cpp`, with descriptions in `utilities/music_parameters_dict.py`. Examples are in `example_inputfiles/`.

## Tests

* **Unit tests** (doctest, in a `-Dunittest=ON` build): `./src/unittest_{bugfix,reconst,cell,minmod,grid}.e`.
  * Run `unittest_bugfix.e` from the repo root: it needs `EOS/`.
* **Integration tests** (the main regression suite). They generate their own small initial states and check invariants (conservation, agreement of code paths), not golden numbers:
  ```bash
  python3 tests/integration/music_integration_tests.py --music ./MUSIChydro [--workdir DIR] [--threads N] [test ...]
  python3 tests/integration/music_integration_tests.py --list
  ```
  * Run a single test by naming it, e.g. `... --music ./MUSIChydro conservation_viscous`.
  * The `gpu_backend*` tests compare the device backends against the CPU code. Backends missing from the build are skipped.
  * `memory_history` needs `--memory-driver build_unittest/src/integration_hydro_memory.e`.
* **Other tests used in CI** (`.github/workflows/BuildTest.yml`): the `RunAutoTest.sh` scripts in `tests/Gubser_flow`, `tests/testIPGlasma2D` and `tests/test_1+1D_with_Akihiko`.
* `tests/` is listed in `.gitignore`, but its files are tracked. Stage changes there with `git add -u` (or `-f` for new files).

## Architecture

### Flow

`main.cpp` → `MUSIC` (`music.cpp`):
1. read the parameters into `InitData DATA` (`data.h`, `read_in_parameters.cpp`);
2. build the EoS (`eos.cpp` dispatches to the `eos_*.cpp` tables) and the hydro source;
3. `Init::InitArena` sets up the initial state (`init.cpp`, selected by `Initial_profile`);
4. `Evolve::EvolveIt` runs the time loop (`evolve.cpp`).

If the fluid reaches the grid edge, `music.cpp` can re-run on a larger grid.

### Fields

`Fields` (`fields.h`) holds one time slice as SoA arrays, with x fastest and 2 ghost cells per side (`getFieldIdx`, `getFieldIdxHalo`). Evolve keeps three slices, `arenaFieldsPrev`/`Curr`/`Next`, and rotates them.

### One step

Each step is a two-stage RK (Heun) step; each stage is `Advance::AdvanceIt` (`advance.cpp`). Per stage:
1. Pressure array.
2. **Face fluxes:** KT central scheme with minmod reconstruction, each face computed once from its 4-cell stencil (`computeFaceFluxes`).
3. **W divergences:** `Diss::computeWDivergences` (`dissipative.cpp`).
4. **Per-cell update:**
   * `FirstRKStepT`: conservative update of τT^{τμ}, J^τ;
   * `Reconst`: Newton root solve of e, n, u from T^{τμ} (`reconst.cpp`);
   * `U_derivative`: gradients, θ, σ, ω (`u_derivative.cpp`);
   * `FirstRKStepW`: relaxation equations for W^{μν}, Π, diffusion q^μ (`dissipative.cpp`, `transport_coeffs.cpp`).
5. In the second stage, `stage2ViscousFixedPoint` makes the viscous correction consistent.

### Diagnostics and output

* `Cell_info::check_conservation_law` (`grid_info.cpp`) writes `global_conservation_laws.dat`. It includes the time-integrated outflow through the boundary faces, which it evaluates itself from the edge cells. Most integration tests check this file.
* The freeze-out surface comes from Cornelius (`freeze.cpp`, `cornelius.cpp`). It is binary, 36 floats per cell without vorticity, one `surface_eps_*_<thread>.dat` file per thread.

### beastMode

* `beastMode` ≥ 2 changes Δτ during the run.
* `beastMode` 3 also enlarges and coarsens the transverse grid (`evolve.cpp`).
* It is incompatible with evolution output and in-memory history, and falls back to 1 when either is on.

### Vacuum skipping

`vacuumSkipEnergy` (`doc/vacuum_skip.md`):
* Interior faces whose stencil is below the threshold get zero flux.
* Frozen cells skip the viscous update.
* The default depends on beastMode: 1e-8 for 0, 1e-6 for 1, 1e-4 for ≥ 2.

### GPU backends (`src/device/`)

The stage kernels are templated on `Real` (float or double) and shared by every backend:
* `device_ideal.h`: KT faces, the cell update and the root solve;
* `device_visc.h`: the viscous cell kernel, W faces and μ_B/T;
* `device_eos.h`, `device_strings.h`, `device_diag.h`.

`DeviceHydro` (`device_hydro.cpp`) drives them from `Evolve`. It fills `StepParams` for each stage and keeps the device and host `Fields` in sync lazily. The backends implement the C interface in `device_api.h`:

| `gpu_backend` | runs on | precision | where |
|---|---|---|---|
| 1 | the CPU (OpenMP) | double | inside `device_hydro.cpp` |
| 2 | Metal | float | `metal_backend.m`; `metal_kernels.metal` and the headers are bundled into a source string compiled at run time |
| 3 | CUDA | double, or float with `gpu_float 1` | `cuda_backend.cu` |

* **The kernels mirror the host functions statement by statement.** Backend 1 reproduces the CPU code bit for bit with clang, and the tests rely on that. Keep the device and host versions of a function in step, with the same statement grouping, because clang contracts multiply-adds within a statement.
* **`StepParams` has a fixed layout** shared across compilers. Some byte offsets are hard-coded in `device_api.h` and checked by `static_assert` in `device_hydro.cpp`. Add new fields at the end.

## Numerical conventions

* Internal units are fm: e is in 1/fm⁴. Convert GeV/fm³ with `Util::hbarc`.
* The NEOS and hotQCD tables store T⁵ on purpose, for interpolation accuracy and the low-e extrapolation (`doc/eos_temperature_lookup.md`).
* Net baryon number must stay non-negative. A negative N_B is a bug, not something to clamp.
* Do not clamp the bulk pressure Π to avoid cavitation. Fix the bulk numerics instead.
* Performance changes are expected to be bit-identical where possible. Verify that `global_conservation_laws.dat` and the freeze-out surface are unchanged.

## Docs

`doc/` holds the design and benchmark reports: GPU porting phases, OpenMP scaling, conservation treatments, vacuum skip. New reports go there as Markdown.
