# Agent guide for MUSIC

MUSIC simulates relativistic heavy-ion collisions with 3+1D second-order viscous hydrodynamics. The core is C++11 with OpenMP and optional Metal and CUDA backends. Follow the existing code style in `.clang-format`. Use the [README](README.md) for input examples and publication citations, and look in [`doc/`](doc/) for previous investigations.

## Physics and numerics

- For changes to the physics, identify the equation, approximation, or reference that motivates the change. State the units and conventions at the interface. In Milne coordinates, check powers of proper time `tau`, especially for rapidity components such as `tau u^eta` and conserved quantities such as `tau T^{tau mu}`. Use `Util::hbarc` for conversions between internal units and GeV-based quantities.
- For changes to evolution, source terms, or regridding, check energy, momentum, and each evolved charge as applicable. Include integrated sources and boundary outflow in the balance. Use `global_conservation_laws.dat` to inspect the residuals. Check relevant limiting cases, such as boost invariance, the ideal-fluid limit, or dilute cells. Investigate unphysical values and document any change to a regulator, floor, or default.
- For discretization changes, use Bjorken or Gubser flow when applicable and check convergence under timestep or grid refinement. For EoS changes, check table seams, interpolation and gradients, the low-density limit, and behavior near the table boundary. For freeze-out changes, check the Cornelius surface and the flux balance across it.
- When touching a shared equation, check its implementation in `src/device/` alongside the host code. Use `gpu_backend 0` as the host reference and `gpu_backend 1` to exercise the device kernels on the CPU. Confirm the runtime log when testing Metal or CUDA so a CPU fallback does not count as GPU coverage. Keep existing exact comparisons where the implementation is intended to agree bit for bit. For comparisons across precisions or hardware, justify tolerances with units, an error scale, and the relevant observables.
- When changing inputs or defaults, check `src/data.h`, `src/read_in_parameters.cpp`, and `utilities/music_parameters_dict.py`. Keep `ReadInParameters::writeParameterFile` in sync with `InitData`. Preserve the replayable `[input]` section of `music_parameters.dat`.

## Build and validate

Run these commands from the repository root. Use separate build directories because `-Dunittest=ON` builds the doctest executables and the memory-history driver in place of `MUSIChydro`. CI installs GSL. GPU builds use `-Dmetal=ON` on macOS or `-Dcuda=ON` with the CUDA toolkit.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2

cmake -S . -B build_unittest -Dunittest=ON
cmake --build build_unittest --parallel 2
# Run from the repository root so the tests can find their fixtures.
test -s EOS/hotQCD/hrg_hotqcd_eos_binary.dat || bash EOS/download_hotQCD.sh
for test in bugfix reconst cell minmod grid; do
    ./build_unittest/src/unittest_${test}.e || exit $?
done
```

The Python integration suite requires NumPy. Use `--list` to find cases relevant to the change. Omit the case names to run the full suite. The memory-history test uses the driver from the unit-test build.

```bash
python3 tests/integration/music_integration_tests.py --list
python3 tests/integration/music_integration_tests.py \
    --music build/src/MUSIChydro --threads 4 \
    --memory-driver build_unittest/src/integration_hydro_memory.e \
    conservation_ideal bjorken_rk_order
```

Use a dedicated run directory for standalone simulations so cleanup of old outputs cannot affect other runs. Set `HYDROPROGRAMPATH` to the absolute directory containing `EOS/` when using tables from a different working directory. Use a fresh directory with `--workdir` to retain integration artifacts. Report skipped cases and missing tables or hardware with the test results.

## Investigation reports

Save a Markdown summary for each investigation under `doc/`. Use the existing topic-based filenames, adding a date when it distinguishes repeated investigations. Include the question, relevant revisions and inputs, reproduction commands, measured results, conclusion, and limits of verification.

For physics results, capture the initial condition and seed, system and collision energy, EoS and table provenance, transport and freeze-out settings, grid and timestep, and the relevant sections of `music_parameters.dat`. For performance results, also record the baseline, compiler and build settings, hardware, backend and precision, thread count, timing scope, and variability across repeats. Compare numerical agreement on the same workload before drawing performance conclusions. Link to source material and retained artifacts where available.

## Commits

Write a short, imperative subject and a brief body summarizing what changed and why. Record relevant validation when available. End the commit message with `Drafted-by: <agent name>` using the name of the AI agent that drafted it, for example `Drafted-by: Codex`.
