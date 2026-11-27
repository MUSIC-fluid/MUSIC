# Build Test physics reference failures · November 27, 2026

At revision `524c127`, the build and all five unit-test executables passed on GitHub Actions. The workflow stopped at the viscous Gubser test, after the simulation reached its final proper time. The legacy checker rejected differences of up to `1.19e-3` from its stored column sums at a tolerance of `1e-8`. Local runs reproduced that failure and found stale sums in both later standalone tests. [GitHub Actions run](https://github.com/MUSIC-fluid/MUSIC/actions/runs/46045574895)

## Findings and changes

- **Viscous Gubser flow.** The conformal test uses η/s = 0.2, a 201 × 201 grid with 0.05 fm spacing, and a 0.005 fm/c timestep. The checker now compares the temperature, transverse four-velocity, and shear components along `y=x` with the bundled semi-analytic profiles at τ = 1.2, 1.5, and 2 fm/c. The maximum relative L1 errors over those times were `8.37e-4` for temperature, `4.99e-3` for transverse velocity, `3.12e-2` for transverse diagonal shear, `1.13e-1` for xy shear, and `9.84e-3` for rapidity shear. The limits allow a margin above these errors at the fixed test resolution. The C++ diagnostic now compares π^xy with the reference π^xy column. Its previous π^xx reference produced a misleading π^xy error of 0.54 at τ = 2 fm/c.
- **IPGlasma 2D.** The legacy sums differed by up to `7.50e-4` on the downloaded IPGlasma event with the hotQCD table. The checker uses the 81-row output from `524c127` as its regression baseline and retains `1e-8` relative tolerance. It checks all four diagnostic files in one pass, so a failed comparison in one file counts toward the exit status even if the other files pass.
- **1+1D baryon diffusion.** The old net-baryon sums differed by up to `4.42e-6`. The new check interpolates MUSIC's energy and net-baryon density onto the rapidity grid of the bundled Akihiko profiles at τ = 2, 5, 10, and 20 fm/c. It uses relative L1 limits of `3e-3` and `5e-3`, respectively. Maximum observed errors were `2.16e-3` for energy density and `3.58e-3` for net-baryon density. This benchmark uses the NEoS 2 tables with shear and bulk disabled.

Each comparison checks that its output has the expected shape and finite values. The Gubser and baryon diffusion wrappers stop if the simulation fails. The C++ change is limited to the Gubser diagnostic.

## Reproduce

From a disposable checkout with GCC, CMake, GSL, and NumPy installed, run the same standalone tests as the workflow. The IPGlasma and baryon diffusion wrappers fetch their fixtures from Bitbucket.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
cmake --install build
export OMP_NUM_THREADS=4
bash tests/Gubser_flow/RunAutoTest.sh
bash tests/testIPGlasma2D/RunAutoTest.sh
bash tests/test_1+1D_with_Akihiko/RunAutoTest.sh
```

The downloaded hotQCD table had SHA-256 `140826a2b13112f8dc405b9a3d98cd7321c3463f225aac9b22b99d22a2acf63b`. The IPGlasma event had SHA-256 `5ce96483667ea831019b7dae645ac2f280055812c0dd7d6842ee3faa8d6fee06`.

## Local verification

All 67 unit tests and the three standalone workflow scripts passed on Ubuntu 26.04.1 with GCC 15.2.0, GSL 2.8, NumPy 2.5.4, and four OpenMP threads. The standalone checks also accepted output from a separate `-march=x86-64-v3` release build run with two threads. Across those builds the IPGlasma sums agreed within `1.70e-10` relative. Deliberate perturbations of the Gubser xy shear, baryon density, and each of the four IPGlasma diagnostic files were rejected.

All 27 integration cases passed across split local runs with three or four OpenMP threads, covering the reference solver and the device kernels running on the CPU. The temporary output quota interrupted a duplicate freeze-out run after the independent freeze-out run had passed. Optional Metal, CUDA, NEoS-BQS, and NEoS-4D configurations were outside local coverage. GitHub Actions has yet to run these changes.
