# Build Test and the incomplete freeze-out exit status · November 27, 2026

Commit `e60244a` made `MUSIChydro` return status 1 when a run reaches its time limit with freeze-out incomplete. Unit tests and Gubser passed on GitHub Actions. IPGlasma failed after step 80 of 80 at τ = 1 fm/c, and the workflow stopped before running the other tests. [Failed Build Test run](https://github.com/MUSIC-fluid/MUSIC/actions/runs/46074057163)

The IPGlasma fixture is a boost-invariant viscous hydro run with the bundled test configuration and the downloaded IPGlasma event and hotQCD EoS. It starts at τ = 0.6 fm/c and evolves for 0.4 fm/c at Δτ = 0.005 fm/c. Its checker uses the eccentricities, inverse Reynolds numbers, momentum anisotropy, and mean transverse momentum diagnostics. Disabling freeze-out in this fixture lets it reach the same time limit with status 0.

The `vacuum_skip` integration case intentionally compares partial freeze-out surfaces over a fixed evolution window. Its surface comparison runs now require status 1 together with the time-limit message, then apply the existing conservation and surface checks. The `threads` case and other complete freeze-out tests continue to require status 0.

## Reproduce

Run in a disposable checkout with NumPy, GSL, and a release build installed at the repository root.

```bash
export OMP_NUM_THREADS=4
bash tests/testIPGlasma2D/RunAutoTest.sh
python3 tests/integration/music_integration_tests.py --music ./MUSIChydro \
    --threads 4 vacuum_skip threads freezeout_balance
```

## Local verification

With GCC 15.2 and four OpenMP threads, the original IPGlasma input reproduced exit status 1 after step 80. The revised input returned 0 at the same step, and the CI wrapper passed all four diagnostic comparisons. The same executable reproduced the old `vacuum_skip` assertion failure on its first partial surface. The revised exit check passed for the host and CPU device-kernel partial-surface runs. The host and device-kernel conservation histories and integrated surface volumes agreed in the comparison. The `threads`, `freezeout_balance`, `beastmode3`, `rerun`, `rerun_2d`, `rerun_source`, and `gpu_backend_freezeout` cases passed with the existing success expectation on the CPU build. [Fixture provenance and baseline](build_test_physics_references_2026-11-27.md)
