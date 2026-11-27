# Build Test unit-test failures · November 27, 2026

The `gpu` branch built on GitHub Actions, but the Unit Tests step failed. The runner used Ubuntu 26.04 and GCC 15.2.0. The failure occurred at revision `33d2bdc`, and commit `524c127` contains the fix. [Failed Build Test run](https://github.com/MUSIC-fluid/MUSIC/actions/runs/45890109288)

## Findings

| Failure | Evidence | Change |
| --- | --- | --- |
| Repository fixtures were missing at runtime | `unittest_bugfix.e` could not open `src/data.h` or `tests/test_source_terms/strings_event_0.dat` from `build_unittest/`. | Run the unit-test executables from the repository root. |
| EoS 9 terminated the test process | The checkout lacked `EOS/hotQCD/hrg_hotqcd_eos_binary.dat`, which the outflow and other device tests require. | Download the table with `EOS/download_hotQCD.sh` before running the suite. |
| Exact floating-point comparisons failed | GCC produced small host-to-device differences, up to `2.32e-9` in the logged viscous comparison. The per-cell and vectorized viscous source comparisons also differed. These results are consistent with different multiply-add contraction. | Set `-ffp-contract=off` for GCC unit-test builds to preserve the exact assertions. |

## Reproduce the verification

Run from a clean checkout at the fix commit with GCC, CMake, and GSL installed.

```bash
bash EOS/download_hotQCD.sh
cmake -S . -B build_unittest -Dunittest=ON
cmake --build build_unittest --parallel 2
for test in bugfix reconst cell minmod grid; do
    ./build_unittest/src/unittest_${test}.e || exit $?
done
```

On x86-64 hardware that supports the v3 instruction set, use `CXX='g++ -march=x86-64-v3'` on the configure command and choose a separate build directory to test a target that supports fused multiply-add instructions.

All 67 test cases passed locally with GCC 15.2.0 in both the default and x86-64-v3 builds. The clean test copy included the downloaded hotQCD table. The optional NEoS tables, GPU hardware, and later workflow steps were outside this verification. No GitHub run of the fix was verified as part of this investigation.
