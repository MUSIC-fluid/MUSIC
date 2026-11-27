#!/usr/bin/env python3

from pathlib import Path

import numpy as np


fixture_dir = Path(__file__).resolve().parent
failed = False

for tau in (2, 5, 10, 20):
    filename = f"1+1D_check_tau_{tau}.dat"
    output = np.loadtxt(filename)
    if output.shape != (695, 5) or not np.isfinite(output).all():
        raise SystemExit(f"{filename} has an unexpected shape or nonfinite values")
    rapidity = output[:, 0]
    np.testing.assert_allclose(rapidity, np.linspace(-6.94, 6.94, 695),
                               rtol=0, atol=1e-8)

    for name, column, limit in (("e", 1, 3e-3), ("rhoB", 2, 5e-3)):
        ref_eta, ref_values = np.loadtxt(
            fixture_dir / f"{name}_baryon_t{tau}.dat", unpack=True)
        # The benchmark profiles use a coarser rapidity grid than MUSIC.
        values = np.interp(ref_eta, rapidity, output[:, column])
        error = np.sum(np.abs(values - ref_values)) / np.sum(np.abs(ref_values))
        print(f"tau={tau:2d} {name:>4} relative L1 error {error:.6e} (limit {limit:.3e})")
        failed |= error > limit

raise SystemExit(int(failed))
