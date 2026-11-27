#!/usr/bin/env python3

from pathlib import Path

import numpy as np


# Relative L1 bounds for the 201 x 201 grid with dtau = 0.005 fm/c.
# The diagonal xy shear component has a smaller norm and larger relative error.
components = ("T", "u^x", "u^y", "pi^xx", "pi^yy", "pi^xy", "pi^etaeta")
limits = np.array([1e-3, 6e-3, 6e-3, 4e-2, 4e-2, 1.2e-1, 1.2e-2])
fixture_dir = Path(__file__).resolve().parent
failed = False

for tau in ("1.2", "1.5", "2"):
    filename = f"Gubser_flow_check_tau_{tau}.dat"
    output = np.loadtxt(filename)
    reference = np.loadtxt(fixture_dir / f"y=x_tau={tau}_SemiAnalytic.dat")
    n = len(reference)
    if output.shape != (n * n, 13) or not np.isfinite(output).all():
        raise SystemExit(f"{filename} has an unexpected shape or nonfinite values")

    diagonal = output[::n + 1]
    np.testing.assert_allclose(diagonal[:, :2], reference[:, :2], rtol=0, atol=1e-8)
    # Skip energy and baryon density columns absent from the reference.
    reference_fields = reference[:, 2:]
    errors = (np.sum(np.abs(diagonal[:, 4:11] - reference_fields), axis=0)
              / np.sum(np.abs(reference_fields), axis=0))
    for name, error, limit in zip(components, errors, limits):
        print(f"tau={tau} {name:>10} relative L1 error {error:.6e} (limit {limit:.3e})")
        failed |= error > limit

raise SystemExit(int(failed))
