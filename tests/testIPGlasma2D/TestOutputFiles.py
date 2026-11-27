#!/usr/bin/env python3

import numpy as np


# Column sums for the bundled CI configuration at revision 524c127,
# using the IPGlasma event and hotQCD table fetched by RunAutoTest.sh.
target_sums = {
    "eccentricities_evo": [
        64.8, 3.7224097661, 1.2120667128, 7.0194025639, 8.3067164154,
        2.0335917408, 11.589587042, 1.4446888516, 4.5391931488,
        2.0349921414, 4.6554529709, 23.760717202, 2.540272379,
        25.30320361, 14.415169907,
    ],
    "inverse_Reynolds_number": [64.8, 24.272490425, 2.9238256434,
                               86.00236429, 23.591850805],
    "momentum_anisotropy": [
        64.8, 2.3596597673, 0.51072557337, 0.424266338776, 0.118299924877,
        0.391517555675, 0.118192168314, 8.886845652, 3.9010497677,
        8.1036695669, 3.8462676151, 8.1476275493, 3.8796610236,
        1.558294238, 1.19685668611, 1.4162378645, 1.21075022705,
        1.42441001359, 1.29411029262,
    ],
    "meanpT_estimators": [64.8, 206216.04687, 50991.755522,
                          4187.2623166, 801.80246801],
}

failed = False
for name, target in target_sums.items():
    filename = f"{name}_eta_-0.5_0.5.dat"
    data = np.loadtxt(filename)
    if data.shape != (81, len(target)) or not np.isfinite(data).all():
        raise SystemExit(f"{filename} has an unexpected shape or nonfinite values")

    sums = np.sum(np.abs(data), axis=0)
    for column, (expected, actual) in enumerate(zip(target, sums), start=1):
        diff = abs(actual - expected) / expected
        if diff > 1e-8:
            print(f"{filename}, column {column}: expected {expected:.12g}, "
                  f"got {actual:.12g} (relative difference {diff:.6e})")
            failed = True

raise SystemExit(int(failed))
