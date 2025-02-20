#!/usr/bin/env python3

import numpy as np

TargetSum = [
    1.0050000000000000e+03,
    3.0711980046800059e+03,
    9.8458481676949802e+01,
    4.3438488751545300e+01,
    0.0000000000000000e+00,
    0.0000000000000000e+00]

data = np.loadtxt("1+1D_RiemannTest_tau_4.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("Diff: ", val_i, ResSum[i])
        Nfailed += 1
