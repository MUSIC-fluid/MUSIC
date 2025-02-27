#!/usr/bin/env python3

import numpy as np

TargetSum = [
    1.2550000000000000e+03,
    9.8862329999999162e+01,
    4.8999999998387217e+01,
    2.5050000000000000e+02]

data = np.loadtxt("1+1D_DiffusionTest_tau_4.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("Diff: ", val_i, ResSum[i])
        Nfailed += 1
