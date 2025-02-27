#!/usr/bin/env python3

import numpy as np

TargetSum = [
    1.0050000000000000e+03,
    3.1085225558300012e+03,
    9.7305773718511347e+01,
    4.4992432148737031e+01,
    3.8278812899518677e+01,
    0.0000000000000000e+00]

data = np.loadtxt("1+1D_RiemannTest_tau_4.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("Diff: ", val_i, ResSum[i])
        Nfailed += 1
