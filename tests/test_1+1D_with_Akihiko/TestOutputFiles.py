#!/usr/bin/env python3

import numpy as np

TargetSum = [
    2.4151200000000008e+03,
    5.4440723046865851e+03,
    8.1389779803540847e+01,
]

filename = "1+1D_check_tau_2.dat"
data = np.loadtxt(filename)
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    diff = abs(val_i - ResSum[i])/val_i
    if diff > 1e-8:
        print("{0}: Diff: {1:.16e} {2:.16e} {3:.6e}".format(filename, val_i,
                                                            ResSum[i], diff))
        Nfailed += 1

TargetSum = [
    2.4151200000000008e+03,
    1.7157548692834894e+03,
    3.2306419625615632e+01,
]

filename = "1+1D_check_tau_5.dat"
data = np.loadtxt(filename)
ResSum = np.sum(abs(data), axis=0)

for i, val_i in enumerate(TargetSum):
    diff = abs(val_i - ResSum[i])/val_i
    if diff > 1e-8:
        print("{0}: Diff: {1:.16e} {2:.16e} {3:.6e}".format(filename, val_i,
                                                            ResSum[i], diff))
        Nfailed += 1

TargetSum = [
    2.4151200000000008e+03,
    7.4088079870732668e+02,
    1.6153535542164590e+01,
]

filename = "1+1D_check_tau_10.dat"
data = np.loadtxt(filename)
ResSum = np.sum(abs(data), axis=0)

for i, val_i in enumerate(TargetSum):
    diff = abs(val_i - ResSum[i])/val_i
    if diff > 1e-8:
        print("{0}: Diff: {1:.16e} {2:.16e} {3:.6e}".format(filename, val_i,
                                                            ResSum[i], diff))
        Nfailed += 1

TargetSum = [
    2.4151200000000008e+03,
    3.3092248507449591e+02,
    8.0851687718083127e+00,
]

filename = "1+1D_check_tau_20.dat"
data = np.loadtxt(filename)
ResSum = np.sum(abs(data), axis=0)

for i, val_i in enumerate(TargetSum):
    diff = abs(val_i - ResSum[i])/val_i
    if diff > 1e-8:
        print("{0}: Diff: {1:.16e} {2:.16e} {3:.6e}".format(filename, val_i,
                                                            ResSum[i], diff))
        Nfailed += 1

exit(Nfailed)
