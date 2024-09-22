#!/usr/bin/env python3

import numpy as np

TargetSum = [
    6.4800000000000011e+01,
    3.7239163071999992e+00,
    1.2122594053000002e+00,
    7.0175462621000015e+00,
    8.3051035398999975e+00,
    2.0319850283000003e+00,
    1.1591788223000000e+01,
    1.4463400459999998e+00,
    4.5466858733999986e+00,
    2.0260406959999999e+00,
    4.6506243910999991e+00,
    2.3769707982000003e+01,
    2.5438350987999989e+00]

filename = "eccentricities_evo_eta_-0.5_0.5.dat"
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
    6.4800000000000011e+01,
    1.2780592846999998e+02,
    1.5114159886800005e+01,
    8.6003927700000020e+01,
    2.3591794753000002e+01]

filename = "inverse_Reynolds_number_eta_-0.5_0.5.dat"
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
    6.4800000000000011e+01,
    2.3607672464999996e+00,
    5.1168141943999979e-01,
    4.2385794246210001e-01,
    1.1816455545524997e-01,
    3.9122687502010006e-01,
    1.1820525990199998e-01,
    8.8910506421999997e+00,
    3.9057346994000017e+00,
    8.1077981776000012e+00,
    3.8510958245000007e+00,
    8.1528788654999982e+00,
    3.8855050015999995e+00,
    1.5626989439171000e+00,
    1.1970320337170000e+00,
    1.4207528854690008e+00,
    1.2122066540400001e+00,
    1.4286712099673002e+00,
    1.2965062316940001e+00]

filename = "momentum_anisotropy_eta_-0.5_0.5.dat"
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
    6.4800000000000011e+01,
    2.0621871678999992e+05,
    5.0993951872999998e+04,
    4.1870226609000010e+03,
    8.0181408006000004e+02]

filename = "meanpT_estimators_eta_-0.5_0.5.dat"
data = np.loadtxt(filename)
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    diff = abs(val_i - ResSum[i])/val_i
    if diff > 1e-8:
        print("{0}: Diff: {1:.16e} {2:.16e} {3:.6e}".format(filename, val_i,
                                                            ResSum[i], diff))
        Nfailed += 1

exit(Nfailed)
