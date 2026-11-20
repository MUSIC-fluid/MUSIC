#!/usr/bin/env python3

import numpy as np

TargetSum = [
    6.4800000000000011e+01,
    3.7224664164999988e+00,
    1.2122756189999999e+00,
    7.0192042334000018e+00,
    8.3066049192999980e+00,
    2.0336951783999995e+00,
    1.1589458342999997e+01,
    1.4452515989999997e+00,
    4.5388483933000003e+00,
    2.0348586049000001e+00,
    4.6548842841000004e+00,
    2.3760602913999989e+01,
    2.5409215311000000e+00,
]

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
    2.4272721439000009e+01,
    2.9243180162000018e+00,
    8.6001725590000007e+01,
    2.3591798932999996e+01,
]

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
    2.3596815059999998e+00,
    5.1106705844000000e-01,
    4.2452473930150003e-01,
    1.1836247222334999e-01,
    3.9175964455200007e-01,
    1.1828083527559999e-01,
    8.8872523881999967e+00,
    3.9026074748000004e+00,
    8.1044593526000011e+00,
    3.8473163456000012e+00,
    8.1479879324000013e+00,
    3.8806242355999996e+00,
    1.5590173586259997e+00,
    1.1972167110970000e+00,
    1.4170718184790998e+00,
    1.2113836224241998e+00,
    1.4250398988447996e+00,
    1.2941183530629998e+00,
]

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
    2.0621687322000001e+05,
    5.0990925959999993e+04,
    4.1872466895000007e+03,
    8.0180949146000012e+02,
]

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
