#!/usr/bin/env python3

import numpy as np

TargetSum = [
    6.4800000000000011e+01,
    3.7245849728000002e+00,
    1.2130409679999996e+00,
    7.0179220035000016e+00,
    8.3047643607999984e+00,
    2.0328048919999993e+00,
    1.1592288020999998e+01,
    1.4475523090999995e+00,
    4.5450341006999988e+00,
    2.0257219332000003e+00,
    4.6492648228000002e+00,
    2.3771217929999995e+01,
    2.5431565151999997e+00,
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
    2.4272356982000009e+01,
    2.9142988158000001e+00,
    8.6001665720000048e+01,
    2.3592283427000002e+01,
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
    2.3594097941999999e+00,
    5.1096212779999994e-01,
    4.2347499017043994e-01,
    1.1778391660247001e-01,
    3.9072130123190002e-01,
    1.1786522862680002e-01,
    8.8864331959000022e+00,
    3.9028621899000000e+00,
    8.1034608426999988e+00,
    3.8480816271000009e+00,
    8.1478879409999987e+00,
    3.8822935825999991e+00,
    1.5601502237092006e+00,
    1.1958904909240000e+00,
    1.4181127808074006e+00,
    1.2104574914840001e+00,
    1.4259424392201003e+00,
    1.2941895140050006e+00,
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
    2.0623359909000000e+05,
    5.0996416666000012e+04,
    4.1874376109000013e+03,
    8.0183229563999998e+02,
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
