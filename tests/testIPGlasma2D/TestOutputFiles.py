#!/usr/bin/env python3

import numpy as np

TargetSum = [
    6.4800000000000011e+01,
    3.7247454603000012e+00,
    1.2131535135999998e+00,
    7.0178815549000015e+00,
    8.3046429019999994e+00,
    2.0327752493999993e+00,
    1.1592388593999999e+01,
    1.4475500940000003e+00,
    4.5452682302999987e+00,
    2.0253628553000005e+00,
    4.6487572964999995e+00,
    2.3771828430000006e+01,
    2.5436323611000002e+00,
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
    2.4272386434000001e+01,
    2.9142964666999993e+00,
    8.6001598059999992e+01,
    2.3592264931999988e+01,
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
    2.3594673760000000e+00,
    5.1107636801000000e-01,
    4.2346686020347996e-01,
    1.1782242248458004e-01,
    3.9071416232980016e-01,
    1.1789995580959997e-01,
    8.8865844449000004e+00,
    3.9030391235999988e+00,
    8.1035027886999984e+00,
    3.8481641586999986e+00,
    8.1479295820999980e+00,
    3.8823736082000009e+00,
    1.5601902816811999e+00,
    1.1958446234299995e+00,
    1.4181293727498994e+00,
    1.2104063096489994e+00,
    1.4259567284557992e+00,
    1.2941383317009996e+00,
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
    2.0623367202999996e+05,
    5.0996332443000007e+04,
    4.1874356487999985e+03,
    8.0183442861000015e+02,
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
