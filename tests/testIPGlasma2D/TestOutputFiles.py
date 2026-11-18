#!/usr/bin/env python3

import numpy as np

TargetSum = [
    6.4800000000000011e+01,
    3.7240086457000001e+00,
    1.2124483670999999e+00,
    7.0176993282000026e+00,
    8.3050525250000007e+00,
    2.0317952437000000e+00,
    1.1591705244999998e+01,
    1.4462400443000001e+00,
    4.5461565184000010e+00,
    2.0265420789000008e+00,
    4.6503499002999993e+00,
    2.3769470552999991e+01,
    2.5443717095999996e+00]

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
    2.4273261630000004e+01,
    2.9141047523000014e+00,
    8.6003949059999982e+01,
    2.3591791416000003e+01]

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
    2.3607848645999998e+00,
    5.1182748999999994e-01,
    4.2386102092070005e-01,
    1.1826998947229002e-01,
    3.9122369147449998e-01,
    1.1826037408460002e-01,
    8.8910668950000016e+00,
    3.9059291730000010e+00,
    8.1077870911000005e+00,
    3.8512778688000018e+00,
    8.1528495306000046e+00,
    3.8855906250999985e+00,
    1.5627874880778005e+00,
    1.1972090483730005e+00,
    1.4208524607749995e+00,
    1.2123835146927999e+00,
    1.4287170725430001e+00,
    1.2966120273800004e+00]

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
    2.0621890782999998e+05,
    5.0993992677000009e+04,
    4.1870207172999990e+03,
    8.0181540491000010e+02]

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
