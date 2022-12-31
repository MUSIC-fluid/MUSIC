#!/usr/bin/env python3

import numpy as np

TargetSum = [
    6.4800000000000011e+01,
    3.7235944538000001e+00,
    1.2120676067999996e+00,
    7.0175700110000001e+00,
    8.3052681264000032e+00,
    2.0320722766999992e+00,
    1.1591267063000000e+01,
    1.4458890885000000e+00,
    4.5460845083999990e+00,
    2.0271981437000002e+00,
    4.6516473991000016e+00,
    2.3768197602999987e+01,
    2.5437955508000005e+00]

data = np.loadtxt("eccentricities_evo_eta_-0.5_0.5.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i])/ResSum[i] > 1e-8:
        print("eccentricities_evo: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    6.4800000000000011e+01,
    1.2780155474000003e+02,
    1.5089618294400003e+01,
    8.6004110200000056e+01,
    2.3591812535000003e+01]

data = np.loadtxt("inverse_Reynolds_number_eta_-0.5_0.5.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i])/ResSum[i] > 1e-8:
        print("inverseReynoldsNumbers: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    6.4800000000000011e+01,
    2.3607696820000004e+00,
    5.1171075872000005e-01,
    4.2389967757479990e-01,
    1.1818478197884000e-01,
    3.9129084032190004e-01,
    1.1822507963059999e-01,
    8.8909357168999978e+00,
    3.9056820904000014e+00,
    8.1077076357000006e+00,
    3.8510225726999985e+00,
    8.1529401965999959e+00,
    3.8854500403999990e+00,
    1.5627119453058989e+00,
    1.1970754529920002e+00,
    1.4208116339400005e+00,
    1.2121877028281998e+00,
    1.4287486738919000e+00,
    1.2965363997980002e+00]

data = np.loadtxt("momentum_anisotropy_eta_-0.5_0.5.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i])/ResSum[i] > 1e-8:
        print("momentum_anisotropy_evo: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    6.4800000000000011e+01,
    2.0622189934999993e+05,
    5.0994775190999993e+04,
    4.1870379982999993e+03,
    8.0181895288000010e+02]

data = np.loadtxt("meanpT_estimators_eta_-0.5_0.5.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i])/ResSum[i] > 1e-8:
        print("meanpTEstimators_evo: Diff: ", val_i, ResSum[i])
        Nfailed += 1

exit(Nfailed)
