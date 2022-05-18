#!/usr/bin/env python3

import numpy as np

TargetSum = [
    6.4800000000000011e+01,
    3.7235941924999998e+00,
    1.2120666356999996e+00,
    7.0175693517999997e+00,
    8.3052687386999988e+00,
    2.0320736821999996e+00,
    1.1591268086999996e+01,
    1.4458899007000001e+00,
    4.5460872620999995e+00,
    2.0271943015999998e+00,
    4.6516470628000004e+00,
    2.3768201454000000e+01,
    2.5437910798000001e+00]

data = np.loadtxt("eccentricities_evo_eta_-0.5_0.5.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("eccentricities_evo: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    6.4800000000000011e+01,
    1.2780155470000000e+02,
    1.5089618188199999e+01,
    8.6004109960000065e+01,
    2.3591812533999999e+01]

data = np.loadtxt("inverse_Reynolds_number_eta_-0.5_0.5.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("inverseReynoldsNumbers: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    6.4800000000000011e+01,
    2.3607697218000010e+00,
    5.1171218192000012e-01,
    4.2389946895520003e-01,
    1.1818463460729002e-01,
    3.9129063052329999e-01,
    1.1822492958919997e-01,
    8.8909360533000008e+00,
    3.9056845554000001e+00,
    8.1077071061000012e+00,
    3.8510229106999989e+00,
    8.1529396719999987e+00,
    3.8854503512999989e+00,
    1.5627131545955000e+00,
    1.1970757158339995e+00,
    1.4208120655639993e+00,
    1.2121881951334001e+00,
    1.4287490303099994e+00,
    1.2965368941610003e+00]

data = np.loadtxt("momentum_anisotropy_eta_-0.5_0.5.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("momentum_anisotropy_evo: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    6.4800000000000011e+01,
    2.0622189801000006e+05,
    5.0994774724999988e+04,
    4.1870380055999985e+03,
    8.0181894598000008e+02]

data = np.loadtxt("meanpT_estimators_eta_-0.5_0.5.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("meanpTEstimators_evo: Diff: ", val_i, ResSum[i])
        Nfailed += 1

exit(Nfailed)
