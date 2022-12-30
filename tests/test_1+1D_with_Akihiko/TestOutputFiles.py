#!/usr/bin/env python3

import numpy as np

TargetSum = [
    2.4151200000000008e+03,
    5.4440724726112649e+03,
    8.1379650394072002e+01]

data = np.loadtxt("1+1D_check_tau_2.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("eccentricities_evo: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    2.4151200000000008e+03,
    1.7157553286169775e+03,
    3.2297206619534670e+01]

data = np.loadtxt("1+1D_check_tau_5.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("inverseReynoldsNumbers: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    2.4151200000000008e+03,
    7.4088124750189286e+02,
    1.6147310790469277e+01]

data = np.loadtxt("1+1D_check_tau_10.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("inverseReynoldsNumbers: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    2.4151200000000008e+03,
    3.3092286816402873e+02,
    8.0815561201914647e+00]

data = np.loadtxt("1+1D_check_tau_20.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("inverseReynoldsNumbers: Diff: ", val_i, ResSum[i])
        Nfailed += 1

exit(Nfailed)
