#!/usr/bin/env python3

import numpy as np

TargetSum = [
    6.4800000000000011e+01,
    3.7235926002000004e+00,
    1.2120665831999999e+00,
    7.0175712151999985e+00,
    8.3052689951000023e+00,
    2.0320704737000002e+00,
    1.1591262477999992e+01,
    1.4458869345000003e+00,
    4.5460791726000007e+00,
    2.0272071558999989e+00,
    4.6516515139999992e+00,
    2.3768187376999997e+01,
    2.5437976524999995e+00]

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
    1.2780150372000001e+02,
    1.5089610773700006e+01,
    8.6004096229999988e+01,
    2.3591812535000003e+01]

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
    2.3607688647000002e+00,
    5.1170939747999988e-01,
    4.2389937602059991e-01,
    1.1818315096435003e-01,
    3.9129061093519996e-01,
    1.1822376212300002e-01,
    8.8909517788999999e+00,
    3.9056880081000007e+00,
    8.1077229196999987e+00,
    3.8510277763999992e+00,
    8.1529560292999985e+00,
    3.8854558958999994e+00,
    1.5627113984728003e+00,
    1.1970795773280003e+00,
    1.4208113737410006e+00,
    1.2121928875283998e+00,
    1.4287487690785001e+00,
    1.2965418171539997e+00]

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
    2.0622194771000001e+05,
    5.0994776798000014e+04,
    4.1870377936000004e+03,
    8.0181897941999978e+02]

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
