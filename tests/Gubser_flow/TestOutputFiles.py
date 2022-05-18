#!/usr/bin/env python3

import numpy as np

TargetSum = [
    1.0150499999999958e+05,
    1.0150500000000000e+05,
    1.1695506944236671e+04,
    0.0000000000000000e+00,
    3.4154486647410031e+03,
    1.7541942582439213e+04,
    1.7541942582439118e+04,
    1.5661044876664125e+03,
    1.5661044876664262e+03,
    2.9834954427416056e+02,
    2.1920881825318106e+03,]

data = np.loadtxt("Gubser_flow_check_tau_1.2.dat")
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("Tau 1.2: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    1.0150499999999958e+05,
    1.0150500000000000e+05,
    6.8042139103484769e+03,
    0.0000000000000000e+00,
    3.0883110192107474e+03,
    2.2538664434060753e+04,
    2.2538664434060909e+04,
    8.8036640789070520e+02,
    8.8036640789070452e+02,
    2.1643887395530089e+02,
    1.0811562366304145e+03,]

data = np.loadtxt("Gubser_flow_check_tau_1.5.dat")
ResSum = np.sum(abs(data), axis=0)

for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("Tau 1.5: Diff: ", val_i, ResSum[i])
        Nfailed += 1

TargetSum = [
    1.0150499999999958e+05,
    1.0150500000000000e+05,
    3.2269746223140432e+03,
    0.0000000000000000e+00,
    2.7196343648173633e+03,
    3.1513260848915426e+04,
    3.1513260848915357e+04,
    4.3990214732176980e+02,
    4.3990214732177088e+02,
    1.4504254181735934e+02,
    4.3268368451460509e+02,]

data = np.loadtxt("Gubser_flow_check_tau_2.dat")
ResSum = np.sum(abs(data), axis=0)

for i, val_i in enumerate(TargetSum):
    if abs(val_i - ResSum[i]) > 1e-8:
        print("Tau 2: Diff: ", val_i, ResSum[i])
        Nfailed += 1

exit(Nfailed)
