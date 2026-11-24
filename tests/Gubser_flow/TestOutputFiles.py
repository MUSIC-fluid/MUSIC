#!/usr/bin/env python3

import numpy as np

TargetSum = [
    1.0150499999999958e+05,
    1.0150500000000000e+05,
    1.1693538506808669e+04,
    0.0000000000000000e+00,
    3.4153034770665554e+03,
    1.7548544785397306e+04,
    1.7548544785397280e+04,
    1.5690231040563344e+03,
    1.5690231040563135e+03,
    2.9951421246331080e+02,
    2.1952922587096359e+03,
]

filename = "Gubser_flow_check_tau_1.2.dat"
data = np.loadtxt(filename)
ResSum = np.sum(abs(data), axis=0)

Nfailed = 0
for i, val_i in enumerate(TargetSum):
    diff = abs(val_i - ResSum[i])/(val_i + 1e-16)
    if diff > 1e-8:
        print("{0}: Diff: {1:.16e} {2:.16e} {3:.6e}".format(filename, val_i,
                                                            ResSum[i], diff))
        Nfailed += 1

TargetSum = [
    1.0150499999999958e+05,
    1.0150500000000000e+05,
    6.8014142654723682e+03,
    0.0000000000000000e+00,
    3.0881322466310285e+03,
    2.2541635202073376e+04,
    2.2541635202073558e+04,
    8.8381017702819122e+02,
    8.8381017702820031e+02,
    2.1796012743174285e+02,
    1.0843355809645861e+03,
]

filename = "Gubser_flow_check_tau_1.5.dat"
data = np.loadtxt(filename)
ResSum = np.sum(abs(data), axis=0)

for i, val_i in enumerate(TargetSum):
    diff = abs(val_i - ResSum[i])/(val_i + 1e-16)
    if diff > 1e-8:
        print("{0}: Diff: {1:.16e} {2:.16e} {3:.6e}".format(filename, val_i,
                                                            ResSum[i], diff))
        Nfailed += 1

TargetSum = [
    1.0150499999999958e+05,
    1.0150500000000000e+05,
    3.2252592379857783e+03,
    0.0000000000000000e+00,
    2.7199245764543934e+03,
    3.1463668916200142e+04,
    3.1463668916200233e+04,
    4.4232754190295941e+02,
    4.4232754190296117e+02,
    1.4589723647530923e+02,
    4.3460915164337143e+02,
]

filename = "Gubser_flow_check_tau_2.dat"
data = np.loadtxt(filename)
ResSum = np.sum(abs(data), axis=0)

for i, val_i in enumerate(TargetSum):
    diff = abs(val_i - ResSum[i])/(val_i + 1e-16)
    if diff > 1e-8:
        print("{0}: Diff: {1:.16e} {2:.16e} {3:.6e}".format(filename, val_i,
                                                            ResSum[i], diff))
        Nfailed += 1

exit(Nfailed)
