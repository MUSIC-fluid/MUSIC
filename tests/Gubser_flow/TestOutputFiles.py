#!/usr/bin/env python3

import numpy as np

TargetSum = [
    1.0150499999999958e+05,
    1.0150500000000000e+05,
    1.1694870233937785e+04,
    0.0000000000000000e+00,
    3.4154097685618144e+03,
    1.7543444293943372e+04,
    1.7543444293943292e+04,
    1.5662404612665234e+03,
    1.5662404612665343e+03,
    2.9851700358226071e+02,
    2.1921051769456840e+03,
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
    6.8035513343971652e+03,
    0.0000000000000000e+00,
    3.0884305781903549e+03,
    2.2527545381182648e+04,
    2.2527545381182856e+04,
    8.8034271772744171e+02,
    8.8034271772744250e+02,
    2.1640913888698677e+02,
    1.0811836156273553e+03,
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
    3.2272517006745034e+03,
    0.0000000000000000e+00,
    2.7204739650200245e+03,
    3.1433678969242406e+04,
    3.1433678969242323e+04,
    4.3920855353045175e+02,
    4.3920855353045766e+02,
    1.4422974619925020e+02,
    4.3272388472231341e+02,
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
