#!/usr/bin/env python3

import numpy as np

TargetSum = [
    1.0150499999999958e+05,
    1.0150500000000000e+05,
    1.1694855392337635e+04,
    0.0000000000000000e+00,
    3.4153889251022010e+03,
    1.7545619295415341e+04,
    1.7545619295415094e+04,
    1.5662522124292821e+03,
    1.5662522124292902e+03,
    2.9854691247734809e+02,
    2.1921159761453432e+03,
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
    6.8034169967035232e+03,
    0.0000000000000000e+00,
    3.0882070167647576e+03,
    2.2546426095969597e+04,
    2.2546426095969706e+04,
    8.8039459213140810e+02,
    8.8039459213140742e+02,
    2.1657204360811133e+02,
    1.0811735486743969e+03,
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
    3.2268089917055659e+03,
    0.0000000000000000e+00,
    2.7195662855431551e+03,
    3.1522342758154817e+04,
    3.1522342758154675e+04,
    4.3947779142881114e+02,
    4.3947779142881808e+02,
    1.4480100682104359e+02,
    4.3264469854643522e+02,
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
