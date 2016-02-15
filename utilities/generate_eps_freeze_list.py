#!/usr/bin/env python

from numpy import *
from os import path

eos_path = "../EOS/s95p-v1.2"

table1 = loadtxt(path.join(eos_path, "s95p-v1.2_dens1.dat"), skiprows=3)
table2 = loadtxt(path.join(eos_path, "s95p-v1.2_par1.dat"), skiprows=3)

temperature_list = linspace(0.1, 0.18, 17)
energy = -interp(-temperature_list, -table2[:, 0], -table1[:, 0])
pressure = -interp(-temperature_list, -table2[:, 0], -table1[:, 1])
entropy = -interp(-temperature_list, -table2[:, 0], -table1[:, 2])
mu_b = -interp(-temperature_list, -table2[:, 0], -table1[:, 3])

f = open(path.join(eos_path, "eps_freeze_list.dat"), "w")
f.write("#eps (GeV/fm^3)   P (GeV/fm^3)      s (fm^-3)    n_b (fm^-3)     ")
f.write("T (MeV)  mu_b (MeV)  mu_s (MeV) \n")
for i in range(len(temperature_list)):
        f.write("%.10f  %.10f  %.10f  %.10f  %.2f  %.5f  %.5f \n"
                % (energy[i], pressure[i], entropy[i], 0.0e0,
                   temperature_list[i]*1000.0, 0.0e0, 0.0e0))
f.close()
