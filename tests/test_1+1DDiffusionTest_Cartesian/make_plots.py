#! /usr/bin/env python

import sys
sys.path.insert(0, '../../utilities')

from numpy import *
import matplotlib.pyplot as plt
plt.style.use('CS_paper')
from os import path
from CSplottools import getPlotElements

plotfontsize = 20
plotLinewidth = 2
plotMarkerSize = 5
nPointMax = 50

# numerical simulation resutls
numeric_data = loadtxt('../../1+1D_DiffusionTest_tau_4.dat')

# analytic resutls
Analytic_rhoB = zeros(len(numeric_data[:, 0]))
idx = (numeric_data[:, 0] > 1) & (numeric_data[:, 0] < 2)
Analytic_rhoB[idx] = 1.


# plot rhoB
fig = plt.figure()
ax = plt.axes([0.13, 0.12, 0.81, 0.83])
plt.plot(numeric_data[:, 0], Analytic_rhoB, color = 'k',
         linestyle = '-', label = r'Analytic')
plt.plot(numeric_data[:, 0], numeric_data[:, 2], color = 'r',
         linestyle = '', marker = 'x', label=r"MUSIC")

hl = plt.legend(loc=0)
hl.draw_frame(False)
plt.xlim(-1, 3)
plt.ylim(0, 1.2)
plt.xlabel(r'$x$ (fm)')
plt.ylabel(r'$\rho_B$ (1/fm$^3$)')
plt.savefig('DiffusionTest_rhoB.pdf')
