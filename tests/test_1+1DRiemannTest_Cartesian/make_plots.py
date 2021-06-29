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

# analytic resutls
Analytic_e = loadtxt("1DRiemannTest_e.txt")
Analytic_v = loadtxt("1DRiemannTest_v.txt")
Analytic_theta = loadtxt("1DRiemannTest_theta.txt")

# numerical simulation resutls
numeric_data = loadtxt('../../1+1D_RiemannTest_tau_4.dat')
nSkip = int(len(numeric_data[:, 0])/nPointMax) + 1

# plot e
fig = plt.figure()
ax = plt.axes([0.13, 0.12, 0.81, 0.83])
plt.plot(Analytic_e[:, 0], Analytic_e[:, 1], color = 'k',
         linestyle = '-', label = r'Analytic')
plt.plot(numeric_data[::nSkip, 0], numeric_data[::nSkip, 1], color = 'r',
         linestyle = '', marker = 'X', label=r"MUSIC")

hl = plt.legend(loc=0)
hl.draw_frame(False)
plt.xlim(-3, 3.5)
plt.ylim(0, 20)
plt.xlabel(r'$x$ (fm)')
plt.ylabel(r'$e$ (GeV/fm$^3$)')
plt.savefig('RiemannTest_e.pdf')

# plot v
fig = plt.figure()
ax = plt.axes([0.13, 0.12, 0.81, 0.83])
plt.plot(Analytic_v[:, 0], Analytic_v[:, 1], color = 'k',
         linestyle = '-', label = r'Analytic')
plt.plot(numeric_data[::nSkip, 0], numeric_data[::nSkip, 2], color = 'r',
         linestyle = '', marker = 'X', label=r"MUSIC")

hl = plt.legend(loc=0)
hl.draw_frame(False)
plt.xlim(-3, 3.5)
plt.ylim(0, 0.6)
plt.xlabel(r'$x$ (fm)')
plt.ylabel(r'$v$')
plt.savefig('RiemannTest_v.pdf')

# plot theta
nSkip = 1
fig = plt.figure()
ax = plt.axes([0.13, 0.12, 0.81, 0.83])
plt.plot(Analytic_theta[:, 0], Analytic_theta[:, 1], color = 'k',
         linestyle = '-', label = r'Analytic')
plt.plot(numeric_data[::nSkip, 0], numeric_data[::nSkip, 3], color = 'r',
         linestyle = '', marker = '+', label=r"MUSIC")

hl = plt.legend(loc=0)
hl.draw_frame(False)
plt.xlim(-3, 3.5)
plt.ylim(-0.4, 0.4)
plt.xlabel(r'$x$ (fm)')
plt.ylabel(r'$\theta$ (1/fm)')
plt.savefig('RiemannTest_theta.pdf')
