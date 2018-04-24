#! /usr/bin/env python3

import sys
sys.path.insert(0, '../../utilities')

from numpy import *
import matplotlib.pyplot as plt
plt.style.use('CS_paper')
from os import path
from CSplottools import getPlotElements

plotfontsize = 20
plotLinewidth = 2
plotMarkerSize = 8
hbarC = 0.19733

# analytic resutls
AK_e_tau1p0 = loadtxt('e_baryon_init.dat')
AK_e_tau2p0 = loadtxt('e_baryon_t2.dat')
AK_e_tau5p0 = loadtxt('e_baryon_t5.dat')
AK_e_tau10p0 = loadtxt('e_baryon_t10.dat')
AK_e_tau20p0 = loadtxt('e_baryon_t20.dat')
AK_rhob_tau1p0 = loadtxt('rhoB_baryon_init.dat')
AK_rhob_tau2p0 = loadtxt('rhoB_baryon_t2.dat')
AK_rhob_tau5p0 = loadtxt('rhoB_baryon_t5.dat')
AK_rhob_tau10p0 = loadtxt('rhoB_baryon_t10.dat')
AK_rhob_tau20p0 = loadtxt('rhoB_baryon_t20.dat')

# numerical simulation resutls
numeric_tau_1p0 = loadtxt('../../1+1D_check_tau_1.dat')
numeric_tau_2p0 = loadtxt('../../1+1D_check_tau_2.dat')
numeric_tau_5p0 = loadtxt('../../1+1D_check_tau_5.dat')
numeric_tau_10p0 = loadtxt('../../1+1D_check_tau_10.dat')
numeric_tau_20p0 = loadtxt('../../1+1D_check_tau_20.dat')

# plot e
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(AK_e_tau1p0[::10, 0], AK_e_tau1p0[::10, 1],
         color = plotColor,
         linestyle = 'none', linewidth = plotLinewidth,
         marker = plotMarker, markersize = plotMarkerSize,
         label = r'$\tau = 1.0$ fm')
plt.plot(numeric_tau_1p0[:, 0], numeric_tau_1p0[:, 1], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(AK_e_tau2p0[::5, 0], AK_e_tau2p0[::5, 1],
         color = plotColor,
         linestyle = 'none', linewidth = plotLinewidth,
         marker = plotMarker, markersize = plotMarkerSize,
         label = r'$\tau = 2.0$ fm')
plt.plot(numeric_tau_2p0[:, 0], numeric_tau_2p0[:, 1], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(AK_e_tau5p0[::5, 0], AK_e_tau5p0[::5, 1],
         color = plotColor,
         linestyle = 'none', linewidth = plotLinewidth,
         marker = plotMarker, markersize = plotMarkerSize,
         label = r'$\tau = 5.0$ fm')
plt.plot(numeric_tau_5p0[:, 0], numeric_tau_5p0[:, 1], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(AK_e_tau10p0[::5, 0], AK_e_tau10p0[::5, 1],
         color = plotColor,
         linestyle = 'none', linewidth = plotLinewidth,
         marker = plotMarker, markersize = plotMarkerSize,
         label = r'$\tau = 10.0$ fm')
plt.plot(numeric_tau_10p0[:, 0], numeric_tau_10p0[:, 1], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(4)
plt.plot(AK_e_tau20p0[::5, 0], AK_e_tau20p0[::5, 1],
         color = plotColor,
         linestyle = 'none', linewidth = plotLinewidth,
         marker = plotMarker, markersize = plotMarkerSize,
         label = r'$\tau = 20.0$ fm')
plt.plot(numeric_tau_20p0[:, 0], numeric_tau_20p0[:, 1], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth)
hl = plt.legend(loc=0, fontsize = 20)
hl.draw_frame(False)
plt.xlim(-6.5, 6.5)
plt.ylim(1e-4, 100)
plt.yscale('log')
plt.xticks(linspace(-6.0, 6.0, 7), color = 'k', size = plotfontsize)
plt.yticks(10**linspace(-4., 2., 7), color = 'k', size = plotfontsize)
plt.xlabel(r'$\eta$', {'fontsize': plotfontsize})
plt.ylabel(r'$e$ (GeV/fm$^{3}$)', fontsize = plotfontsize)
plt.savefig('1+1D_e.pdf', format='pdf')

# plot rhob
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(AK_rhob_tau1p0[::10, 0], AK_rhob_tau1p0[::10, 1],
         color = plotColor,
         linestyle = 'none', linewidth = plotLinewidth,
         marker = plotMarker, markersize = plotMarkerSize,
         label = r'$\tau = 1.0$ fm')
plt.plot(numeric_tau_1p0[:, 0], numeric_tau_1p0[:, 2], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(AK_rhob_tau2p0[::5, 0], AK_rhob_tau2p0[::5, 1],
         color = plotColor,
         linestyle = 'none', linewidth = plotLinewidth,
         marker = plotMarker, markersize = plotMarkerSize,
         label = r'$\tau = 2.0$ fm')
plt.plot(numeric_tau_2p0[:, 0], numeric_tau_2p0[:, 2], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(AK_rhob_tau5p0[::5, 0], AK_rhob_tau5p0[::5, 1],
         color = plotColor,
         linestyle = 'none', linewidth = plotLinewidth,
         marker = plotMarker, markersize = plotMarkerSize,
         label = r'$\tau = 5.0$ fm')
plt.plot(numeric_tau_5p0[:, 0], numeric_tau_5p0[:, 2], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(AK_rhob_tau10p0[::5, 0], AK_rhob_tau10p0[::5, 1],
         color = plotColor,
         linestyle = 'none', linewidth = plotLinewidth,
         marker = plotMarker, markersize = plotMarkerSize,
         label = r'$\tau = 10.0$ fm')
plt.plot(numeric_tau_10p0[:, 0], numeric_tau_10p0[:, 2], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(4)
plt.plot(AK_rhob_tau20p0[::5, 0], AK_rhob_tau20p0[::5, 1],
         color = plotColor,
         linestyle = 'none', linewidth = plotLinewidth,
         marker = plotMarker, markersize = plotMarkerSize,
         label = r'$\tau = 20.0$ fm')
plt.plot(numeric_tau_20p0[:, 0], numeric_tau_20p0[:, 2], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth)
         
hl = plt.legend(loc=0, fontsize = 20)
hl.draw_frame(False)
plt.xlim(-6.5, 6.5)
plt.ylim(1e-7, 1.0)
plt.yscale('log')
plt.xticks(linspace(-6.0, 6.0, 7), color = 'k', size = plotfontsize)
plt.yticks(10**linspace(-7.0, 0.0, 8), color = 'k', size = plotfontsize)
plt.xlabel(r'$\eta$', {'fontsize': plotfontsize})
plt.ylabel(r'$n_B$ (1/fm$^3$)', fontsize = plotfontsize)
plt.savefig('1+1D_rhob.pdf', format='pdf')

