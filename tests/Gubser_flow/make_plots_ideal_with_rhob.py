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
plotMarkerSize = 8
hbarC = 0.19733

# analytic resutls
Analytic_tau1p0 = loadtxt('y=0_tau=1.00_ideal.dat')
Analytic_tau1p2 = loadtxt('y=0_tau=1.20_ideal.dat')
Analytic_tau1p5 = loadtxt('y=0_tau=1.50_ideal.dat')
Analytic_tau2p0 = loadtxt('y=0_tau=2.00_ideal.dat')
Analytic_tau3p0 = loadtxt('y=0_tau=3.00_ideal.dat')

# numerical simulation resutls
numeric_tau_1p0 = loadtxt('../../Gubser_flow_check_tau_1.dat')
numeric_tau_1p2 = loadtxt('../../Gubser_flow_check_tau_1.2.dat')
numeric_tau_1p5 = loadtxt('../../Gubser_flow_check_tau_1.5.dat')
numeric_tau_2p0 = loadtxt('../../Gubser_flow_check_tau_2.dat')
numeric_tau_3p0 = loadtxt('../../Gubser_flow_check_tau_3.dat')
idx = fabs(numeric_tau_1p2[:, 1]) < 1e-6

# plot e
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(Analytic_tau1p0[idx, 0], Analytic_tau1p0[idx, 2]*hbarC,
         color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4,
         label = r'$\tau = 1.0$ fm')
plt.plot(numeric_tau_1p0[idx, 0], numeric_tau_1p0[idx, 2], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(Analytic_tau1p2[idx, 0], Analytic_tau1p2[idx, 2]*hbarC,
         color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4,
         label = r'$\tau = 1.2$ fm')
plt.plot(numeric_tau_1p2[idx, 0], numeric_tau_1p2[idx, 2], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(Analytic_tau1p5[idx, 0], Analytic_tau1p5[idx, 2]*hbarC,
         color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4,
         label = r'$\tau = 1.5$ fm')
plt.plot(numeric_tau_1p5[idx, 0], numeric_tau_1p5[idx, 2], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(Analytic_tau2p0[idx, 0], Analytic_tau2p0[idx, 2]*hbarC,
         color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4,
         label = r'$\tau = 2.0$ fm')
plt.plot(numeric_tau_2p0[idx, 0], numeric_tau_2p0[idx, 2], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(4)
plt.plot(Analytic_tau3p0[idx, 0], Analytic_tau3p0[idx, 2]*hbarC,
         color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4,
         label = r'$\tau = 2.0$ fm')
plt.plot(numeric_tau_3p0[idx, 0], numeric_tau_3p0[idx, 2], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth)
         
hl = plt.legend(loc=(2), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plt.ylim(1e-4, 1.0)
plt.yscale('log')
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
plt.yticks(10**linspace(-4, 0, 5), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$e$ (GeV/fm$^{3}$)', fontsize = plotfontsize)
plt.savefig('Gubser_e.pdf', format='pdf')


# plot rhob
fig = plt.figure()
ax = plt.axes([0.13, 0.12, 0.82, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(Analytic_tau1p0[idx, 0], Analytic_tau1p0[idx, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4)
plt.plot(numeric_tau_1p0[idx, 0], numeric_tau_1p0[idx, 3], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth,
         label = r'$\tau = 1.0$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(Analytic_tau1p2[idx, 0], Analytic_tau1p2[idx, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4)
plt.plot(numeric_tau_1p2[idx, 0], numeric_tau_1p2[idx, 3], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth,
         label = r'$\tau = 1.2$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(Analytic_tau1p5[idx, 0], Analytic_tau1p5[idx, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4)
plt.plot(numeric_tau_1p5[idx, 0], numeric_tau_1p5[idx, 3], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth,
         label = r'$\tau = 1.5$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(Analytic_tau2p0[idx, 0], Analytic_tau2p0[idx, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4)
plt.plot(numeric_tau_2p0[idx, 0], numeric_tau_2p0[idx, 3], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth,
         label = r'$\tau = 2.0$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(4)
plt.plot(Analytic_tau3p0[idx, 0], Analytic_tau3p0[idx, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4)
plt.plot(numeric_tau_3p0[idx, 0], numeric_tau_3p0[idx, 3], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth,
         label = r'$\tau = 3.0$ fm')
         
hl = plt.legend(loc=0, ncol = 2)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plt.ylim(2e-3, 5.0)
plt.yscale('log')
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
plt.yticks(10**linspace(-2, 0, 3), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$n_B$ (1/fm$^{3}$)', fontsize = plotfontsize)
plt.savefig('Gubser_rhob.pdf', format='pdf')

# plot ux
fig = plt.figure()
ax = plt.axes([0.13, 0.12, 0.82, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(Analytic_tau1p0[idx, 0], Analytic_tau1p0[idx, 4], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4)
plt.plot(numeric_tau_1p0[idx, 0], numeric_tau_1p0[idx, 5], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth,
         label = r'$\tau = 1.0$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(Analytic_tau1p2[idx, 0], Analytic_tau1p2[idx, 4], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4)
plt.plot(numeric_tau_1p2[idx, 0], numeric_tau_1p2[idx, 5], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth,
         label = r'$\tau = 1.2$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(Analytic_tau1p5[idx, 0], Analytic_tau1p5[idx, 4], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4)
plt.plot(numeric_tau_1p5[idx, 0], numeric_tau_1p5[idx, 5], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth,
         label = r'$\tau = 1.5$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(Analytic_tau2p0[idx, 0], Analytic_tau2p0[idx, 4], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4)
plt.plot(numeric_tau_2p0[idx, 0], numeric_tau_2p0[idx, 5], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth,
         label = r'$\tau = 2.0$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(4)
plt.plot(Analytic_tau3p0[idx, 0], Analytic_tau3p0[idx, 4], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth, alpha = 0.4)
plt.plot(numeric_tau_3p0[idx, 0], numeric_tau_3p0[idx, 5], color = plotColor,
         linestyle = '--', linewidth = plotLinewidth,
         label = r'$\tau = 3.0$ fm')
         
hl = plt.legend(loc=0, fontsize = 20)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plt.ylim(-3.5, 3.5)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
plt.yticks(linspace(-3.0, 3.0, 7), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$u^x$', fontsize = plotfontsize)
plt.savefig('Gubser_ux.pdf', format='pdf')

