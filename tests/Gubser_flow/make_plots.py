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

# analytic resutls
Analytic_tau1p2 = loadtxt('y=0_tau=1.2_SemiAnalytic.dat')
Analytic_tau1p5 = loadtxt('y=0_tau=1.5_SemiAnalytic.dat')
Analytic_tau2p0 = loadtxt('y=0_tau=2_SemiAnalytic.dat')

# numerical simulation resutls
numeric_tau_1p0 = loadtxt('../../Gubser_flow_check_tau_1.dat')
numeric_tau_1p2 = loadtxt('../../Gubser_flow_check_tau_1.2.dat')
numeric_tau_1p5 = loadtxt('../../Gubser_flow_check_tau_1.5.dat')
numeric_tau_2p0 = loadtxt('../../Gubser_flow_check_tau_2.dat')
idx = fabs(numeric_tau_1p2[:, 1]) < 1e-6

# plot ux
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(numeric_tau_1p0[idx, 0], numeric_tau_1p0[idx, 5], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.0$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(numeric_tau_1p2[idx, 0], numeric_tau_1p2[idx, 5], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p2[:, 0], Analytic_tau1p2[:, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.2$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(numeric_tau_1p5[idx, 0], numeric_tau_1p5[idx, 5], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p5[:, 0], Analytic_tau1p5[:, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.5$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(numeric_tau_2p0[idx, 0], numeric_tau_2p0[idx, 5], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau2p0[:, 0], Analytic_tau2p0[:, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 2.0$ fm')
         
hl = plt.legend(loc=(2), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.ylim(-2.2, 2.2)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
plt.yticks(linspace(-2.0, 2.0, 5), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$u^x$', fontsize = plotfontsize)
plt.savefig('Gubser_ux_y=0.pdf', format='pdf')

# plot pixx
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(numeric_tau_1p0[idx, 0], numeric_tau_1p0[idx, 7], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.0$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(numeric_tau_1p2[idx, 0], numeric_tau_1p2[idx, 7], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p2[:, 0], Analytic_tau1p2[:, 5], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.2$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(numeric_tau_1p5[idx, 0], numeric_tau_1p5[idx, 7], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p5[:, 0], Analytic_tau1p5[:, 5], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.5$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(numeric_tau_2p0[idx, 0], numeric_tau_2p0[idx, 7], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau2p0[:, 0], Analytic_tau2p0[:, 5], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 2.0$ fm')
         
hl = plt.legend(loc=(3), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plt.ylim(-0.45, 0.0)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
plt.yticks(linspace(-0.4, 0.0, 5), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$\pi^{xx}$ (GeV/fm$^3$)', fontsize = plotfontsize)
plt.savefig('Gubser_pixx_y=0.pdf', format='pdf')


# analytic resutls
Analytic_tau1p2 = loadtxt('y=x_tau=1.2_SemiAnalytic.dat')
Analytic_tau1p5 = loadtxt('y=x_tau=1.5_SemiAnalytic.dat')
Analytic_tau2p0 = loadtxt('y=x_tau=2_SemiAnalytic.dat')

idx = fabs(numeric_tau_1p0[:, 0] - numeric_tau_1p0[:, 1]) < 1e-6

# plot ux
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(numeric_tau_1p0[idx, 0], numeric_tau_1p0[idx, 4], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.0$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(numeric_tau_1p2[idx, 0], numeric_tau_1p2[idx, 4], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p2[:, 0], Analytic_tau1p2[:, 2], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.2$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(numeric_tau_1p5[idx, 0], numeric_tau_1p5[idx, 4], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p5[:, 0], Analytic_tau1p5[:, 2], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.5$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(numeric_tau_2p0[idx, 0], numeric_tau_2p0[idx, 4], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau2p0[:, 0], Analytic_tau2p0[:, 2], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 2.0$ fm')
         
hl = plt.legend(loc=(2), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.ylim(0.0, 0.25)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
plt.yticks(linspace(0.0, 0.25, 6), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$T$ (GeV)', fontsize = plotfontsize)
plt.savefig('Gubser_T_y=x.pdf', format='pdf')


# plot ux
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(numeric_tau_1p0[idx, 0], numeric_tau_1p0[idx, 5], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.0$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(numeric_tau_1p2[idx, 0], numeric_tau_1p2[idx, 5], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p2[:, 0], Analytic_tau1p2[:, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.2$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(numeric_tau_1p5[idx, 0], numeric_tau_1p5[idx, 5], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p5[:, 0], Analytic_tau1p5[:, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.5$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(numeric_tau_2p0[idx, 0], numeric_tau_2p0[idx, 5], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau2p0[:, 0], Analytic_tau2p0[:, 3], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 2.0$ fm')
         
hl = plt.legend(loc=(2), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.ylim(-2.0, 2.0)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
plt.yticks(linspace(-2.0, 2.0, 5), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$u^x$', fontsize = plotfontsize)
plt.savefig('Gubser_ux_y=x.pdf', format='pdf')

# plot pixx
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(numeric_tau_1p0[idx, 0], numeric_tau_1p0[idx, 7], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.0$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(numeric_tau_1p2[idx, 0], numeric_tau_1p2[idx, 7], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p2[:, 0], Analytic_tau1p2[:, 5], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.2$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(numeric_tau_1p5[idx, 0], numeric_tau_1p5[idx, 7], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p5[:, 0], Analytic_tau1p5[:, 5], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.5$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(numeric_tau_2p0[idx, 0], numeric_tau_2p0[idx, 7], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau2p0[:, 0], Analytic_tau2p0[:, 5], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 2.0$ fm')
         
hl = plt.legend(loc=(3), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plt.ylim(-0.32, 0.0)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
plt.yticks(linspace(-0.3, 0.0, 4), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$\pi^{xx}$ (GeV/fm$^3$)', fontsize = plotfontsize)
plt.savefig('Gubser_pixx_y=x.pdf', format='pdf')

# plot pixy
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0

plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.plot(numeric_tau_1p0[idx, 0], numeric_tau_1p0[idx, 9], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.0$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(1)
plt.plot(numeric_tau_1p2[idx, 0], numeric_tau_1p2[idx, 9], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p2[:, 0], Analytic_tau1p2[:, 7], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.2$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(2)
plt.plot(numeric_tau_1p5[idx, 0], numeric_tau_1p5[idx, 9], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau1p5[:, 0], Analytic_tau1p5[:, 7], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 1.5$ fm')
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(3)
plt.plot(numeric_tau_2p0[idx, 0], numeric_tau_2p0[idx, 9], color = 'k',
         linestyle = '', linewidth = plotLinewidth,
         marker = 'o', markersize = plotMarkerSize)
plt.plot(Analytic_tau2p0[:, 0], Analytic_tau2p0[:, 7], color = plotColor,
         linestyle = '-', linewidth = plotLinewidth,
         label = r'$\tau = 2.0$ fm')
         
hl = plt.legend(loc=(3), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plt.ylim(-0.11, 0.01)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
plt.yticks(linspace(-0.1, 0.0, 6), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$\pi^{xx}$ (GeV/fm$^3$)', fontsize = plotfontsize)
plt.savefig('Gubser_pixy_y=x.pdf', format='pdf')
