#! /usr/bin/env python

import sys
sys.path.insert(0, '../../utilities')

from numpy import *
import matplotlib.pyplot as plt
from os import path
from CSplottools import getPlotElements

plotfontsize = 20
plotLinewidth = 2
plotMarkerSize = 8

tau_list = linspace(1, 1.1, 11)
tau_list = [1.02,]
filename = "../../Check_velocity_shear_tensor_tau_%g.dat"
# numerical simulation resutls

# plot sigma^{00}
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0
for itau, tau in enumerate(tau_list):
    data = loadtxt(filename % tau)
    plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(itau)
    plt.plot(data[:, 0], data[:, 6], color = plotColor,
             linestyle = '-', linewidth = plotLinewidth,
             label = r'$\tau = %g$ fm' % tau)
hl = plt.legend(loc=(2), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
#plt.ylim(-2.0, 2.0)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
#plt.yticks(linspace(-2.0, 2.0, 5), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$\sigma^{00}$', fontsize = plotfontsize)
plt.savefig('check_sigma00_x=y.pdf', format='pdf')
plt.close()

# plot sigma^{12}
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0
for itau, tau in enumerate(tau_list):
    data = loadtxt(filename % tau)
    plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
    plt.plot(data[:, 0], data[:, 11], color = plotColor,
             linestyle = '-', linewidth = plotLinewidth,
             label = r'$\tau = %g$ fm' % tau)
hl = plt.legend(loc=(2), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
#plt.ylim(-2.0, 2.0)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
#plt.yticks(linspace(-2.0, 2.0, 5), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$\sigma^{12}$', fontsize = plotfontsize)
plt.savefig('check_sigma12_x=y.pdf', format='pdf')
plt.close()

# plot sigma^{12}
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0
for itau, tau in enumerate(tau_list):
    data = loadtxt(filename % tau)
    plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
    plt.plot(data[:, 0], data[:, 17] - data[:, 3], color = plotColor,
             linestyle = '-', linewidth = plotLinewidth,
             label = r'$\tau = %g$ fm' % tau)
hl = plt.legend(loc=(2), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
#plt.ylim(-2.0, 2.0)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
#plt.yticks(linspace(-2.0, 2.0, 5), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'$\sigma^{12}$', fontsize = plotfontsize)
plt.savefig('check_ux_x=y.pdf', format='pdf')
plt.show()
plt.close()


# plot sigma^mu_mu
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
for itau, tau in enumerate(tau_list):
    data = loadtxt(filename % tau)
    trace = data[:, 2]**2. - data[:, 3]**2. - data[:, 4]**2. - data[:, 5]**2.
    plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(itau)
    plt.plot(data[:, 0], trace, color = plotColor,
             linestyle = '-', linewidth = plotLinewidth,
             label = r'$\tau = %g$ fm' % tau)
hl = plt.legend(loc=(2), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
plt.ylim(0.99999, 1.00001)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
#plt.yticks(linspace(-2.0, 2.0, 5), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'trace', fontsize = plotfontsize)
plt.savefig('check_usq_x=y.pdf', format='pdf')
plt.close()


# plot sigma^mu_mu
fig = plt.figure()
ax = plt.axes([0.14, 0.12, 0.81, 0.83])
iplot = 0
for itau, tau in enumerate(tau_list):
    data = loadtxt(filename % tau)
    trace = data[:, 6] - data[:, 10] - data[:, 13] - data[:, 15]
    plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(itau)
    plt.plot(data[:, 0], trace, color = plotColor,
             linestyle = '-', linewidth = plotLinewidth,
             label = r'$\tau = %g$ fm' % tau)
hl = plt.legend(loc=(2), fontsize = 17)
hl.draw_frame(False)
plt.xlim(-5.0, 5.0)
plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(0)
#plt.ylim(-1e-5, 1e-5)
plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
#plt.yticks(linspace(-2.0, 2.0, 5), color = 'k', size = plotfontsize)
plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
plt.ylabel(r'trace', fontsize = plotfontsize)
plt.savefig('check_sigma_trace_x=y.pdf', format='pdf')
plt.close()

# plot u_mu sigma^{mu i}
for i in range(4):
    fig = plt.figure()
    ax = plt.axes([0.14, 0.12, 0.81, 0.83])
    for itau, tau in enumerate(tau_list):
        data = loadtxt(filename % tau)
        if i == 0:
            u_dot_sigma = (data[:, 2]*data[:, 6] - data[:, 3]*data[:, 7]
                           - data[:, 4]*data[:, 8] - data[:, 5]*data[:, 9])
        elif i == 1:
            u_dot_sigma = (data[:, 2]*data[:, 7] - data[:, 3]*data[:, 10]
                           - data[:, 4]*data[:, 11] - data[:, 5]*data[:, 12])
        elif i == 2:
            u_dot_sigma = (data[:, 2]*data[:, 8] - data[:, 3]*data[:, 11]
                           - data[:, 4]*data[:, 13] - data[:, 5]*data[:, 14])
        elif i == 3:
            u_dot_sigma = (data[:, 2]*data[:, 9] - data[:, 3]*data[:, 12]
                           - data[:, 4]*data[:, 14] - data[:, 5]*data[:, 15])
        plotlinestyle, plotMarker, plotColor, plotshadowColor = getPlotElements(itau)
        plt.plot(data[:, 0], u_dot_sigma, color = plotColor,
                 linestyle = '-', linewidth = plotLinewidth,
                 label = r'$\tau = %g$ fm' % tau)
    hl = plt.legend(loc=(2), fontsize = 17)
    hl.draw_frame(False)
    plt.xlim(-5.0, 5.0)
    #plt.ylim(-1e-5, 1e-5)
    plt.xticks(linspace(-5.0, 5.0, 5), color = 'k', size = plotfontsize)
    #plt.yticks(linspace(-2.0, 2.0, 5), color = 'k', size = plotfontsize)
    plt.xlabel(r'$x$ (fm)', {'fontsize': plotfontsize})
    plt.ylabel(r'$u^\mu \sigma_{\mu\nu}$', fontsize = plotfontsize)
    plt.savefig('check_sigma_transversality_%d.pdf' % i, format='pdf')
    plt.close()
