#! /usr/bin/env python

from numpy import *

tau_list = [1.0, 1.2, 1.5, 2.0]
for tau in tau_list:
    q = 1.0
    e_0 = 1.0
    eps = 1e-10
    x = linspace(-5., 5., 201)
    y = linspace(-5., 5., 201)
    X, Y = meshgrid(x, y)
    
    r = sqrt(X**2. + Y**2.)
    e = (e_0/(tau**(4./3.))
         *((2.*q)**(8./3.))/((1. + 2.*(q**2.)*(tau**2. + r**2.)
                              + q**4.*(tau**2. - r**2.)**2.)**(4./3.)))
    kappa = arctanh((2.*q**2.*tau*r)/(1. + q**2.*tau**2. + q**2.*r**2.))
    ux = sinh(kappa)*X/(r + eps)
    uy = sinh(kappa)*Y/(r + eps)
    
    filename = "y=0_tau=%3.1f_ideal.dat" % tau
    f = open(filename, "w")
    for i in range(201):
        for j in range(201):
            f.write("%.8e  %.8e  %.8e  %.8e  %.8e\n"
                    % (x[i], y[j], e[j, i], ux[j, i], uy[j, i]))
    f.close()

    # check
    #import matplotlib.pyplot as plt
    #plt.plot(x, e[100, :])
    #plt.plot(x, ux[100, :])
    #plt.show()
