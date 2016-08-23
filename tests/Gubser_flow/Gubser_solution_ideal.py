#! /usr/bin/env python

from numpy import *

dim = 201
tau_list = [0.98, 1.0, 1.2, 1.5, 2.0, 3.0]
for tau in tau_list:
    q = 1.0
    e_0 = 1.0
    nb_0 = 0.5
    eps = 1e-10
    x = linspace(-5., 5., dim)
    y = linspace(-5., 5., dim)
    X, Y = meshgrid(x, y)
    
    r = sqrt(X**2. + Y**2.)
    e = (e_0/(tau**(4./3.))
         *((2.*q)**(8./3.))/((1. + 2.*(q**2.)*(tau**2. + r**2.)
                              + q**4.*(tau**2. - r**2.)**2.)**(4./3.)))
    nb = (nb_0/(tau**3.)
          *(4.*(q*tau)**2./(1. + 2.*(q**2.)*(tau**2. + r**2.)
                            + q**4.*(tau**2. - r**2.)**2.)))
    kappa = arctanh((2.*q**2.*tau*r)/(1. + q**2.*tau**2. + q**2.*r**2.))
    ux = sinh(kappa)*X/(r + eps)
    uy = sinh(kappa)*Y/(r + eps)
    
    filename = "y=0_tau=%4.2f_ideal.dat" % tau
    f = open(filename, "w")
    for i in range(dim):
        for j in range(dim):
            f.write("%.8e  %.8e  %.8e  %.8e  %.8e  %.8e\n"
                    % (x[i], y[j], e[j, i], nb[j, i], ux[j, i], uy[j, i]))
    f.close()

    # check
    #import matplotlib.pyplot as plt
    #plt.plot(x, e[100, :])
    #plt.plot(x, ux[100, :])
    #plt.show()
