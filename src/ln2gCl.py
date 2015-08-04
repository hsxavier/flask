#! /usr/bin/env python

import numpy as np
import sys
from scipy import interpolate as scint

#lsup   = 1000
# ATTENTION!! alpha are hard-coded here for density power spectrum:
alpha1 = 1
alpha2 = 1

# Get input
clin   = sys.argv[1]
xiout  = sys.argv[2]
l, cl  = np.loadtxt(clin, unpack=True)

# Set monopole to zero:
l      = np.insert(l,  0, 0.0)
cl     = np.insert(cl, 0, 0.0)
# Interpolate dipole:
c1     = scint.splev(1, scint.splrep(l, cl))
l      = np.insert(l, 1, 1)
cl     = np.insert(cl, 1, c1)

# Compute Xi (correlation function):
print "Get Xi..."
Ntheta = 2*l.size 
theta  = np.linspace(np.pi/Ntheta, ((Ntheta-1)*np.pi)/Ntheta, Ntheta)
mu     = np.cos(theta)
cl     = (2.0*l+1.0)*cl/(4.0*np.pi)#*np.exp(-l*(l+1)/(lsup*lsup))
xi     = np.polynomial.legendre.legval(mu, cl)

# Get Xi for associated gaussian variables:
print "Take the log..."
xi     = np.log(1.0 + xi/alpha1/alpha2)

# Compute Cls for associated gaussian variables:
print "Go back to Cls..."
cl     = np.polynomial.legendre.legfit(mu, xi, l.size-1, w=None)
cl     = 4.0*np.pi/(2.0*l+1.0)*cl

# Export these Cls to file:
np.savetxt(xiout, np.transpose([l, cl]))
