#! /usr/bin/env python

import numpy as np
import sys
from scipy import interpolate as scint

lsup=1000


# Get input
clin   = sys.argv[1]
xiout  = sys.argv[2]
l, cl  = np.loadtxt(clin, unpack=True)

# Set monopole to zero:
l  = np.insert(l,  0, 0.0)
cl = np.insert(cl, 0, 0.0)
# Interpolate dipole:
c1 = scint.splev(1, scint.splrep(l, cl))
l  = np.insert(l, 1, 1)
cl = np.insert(cl, 1, c1)

# Compute Xi (correlation function)
Ntheta = 2*l.size 
theta  = np.linspace(np.pi/Ntheta, ((Ntheta-1)*np.pi)/Ntheta, Ntheta)
cl     = (2.0*l+1.0)*cl/(4.0*np.pi)*np.exp(-l*(l+1)/(lsup*lsup))
xi     = np.polynomial.legendre.legval(np.cos(theta), cl)

# Export Xi to file:
np.savetxt(xiout, np.transpose([np.degrees(theta), xi]))
