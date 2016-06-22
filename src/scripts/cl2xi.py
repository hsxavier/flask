#! /usr/bin/env python

"""
USAGE:   cl2xi.py <CL_IN> <XI_OUT>
EXAMPLE: cl2xi.py Cl-f1z1f1z1.dat Xi-f1z1f1z1.dat

This script take an angular power spectrum C(l) and computes the 
angular correlation function Xi(theta).

ATTENTION:
- It assumes that the input C(l) starts at l=2 and then 
  sets the monopole to zero and interpolate the dipole.
- It instroduces a hard-coded Gaussian suppression of power to 
  avoid oscillations in Xi caused by the hard bandlimit. 
  Check the lsup variable.

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 28/jul/2015 
"""

import numpy as np
import sys
from scipy import interpolate as scint

# Docstring output:
if len(sys.argv) != 1 + 2: 
    print(__doc__)
    sys.exit(1)


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
