#! /usr/bin/env python

"""
USAGE:   ln2gCl.py <input Cl^{ij}> <redshift of kappa i> <redshift of kappa j> <output Gaussian Cl^{ij}>
If dealing with galaxy density i, set <redshift of kappa i> = -1
EXAMPLE: ln2gCl.py Cl-f1z5f2z7.dat -1 0.65 gCl-f1z5f2z7.dat
OUTPUT:  C(l) for the associated Gaussian fields written to <output Gaussian Cl^{ij}>

This script takes an input angular power spectrum C(l) between lognormal fields i and j 
[an ASCII file with two columns: l and C(l)], assumes that it starts at l=2 (monopole 
is set to zero and dipole is interpolated), does not apply any suppression to C(l), 
compute the angular correlation function Xi, compute Xi of the associated Gaussian 
variables (assuming the lognormal field has mean zero and shift calculated from redshift 
by one of the formulas below in the code or is equal to one for the field i if 
<redshift of kappa i> is set to -1) and goes back to C(l). The current formula for 
shift(z) is the one obtained from FLASK density line of sight integration.  

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 06/aug/2015.
"""

import numpy as np
import sys
from scipy import interpolate as scint

# Convergence kappa shift formula from Hilbert, Hartlap & Schneider (2011)
def HilbertShift(z):
    return 0.008*z + 0.029*(z**2) - 0.0079*(z**3) + 0.00065*(z**4) 
# Convergence kappa shift formula from Xavier et al. 2015.
def XavierShift(z):
    a0 = 0.2;
    s0 = 0.568591;
    return a0*(((z*s0)**2 + z*s0 + 1)/(z*s0 + 1) - 1)


# Get input
clin   = sys.argv[1]
z1     = float(sys.argv[2])
z2     = float(sys.argv[3])
gclout = sys.argv[4]
l, cl  = np.loadtxt(clin, unpack=True)


# Set monopole to zero:            # Assuming the input Cl start at l=2.
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
cl     = (2.0*l+1.0)*cl/(4.0*np.pi)
xi     = np.polynomial.legendre.legval(mu, cl)


# Get Xi for associated gaussian variables:
print "Take the log..."
if z1 < 0: 
    alpha1 = 1
else:
    alpha1 = XavierShift(z1)
if z2 < 0: 
    alpha2 = 1
else:
    alpha2 = XavierShift(z2)
xi     = np.log(1.0 + xi/alpha1/alpha2)


# Compute Cls for associated gaussian variables:
print "Go back to Cls..."
cl     = np.polynomial.legendre.legfit(mu, xi, l.size-1, w=None)
cl     = 4.0*np.pi/(2.0*l+1.0)*cl


# Export these Cls to file:
np.savetxt(gclout, np.transpose([l, cl]))
