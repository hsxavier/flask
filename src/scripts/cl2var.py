#! /usr/bin/env python

"""
USAGE:   cl2var.py <INPUT_CL>
EXAMPLE: cl2var.py Cl-f1z2f1z2.dat

This script computes the expected variance in a pixel (with a real space 
window function equal to a Dirac delta function) from a angular power 
spectrum C(l). It assumes the input C(l) starts at l=2, and it can interpolate 
and use the dipole (with the monopole set to zero) or not (current version 
does not include the dipole). It can include an exponential suppression to 
the C(l) or not (current version does not include the suppression). These 
options are hard-coded. 

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 04/aug/2015.
"""

import sys
import numpy as np
from scipy import interpolate as scint

# Docstring output:
if len(sys.argv) != 1 + 1: 
    print(__doc__)
    sys.exit(1)


# Get input parameters:
Clfile   = sys.argv[1]
#lsup     = 4000
#supindex = 6

# Load Cl:
l, cl = np.loadtxt(Clfile, unpack=True)

# Set monopole to zero:
#l  = np.insert(l,  0, 0.0)
#cl = np.insert(cl, 0, 0.0)
# Interpolate dipole:
#c1 = scint.splev(1, scint.splrep(l, cl))
#l  = np.insert(l, 1, 1)
#cl = np.insert(cl, 1, c1)

# Compute the variance in the pixel (Dirac delta function) according to the input C(l):
#variance = np.sum((2.0*l+1.0)/(4.0*np.pi)*cl*np.exp(-(l/lsup)**supindex))
variance = np.sum((2.0*l+1.0)/(4.0*np.pi)*cl)

print variance
