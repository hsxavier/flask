#! /usr/bin/env python

import sys
import numpy as np
from scipy import interpolate as scint

Clfile   = sys.argv[1]
#lsup     = 4000
#supindex = 6

l, cl = np.loadtxt(Clfile, unpack=True)

# Set monopole to zero:
#l  = np.insert(l,  0, 0.0)
#cl = np.insert(cl, 0, 0.0)
# Interpolate dipole:
#c1 = scint.splev(1, scint.splrep(l, cl))
#l  = np.insert(l, 1, 1)
#cl = np.insert(cl, 1, c1)


#variance = np.sum((2.0*l+1.0)/(4.0*np.pi)*cl*np.exp(-(l/lsup)**supindex))
variance = np.sum((2.0*l+1.0)/(4.0*np.pi)*cl)

print variance
