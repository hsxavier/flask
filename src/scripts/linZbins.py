#! /usr/bin/env python

"""
USAGE: linZbins.py <zmin> <zmax> <nbins>

Generates a CSV list of <nbins> central values of bins that covers from 
<zmin> to <zmax>, where <zmin> is the lower limit of the first bin and 
<zmax> is the upper limit of the last bin. Their spacing is constant.

RETURNS:
bin half-width
list of central values

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 26/oct/2015.
"""

import sys
import numpy as np

# Docstring output:
if len(sys.argv) != 1 + 3: 
    print(__doc__)
    sys.exit(1)

# Internal definitions:

# Load shift file:
# ATTENTION!! Values not extrapolated: extremes are used!
#fileZ, fileSh = np.loadtxt(shiftfile, unpack=True)


# Get input:
zmin  = float(sys.argv[1])
zmax  = float(sys.argv[2])
nbins = int  (sys.argv[3])

# Compute stuff:
halfDz = (zmax-zmin)/nbins/2.0
central = np.linspace(zmin+halfDz, zmax-halfDz, nbins)

# Output:
print (halfDz)
for zi in central[:-1]: 
    sys.stdout.write(str(zi))
    sys.stdout.write(", ")
print (central[-1])
