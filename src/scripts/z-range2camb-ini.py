#!/usr/bin/env python

"""
Prints the CAMB sources .ini file window information to the screen according to 
the redshift range and number of bins provided. This output can be pasted in a 
CAMB sources .ini file. The bins are equally spaced and the width is set to make 
them contiguous (assuming the bin shape is top-hat). Currently this only accepts 
counts entries.

USAGE:   z-range2camb-ini.py <zmin> <zmax> <Nbins>
EXAMPLE: z-range2camb-ini.py 0.55 0.65 40


Like linZbins.py, zmin is the lower bound of the first bin and zmax is the upper 
bound of the last bin.

Written by: Henrique S. Xavier, hsxavier@if.usp.br, on 22-jun-2016.
"""

import sys
import numpy as np

# Docstring output:
if len(sys.argv) != 4: 
    print(__doc__)
    sys.exit(1)

# Input:
zmin  = float(sys.argv[1])
zmax  = float(sys.argv[2])
Nbins = int(sys.argv[3])

# Compute stuff:
halfDz = (zmax-zmin)/Nbins/2.0
central = np.linspace(zmin+halfDz, zmax-halfDz, Nbins)

# Output:
textDz = str(halfDz)
print "num_redshiftwindows =", Nbins 
for i in range(0, Nbins):
    textN = str(i+1)
    textZ = str(central[i])
    print ""
    print "redshift("+textN+") = "+textZ
    print "redshift_kind("+textN+") = counts"
    print "redshift_bias("+textN+") = 1.0"
    print "redshift_sigma("+textN+") = "+textDz
    print "redshift_dlog10Ndm("+textN+") = 0"
    
