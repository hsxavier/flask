#! /usr/bin/env python

"""
USAGE:   prepCambInput.py <CAMB_COV_FILE> <FIELDS_INFO_FILE> <FLASK_CL_PREFIX>
EXAMPLE: prepCambInput.py test08_scalCovCls.dat test08/fields-info.dat Cl-
OUTPUT:  <FLASK_CL_PREFIX>f1z1f1z1.dat, <FLASK_CL_PREFIX>f1z1f1z2.dat, ...

This script takes a CAMB angular power spectra output file with repeated columns and 
the FLASK FIELDS_INFO file <FIELDS_INFO_FILE> prepared by the camb2info.py script 
and writes a separate file for each of the power spectra in the <CAMB_COV_FILE> file 
containing two columns, l and C(l). Moreover, it removes the pi and 
l factors to return the pure C(l)s. 

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 04/aug/2015.
"""


import numpy as np
import sys

# Docstring output:
if len(sys.argv) != 1 + 3: 
    print(__doc__)
    sys.exit(1)


# Function to get position of Cl in CAMBsources Cov file:
def CovPosition(f1, f2, Nfields):
    return (Nfields+3)*(3+f1-1) + (3+f2)

# Get input:
cambfile  = sys.argv[1]
infofile  = sys.argv[2]
outprefix = sys.argv[3]
camb    = np.loadtxt(cambfile, unpack=True)
info    = np.loadtxt(infofile)
Nfields = info.shape[0]
print 'Nfields =', Nfields

# Prepare factors:
l   = camb[0]
fac = (2.0*np.pi)/l/(l+1)

# LOOP over fields x fields:
for f1 in range(1,Nfields+1):
    for f2 in range(f1,Nfields+1):
        # Prepare Cls:
        Cl = np.transpose( [l, fac*camb[CovPosition(f1,f2,Nfields)]] )
        # Export:
        outfile = outprefix +'f'+str(int(info[f1-1][0]))+'z'+str(int(info[f1-1][1]))+'f'+str(int(info[f2-1][0]))+'z'+str(int(info[f2-1][1]))+'.dat'
        print "Writing file "+outfile
        np.savetxt(outfile, Cl, fmt=['%d','%e'])

print "Done."
