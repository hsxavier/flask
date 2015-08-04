#! /usr/bin/env python

import numpy as np
import sys

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
        np.savetxt(outfile, Cl)

print "Done."
