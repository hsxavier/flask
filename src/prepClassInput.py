#! /usr/bin/env python

import numpy as np
import sys
import math as m

# Get input
classfile = sys.argv[1]
outprefix = sys.argv[2]

# Get file header:
fp = open(classfile, "r")
for i, line in enumerate(fp):
    if i == 6:
        header0 = line
    elif i > 6:
        break
fp.close()
header = header0.split()
header.pop(0)

# Find out if there is more than one field type:
if "dens[" in header0 and "lens[" in header0:
    NfieldTypes=2
else:
    NfieldTypes=1

# Generate output filenames:
outfiles=[]
for name in header:
    ini=name.find(":")+1;
    name = name[ini:]
    if NfieldTypes==1:
        name = name.replace("dens[","f1z")
        name = name.replace("lens[","f1z")
    else:
        name = name.replace("dens[","f1z")
        name = name.replace("lens[","f2z")
    name = name.replace("]-","")
    name = name.replace("]","")
    name = outprefix+name+".dat"
    outfiles.append(name)
    
# Load data:
classout = np.loadtxt(classfile, unpack=True)
ncols    = classout.shape[0]
length   = classout.shape[1]

# First column must be ell:
l = classout[0]

# Precompute basic factor:

fac0 = (2.0*m.pi)/l/(l+1)
fac1 = (-1.0*l*(l+1))/2.0

# LOOP over columns 
for i in range(1, ncols):

    # Discover if column is gal-gal, gal-lens ou lens-lens and compute appropriate factors:
    factor = fac0
    first = header[i].find("lens")
    if first != -1:
        # One lens so far, use fac1
        factor = factor*fac1
        second = header[i].find("lens", first+1)
        if second !=-1:
            # Two lens, use another fac1 
            factor = factor*fac1
    
    # Multiply column by factor
    classout[i] = classout[i]*factor
    
    # Export
    print "Writing file "+outfiles[i]
    np.savetxt(outfiles[i], np.transpose([l,classout[i]]))

print "Done."
