#! /usr/bin/env python

"""
USAGE:   prepClassInput.py <CLASS_CL_FILE> <FLASK_CL_PREFIX>
EXAMPLE: prepClassInput.py prec17_20_cl.dat Cl-
OUTPUT:  <FLASK_CL_PREFIX>f1z1f1z1.dat, <FLASK_CL_PREFIX>f1z1f1z2.dat, ...

This script takes a CLASS angular power spectra output file with header 
<CLASS_CL_FILE> and writes a separate file for each of the power spectra in 
this file containing two columns, l and C(l). Moreover, it removes the pi and 
l factors to return the pure C(l)s and transforms the lensing potential C(l)s 
[and cross C(l)s] to ones related to the convergence by applying l factors. 

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 10/sep/2015.
"""

import numpy as np
import sys
import math as m
import re

# Docstring output:
if len(sys.argv) != 1 + 2: 
    print(__doc__)
    sys.exit(1)


# Get input
classfile = sys.argv[1]
outprefix = sys.argv[2]



#########################
### Header operations ###
#########################

### Get file header (assumed to be last commented line):
fp = open(classfile, "r")
for line in fp:
    if line[0]=='#':
        header0=line
    else:
        break
fp.close()


### Correct for missing space between # and column names:
if header0[1]!=' ':
    header0 = header0[0]+' '+header0[1:]

### Add one to bin numbering if they start at zero (find such situation by looking for a [0] bin):
def addOne(s):
    if s.isdigit():
        return '['+str(int(s)+1)+']'
    else:
        return s
# Look for a zeroth bin, and do the dance:
if '[0]' in header0:
    header0 = ''.join([addOne(s) for s in re.split('\[|\]',header0)])


# Find out if there is more than one field type:
if "dens[" in header0 and "lens[" in header0 and "temp[" in header0:
    FieldTypes=3
elif "dens[" in header0 and "lens[" in header0:
    FieldTypes=2
elif "dens[" in header0 and "temp[" in header0:
    FieldTypes=4
elif "lens[" in header0 and "temp[" in header0:
    FieldTypes=6
else:
    FieldTypes=1

# Select column names from header:
header = header0.split()
header.pop(0)
# Generate output filenames:
# Let's adopt the field type priority in numbering: temperature, density and lensing.
outfiles=[]
for name in header:
    ini=name.find(":")+1;
    name = name[ini:]
    if FieldTypes==1:
        # Only one field type:
        name = name.replace("dens[","f1z")
        name = name.replace("lens[","f1z")
        name = name.replace("temp[","f1z")
    elif FieldTypes==2:
        # Density and convergence:
        name = name.replace("dens[","f1z")
        name = name.replace("lens[","f2z")
    elif FieldTypes==4:
        # Temperature and density:
        name = name.replace("temp[","f1z")
        name = name.replace("dens[","f2z")
    elif FieldTypes==6:
        # Temperature and convergence:
        name = name.replace("temp[","f1z")
        name = name.replace("lens[","f2z") 
    elif FieldTypes==3:
        # Temperature, density and convergence:
        name = name.replace("temp[","f1z")
        name = name.replace("dens[","f2z")
        name = name.replace("lens[","f3z")
    name = name.replace("]-","")
    name = name.replace("]","")
    name = outprefix+name+".dat"
    outfiles.append(name)



#######################
### Data operations ###
#######################    

# Load data:
classout = np.loadtxt(classfile, unpack=True)
ncols    = classout.shape[0]
length   = classout.shape[1]

# First column must be ell:
l = classout[0]

# Precompute basic factor:

fac0 = (2.0*m.pi)/l/(l+1)
fac1 = l*(l+1)/2.0   # The sign was supposed to be negative according to Hu 2000 (PRD 62:043007), 
                     # but that would lead to negative kappa-delta cross-Cls, which according to 
                     # numerical integration of delta auto-Cls is wrong. This is probably caused 
                     # by differences in definitions.

# LOOP over columns 
for i in range(1, ncols):
    factor = fac0

    # Discover if column is gal-gal, gal-lens ou lens-lens and compute appropriate factors:
    # (Factor-wise, temperature is treated as galaxy)
    first = header[i].find("lens") 
    if first != -1:
        # One lens so far, use fac1
        factor = factor*fac1
        second = header[i].find("lens", first+1)
        if second !=-1:
            # Two lens, use another fac1 
            factor = factor*fac1

    # Find if columns are CMB lensing potentials and convert them to convergence:
    first = header[i].find("phi") 
    if first != -1:
        # One lens so far, use fac1
        factor = factor*fac1
        second = header[i].find("phi", first+1)
        if second !=-1:
            # Two lens, use another fac1 
            factor = factor*fac1
    
    # Multiply column by factor
    classout[i] = classout[i]*factor
    
    # Export
    print "Writing file "+outfiles[i]
    np.savetxt(outfiles[i], np.transpose([l,classout[i]]), fmt=['%d','%e'])

print "Done."
