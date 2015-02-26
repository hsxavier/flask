#! /usr/bin/env python

import sys
import numpy as np

# Internal definitions:
GalMean    = 0.0
KappaMean  = 0.0
GalShift   = 1.0
FixKappa   = 0
KappaShift = 1.0
GalType    = 1
KappaType  = 2
shiftfile  = "/home/skems/pkphotoz/prog/corrlnfields/data/Hilbert2011_shifts.dat"

# Load shift file:
# ATTENTION!! Values not extrapolated: extremes are used!
fileZ, fileSh = np.loadtxt(shiftfile, unpack=True)


# Get input:
classinput = sys.argv[1]
fieldinfo  = sys.argv[2]

# Read input file:
fin = open(classinput, 'r')
lines = fin.readlines()

HasGal   = 0
HasKappa = 0
# LOOP over lines:
for line in lines:
    # Avoid comments:
    if line[0]!='#':
        # Find an input: 
        if '=' in line:
            eqpos = line.index('=')
            # Count number of fields:
            if 'output' == line[:eqpos] or 'output'==line.split()[0]:
                words = line[eqpos+1:].split(',')
                nf    = len(words)
                if 'nCl' in line[eqpos+1:]:
                    HasGal   = 1
                if 'sCl' in line[eqpos+1:]:
                    HasKappa = 1
                if HasGal+HasKappa!=nf:
                    print "ERROR: unkown field settings in output."
                    sys.exit()
            # Count number of redshift bins and get their means:
            if 'selection_mean' == line[:eqpos] or 'selection_mean'==line.split()[0]:
                words = line[eqpos+1:].split(',')
                nz    = len(words)
                meanz = [float(z) for z in words]
            # Get the bins width:
            if 'selection_width' == line[:eqpos] or 'selection_width'==line.split()[0]:
                words      = line[eqpos+1:].split(',')
                halfwidthz = [float(dz) for dz in words]
                if len(halfwidthz) == 1:
                    halfwidithz = halfwidthz * nz;
                elif len(halfwidthz) != nz:
                    print "ERROR: redshift bins are not well defined."
                    sys.exit()

fin = open(classinput, 'r')
lines = fin.readlines()

HasGal   = 0
HasKappa = 0
# LOOP over lines:
for line in lines:
    # Avoid comments:
    if line[0]!='#':
        # Find an input: 
        if '=' in line:
            eqpos = line.index('=')
            # Count number of fields:
            if 'output' == line[:eqpos] or 'output'==line.split()[0]:
                words = line[eqpos+1:].split(',')
                nf    = len(words)
                if 'nCl' in line[eqpos+1:]:
                    HasGal   = 1
                if 'sCl' in line[eqpos+1:]:
                    HasKappa = 1
                if HasGal+HasKappa!=nf:
                    print "ERROR: unkown field settings in output."
                    sys.exit()
            # Count number of redshift bins and get their means:
            if 'selection_mean' == line[:eqpos] or 'selection_mean'==line.split()[0]:
                words = line[eqpos+1:].split(',')
                nz    = len(words)
                meanz = [float(z) for z in words]
            # Get the bins width:
            if 'selection_width' == line[:eqpos] or 'selection_width'==line.split()[0]:
                words      = line[eqpos+1:].split(',')
                halfwidthz = [float(dz) for dz in words]
                if len(halfwidthz) == 1:
                    halfwidthz = halfwidthz * nz;

# Functions for output:
def mean(f):
    m = [GalMean, KappaMean]
    if nf==2:
        return m[f-1]
    elif nf==1:
        if HasGal==1: 
            return GalMean
        if HasKappa==1:
            return KappaMean
# Functions for output:
def shift(f, z):
    s = [GalShift, KappaShift]
    if nf==2:
        if f==1:
            return GalShift
        if f==2:
            if FixKappa==1:
                return KappaShift
            if FixKappa==0:
                return np.interp(z, fileZ, fileSh)
    elif nf==1:
        if HasGal==1: 
            return GalShift
        if HasKappa==1:
            if FixKappa==1:
                return KappaShift
            if FixKappa==0:
                return np.interp(z, fileZ, fileSh)
# Functions for output:
def ftype(f):
    t = [GalType, KappaType]
    if nf==2:
        return t[f-1]
    elif nf==1:
        if HasGal==1: 
            return GalType
        if HasKappa==1:
            return KappaType


# Output:
fout = open(fieldinfo, 'w')

print >>fout, "# Field number, z bin number, mean, shift, field type, zmin, zmax"
print >>fout, "# Types: 1-galaxies 2-shear\n"

for f in range(1,nf+1):
    for z in range(1,nz+1):
        print >>fout,'{0:6d} {1:6d}   {2:.4f}   {3:.4f} {4:6d}   {5:.4f}   {6:.4f}'.format(f, z, mean(f), shift(f, meanz[z-1]), ftype(f), meanz[z-1]-halfwidthz[z-1], meanz[z-1]+halfwidthz[z-1])

