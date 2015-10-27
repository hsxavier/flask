#! /usr/bin/env python

"""
USAGE:   camb2info.py <CAMB_INPUT> <FIELDS_INFO_FILENAME> 
EXAMPLE: camb2info.py lcdm_10bins.ini lcdm_10bins/fields-info.dat

This script takes a CAMBsources input file and creates a file used by 
FLASK describing the simulated fields and redshift bins (basically a 
table of field and redshift bin IDs, mean values, shift parameters and 
redshift ranges).

There are three methods for computing the shift parameter for convergence:
A table read from a file which is interpolated, the formula from 
Hilbert, Hartlap & Schneider (2011), and a formula computed from FLASK 
density line of sight integration. The latter is currently used (the other 
ones are commented).

Written by Henrique S. Xavier, hsxavier@if.usp.br, 05/aug/2015.
"""

import sys
import numpy as np
import re

# Internal definitions:
GalMean    = 0.0
KappaMean  = 0.0
GalShift   = 1.0
FixKappa   = 0
KappaShift = 1.0
GalType    = 1
KappaType  = 2
 
#shiftfile  = "/home/skems/pkphotoz/prog/corrlnfields/data/Hilbert2011_shifts.dat"
#shiftfile  = "/home/skems/pkphotoz/prog/corrlnfields/data/k0_empty_LCDM_Om30.dat"

# Load shift file:
# ATTENTION!! Values not extrapolated: extremes are used!
#fileZ, fileSh = np.loadtxt(shiftfile, unpack=True)


# Get input:
cambinput = sys.argv[1]
fieldinfo  = sys.argv[2]

# Read input file:
fin = open(cambinput, 'r')
lines = fin.readlines()

# LOOP over lines:
for line in lines:
    # Avoid comments:
    if line[0]!='#':
        # Find an input: 
        if '=' in line:
            eqpos        = line.index('=')
            parameter    = line[:eqpos].split()[0]
            # Find out number of redshift windows:
            if 'num_redshiftwindows' == parameter:
                Nwindows = int(line[eqpos+1:].split()[0])
                zmean    = np.zeros(Nwindows)
                zwidth   = np.zeros(Nwindows)
                bias     = np.zeros(Nwindows)
                ftype    = np.zeros(Nwindows)
                zindex   = np.zeros(Nwindows)
                findex   = np.zeros(Nwindows)
            # Get mean redshifts:
            if re.match(re.compile('redshift\(.*\)'), parameter):   
                index           = int(parameter.split('(')[1].split(')')[0])
                zmean[index-1]  = float(line[eqpos+1:].split()[0]) 
            # Get redshift widths:
            if re.match(re.compile('redshift_sigma\(.*\)'), parameter):   
                index           = int(parameter.split('(')[1].split(')')[0])
                zwidth[index-1] = float(line[eqpos+1:].split()[0]) 
            # Get bias:
            if re.match(re.compile('redshift_bias\(.*\)'), parameter):   
                index           = int(parameter.split('(')[1].split(')')[0])
                bias[index-1]   = float(line[eqpos+1:].split()[0])
            # Get field type:
            if re.match(re.compile('redshift_kind\(.*\)'), parameter):   
                index           = int(parameter.split('(')[1].split(')')[0])
                if line[eqpos+1:].split()[0] == 'counts':
                    ftype[index-1] = GalType
                if line[eqpos+1:].split()[0] == 'lensing':
                    ftype[index-1] = KappaType
# End of scan over CAMB sources .ini file.


# Assign new index to new value in list or the same previously assigned index if value is repeated:
# (indexes must be previously allocated to zero vector)
def Values2Index(values, indexes):
    # Recursion initial condition:
    lasti = values.size-1
    if lasti==0: indexes[lasti]=1
    # Recursion:
    else:
        Values2Index(values[:lasti], indexes[:lasti])
        # If new value is repeated, use previously assigned index:
        for i in range(lasti): 
            if values[i]==values[lasti]: indexes[lasti] = indexes[i]
        # Else use new index:
        if indexes[lasti]==0: indexes[lasti] = np.amax(indexes[:lasti])+1
    return

# Define indexes for redshift bins and fields:
Values2Index(zmean, zindex)
Values2Index(bias,  findex)
        

# Convergence kappa shift formula from Hilbert, Hartlap & Schneider (2011)
def HilbertShift(z):
    return 0.008*z + 0.029*(z**2) - 0.0079*(z**3) + 0.00065*(z**4) 

def XavierShift(z):
    a0 = 0.2;
    s0 = 0.568591;
    return a0*(((z*s0)**2 + z*s0 + 1)/(z*s0 + 1) - 1)

# Functions for output:
def mean(ftype):
    if ftype==GalType:   return GalMean
    if ftype==KappaType: return KappaMean

# Functions for output:
def shift(ftype, z):
    if ftype==GalType:  return GalShift
    if ftype==KappaType:
        if FixKappa==1: return KappaShift
        if FixKappa==0:
                #return np.interp(z, fileZ, fileSh)
                #return HilbertShift(z)
                return XavierShift(z)


# Output:
fout = open(fieldinfo, 'w')

print >>fout, "# Field number, z bin number, mean, shift, field type, zmin, zmax"
print >>fout, "# Types: 1-galaxies 2-shear\n"

for i in range(Nwindows):
        print >>fout,'{0:6d} {1:6d}   {2:.4f}   {3:.4f} {4:6d}   {5:.4f}   {6:.4f}'.format(int(findex[i]), int(zindex[i]), mean(ftype[i]), shift(ftype[i], zmean[i]), int(ftype[i]), zmean[i]-zwidth[i], zmean[i]+zwidth[i])

