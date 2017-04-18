#! /usr/bin/env python

"""
This script takes a bunch of CAMB Cl output 'scalCov' files <camb_files> 
which must have exactly two redshift bins and must be called 
<something>f<fname1>z<zbin1>f<fname2>z<zbin2>_scalCovCls.dat and creates 
a separate file for each Cl in each file in <camb_files> named 
<prefix>f<fname_n>z<zbin_i>f<fname_m>z<zbin_j>.dat, where <zbin_i> and 
<zbin_j> can only be either <zbin1> or <zbin2> and <fname_n> and <fname_m> 
can only be either <fname1> or <fname2>. Certain factors are applied to 
the Cls such that they represent the pure Cls [without (2l+1) and 4*pi 
factors] of density and convergence.

USAGE:   prepCambBreak.py <class_files> <prefix>
EXAMPLE: prepCambBreak.py break-*_scalCovCls.dat Cl-

NOTE: This code does not overwrite previous files with the same output name.

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 18/apr/2017.
"""

import numpy as np
import sys
import math as m
import os.path
import re

# Docstring output:
if len(sys.argv) < 1 + 2: 
    print(__doc__)
    sys.exit(1)



#####################
### Prelimiraries ###
#####################

# Hard-coded parameters:
Nfields  = 2
CovNCols = 26

# Function to get position of Cl in CAMBsources Cov file:
def CovPosition(f1, f2, Nfields):
    return (Nfields+3)*(3+f1-1) + (3+f2)



#################
### Main code ###
#################

# Get input (list of input files and output prefix):
camblist  = sys.argv[1:-1]
outprefix = sys.argv[-1]

# LOOP over input files:
for cambfile in camblist:

    # Get Field names:
    Fname  = re.findall('f[0-9]*z[0-9]*', cambfile)

    # Load Cls:
    camb = np.loadtxt(cambfile, unpack=True)
    if len(camb) != 26:
        print "ERROR! Current version only accepts files with two fields."
        sys.exit(1)
    l    = camb[0]
    fac  = (2.0*np.pi)/l/(l+1)
    # LOOP over Cls in one file:
    for f1 in range(1, Nfields+1):
        for f2 in range(f1, Nfields+1):
            outfile = outprefix+Fname[f1-1]+Fname[f2-1]+'.dat'
            if os.path.isfile(outfile) is False:
                Cl = np.transpose( [l, fac*camb[CovPosition(f1,f2,Nfields)]] )
                print 'Writing file '+outfile
                np.savetxt(outfile, Cl, fmt=['%d','%e'])
print 'Done.'
