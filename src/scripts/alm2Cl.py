#! /usr/bin/env python

"""
USAGE:   alm2Cl.py <INPUT_ALM_FILE>
EXAMPLE: alm2CL.py test_alm.dat
OUTPUT:  l, Cl(1,1), Cl(1,2), ... (on screen)

This script takes a FLASK alm output table (ASCII .dat file) 
and computes the Cls from it by computing the covariance of the 
alms. It has many important hard-coded parameters (see code). 
It writes the output to stdout. This computation is now made 
by FLASK itself, so this script is obsolete.

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 06/may/2015.
"""

import numpy as np
import sys

# Get input parameters:
inputfile = sys.argv[1]


field   = np.loadtxt(inputfile, unpack=True)
fRe     = field[::2]
fIm     = field[1::2]
Nfields = fRe.shape[0]

# Hard-coded values:
nrealiz = 5000*10    # Number of realizations of each l (note that each different 
                     # m count as one realization and all ls must have the same number of 
                     # realizations.
lmin    = 10         # minimum l for which the alms were computed.
lmax    = 109        # maximum l for which the alms were computed.



# Header
print "# l",
for f1 in range(1, Nfields+1):
    for f2 in range(f1, Nfields+1):
        print "f"+`f1`+"-f"+`f2`,
print ""
print ""

for l in range(0, lmax-lmin+1):

    covmatRe = np.cov(np.transpose(fRe)[nrealiz*l:nrealiz*(l+1)-1], rowvar=0)
    covmatIm = np.cov(np.transpose(fIm)[nrealiz*l:nrealiz*(l+1)-1], rowvar=0)
    covmat   = covmatRe + covmatIm 
    
    print l+lmin,
    for f1 in range(0, Nfields):
        for f2 in range(f1, Nfields):
            print covmat[f1,f2],
    print ""
