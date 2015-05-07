#! /usr/bin/env python

import numpy as np
import sys

inputfile = sys.argv[1]

field   = np.loadtxt(inputfile, unpack=True)
fRe     = field[::2]
fIm     = field[1::2]
Nfields = fRe.shape[0]

nrealiz = 5000*10
lmin    = 10
lmax    = 109



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
