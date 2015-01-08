#! /usr/bin/env python

import numpy as np
import sys

inputfile = sys.argv[1]

f1Re, f1Im, f2Re, f2Im = np.loadtxt(inputfile, unpack=True)
nrealiz=10000*11
#nrealiz=3*11
lmin=10
lmax=50

#print("l f1f1 f2f2 f1f2")
for l in range(0, lmax-lmin+1):

    #covmat = np.cov(f1Re[nrealiz*l*(l+1)/2:nrealiz*(l+1)*(l+2)/2-1],f2Re[nrealiz*l*(l+1)/2:nrealiz*(l+1)*(l+2)/2-1]) + np.cov(f1Im[nrealiz*l*(l+1)/2:nrealiz*(l+1)*(l+2)/2-1],f2Im[nrealiz*l*(l+1)/2:nrealiz*(l+1)*(l+2)/2-1])
    covmat = np.cov(f1Re[nrealiz*l:nrealiz*(l+1)-1],f2Re[nrealiz*l:nrealiz*(l+1)-1]) + np.cov(f1Im[nrealiz*l:nrealiz*(l+1)-1],f2Im[nrealiz*l:nrealiz*(l+1)-1])
    print l+lmin, covmat[0,0], covmat[1,1], covmat[0,1]
