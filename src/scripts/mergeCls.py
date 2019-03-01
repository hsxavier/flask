#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script takes several Cl <CL_FILES> that have first column L and second column 
Cl and place them in a single file <OUTFILE> with a common L column and each Cl as 
a different column.

USAGE:   mergeCls.py <CL_FILES> <OUTFILE>
EXAMPLE: mergeCls.py ../../data/exampleCl-f* mergedCls.dat

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 01/mar/2019
"""

import numpy as np
import sys

# Docstring output:
if len(sys.argv) < 1 + 2: 
    print(__doc__)
    sys.exit(1)

# Get input:
ClFiles = sys.argv[1:-1]
OutFile = sys.argv[-1]

# Load data:
lCl = [np.loadtxt(f, unpack=True) for f in ClFiles]
if np.all([np.all(t[0]==lCl[0][0]) for t in lCl])==False:
    print "All input Cls must cover the same L range."
    sys.exit(1)

# Create new table:
mergedCls = np.transpose(np.insert(np.array(lCl)[:,1,:], 0, lCl[0][0], axis=0))
ClNames   = [f.split('/')[-1] for f in ClFiles]
header    = 'l '+' '.join(ClNames)

# Save to file:
fmtArr = ['%.18e']*len(ClNames)
fmtArr.insert(0,'%d')
np.savetxt(OutFile, mergedCls, fmt=fmtArr, header=header)
