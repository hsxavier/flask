#! /usr/bin/env python

import sys
import numpy as np

# Get input:
Clfile     = sys.argv[1]
factor     = sys.argv[2]
outfile    = sys.argv[3]

# Load Cl file:
l, Cl = np.loadtxt(Clfile, unpack=True)

Cl = Cl * float(factor)

np.savetxt(outfile, np.array([l,Cl]).transpose())
