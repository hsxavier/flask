#! /usr/bin/env python

# This script generates a Healpix map of size <NSIDE> with 
# all values set to one and exports to <OUTFILE>.
# USAGE: ./GenUnitMap.py <NSIDE> <OUTFILE> 

import numpy as np
import healpy as hp
import sys

# Get input:
Nside   = int(sys.argv[1])
outfile =     sys.argv[2]

# Generate map:
Npixels = 12*Nside*Nside
m = np.ones(Npixels)

# Wrtie to file:
hp.write_map(outfile, m)
