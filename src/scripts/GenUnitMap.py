#! /usr/bin/env python

"""
USAGE:   GenUnitMap.py <NSIDE> <OUTFILE> 
EXAMPLE: GenUnitMap.py 1024 full-sky_1024.fits

This script generates a Healpix map of size <NSIDE> with 
all values set to one and exports to <OUTFILE>.
It was used when FLASK did not accept no angular selection function.
Now this script is not that useful.

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 23/apr/2015.
"""

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
