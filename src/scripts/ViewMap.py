#! /usr/bin/env python

"""
This script plots in Mollweide projection the Healpix maps saved as fits files 
and passed to it as arguments. 

USAGE:   ViewMap.py <HEALPIX_MAP_1> <HEALPIX_MAP_2> ...
EXAMPLE: ViewMap.py euclid_footprint.fits

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 23/apr/2015.
"""

import healpy as hp
import matplotlib.pyplot as plt
import sys

# Docstring output:
if len(sys.argv) < 1 + 1: 
    print(__doc__)
    sys.exit(1)


# Get input:

for infile in sys.argv[1:]:
    m = hp.read_map(infile)
    hp.mollview(m, title=infile)

plt.show()
