#! /usr/bin/env python

"""
USAGE:   ViewMap.py <HEALPIX_MAP>
EXAMPLE: ViewMap.py euclid_footprint.fits

This script plots in Mollweide projection a Healpix map saved as a fits file. 

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 23/apr/2015.
"""

import healpy as hp
import matplotlib.pyplot as plt
import sys

# Docstring output:
if len(sys.argv) != 1 + 1: 
    print(__doc__)
    sys.exit(1)


# Get input:
infile = sys.argv[1]
m      = hp.read_map(infile)

hp.mollview(m)
plt.show()
