#! /usr/bin/env python

"""
USAGE: ./ChangeMapResolution.py <IN> <NSIDE> <OUT>
EXAMPLE: ./ChangeMapResolution.py euclid_footprint_4096.fits 1024 euclid_footprint_1024.fits
 
This script gets a Healpix map in FITS file <IN>, changes its 
resolution (Nside) according to <NSIDE> and write the new map 
to <OUT>.

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 23/apr/2015.
"""

import healpy as hp
import sys

# Docstring output:
if len(sys.argv) != 1 + 3: 
    print(__doc__)
    sys.exit(1)

# Get input:
mapIn    =     sys.argv[1]
newNside = int(sys.argv[2])
mapOut   =     sys.argv[3]

# Operate on map:
mapx = hp.read_map(mapIn)
mapy = hp.ud_grade(mapx, nside_out=newNside)
print "NEW NSIDE =", newNside

# Write it out:
hp.write_map(mapOut, mapy)
print "New map written to", mapOut


