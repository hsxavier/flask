#! /usr/bin/env python

import healpy as hp
import matplotlib.pyplot as plt
import sys

# Get input:
infile = sys.argv[1]
m      = hp.read_map(infile)

hp.mollview(m)
plt.show()
