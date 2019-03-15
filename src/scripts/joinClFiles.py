#! /usr/bin/env python

"""
This script takes a bunch of individual Cl files <CLS> (with columns L and Cl, 
in this order) and concatenates all of them into a single Cl table <OUTFILE>. 
Note that it is required the Cls to have the same Ls. The header will be based 
on the Cl filenames.
 

USAGE:   joinClFiles.py <CLS> <OUTFILE>
EXAMPLE: joinClFiles.py sims/Cl05/fidCl-f1z* test.dat

Written by Henrique S. Xavier, hsxavier@if.usp.br, 14/mar/2019.
"""

import sys
import numpy as np
import re

# Docstring output:
if len(sys.argv) < 1 + 2: 
    print(__doc__)
    sys.exit(1)

# Get input:
ClFiles = sys.argv[1:-1]
OutFile = sys.argv[-1]

# Transforms filenames into labels:
def File2Label(f):
    return f.split('/')[-1].split('.dat')[0]
Labels = ['l']+[File2Label(f) for f in ClFiles]

# Loads Cl files:
Cls = [np.loadtxt(f, unpack=True) for f in ClFiles]
if np.all([np.all(Cls[0][0]==cl[0]) for cl in Cls])==False:
    print "!! ERROR: All Cl files must cover the same Ls. !!"
    sys.exit()

# Creates a single table with all Cls:
Table = np.transpose(np.insert(np.array(Cls)[:,1,:],0,Cls[0][0], axis=0))

# Export:
Header = ' '.join(Labels)
fmt    = ['%d']+['%.12e']*len(ClFiles)
np.savetxt(OutFile,Table,fmt=fmt,header=Header)
