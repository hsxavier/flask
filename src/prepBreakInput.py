#! /usr/bin/env python

"""
USAGE:   prepBreakInput.py <class_files> <prefix>
EXAMPLE: prepBreakInput.py break-*_cl.dat Cl-

This script takes a bunch of CLASS Cl output files <class_files> which must have exactly 
two redshift bins and must be called <something>-<zbin1>-<zbin2>_cl.dat and creates 
a separate file for each Cl in each file in <class_files> named 
<prefix>f<type_i>z<zbin_i>f<type_j>z<zbin_j>.dat, where <zbin_i> and <zbin_j> can 
only be either <zbin1> or <zbin2>. Certain factors are applied to the Cls such 
that they represent the pure Cls [without (2l+1) and 4*pi factors] of density and 
convergence (instead of lensing potentials).

NOTE: This code does not overwrite previous files with the same output name.

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 26/oct/2015.
"""

import numpy as np
import sys
import math as m
import os.path

# Get script parameters (list of input files and output prefix):
classlist = sys.argv[1:-1]
outprefix = sys.argv[-1]

# LOOP over input files:
for classfile in classlist:

    # Get bin numbers from filename:
    bins=classfile.split("_")[-2].split("-")[-2:]

    # Get file header:
    fp = open(classfile, "r")
    for i, line in enumerate(fp):
        if i == 6:
            header0 = line
        elif i > 6:
            break
    fp.close()
    header = header0.split()
    header.pop(0)

    # Find out if there is more than one field type:
    if "dens[" in header0 and "lens[" in header0:
        NfieldTypes=2
    else:
        NfieldTypes=1
    # Generate output filenames:
    outfiles=[]
    for name in header:
        ini=name.find(":")+1
        name = name[ini:]
        # Find out if it is dens or lens type Cl:
        if NfieldTypes==1:
            # If there is just one type, call the field '1':
            name = name.replace("dens[","f1z")
            name = name.replace("lens[","f1z")
        else:
            # If there are two types, dens is '1' and lens is '2':
            name = name.replace("dens[","f1z")
            name = name.replace("lens[","f2z")
        # Give z bin name according to filename:
        name = name.replace("1]-",bins[0])
        name = name.replace("1]" ,bins[0])
        name = name.replace("2]-",bins[1])
        name = name.replace("2]" ,bins[1])
        name = outprefix+name+".dat"
        outfiles.append(name)
   
            
    # Load data:
    classout = np.loadtxt(classfile, unpack=True)
    ncols    = classout.shape[0]
    length   = classout.shape[1]

    # First column must be ell:
    l = classout[0]

    # Precompute basic factor:
    fac0 = (2.0*m.pi)/l/(l+1)
    fac1 = l*(l+1)/2.0   # The sign was supposed to be negative according to Hu 2000 (PRD 62:043007), 
    # but that would lead to negative kappa-delta cross-Cls, which according to 
    # numerical integration of delta auto-Cls is wrong. This is probably caused 
    # by differences in definitions.

    # LOOP over columns 
    for i in range(1, ncols):

        # Discover if column is gal-gal, gal-lens ou lens-lens and compute appropriate factors:
        factor = fac0
        first = header[i].find("lens")
        if first != -1:
            # One lens so far, use fac1
            factor = factor*fac1
            second = header[i].find("lens", first+1)
            if second !=-1:
                # Two lens, use another fac1 
                factor = factor*fac1

        # Multiply column by factor
        classout[i] = classout[i]*factor

        # Export files that do not exist:
        if os.path.isfile(outfiles[i]) is False:
            print "Writing file "+outfiles[i]
            np.savetxt(outfiles[i], np.transpose([l,classout[i]]))
