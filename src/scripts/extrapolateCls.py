#! /usr/bin/env python

"""
USAGE:   extrapolateCls.dat <FIELDS_INFO_FILE> <CL_IN_PREFIX>
EXAMPLE: extrapolateCls.dat dens-info.dat ClIn/Cl-

This script extrapolates the Cls for one extra redshift bin that goes from 
z=0 to the zmin of the lowest redshift in the input. 

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 20/nov/2015.
"""

import numpy as np
import os.path
import sys


# Loads a two column file and only returns the second column:
def loadCl(file1, file2):
    if os.path.isfile(file1):
        #print "Will load "+file1
        l, cl = np.loadtxt(file1, unpack=True)
    elif os.path.isfile(file2):
        #print "Will load "+file2
        l, cl = np.loadtxt(file2, unpack=True)
    else:
        print "ERROR! Could not find file "+file1+" nor "+file2
        sys.exit(1)
    return cl


# Creates string that specifies the Field (field and z bin) given by f_z_:
def makeFname(fName, zName):
    return 'f'+str(int(fName))+'z'+str(int(zName))


##################################################################################
#                                 Main code                                      #
##################################################################################
 
# Get input:
infofile = sys.argv[1]
clprefix = sys.argv[2]


# Load fields info file:
fName, zName, mean, shift, fType, zMin, zMax = np.loadtxt(infofile, unpack=True)
zMean   = (zMax+zMin) / 2.0  # Bins mean redshift
z0      = zMin[0] / 2.0      # New bin mean redshift
z0hw    = zMin[0] / 2.0      # New bin half-width
Nfields = len(fName)

if len(set(fName)) != 1: 
    print "ERROR! Current version expects only one field."
    sys.exit(1)
if set(fType) != set([1]):
    print "ERROR! Current version expects only density fields (fType=1)." 
    sys.exit(1)
if not np.array_equal(zMin, np.sort(zMin)):
    print "ERROR! Current version needs z bins (z min) in increasing order." 
    sys.exit(1) 
if not np.array_equal(zMax, np.sort(zMax)):
    print "ERROR! Current version needs z bins (z max) in increasing order." 
    sys.exit(1)


# Create matrix of Cl filenames:
fString      = np.char.array(map(makeFname, fName, zName))
ClFileMatrix = clprefix + (fString[:,None]+fString) + '.dat'
ClFileList   = ClFileMatrix.reshape(Nfields*Nfields) 
ClFileTransp = np.transpose(ClFileMatrix).reshape(Nfields*Nfields)

# Load Cls in covariance matrices for each ell:
print "Loading Cls..."
CovByl = np.transpose(np.array(map(loadCl, ClFileList, ClFileTransp)))
Nl     = CovByl.shape[0]
CovByl = CovByl.reshape((Nl,Nfields,Nfields))

# Preparing to extrapolate:
print "Extrapolating Cls..."
Cl0    = np.zeros((Nfields+1, Nl))
dzFrac = (z0-zMean[0])/(zMean[0]-zMean[1])
# Extrapolate each l:
for l in range(0, Nl):
    # Extrapolate main diagonal (auto-Cl is sharply peaked inside the bin, therefore the extra factor):
    Cl0[0,l] = ((CovByl[l,0,0]-CovByl[l,1,1])*dzFrac + CovByl[l,0,0]) * (zMax[0]-zMin[0])/(2.0*z0hw) 

    # Extrapolate remaining diagonals but last two (which have only one term):
    for z in range(1, Nfields-1):
        Cl0[z,l] = (CovByl[l,0,z]-CovByl[l,1,1+z])*dzFrac + CovByl[l,0,z]

    # Extrapolate last terms using the column and not the diagonal:
    Cl0[Nfields-1,l] = (CovByl[l,0,Nfields-2]-CovByl[l,1,Nfields-2])*dzFrac + CovByl[l,0,Nfields-2]
    Cl0[Nfields  ,l] = (CovByl[l,0,Nfields-1]-CovByl[l,1,Nfields-1])*dzFrac + CovByl[l,0,Nfields-1]

# Export Cls:
print "Exporting files:"
lList    = np.loadtxt(ClFileList[0], usecols=(0,))
filename = clprefix+"f"+str(int(fName[0]))+"z0f"+str(int(fName[0]))+"z0.dat"
print filename
np.savetxt(filename,np.transpose([lList,Cl0[0]]))
for z in range (1, Nfields+1):
    filename = clprefix+"f"+str(int(fName[0]))+"z0f"+str(int(fName[z-1]))+"z"+str(int(zName[z-1]))+".dat"
    print filename
    np.savetxt(filename,np.transpose([lList,Cl0[z]]))

