#! /usr/bin/env python

"""
USAGE:   ini2info.dat <CLASS_INPUT> <FIELDS_INFO_FILE>
EXAMPLE: ini2info.dat prec17_20_dens.ini prec17_20_dens/fields-info.dat

NOTE: Current version expects that input Cls are a single density field Cls.
      This is relevant because density Cl(z,z') are sharply peaked at z=z', 
      which we use to extrapolate the Cls. Furthermore, we assume that z=z' 
      happens at the covariance matrix diagonal only, so only one density 
      field can appear and the redshifts must be ordered.
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
Nfields = len(fName)
if len(set(fName)) != 1: 
    print "ERROR! Current version expects only one field."
    sys.exit(1)
if set(fType) != set([1]):
    print "ERROR! Current version expects only density fields (fType=1)." 
    sys.exit(1)
# !!! Check if redshifts are correctly organized.

# Create matrix of Cl filenames:
fString      = np.char.array(map(makeFname, fName, zName))
ClFileMatrix = clprefix + (fString[:,None]+fString) + '.dat'
ClFileList   = ClFileMatrix.reshape(Nfields*Nfields) 
ClFileTransp = np.transpose(ClFileMatrix).reshape(Nfields*Nfields)

# Load Cls in covariance matrices for each ell:
CovByl = np.transpose(np.array(map(loadCl, ClFileList, ClFileTransp)))
Nl     = CovByl.shape[0]
CovByl = CovByl.reshape((Nl,Nfields,Nfields))

print np.diagonal(CovByl[0], 4)


#print ClFileList
#print ClFileTransp

#isfile2 = np.apply_over_axes(os.path.isfile, ClFileMatrix, )

#print isfile2(ClFileMatrix)
#map(loadCl, ClFileMatrix)
