#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as pl

# Get input files:
dist     = sys.argv[1]
InfoFile = sys.argv[2] 
mapFile  = sys.argv[3]

# Load fields' informations:
fname, zname, mean, shift, ftype, zmin, zmax = np.loadtxt(InfoFile, unpack=True)
# Loads map file (first two columns are angular coordinates):
Map = np.loadtxt(mapFile, unpack=True)

for i in range(0, len(fname)):
    obsMean = np.mean(Map[2+i])
    obsDev  = np.std(Map[2+i])/np.sqrt(len(Map[2+i]))
    obsMin  = np.min(Map[2+i])
    if dist!='HOMOGENEOUS': 
        nSigma  = (obsMean-mean[i])/obsDev
        if np.abs(nSigma)>3:
            name = 'f{}z{}:'.format(int(fname[i]),int(zname[i]))
            print name
            print 'Deviation from expected mean:', nSigma
            pl.hist(Map[2+i]-mean[i], bins=100, label=name)
            pl.legend()
            pl.show()
    else:
        if np.all(Map[2+i]==mean[i])==False:
            print '!! Field is not homogeneous !!'
            print 'Mean Obs:', obsMean, 'Mean Exp:', mean[i]
    if dist=='LOGNORMAL' and obsMin<-shift[i]:
        print '!! Minimum value beyond expected !!'
        print 'Min Obs:', obsMin, 'Min Exp:', -shift[i]
