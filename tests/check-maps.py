#!/usr/bin/env python

import sys
import numpy as np
import healpy as hp
import matplotlib.pyplot as pl

# Get input files:
InfoFile  = sys.argv[1]
Nside     = int(sys.argv[2])
MaskFile  = sys.argv[3]
nzPrefix  = sys.argv[4]
SelScale  = float(sys.argv[5])
StarFile  = sys.argv[6]
mapPrefix = sys.argv[7]
separable = int(sys.argv[8])
Npix = 12*Nside**2



#############################
### Loads all information ###
#############################

# Load fields' informations:
fname, zname, mean, shift, ftype, zmin, zmax = np.loadtxt(InfoFile, unpack=True)
OnlyLens = np.all(ftype==2)

# Loads mask:
def GetMask(File,Type):
    if Type==1 or OnlyLens==True:
        return hp.read_map(File, verbose=False)
    else:
        return np.ones(Npix)
if separable==0:
    mask = [GetMask(MaskFile+'f'+str(int(fname[i]))+'z'+str(int(zname[i]))+'.fits', ftype[i]) for i in range(0,len(fname))]
if separable==1:
    mask = [hp.read_map(MaskFile, verbose=False) for i in range(0,len(fname))]
if separable==2:
    mask = [GetMask(MaskFile+'f'+str(int(fname[i]))+'.fits', ftype[i]) for i in range(0,len(fname))]

# Loads radial selection function and integrate dn/dz over redshift intervals:
def GetNz(File,Type):
   if Type==1:
       return np.loadtxt(File)
   else:
       return -1
if separable!=0:
    nzRaw = [GetNz(nzPrefix+'f'+str(int(fname[i]))+'.dat',ftype[i]) for i in range(0,len(fname))]
    # Integration:
    nsteps=20
    def nzInt(i):
        if ftype[i]==1:
            y = np.trapz(np.interp(np.linspace(zmin[i],zmax[i],nsteps),nzRaw[i][:,0],nzRaw[i][:,1]),x=np.linspace(zmin[i],zmax[i],nsteps)) 
            return y*4*np.pi*(180*60.0/np.pi)**2/Npix
        else:
            return -1
    GalPerPix = [nzInt(i) for i in range(0,len(fname))]

# Loads star mask:
starmask = hp.read_map(StarFile,verbose=False)

# Loads maps:
Map = [hp.read_map(mapPrefix+'f'+str(int(fname[i]))+'z'+str(int(zname[i]))+'.fits', verbose=False) for i in range(0,len(fname))]



####################
### Do the tests ###
####################

fullmask = [m*starmask for m in mask]

# Test that maps as properly masked:
for i in range(0, len(fname)):
    if (ftype[i]==1 or separable==1):
        if (np.all(Map[i][fullmask[i]==0]==0) or np.all(Map[i][fullmask[i]==0]==hp.UNSEEN))==False:
            print '!! ERROR: map f{}z{} has regions that should be masked and are not !!'.format(fname[i],zname[i])

# Test that galaxies are sampled properly:
for i in range(0, len(fname)):
    if (ftype[i]==1):
        residual = Map[i] - SelScale*GalPerPix[i]*fullmask[i]
        #hp.mollview(residual)
        #pl.show()
        sample = residual[fullmask[i]>0]
        if np.abs(np.mean(sample)/(np.std(sample)/np.sqrt(len(sample)-1)))>7:
            print 'Found excess of difference in the mean: f{}z{}'.format(int(fname[i]),int(zname[i]))
            print 'Observed mean dev:  {:.3f}'.format(np.std(sample)/np.sqrt(len(sample)-1))
            print 'Expected density:  {:.3f}'.format(SelScale*GalPerPix[i]*np.mean(fullmask[i][fullmask[i]>0]))
            print 'Observed density:  {:.3f}'.format(np.mean(Map[i][fullmask[i]>0]))
            print 'Num. of sigma:     {:.3f}'.format(np.mean(sample)/(np.std(sample)/np.sqrt(len(sample)-1)))
            pl.hist(sample, bins=np.arange(-np.floor(np.max(SelScale*GalPerPix[i]*fullmask[i])), \
                                             np.ceil(np.max(Map[i]))))
            pl.show()

