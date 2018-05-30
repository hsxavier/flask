#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as pl

# Get input files:
headerFile = sys.argv[1]
meanFile   = sys.argv[2]
devFile    = sys.argv[3]
Nsample    = int(sys.argv[4])
theoFiles  = sys.argv[5:]


# Function for getting field and redshift numbers:
def getFZ(clFile):
    fileLen = len(clFile)
    reverse = clFile[::-1]
    z2end = fileLen - 4
    z2ini = fileLen - reverse.find('z')
    z2 = int(clFile[z2ini:z2end])
    f2end = z2ini - 1
    f2ini = fileLen - reverse.find('f',fileLen-f2end)
    f2 = int(clFile[f2ini:f2end])
    z1end = f2ini - 1
    z1ini = fileLen - reverse.find('z',fileLen-z1end)
    z1 = int(clFile[z1ini:z1end])
    f1end = z1ini - 1
    f1ini = fileLen - reverse.find('f',fileLen-f1end)
    f1 = int(clFile[f1ini:f1end])
    return [f1,z1,f2,z2]

# Sort theoretical Cl files according to summarizeData.py output order:
#FZs = [ getFZ(theoFile) for theoFile in theoFiles ]
#[FZs[i].append(theoFiles[i]) for i in range(0,len(theoFiles))]
#sortedTheoFiles = np.array(sorted(FZs, key = lambda x: (x[0], x[1], x[2], x[3])))[:,4]
sortedTheoFiles = np.array(theoFiles)

# Sort header of recovered Cls:
sortedHeader = [ a[0] for a in sorted(zip( range(0,len(theoFiles)), np.loadtxt(headerFile, dtype=str)[1:] ), key = lambda x: x[1]) ]

# Load theoretical Cls:
theoCls = []
theols  = []
for theoFile in sortedTheoFiles:
    l, cl = np.loadtxt(theoFile, unpack=True)
    theoCls.append(cl)
    if len(theols)==0:
        theols=l
    else:
        if np.all(theols==l)==False:
            print "ERROR: Theoretical Cls do not have all the same ells"
            sys.exit(1)

# Load mean recovered Cls and standard deviation:
meanClsRaw = np.loadtxt(meanFile, unpack=True)
devClsRaw  = np.loadtxt(devFile,  unpack=True)

# Find the same ell range for theoretical predictions and mean recovery: 
minEll  = np.max( [np.min(meanClsRaw[0]), np.min(theols)] )
maxEll  = np.min( [np.max(meanClsRaw[0]), np.max(theols)] )
meanIni = np.where(meanClsRaw[0]==minEll)[0][0]
meanFin = np.where(meanClsRaw[0]==maxEll)[0][0]
theoIni = np.where(theols==minEll)[0][0]
theoFin = np.where(theols==maxEll)[0][0]

# Separate ls and Cls and sort Cls:
meanl   = meanClsRaw[0,meanIni:meanFin+1] 
meanCls = (meanClsRaw[1:,meanIni:meanFin+1])[sortedHeader]
devCls  = (devClsRaw[1:,meanIni:meanFin+1])[sortedHeader]
theoCls = np.array(theoCls)[:,theoIni:theoFin+1]

# Compute deviations of mean recovered Cls from original Cls, normalized by std. dev of the mean :
nSigmaCls = (meanCls - theoCls) / devCls * np.sqrt(Nsample)

YesPlot = False
for i in range(0, len(theoCls)):
    nSfinal = np.mean(nSigmaCls[i])/np.std(nSigmaCls[i])*np.sqrt(len(nSigmaCls[i])-1)
    if np.abs(nSfinal)>5:
        name = sortedTheoFiles[i]
        print name, '-- Avg. dev. (sigma units):', nSfinal
        YesPlot = True
        pl.errorbar(meanl, nSigmaCls[i], yerr=1, fmt='o', label=name)
if YesPlot==True:
    pl.legend(loc='lower left')
    pl.show()

# Plots of the mean of recovered Cls and original Cls (regularized ones):
if 1==0:
    for i in range(0, len(meanCls)):
        pl.plot(meanCls[i],'r+')
        pl.plot(theoCls[i],'b-')
        pl.xscale('log')
        pl.yscale('log')
        pl.show()

# Plots of difference between mean of recovered Cls and original Cls (regularized ones):
if 1==0:
    for i in range(0,5):
        print sortedTheoFiles[i]
        pl.errorbar(meanCls[0,meanIni:meanFin+1], (meanCls[1:,meanIni:meanFin+1]-np.array(theoCls)[: ,theoIni:theoFin+1])[i]/((devCls[1:,meanIni:meanFin+1])[i] / np.sqrt(Nsample)), yerr=1, fmt='.', label=sortedTheoFiles[i])
    for i in range(5,10):
        print sortedTheoFiles[i]
        pl.errorbar(meanCls[0,meanIni:meanFin+1], (meanCls[1:,meanIni:meanFin+1]-np.array(theoCls)[: ,theoIni:theoFin+1])[i]/((devCls[1:,meanIni:meanFin+1])[i] / np.sqrt(Nsample)), yerr=1, fmt='o', label=sortedTheoFiles[i])
    pl.plot(meanCls[0,meanIni:meanFin+1], np.zeros(len(meanCls[0,meanIni:meanFin+1])))
    pl.legend(loc='lower left')
    pl.show()

# Plot histogram of deviations:
if 1==0:
    print np.shape(nSigmaCls)
    pl.hist(nSigmaCls.flatten(),bins=50)
    pl.show()


