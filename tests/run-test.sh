#!/usr/bin/env bash

zerowarnings="Zero warnings.................................... OK"

# Run test 01:
testN=$1
echo "** TEST $testN **"

# Check Cl processing and output:
echo     "Outputting Cl preparation files..."
../bin/flask test${testN}.config FLIST_OUT: tout${testN}_flist.txt SMOOTH_CL_PREFIX: tout${testN}_sCl- XIOUT_PREFIX: tout${testN}_xi- GXIOUT_PREFIX: tout${testN}_gxi- GCLOUT_PREFIX: tout${testN}_gCl- EXIT_AT: GCLOUT_PREFIX 1>tout${testN}.log 2>>tout${testN}.log

warnings=`grep "Total number of warnings:" tout${testN}.log | cut -d: -f2`
if [ $warnings -eq 0 ]; then
    echo $zerowarnings
else
    echo "!! Run has $warnings warnings or crashed !!"
    echo "Aborting test..."
    exit 1
fi


# Check map distribution:
echo "Checking full-sky noisless map distribution..."
declare -a arr0=(LOGNORMAL GAUSSIAN HOMOGENEOUS)
for dist in "${arr0[@]}"; do
    echo $dist
    ../bin/flask test${testN}.config DIST: $dist EXIT_AT: MAP_OUT MAP_OUT: tout${testN}_map.dat 1>tout${testN}.log 2>>tout${testN}.log
    warnings=`grep "Total number of warnings:" tout${testN}.log | cut -d: -f2`
    if [ $warnings -eq 0 ]; then
	echo $zerowarnings
    else
	echo "!! Run has $warnings warnings or crashed!!"
	echo "Aborting test..."
	exit 1
    fi
    info=`grep "FIELDS_INFO:" test${testN}.config | awk '{print $2}'`
    ./check-dist.py $dist $info tout${testN}_map.dat
done

# Checking recovery of input Cls:
declare -a arr1=(LOGNORMAL GAUSSIAN)
nSamples=50
for dist in "${arr1[@]}"; do
    echo $dist
    echo "Creating regularized Cls..."
    ../bin/flask test${testN}.config DIST: $dist EXIT_AT: REG_CL_PREFIX REG_CL_PREFIX: tout${testN}_regCl- 1>tout${testN}.log 2>>tout${testN}.log
    warnings=`grep "Total number of warnings:" tout${testN}.log | cut -d: -f2`
    if [ $warnings -eq 0 ]; then
	echo $zerowarnings
    else
	echo "!! Run has $warnings warnings or crashed!!"
	echo "Aborting test..."
	exit 1
    fi
    echo     "Checking recovery from full-sky noiseless maps..."
    ./many_recov-cls.sh ${testN} $dist 1 $nSamples
    ../src/scripts/summarizeData.py tout${testN}_mean-recCl.dat tout${testN}_dev-recCl.dat tout${testN}_recCl-*.dat
    head -n 1 tout${testN}_recCl-0001.dat | sed 's/#//g' > tout${testN}_header-recCl.dat
    rm -f tout${testN}_recCl-*
    ./compare-cls.py tout${testN}_header-recCl.dat tout${testN}_mean-recCl.dat tout${testN}_dev-recCl.dat $nSamples tout${testN}_regCl-f*.dat
done    

# Check selection function operations:
echo "Checking maps sampling (including selection functions)..."
declare -a arr2=(0 1)
for dist in "${arr1[@]}"; do
    echo $dist
    for poisson in "${arr2[@]}"; do
	echo "POISSON: $poisson"
	../bin/flask test${testN}.config DIST: $dist POISSON: $poisson EXIT_AT: MAPWERFITS_PREFIX MAPWERFITS_PREFIX: tout${testN}_map- 1>tout${testN}.log 2>>tout${testN}.log
	warnings=`grep "Total number of warnings:" tout${testN}.log | cut -d: -f2`
	if [ $warnings -eq 0 ]; then
	    echo $zerowarnings
	else
	    echo "!! Run has $warnings warnings or crashed!!"
	    echo "Aborting test..."
	    exit 1
	fi
	info=`grep "FIELDS_INFO:" test${testN}.config | awk '{print $2}'`
	nside=`grep "NSIDE:" test${testN}.config | awk '{print $2}'`
	mask=`grep "SELEC_PREFIX:" test${testN}.config | awk '{print $2}'`
	nzFile=`grep "SELEC_Z_PREFIX:" test${testN}.config | awk '{print $2}'`
	selScale=`grep "SELEC_SCALE:" test${testN}.config | awk '{print $2}'`
	starmask=`grep "STARMASK:" test${testN}.config | awk '{print $2}'`
	separable=`grep "SELEC_SEPARABLE:" test${testN}.config | awk '{print $2}'`
	./check-maps.py $info $nside $mask $nzFile $selScale $starmask tout${testN}_map- $separable
    done
done

dist=HOMOGENEOUS
poisson=1
echo $dist
echo "POISSON: $poisson"
../bin/flask test${testN}.config DIST: $dist POISSON: $poisson EXIT_AT: MAPWERFITS_PREFIX MAPWERFITS_PREFIX: tout${testN}_map- 1>tout${testN}.log 2>>tout${testN}.log
warnings=`grep "Total number of warnings:" tout${testN}.log | cut -d: -f2`
if [ $warnings -eq 0 ]; then
    echo $zerowarnings
else
    echo "!! Run has $warnings warnings or crashed!!"
    echo "Aborting test..."
    exit 1
fi
info=`grep "FIELDS_INFO:" test${testN}.config | awk '{print $2}'`
nside=`grep "NSIDE:" test${testN}.config | awk '{print $2}'`
mask=`grep "SELEC_PREFIX:" test${testN}.config | awk '{print $2}'`
nzFile=`grep "SELEC_Z_PREFIX:" test${testN}.config | awk '{print $2}'`
selScale=`grep "SELEC_SCALE:" test${testN}.config | awk '{print $2}'`
starmask=`grep "STARMASK:" test${testN}.config | awk '{print $2}'`
separable=`grep "SELEC_SEPARABLE:" test${testN}.config | awk '{print $2}'`
./check-maps.py $info $nside $mask $nzFile $selScale $starmask tout${testN}_map- $separable

	


# Remove test files
#rm -f recCl-*.dat 
