#!/usr/bin/env bash

# 
# USAGE: breakClassIni.sh <class_file.ini>
# This script takes a CLASS input file and creates in the current directory 
# one input file for each combination of two redshift bins. Each input file 
# is called "break-<zbin1>-<zbin2>.ini". This avoids CLASS crashing issues 
# while increating the number of Cls to be computed by slightly less than 
# 3 times (with the advantage that the computation can be distributed in a 
# cluster).
# 
# Written by: Henrique S. Xavier
# hsxavier@if.usp.br
# 2015-10-17
#

# Get parameters from command line:
inifile=$1
outpath=$2

# Get prefix for output:
ndots=`echo $inifile | grep "\." -o | wc -l`
prefix=`echo $inifile | cut -d. -f$ndots`

# Get list of redshifts:
zCSVlist=`grep selection_mean $inifile | cut -d= -f2`
nbins=`echo $zCSVlist | grep "," -o | wc -l`
nbins=`echo "$nbins + 1" | bc`

# Get list of bin widths:
dzCSVlist=`grep selection_width $inifile | cut -d= -f2`
dbins=`echo $dzCSVlist | grep "," -o | wc -l`
dbins=`echo "$dbins + 1" | bc`


# Loop over pairs of bins:
bin1=1
while [ $bin1 -lt $nbins ]; do
    bin2=`echo "$bin1 + 1" | bc`
    while [ $bin2 -le $nbins ]; do
	# Get pair of bins and bin widths:
	newbins=`echo $zCSVlist | cut -d, -f$bin1,$bin2`
	if [ $dbins -gt 1 ]; then
	    # If bins have different widths:
	    newdz=`echo $dzCSVlist | cut -d, -f$bin1,$bin2`
	else
	    # If all bins have the same width:
	    newdz=$dzCSVlist
	fi
	# Create .ini file with just two bins:
	sed -e "s/selection_mean.*/selection_mean = $newbins/g" -e "s/selection_width.*/selection_width = $newdz/g" \
	    -e "s/non_diagonal.*/non_diagonal = 1/g" -e "s/root.*/root = ${prefix}-${bin1}-${bin2}_/g" \
	    $inifile > ${outpath}${prefix}-${bin1}-${bin2}.ini
	echo $newbins
	bin2=`echo "$bin2 + 1" | bc`
    done
    bin1=`echo "$bin1 + 1" | bc`
done
