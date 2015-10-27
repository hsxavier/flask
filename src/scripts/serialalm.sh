#!/usr/bin/env bash

#
# USAGE:   serialalm.sh <config ><first> <last>
# EXAMPLE: serialalm.sh sree.config 10 100
#
# This script runs FLASK many times, each time with a different RNDSEED and RECOVALM_OUT.
# It can be used as a propotype for running FLASK many times with different seeds.
#
# Written by: Henrique S. Xavier, hsxavier@if.usp.br, 06/jul/2015. 
#
 
config=$1
first=$2
last=$3

i=$first
while [ $i -le $last ]; do
    /home/skems/pkphotoz/prog/corrlnfields/bin/corrlnfields $config RNDSEED: $i RECOVALM_OUT: recov-alm-$i.dat > /dev/null
    i=`expr $i + 1`
done
