#!/usr/bin/env bash

testN=$1
dist=$2
min=$3
max=$4
counter=$min

config1=test${testN}.config
outprefix1=tout${testN}_recCl-

offset=7000

while [ $counter -le $max ]; do
    if [ $counter -lt 10 ]; then
	tag=000$counter
    elif [ $counter -lt 100 ]; then
	tag=00$counter
    elif [ $counter -lt 1000 ]; then
	tag=0$counter
    else
	tag=$counter
    fi

    seed=`echo "$offset + $counter" | bc`
    ../bin/flask $config1 DIST: $dist RNDSEED: $seed EXIT_AT: RECOVCLS_OUT RECOVCLS_OUT: $outprefix1$tag.dat > $outprefix1$tag.log

    counter=`echo "$counter + 1" | bc`
done

