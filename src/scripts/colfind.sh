#!/usr/bin/env bash

#
# USAGE:   colfind.sh <COL_NAME> <FILE_W_HEADER>
# EXAMPLE: colfind.sh 
#
# This script take the header of <FILE_W_HEADER> (which is assumed to 
# be on the first line and to start with some comment symbol like '#' 
# followed by an empty space) and looks for <COL_NAME> (the collumn 
# names must be separated by spaces). It prints on stdout the number of 
# the column corresponding to <COL_NAME>.
#
# Written by: Henrique S. Xavier, hsxavier@if.usp.br, 21/aug/2015.
#

name=$1
file=$2

header=`head -1 $file`
arr=($header)

col=0
for word in ${arr[@]}; do
    [[ $word == "$name" ]] && echo $col && break
    ((++col))
done
