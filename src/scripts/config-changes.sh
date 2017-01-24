#!/usr/bin/env bash

#
# This script takes a valid FLASK configuration (input) file, i.e. one that
# contains all necessary keywords and no extra ones (tipically the example.config 
# file) and compares with another one to check which keywords must be added or 
# removed from the latter to make it compatible with the corrent version of 
# FLASK (assuming the former one is). The idea is to aid adaptation to different
# FLASK versions that might contain a different set of keywords. 
#
# USAGE:   config-changes.sh <correct config file> <other config file>
# EXAMPLE: config-changes.sh ~/prog/flask/example.config ~/sims/sim01.config
#
# Written by Henrique S. Xavier, hsxavier@if.usp.br, on 24/jan/2017.
#

# Get input:
new=$1
old=$2

# Get keywords from .config files:
grep ".*:" $new -o | sort > temp-new.txt
grep ".*:" $old -o | sort > temp-old.txt

# Check for changes:
toadd=`diff temp-old.txt temp-new.txt | grep  "^>" | awk '{print $2}'`
torem=`diff temp-old.txt temp-new.txt | grep  "^<" | awk '{print $2}'`

# Print changes:
echo ""
echo "** Add these keywords **"
echo ""
for keyword in $toadd; do
    grep "$keyword" $new
done
echo ""
echo "** Remove these keywords **"
echo ""
for keyword in $torem; do
    echo "$keyword"
done
echo ""

# Remove temporary files:
rm -f temp-new.txt
rm -f temp-old.txt
