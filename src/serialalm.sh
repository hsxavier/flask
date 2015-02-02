# This script runs corrlnfields many times, each time with a different RNDSEED and LNALM_OUT.
# USAGE: serialalm.sh <first> <last>

first=$1
last=$2

i=$first
while [ $i -le $last ]; do
    ../bin/corrlnfields ../sree.config RNDSEED: $i RECOVALM_OUT: recov-alm-$i.dat > /dev/null #AUXALM_OUT: aux-alm-$i.dat > /dev/null
    i=`expr $i + 1`
done
