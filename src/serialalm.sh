# This script runs corrlnfields many times, each time with a different RNDSEED and LNALM_OUT.
# USAGE: serialalm.sh <config ><first> <last>

config=$1
first=$2
last=$3

i=$first
while [ $i -le $last ]; do
    /home/skems/pkphotoz/prog/corrlnfields/bin/corrlnfields $config RNDSEED: $i RECOVALM_OUT: recov-alm-$i.dat > /dev/null
    i=`expr $i + 1`
done
