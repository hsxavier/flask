# Runs corrlnfields many times, in parallel.
# USAGE:   manyalm.sh <config file> <RANDoffset> <Ntimes> <Nprocessors>
# EXAMPLE: manyalm.sh sree.config 0 1000 5

config=$1
offset=$2
total=$3
nproc=$4
batchsize=`expr $total / $nproc`
lastp=`expr $nproc - 1`
#echo "batch: $batchsize lastp: $lastp"

# Loop over processors
p=0
while [ $p -lt $lastp ]; do
    first=`expr $p \* $batchsize + 1 + $offset`
    last=`expr $p \* $batchsize + $batchsize + $offset`
    echo "p$p: /home/skems/pkphotoz/prog/corrlnfields/src/serialalm.sh $config $first $last &"
    /home/skems/pkphotoz/prog/corrlnfields/src/serialalm.sh $config $first $last &
    p=`expr $p + 1`
done

# Last processor
first=`expr $p \* $batchsize + 1 + $offset`
last=`expr $total + $offset`
echo "p$p: /home/skems/pkphotoz/prog/corrlnfields/src/serialalm.sh $config $first $last &"
/home/skems/pkphotoz/prog/corrlnfields/src/serialalm.sh $config $first $last &
