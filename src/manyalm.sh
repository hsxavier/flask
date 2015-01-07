# Runs corrlnfields many times, in parallel.
# USAGE:   manyalm.sh <RANDoffset> <Ntimes> <Nprocessors>
# EXAMPLE: manyalm.sh 0 1000 5

offset=$1
total=$2
nproc=$3
batchsize=`expr $total / $nproc`
lastp=`expr $nproc - 1`
#echo "batch: $batchsize lastp: $lastp"

# Loop over processors
p=0
while [ $p -lt $lastp ]; do
    first=`expr $p \* $batchsize + 1 + $offset`
    last=`expr $p \* $batchsize + $batchsize + $offset`
    echo "p$p: ../bin/serialalm.sh $first $last &"
    ../bin/serialalm.sh $first $last &
    p=`expr $p + 1`
done

# Last processor
first=`expr $p \* $batchsize + 1 + $offset`
last=`expr $total + $offset`
echo "p$p: ../bin/serialalm.sh $first $last &"
../bin/serialalm.sh $first $last &
