name=$1
file=$2

header=`head -1 $file`
arr=($header)

col=0
for word in ${arr[@]}; do
    [[ $word == "$name" ]] && echo $col && break
    ((++col))
done
