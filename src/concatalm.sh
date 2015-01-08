
prefix=recov-alm 
# Join multiple realizations into one file:
echo "Joining $prefix..."
grep -vh , $prefix-*.dat | grep -vh '^$' > ${prefix}_unsort.dat
echo "Sorting $prefix..."
sort -n ${prefix}_unsort.dat > ${prefix}_wlm.dat
echo "Removing lm's from $prefix..."
awk '{print $3, $4, $5, $6}' ${prefix}_wlm.dat > $prefix.dat
# Removes intermediary files:
echo "Erasing temp files from $prefix..."
rm $prefix-*.dat -f
rm ${prefix}_*.dat -f
# Run python script to compute C(l)s:
echo "Cling for $prefix..."
../bin/alm2Cl.py $prefix.dat > ${prefix}_cl.dat
rm $prefix.dat -f

prefix=aux-alm
# Join multiple realizations into one file:
echo "Joining $prefix..."
grep -vh , $prefix-*.dat | grep -vh '^$' > ${prefix}_unsort.dat
echo "Sorting $prefix..."
sort -n ${prefix}_unsort.dat > ${prefix}_wlm.dat
echo "Removing lm's from $prefix..."
awk '{print $3, $4, $5, $6}' ${prefix}_wlm.dat > $prefix.dat
# Removes intermediary files:
echo "Erasing temp files from $prefix..."
rm $prefix-*.dat -f
rm ${prefix}_*.dat -f
# Run python script to compute C(l)s:
echo "Cling for $prefix..."
../bin/alm2Cl.py $prefix.dat > ${prefix}_cl.dat
rm $prefix.dat -f
