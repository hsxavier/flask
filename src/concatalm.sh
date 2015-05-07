alm2cl=$1 # 1 - Calculates Cl and erases alms; 0 - Keep alms and do not calculate Cls.

prefix=recov-alm 
# Join multiple realizations into one file:
echo "Joining $prefix..."
grep -vh , $prefix-*.dat | grep -vh '^$' > ${prefix}_unsort.dat
echo "Sorting $prefix..."
sort -n ${prefix}_unsort.dat > ${prefix}_wlm.dat
echo "Removing lm's from $prefix..."
awk '{for(i=3;i<NF;i++) printf"%s",$i OFS;if(NF)printf"%s",$NF;printf ORS}' ${prefix}_wlm.dat > $prefix.dat
# Removes intermediary files:
echo "Erasing temp files from $prefix..."
rm $prefix-*.dat -f
rm ${prefix}_*.dat -f
# Run python script to compute C(l)s:
if [ $alm2cl -eq 1 ]; then
    echo "Cling for $prefix..."
    ../bin/alm2Cl.py $prefix.dat > ${prefix}_cl.dat
    rm $prefix.dat -f
fi
