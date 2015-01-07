grep -vh , first-alm-*.dat | grep -vh '^$' > all-alm.dat
sort -n all-alm.dat > recov-alm.dat
rm all-alm.dat -f
rm first-alm-* -f
