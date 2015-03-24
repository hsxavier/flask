
# Define variations:
facmin=1.05
facmax=-1.05
facvar=0.05
lambmin=0.001
lambmax=0.050
lambvar=0.001

rm -f eigentab.dat

# LOOP over Correlation:
faccont=1
facnew=$facmin
while [ $faccont -eq 1 ]; do
    
    # Scale cross-Cl:
    ./scaleCl.py bkp_Cl-f1z1f1z2.dat $facnew Cl-f1z1f1z2.dat
    facnew=`echo "scale=10; $facnew - $facvar" | bc`
    faccont=`echo "$facnew > $facmax" | bc`

    # LOOP over shifts:
    cont=1
    lambnew=$lambmin
    while [ $cont -eq 1 ]; do
	# Changing one lambda:
	sed "s/0.0097/0$lambnew/g" info-temp.dat > fields-info.dat
	lambnew=`echo "scale=10; $lambnew + $lambvar" | bc`
	cont=`echo "$lambnew < $lambmax" | bc`
	
	# Run code and save result:
	echo "fac: $facnew shift: $lambnew"
	./corrlnfields input.config > dump.txt
	grep -v "#" debug.dat >> eigentab.dat
    done
done
