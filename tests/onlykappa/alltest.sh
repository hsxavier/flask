echo "rm -f ./Cl-*.dat"
rm -f ./Cl-*.dat
echo "rm -f ./Xi-*.dat"
rm -f ./Xi-*.dat
echo "rm -f ./Cl-*.dat"
rm -f ./gCl-*.dat
echo "rm -f ./Xi-*.dat"
rm -f ./gXi-*.dat
echo "rm -f ./covl-*.dat"
rm -f ./covl-*.dat
echo "rm -f ./regl-*.dat"
rm -f ./regl-*.dat
echo "rm -f ./fields-info.dat"
rm -f ./fields-info.dat

echo "cd /home/skems/pkphotoz/prog/class/"
cd /home/skems/pkphotoz/prog/class/
echo "./class /home/skems/pkphotoz/prog/corrlnfields/tests/onlykappa/test.ini"
./class /home/skems/pkphotoz/prog/corrlnfields/tests/onlykappa/test.ini
echo "cd  /home/skems/pkphotoz/prog/corrlnfields/tests/onlykappa/"
cd  /home/skems/pkphotoz/prog/corrlnfields/tests/onlykappa/
echo "./prepClassInput.py /home/skems/pkphotoz/prog/class/output/test_cl.dat Cl-"
./prepClassInput.py /home/skems/pkphotoz/prog/class/output/test_cl.dat Cl-
echo "./ini2info.py test.ini fields-info.dat"
./ini2info.py test.ini fields-info.dat
echo "./corrlnfields input.config"
./corrlnfields input.config
