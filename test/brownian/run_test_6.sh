#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

#####################

for f in out/brownian_mc_????.nc ; do
    f1=${f/_mc_/_mc_size_mass_}
    f2=${f1/.nc/.txt}
    ../../extract_summary_aero_size_mass $f $f2
    f3=${f2/_mc_size_mass_/_mc_size_mass_error_}
    ../../numeric_diff $f2 out/brownian_sect_size_mass.txt 0 0 0 0 2 0 > $f3
done
../../numeric_average out/brownian_mc_average_mass_error.txt out/brownian_mc_size_mass_error_????.txt

#####################

../../numeric_average out/brownian_mc_size_mass_average.txt out/brownian_mc_size_mass_????.txt

../../extract_summary_aero_size_mass out/brownian_sect_0001.nc out/brownian_sect_size_mass.txt

../../numeric_diff out/brownian_mc_size_mass_average.txt out/brownian_sect_size_mass.txt 0 0 0 0 2 0 >> out/brownian_mc_average_mass_error.txt

#####################

echo "**********************************************************"
echo "Testing for 1/sqrt(n) falloff in error from MC averaging"
echo "The second line below should be a factor of 3 smaller than the first line"
cat out/brownian_mc_average_mass_error.txt

../../numeric_diff expected_error_mass.txt out/brownian_mc_average_mass_error.txt 0 1.0 1 1 2 0
status1=$?
../../numeric_diff expected_error_mass.txt out/brownian_mc_average_mass_error.txt 0 1.0 2 2 2 0
status2=$?
exit $(( $status1 && $status2 ))
