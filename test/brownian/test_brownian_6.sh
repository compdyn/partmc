#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

for f in out/brownian_part_????_00000001.nc ; do
    ../../extract_aero_size --mass --dmin 1e-10 --dmax 1e-4 --nbin 220 ${f/_00000001.nc/}
    ../../numeric_diff ${f/00000001.nc/size_num.txt} out/brownian_sect_aero_size_mass.txt 0 0 0 0 2 0 > ${f/00000001.nc/size_mass_error.txt}
done
../../numeric_average out/brownian_part_average_mass_error.txt out/brownian_part_????_size_mass_error.txt

../../numeric_average out/brownian_part_size_mass_average.txt out/brownian_part_????_size_mass.txt
../../extract_sectional_aero_size --mass out/brownian_sect
../../numeric_diff out/brownian_part_size_mass_average.txt out/brownian_sect_aero_size_mass.txt 0 0 0 0 2 0 >> out/brownian_part_average_mass_error.txt

# **********************************************************
# Testing for 1/sqrt(n) falloff in error from MC averaging
# The second line below should be a factor of 3 smaller than the first line
cat out/brownian_part_average_mass_error.txt

../../numeric_diff expected_error_mass.txt out/brownian_part_average_mass_error.txt 0 0.7 1 1 2 0
../../numeric_diff expected_error_mass.txt out/brownian_part_average_mass_error.txt 0 0.7 2 2 2 0
