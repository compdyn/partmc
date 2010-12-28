#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

#####################

for f in out/brownian_part_????_00000001.nc ; do
    f1=${f/_00000001.nc/}
    f2=${f1/_part_/_part_size_num_}.txt
    ../../extract_aero_size_num 1e-10 1e-4 220 ${f1}_ $f2
    f3=${f2/_part_size_num_/_part_size_num_error_}
    ../../numeric_diff $f2 out/brownian_sect_size_num.txt 0 0 0 0 2 0 > $f3
done
../../numeric_average out/brownian_part_average_num_error.txt out/brownian_part_size_num_error_????.txt

#####################

../../numeric_average out/brownian_part_size_num_average.txt out/brownian_part_size_num_????.txt

../../extract_sectional_aero_size_num out/brownian_sect_ out/brownian_sect_size_num.txt

../../numeric_diff out/brownian_part_size_num_average.txt out/brownian_sect_size_num.txt 0 0 0 0 2 0 >> out/brownian_part_average_num_error.txt

#####################

# **********************************************************
# Testing for 1/sqrt(n) falloff in error from MC averaging
# The second line below should be a factor of 3 smaller than the first line
cat out/brownian_part_average_num_error.txt

../../numeric_diff expected_error_num.txt out/brownian_part_average_num_error.txt 0 0.5 1 1 2 0
status1=$?
../../numeric_diff expected_error_num.txt out/brownian_part_average_num_error.txt 0 0.5 2 2 2 0
status2=$?
exit $(( $status1 && $status2 ))
