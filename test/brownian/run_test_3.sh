#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

for f in out/brownian_part_????_00000001.nc ; do
    f1=${f/_00000001.nc/}
    f2=${f1/_part_/_part_size_num_}.txt
    ../../extract_aero_size_num 1e-10 1e-4 220 ${f1}_ $f2
done
../../numeric_average out/brownian_part_size_num_average.txt out/brownian_part_size_num_????.txt

../../extract_sectional_aero_size_num out/brownian_sect_ out/brownian_sect_size_num.txt

../../numeric_diff out/brownian_part_size_num_average.txt out/brownian_sect_size_num.txt 0 0.1 0 0 2 0
