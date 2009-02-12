#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

for f in out/brownian_mc_????.nc ; do
    f1=${f/_mc_/_mc_size_num_}
    ../../extract_summary_aero_size_num $f ${f1/.nc/.txt};
done
../../numeric_average out/brownian_mc_size_num_average.txt out/brownian_mc_size_num_????.txt

../../extract_summary_aero_size_num out/brownian_sect_0001.nc out/brownian_sect_size_num.txt

../../numeric_diff out/brownian_mc_size_num_average.txt out/brownian_sect_size_num.txt 0 0.2 0 0 2 0
exit $?
