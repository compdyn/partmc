#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../bin_average_size 1e-10 1e-4 24 wet average number out/average_comp_0001_00000001.nc out/average_compsizenum

../../extract_aero_size_num 1e-10 1e-4 24 out/average_compsizenum_0001_ out/average_compsizenum_size_num.txt
../../numeric_diff out/average_size_num.txt out/average_compsizenum_size_num.txt 0 1e-12 0 0 2 0
