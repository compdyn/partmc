#!/bin/bash

# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part.spec
../../bin_average_comp 1e-10 1e-4 24 wet out/average_0001_00000001.nc out/average_comp

../../extract_aero_size_num 1e-10 1e-4 24 out/average_0001_ out/average_size_num.txt
../../extract_aero_size_num 1e-10 1e-4 24 out/average_comp_0001_ out/average_comp_size_num.txt
../../numeric_diff out/average_size_num.txt out/average_comp_size_num.txt 0 1e-12 0 0 2 0
