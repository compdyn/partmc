#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../bin_average_size 1e-10 1e-4 24 wet average number out/average_0001_00000001.nc out/average_sizenum

../../extract_aero_size --num --dmin 1e-10 --dmax 1e-4 --nbin 24 out/average_sizenum_0001
../../numeric_diff --by col --rel-tol 1e-12 out/average_0001_aero_size_num.txt out/average_sizenum_0001_aero_size_num.txt
