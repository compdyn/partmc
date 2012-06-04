#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_size --num --dmin 1e-9 --dmax 1e-5 --nbin 100 out/part_naumann_free
../../extract_sectional_aero_size --num out/sect_naumann_free

../../numeric_diff --by col --rel-tol 0.7 out/sect_naumann_free_aero_size_num.txt out/part_naumann_free_aero_size_num.txt
