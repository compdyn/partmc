#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part_flat_specified.spec

../../extract_aero_size --num --dmin 1e-8 --dmax 1e-3 --nbin 160 out/weighting_flat_specified_0001
../../numeric_diff --by col --rel-tol 0.1 out/weighting_exact_aero_size_num.txt out/weighting_flat_specified_0001_aero_size_num.txt
