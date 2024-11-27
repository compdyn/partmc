#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part_nummass.spec

../../extract_aero_size --num --dmin 1e-8 --dmax 1e-3 --nbin 160 out/weighting_nummass_0001
../../extract_sectional_aero_size --num out/weighting_exact
../../numeric_diff --by col --rel-tol 0.1 out/weighting_exact_aero_size_num.txt out/weighting_nummass_0001_aero_size_num.txt
