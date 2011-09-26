#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_size --mass --dmin 1e-10 --dmax 1e-4 --nbin 24 out/average_0001
../../extract_aero_size --mass --dmin 1e-10 --dmax 1e-4 --nbin 24 out/average_comp_0001
../../numeric_diff --by col --rel-tol 1e-12 out/average_0001_aero_size_mass.txt out/average_comp_0001_aero_size_mass.txt
