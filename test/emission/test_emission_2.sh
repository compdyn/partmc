#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_size --mass --dmin 1e-8 --dmax 1e-3 --nbin 160 out/emission_part_0001
../../extract_sectional_aero_size_mass out/emission_exact_ out/emission_exact_size_mass.txt

../../numeric_diff out/emission_part_0001_size_mass.txt out/emission_exact_size_mass.txt 0 5e-2 0 0 2 0
