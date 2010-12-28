#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_size_mass 1e-8 1e-3 160 out/golovin_part_0001_ out/golovin_part_size_mass.txt
../../extract_sectional_aero_size_mass out/golovin_exact_ out/golovin_exact_size_mass.txt

../../numeric_diff out/golovin_part_size_mass.txt out/golovin_exact_size_mass.txt 0 1e-1
exit $?
