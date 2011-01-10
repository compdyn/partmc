#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_size_mass 1e-10 1e-4 24 out/average_0001_ out/average_size_mass.txt
../../extract_aero_size_mass 1e-10 1e-4 24 out/average_comp_0001_ out/average_comp_size_mass.txt
../../numeric_diff out/average_size_mass.txt out/average_comp_size_mass.txt 0 1e-12 0 0 2 0
