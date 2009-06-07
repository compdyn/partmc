#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_state_aero_size_mass 1e-8 1e-3 160 out/mosaic_0001_ out/mosaic_aero_size_mass.txt"
../../extract_state_aero_size_mass 1e-8 1e-3 160 out/mosaic_0001_ out/mosaic_aero_size_mass.txt
echo "../../numeric_diff true_aero_size_mass.txt out/mosaic_aero_size_mass.txt 0 1e-8 0 0 0 0"
../../numeric_diff true_aero_size_mass.txt out/mosaic_aero_size_mass.txt 0 1e-8 0 0 0 0
exit $?
