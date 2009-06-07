#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_state_aero_species out/mosaic_0001_ out/mosaic_aero_species.txt"
../../extract_state_aero_species out/mosaic_0001_ out/mosaic_aero_species.txt
echo "../../numeric_diff true_aero_species.txt out/mosaic_aero_species.txt 0 1e-8 0 0 0 0"
../../numeric_diff true_aero_species.txt out/mosaic_aero_species.txt 0 1e-8 0 0 0 0
exit $?
