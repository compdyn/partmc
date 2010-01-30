#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../partmc run_part.spec"
../../partmc run_part.spec
echo "../../extract_aero_species out/mosaic_0001_ out/mosaic_aero_species.txt"
../../extract_aero_species out/mosaic_0001_ out/mosaic_aero_species.txt
echo "../../numeric_diff true_aero_species.txt out/mosaic_aero_species.txt 0 0.1 0 0 2 0"
../../numeric_diff true_aero_species.txt out/mosaic_aero_species.txt 0 0.1 0 0 2 0
exit $?
