#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_aero_species out/condense_0001_ out/condense_aero_species.txt"
../../extract_aero_species out/condense_0001_ out/condense_aero_species.txt

echo "../../numeric_diff true_aero_species.txt out/condense_aero_species.txt 0 1e-8"
../../numeric_diff true_aero_species.txt out/condense_aero_species.txt 0 1e-6
exit $?
