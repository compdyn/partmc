#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_aero_species out/nucleate_part_0001_ out/aero_species.txt"
../../extract_aero_species out/nucleate_part_0001_ out/aero_species.txt
echo "../../numeric_diff out/aero_species.txt out/nucleate_ode_aero_mass.txt 0 0.1 0 0 2 0"
../../numeric_diff out/aero_species.txt out/nucleate_ode_aero_mass.txt 0 0.1 0 0 2 0
exit $?
