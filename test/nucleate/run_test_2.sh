#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_species out/nucleate_part_0001_ out/aero_species.txt
../../numeric_diff out/aero_species.txt out/nucleate_ode_aero_mass.txt 0 0.15 0 0 2 0
exit $?
