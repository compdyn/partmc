#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part.spec
../../extract_aero_species out/mosaic_0001_ out/mosaic_aero_species.txt
../../numeric_diff true_aero_species.txt out/mosaic_aero_species.txt 0 0.1 0 0 2 0
exit $?
