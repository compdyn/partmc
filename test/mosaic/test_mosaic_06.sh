#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_part_restarted.spec

../../extract_aero_species out/mosaic_restarted_0001_ out/mosaic_aero_species_restarted.txt
tail -n +13 out/mosaic_aero_species.txt > out/mosaic_aero_species_tail.txt

../../numeric_diff out/mosaic_aero_species_restarted.txt out/mosaic_aero_species_tail.txt 0 0.1 0 0 2 0
