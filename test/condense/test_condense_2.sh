#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_species out/condense_0001_ out/condense_aero_species.txt

../../numeric_diff ref_aero_species.txt out/condense_aero_species.txt 0 1e-6
