#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_species out/average_compsizevol_0001_ out/average_compsizevol_species.txt
../../numeric_diff out/average_species.txt out/average_compsizevol_species.txt 0 1e-12 0 0 2 0
