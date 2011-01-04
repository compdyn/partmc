#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_total out/mosaic_0001_ out/mosaic_aero_total.txt
../../numeric_diff true_aero_total.txt out/mosaic_aero_total.txt 0 0.3 0 0 2 2
