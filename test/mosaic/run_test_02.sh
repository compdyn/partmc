#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_aero_total out/mosaic_0001_ out/mosaic_aero_total.txt"
../../extract_aero_total out/mosaic_0001_ out/mosaic_aero_total.txt
echo "../../numeric_diff true_aero_total.txt out/mosaic_aero_total.txt 0 0.2 0 0 2 2"
../../numeric_diff true_aero_total.txt out/mosaic_aero_total.txt 0 0.2 0 0 2 2
exit $?
