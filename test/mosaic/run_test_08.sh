#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_total out/mosaic_restarted_0001_ out/mosaic_aero_total_restarted.txt
tail -n +13 out/mosaic_aero_total.txt > out/mosaic_aero_total_tail.txt

../../numeric_diff out/mosaic_aero_total_restarted.txt out/mosaic_aero_total_tail.txt 0 0.1 0 0 3 3
exit $?
