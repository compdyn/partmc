#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_aero_total out/mosaic_restarted_0001_ out/mosaic_aero_total_restarted.txt"
../../extract_aero_total out/mosaic_restarted_0001_ out/mosaic_aero_total_restarted.txt
echo "tail -n +13 out/mosaic_aero_total.txt > out/mosaic_aero_total_tail.txt"
tail -n +13 out/mosaic_aero_total.txt > out/mosaic_aero_total_tail.txt

echo "../../numeric_diff out/mosaic_aero_total_restarted.txt out/mosaic_aero_total_tail.txt 0 0.1 0 0 3 3"
../../numeric_diff out/mosaic_aero_total_restarted.txt out/mosaic_aero_total_tail.txt 0 0.1 0 0 3 3
exit $?
