#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_gas out/mosaic_restarted_0001_ out/mosaic_gas_restarted.txt
tail -n +13 out/mosaic_gas.txt > out/mosaic_gas_tail.txt

../../numeric_diff out/mosaic_gas_restarted.txt out/mosaic_gas_tail.txt 0 0.01 0 0 2 0
exit $?
