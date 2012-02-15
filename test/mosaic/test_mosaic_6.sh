#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_gas out/mosaic_restarted_0001
tail -n +13 out/mosaic_0001_gas.txt > out/mosaic_0001_gas_tail.txt

../../numeric_diff --by elem --min-col 2 --rel-tol 1e-4 out/mosaic_0001_gas_tail.txt out/mosaic_restarted_0001_gas.txt
