#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_gas out/mosaic_restarted_0001_ out/mosaic_gas_restarted.txt"
../../extract_gas out/mosaic_restarted_0001_ out/mosaic_gas_restarted.txt
echo "tail -n +13 out/mosaic_gas.txt > out/mosaic_gas_tail.txt"
tail -n +13 out/mosaic_gas.txt > out/mosaic_gas_tail.txt

echo "../../numeric_diff out/mosaic_gas_restarted.txt out/mosaic_gas_tail.txt 0 0.01 0 0 2 0"
../../numeric_diff out/mosaic_gas_restarted.txt out/mosaic_gas_tail.txt 0 0.01 0 0 2 0
exit $?
