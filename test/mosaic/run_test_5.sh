#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_state_env out/mosaic_0001_ out/mosaic_env.txt"
../../extract_state_env out/mosaic_0001_ out/mosaic_env.txt
echo "../../numeric_diff true_env.txt out/mosaic_env.txt 0 1e-8 0 0 0 0"
../../numeric_diff true_env.txt out/mosaic_env.txt 0 1e-8 0 0 0 0
exit $?
