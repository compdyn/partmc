#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_env out/mosaic_0001_ out/mosaic_env.txt
../../numeric_diff true_env.txt out/mosaic_env.txt 0 1e-8 0 0 0 0
exit $?
