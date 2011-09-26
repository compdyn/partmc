#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_env out/mosaic_restarted_0001
tail -n +13 out/mosaic_0001_env.txt > out/mosaic_0001_env_tail.txt

../../numeric_diff out/mosaic_restarted_0001_env.txt out/mosaic_0001_env_tail.txt 0 1e-8 0 0 2 0
