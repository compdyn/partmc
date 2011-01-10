#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_env out/mosaic_restarted_0001_ out/mosaic_env_restarted.txt
tail -n +13 out/mosaic_env.txt > out/mosaic_env_tail.txt

../../numeric_diff out/mosaic_env_restarted.txt out/mosaic_env_tail.txt 0 1e-8 0 0 2 0
