#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_part_restarted.spec

../../extract_aero_time out/mosaic_restarted_0001
tail -n +13 out/mosaic_0001_time.txt > out/mosaic_0001_time_tail.txt

../../numeric_diff out/mosaic_restarted_0001_time.txt out/mosaic_0001_time_tail.txt 0 0.1 0 0 2 0
