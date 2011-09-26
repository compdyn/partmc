#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part.spec
../../extract_aero_time out/mosaic_0001
../../numeric_diff ref_mosaic_0001_aero_time.txt out/mosaic_0001_aero_time.txt 0 0.1 0 0 4 0
