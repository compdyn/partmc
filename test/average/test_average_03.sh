#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_time out/average_0001
../../extract_aero_time out/average_comp_0001
../../numeric_diff out/average_0001_aero_time.txt out/average_comp_0001_aero_time.txt 0 1e-12 0 0 2 0
