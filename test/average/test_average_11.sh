#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_time out/average_sizevol_0001
../../numeric_diff out/average_0001_time.txt out/average_sizevol_0001_time.txt 0 0.1 0 0 2 0
