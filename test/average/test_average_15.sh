#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_time out/average_compsizenum_0001
../../numeric_diff --rel-tol 1e-10 out/average_0001_aero_time.txt out/average_compsizenum_0001_aero_time.txt
