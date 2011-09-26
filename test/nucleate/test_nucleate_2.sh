#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_gas out/nucleate_part_0001
../../numeric_diff --by col --rel-tol 0.05  out/nucleate_ode_gas.txt out/nucleate_part_0001_gas.txt
