#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_constant_part.spec
../../partmc run_constant_exact.spec

../../extract_aero_time out/loss_part_constant_0001
../../extract_sectional_aero_time out/loss_exact_constant

../../numeric_diff --by col --rel-tol 0.05 out/loss_exact_constant_aero_time.txt out/loss_part_constant_0001_aero_time.txt
