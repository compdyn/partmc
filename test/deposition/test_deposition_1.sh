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
../../extract_aero_time out/deposition_0001
../../test_deposition_exact run_part.spec

../../numeric_diff --by col -c 4 -C 23 --rel-tol 0.1 out/deposition_0001_aero_time.txt out/deposition_data.txt
