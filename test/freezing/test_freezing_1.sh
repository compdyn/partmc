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
#./extract_data.py
../../test_freezing_extract

#mkdir -p theoretical
#./theoretical_freezing.py
../../test_freezing_theoretical


../../numeric_diff --rel-tol 0.02 out/freezing_part_data.txt out/freezing_theoretical_data.txt
