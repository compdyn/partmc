#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_part.spec
../../partmc run_exact.spec

../../extract_aero_size_num 1e-8 1e-3 160 out/golovin_part_0001_ out/golovin_part_size_num.txt
../../extract_sectional_aero_size_num out/golovin_exact_ out/golovin_exact_size_num.txt

../../numeric_diff out/golovin_part_size_num.txt out/golovin_exact_size_num.txt 0 5e-2
exit $?
