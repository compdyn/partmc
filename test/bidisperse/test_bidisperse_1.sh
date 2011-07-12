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
../../test_bidisperse_extract

../../test_bidisperse_ode

# extract size distributions for plotting
../../extract_aero_size_num 1e-5 1e-3 255 out/bidisperse_part_0001_ out/bidisperse_part_aero_size_num.txt

../../numeric_diff out/bidisperse_part_data.txt out/bidisperse_ode_data.txt 0 1e-5 0 0 1 1
