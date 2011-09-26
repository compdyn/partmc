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
../../partmc run_sect.spec

../../extract_aero_size --num --dmin 1e-10 --dmax 1e-4 --nbin 220 out/brownian_part_0001
../../extract_sectional_aero_size_num out/brownian_sect_ out/brownian_sect_size_num.txt
../../numeric_diff out/brownian_part_0001_size_num.txt out/brownian_sect_size_num.txt 0 0.3 0 0 2 0
