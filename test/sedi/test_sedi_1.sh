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

../../extract_aero_size --num --dmin 1e-8 --dmax 1e-2 --nbin 100 out/sedi_part_0001
../../extract_sectional_aero_size_num out/sedi_sect_ out/sedi_sect_size_num.txt

../../numeric_diff out/sedi_part_0001_size_num.txt out/sedi_sect_size_num.txt 0 0.2 0 0 2 0
