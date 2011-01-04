#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_size_mass 1e-10 1e-4 220 out/brownian_part_0001_ out/brownian_part_size_mass.txt
../../extract_sectional_aero_size_mass out/brownian_sect_ out/brownian_sect_size_mass.txt

../../numeric_diff out/brownian_part_size_mass.txt out/brownian_sect_size_mass.txt 0 0.8 0 0 2 0
