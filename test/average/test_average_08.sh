#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_particles -o out/average_sizenum_particles_unsorted.txt out/average_sizenum_0001_00000001.nc
sort out/average_sizenum_particles_unsorted.txt > out/average_sizenum_particles.txt
../../numeric_diff --by elem --min-col 3 --max-col 4 --rel-tol 3 out/average_particles.txt out/average_sizenum_particles.txt
