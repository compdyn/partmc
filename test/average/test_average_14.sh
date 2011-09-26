#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_size --mass --dmin 1e-10 --dmax 1e-4 --nbin 24 out/average_compsizenum_0001
../../numeric_diff out/average_0001_size_mass.txt out/average_compsizenum_0001_size_mass.txt 0 0.1 0 0 2 0
