#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_size --mass --dmin 1e-9 --dmax 1e-3 --nbin 160 out/loss_part_chamber_0001
../../extract_sectional_aero_size --mass out/loss_exact_chamber

../../numeric_diff --by col --rel-tol 0.1 out/loss_exact_chamber_aero_size_mass.txt out/loss_part_chamber_0001_aero_size_mass.txt
