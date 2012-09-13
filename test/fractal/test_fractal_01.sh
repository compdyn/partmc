#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../test_fractal_radii_conversion

../../partmc run_part_brown_free_df_3.spec
../../partmc run_sect_brown_free_df_3.spec

../../extract_aero_size --num --dmin 1e-9 --dmax 1e-7 --nbin 100 out/part_brown_free_df_3_0001
../../extract_sectional_aero_size --num out/sect_brown_free_df_3

../../numeric_diff --by col --rel-tol 0.1 out/sect_brown_free_df_3_aero_size_num.txt out/part_brown_free_df_3_0001_aero_size_num.txt
