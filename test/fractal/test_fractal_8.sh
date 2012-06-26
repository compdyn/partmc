#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../test_fractal_dimless_time --free out/part_vemury_free_df_3_0001

../../numeric_diff --by col --rel-tol 0.1 ref_vemury_free_df_3_dimless_time_regrid.txt out/part_vemury_free_df_3_0001_dimless_time.txt
