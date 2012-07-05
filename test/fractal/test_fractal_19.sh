#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out_dimless_t
mkdir -p out_dimless_t/restart

../../partmc run_part_vemury_free_df_2_upto1000s.spec
../../partmc run_part_vemury_free_df_2_restart.spec

../../test_fractal_dimless_time --free out_dimless_t/part_vemury_free_df_2_0001
../../test_fractal_dimless_time --free out_dimless_t/restart/part_vemury_free_df_2_0001
../../test_fractal_merge_files part_vemury_free_df_2_0001

../../numeric_diff --by col --rel-tol 0.1 ref_free_df_2_dimless_time_regrid.txt out_dimless_t/part_vemury_free_df_2_0001_dimless_t_series.txt
