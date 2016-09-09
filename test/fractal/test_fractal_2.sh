#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part_brown_free_df_2_4_upto1000s.spec
../../partmc run_part_brown_free_df_2_4_restart.spec

../../test_fractal_dimless_time --free out/part_brown_free_df_2_4_0001
../../test_fractal_dimless_time --free out/restart_part_brown_free_df_2_4_0001
tail -n +2 out/restart_part_brown_free_df_2_4_0001_dimless_time.txt > out/restart_part_brown_free_df_2_4_0001_dimless_time_tailed.txt
cat out/part_brown_free_df_2_4_0001_dimless_time.txt out/restart_part_brown_free_df_2_4_0001_dimless_time_tailed.txt > out/part_brown_free_df_2_4_0001_dimless_t_series.txt

../../numeric_diff --by col --rel-tol 0.1 ref_free_df_2_4_dimless_time_regrid.txt out/part_brown_free_df_2_4_0001_dimless_t_series.txt
