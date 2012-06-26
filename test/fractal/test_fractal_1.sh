#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part_vemury_free_df_3.spec
../../partmc run_sect_vemury_free_df_3.spec

../../extract_aero_size --num --dmin 1e-9 --dmax 1e-7 --nbin 100 out/part_vemury_free_df_3_0001
../../extract_sectional_aero_size --num out/sect_vemury_free_df_3

../../numeric_diff --by col --rel-tol 0.2 out/sect_vemury_free_df_3_aero_size_num.txt out/part_vemury_free_df_3_0001_aero_size_num.txt
