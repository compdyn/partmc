#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out_SPSD

../../partmc run_part_brown_cont_df_3.spec

../../test_fractal_self_preserve --dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 out_SPSD/part_brown_cont_df_3_0001
../../test_fractal_self_preserve --dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 out_SPSD/part_brown_cont_df_3_0002
../../test_fractal_self_preserve --dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 out_SPSD/part_brown_cont_df_3_0003
../../test_fractal_self_preserve --dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 out_SPSD/part_brown_cont_df_3_0004
../../test_fractal_self_preserve --dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 out_SPSD/part_brown_cont_df_3_0005
../../test_fractal_self_preserve --dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 out_SPSD/part_brown_cont_df_3_0006
../../test_fractal_self_preserve --dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 out_SPSD/part_brown_cont_df_3_0007
../../test_fractal_self_preserve --dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 out_SPSD/part_brown_cont_df_3_0008
../../test_fractal_self_preserve --dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 out_SPSD/part_brown_cont_df_3_0009
../../test_fractal_self_preserve --dimless_vol_min 1e-3 --dimless_vol_max 10 --nbin 100 out_SPSD/part_brown_cont_df_3_0010

../../numeric_diff --by col --rel-tol 0.1 ref_cont_df_3_self_preserve_regrid.txt out_SPSD/part_brown_cont_df_3_0001_self_preserve.txt
