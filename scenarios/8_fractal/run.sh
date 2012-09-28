#!/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

mkdir -p out
mkdir -p out_dimless_t
mkdir -p out_dimless_t/restart

../../build/partmc run_part_brown_free_df_3.spec
../../build/partmc run_sect_brown_free_df_3.spec

../../build/partmc run_part_brown_free_df_2_6.spec
../../build/partmc run_part_brown_free_df_2_2.spec
../../build/partmc run_part_brown_free_df_2.spec
../../build/partmc run_part_brown_cont_df_3.spec
../../build/partmc run_part_brown_cont_df_2_2.spec
../../build/partmc run_part_brown_cont_df_1_8.spec

../../build/partmc run_part_brown_free_df_3_upto1000s.spec
../../build/partmc run_part_brown_free_df_3_restart.spec
../../build/partmc run_part_brown_free_df_2_8_upto1000s.spec
../../build/partmc run_part_brown_free_df_2_8_restart.spec
../../build/partmc run_part_brown_free_df_2_4_upto1000s.spec
../../build/partmc run_part_brown_free_df_2_4_restart.spec
../../build/partmc run_part_brown_free_df_2_upto1000s.spec
../../build/partmc run_part_brown_free_df_2_restart.spec
../../build/partmc run_part_brown_cont_df_3_upto1000s.spec
../../build/partmc run_part_brown_cont_df_3_restart.spec
../../build/partmc run_part_brown_cont_df_1_8_upto1000s.spec
../../build/partmc run_part_brown_cont_df_1_8_restart.spec
../../build/partmc run_part_brown_cont_df_1_upto1000s.spec
../../build/partmc run_part_brown_cont_df_1_restart.spec

# Now run ./process.sh to process the data
