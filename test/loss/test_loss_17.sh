#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_drydep_sect_dil.spec
../../partmc run_drydep_part.spec

../../extract_sectional_aero_time out/loss_drydep_sect_dil
../../extract_aero_time out/loss_part_drydep_0001

../../numeric_diff --by col --rel-tol 0.05 out/loss_drydep_sect_dil_aero_time.txt out/loss_part_drydep_0001_aero_time.txt
