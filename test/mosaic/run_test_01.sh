#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_mc.spec
../../extract_summary_aero_size_num out/mosaic_0001.nc out/mosaic_aero_size_num_summary.txt
../../extract_state_aero_size_num 1e-8 1e-3 160 out/mosaic_state_0001_ out/mosaic_aero_size_num_state.txt
../../numeric_diff out/mosaic_aero_size_num_summary.txt out/mosaic_aero_size_num_state.txt 0 1e-8 0 0 0 0
exit $?
