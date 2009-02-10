#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_aero_size_num out/emission_mc_0001.nc out/emission_mc_size_num_summary.txt
../../extract_state_aero_size_num 1e-8 1e-3 160 out/emission_mc_state_0001_ out/emission_mc_size_num_state.txt
../../numeric_diff out/emission_mc_size_num_summary.txt out/emission_mc_size_num_state.txt 0 1e-8
exit $?
