#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_aero_total out/emission_mc_0001.nc out/emission_mc_total_summary.txt
../../extract_state_aero_total out/emission_mc_state_0001_ out/emission_mc_total_state.txt
../../numeric_diff out/emission_mc_total_summary.txt out/emission_mc_total_state.txt 0 1e-10 0 0 2 2
exit $?
