#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_state_aero_total out/emission_mc_0001_ out/emission_mc_total.txt"
../../extract_state_aero_total out/emission_mc_0001_ out/emission_mc_total.txt
echo "../../extract_summary_aero_total out/emission_exact_0001.nc out/emission_exact_total.txt"
../../extract_summary_aero_total out/emission_exact_0001.nc out/emission_exact_total.txt
echo "../../numeric_diff out/emission_mc_total.txt out/emission_exact_total.txt 0 5e-2 0 0 3 3"
../../numeric_diff out/emission_mc_total.txt out/emission_exact_total.txt 0 5e-2 0 0 3 3
exit $?
