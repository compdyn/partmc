#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_state_aero_size_mass 1e-8 1e-3 160 out/emission_mc_0001_ out/emission_mc_size_mass.txt"
../../extract_state_aero_size_mass 1e-8 1e-3 160 out/emission_mc_0001_ out/emission_mc_size_mass.txt
echo "../../extract_summary_aero_size_mass out/emission_exact_0001.nc out/emission_exact_size_mass.txt"
../../extract_summary_aero_size_mass out/emission_exact_0001.nc out/emission_exact_size_mass.txt
echo "../../numeric_diff out/emission_mc_size_mass.txt out/emission_exact_size_mass.txt 0 5e-2 0 0 2 0"
../../numeric_diff out/emission_mc_size_mass.txt out/emission_exact_size_mass.txt 0 5e-2 0 0 2 0
exit $?
