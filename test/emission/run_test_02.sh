#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_aero_size_mass out/emission_mc_0001.nc out/emission_mc_size_mass.txt
../../extract_summary_aero_size_mass out/emission_exact_0001.nc out/emission_exact_size_mass.txt
../../numeric_diff out/emission_mc_size_mass.txt out/emission_exact_size_mass.txt 0 5e-2 0 0 2 0
exit $?
