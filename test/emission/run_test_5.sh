#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_aero_species out/emission_mc_0001.nc out/emission_mc_species.txt
../../extract_summary_aero_species out/emission_exact_0001.nc out/emission_exact_species.txt
../../numeric_diff out/emission_mc_species.txt out/emission_exact_species.txt 0 1e-2 0 0 2 0
exit $?
