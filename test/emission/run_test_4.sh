#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_aero_species out/emission_mc_0001.nc out/emission_mc_species_summary.txt
../../extract_state_aero_species out/emission_mc_state_0001_ out/emission_mc_species_state.txt
../../numeric_diff out/emission_mc_species_summary.txt out/emission_mc_species_state.txt 1e-5 1e-8
exit $?
