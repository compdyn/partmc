#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_aero_size_mass out/mosaic_0001.nc out/mosaic_aero_size_mass_summary.txt
../../extract_state_aero_size_mass 1e-8 1e-3 160 out/mosaic_state_0001_ out/mosaic_aero_size_mass_state.txt
../../numeric_diff out/mosaic_aero_size_mass_summary.txt out/mosaic_aero_size_mass_state.txt 0 1e-8 0 0 0 0
exit $?
