#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_aero_total out/emission_mc_0001.nc out/emission_mc_total.txt
../../extract_summary_aero_total out/emission_exact_0001.nc out/emission_exact_total.txt
../../numeric_diff out/emission_mc_total.txt out/emission_exact_total.txt 0 3e-2 0 0 3 3
exit $?
