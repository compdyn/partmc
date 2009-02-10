#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_mc.spec
../../partmc run_exact.spec
../../extract_summary_aero_size_num out/emission_mc_0001.nc out/emission_mc_size_num.txt
../../extract_summary_aero_size_num out/emission_exact_0001.nc out/emission_exact_size_num.txt
../../numeric_diff out/emission_mc_size_num.txt out/emission_exact_size_num.txt 0 5e-2 0 0 2 0
exit $?
