#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_mc.spec
../../partmc run_exact.spec
../../extract_summary_aero num out/emission_mc_0001.nc out/emission_mc.txt
../../extract_summary_aero num out/emission_exact_0001.nc out/emission_exact.txt
../../numeric_diff out/emission_mc.txt out/emission_exact.txt 0 1e-2
exit $?
