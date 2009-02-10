#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_mc.spec
../../partmc run_exact.spec

../../extract_summary_aero_size_num out/golovin_mc_0001.nc out/golovin_mc_num.txt
../../extract_summary_aero_size_num out/golovin_exact_0001.nc out/golovin_exact_num.txt

../../numeric_diff out/golovin_mc_num.txt out/golovin_exact_num.txt 0 2e-2 0 0 2 0
exit $?
