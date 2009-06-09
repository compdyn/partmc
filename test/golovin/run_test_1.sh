#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../partmc run_mc.spec"
../../partmc run_mc.spec
echo "../../partmc run_exact.spec"
../../partmc run_exact.spec

echo "../../extract_state_aero_size_num 1e-8 1e-3 160 out/golovin_mc_0001_ out/golovin_mc_num.txt"
../../extract_state_aero_size_num 1e-8 1e-3 160 out/golovin_mc_0001_ out/golovin_mc_num.txt
echo "../../extract_sectional_aero_size_num out/golovin_exact_ out/golovin_exact_num.txt"
../../extract_sectional_aero_size_num out/golovin_exact_ out/golovin_exact_num.txt

echo "../../numeric_diff out/golovin_mc_num.txt out/golovin_exact_num.txt 0 2e-2"
../../numeric_diff out/golovin_mc_num.txt out/golovin_exact_num.txt 0 2e-2
exit $?
