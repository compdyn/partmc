#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../partmc run_mc.spec"
../../partmc run_mc.spec
echo "../../partmc run_sect.spec"
../../partmc run_sect.spec

echo "../../extract_aero_size_num 1e-10 1e-4 220 out/brownian_mc_0001_ out/brownian_mc_size_num.txt"
../../extract_aero_size_num 1e-10 1e-4 220 out/brownian_mc_0001_ out/brownian_mc_size_num.txt
echo "../../extract_sectional_aero_size_num out/brownian_sect_ out/brownian_sect_size_num.txt"
../../extract_sectional_aero_size_num out/brownian_sect_ out/brownian_sect_size_num.txt

echo "../../numeric_diff out/brownian_mc_size_num.txt out/brownian_sect_size_num.txt 0 0.5 0 0 2 0"
../../numeric_diff out/brownian_mc_size_num.txt out/brownian_sect_size_num.txt 0 0.5 0 0 2 0
exit $?
