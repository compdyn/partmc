#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../partmc run_mc.spec"
../../partmc run_mc.spec
echo "../../partmc run_sect.spec"
../../partmc run_sect.spec

echo "../../extract_state_aero_size_num 1e-8 1e-2 220 out/sedi_mc_0001_ out/sedi_mc_size_num.txt"
../../extract_state_aero_size_num 1e-8 1e-2 220 out/sedi_mc_0001_ out/sedi_mc_size_num.txt
echo "../../extract_sectional_aero_size_num out/sedi_sect_ out/sedi_sect_size_num.txt"
../../extract_sectional_aero_size_num out/sedi_sect_ out/sedi_sect_size_num.txt

echo "../../numeric_diff out/sedi_mc_size_num.txt out/sedi_sect_size_num.txt 0 5e-2 0 0 2 0"
../../numeric_diff out/sedi_mc_size_num.txt out/sedi_sect_size_num.txt 0 5e-2 0 0 2 0
exit $?
