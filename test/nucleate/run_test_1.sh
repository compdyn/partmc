#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../partmc run_part.spec"
../../partmc run_part.spec

echo "../../extract_aero_species out/nucleate_part_0001_ out/aero_species.txt"
../../extract_aero_species out/nucleate_part_0001_ out/aero_species.txt
echo "../../extract_gas out/nucleate_part_0001_ out/gas.txt"
../../extract_gas out/nucleate_part_0001_ out/gas.txt

#echo "../../numeric_diff out/emission_part_size_num.txt out/emission_exact_size_num.txt 0 5e-2 0 0 2 0"
#../../numeric_diff out/emission_part_size_num.txt out/emission_exact_size_num.txt 0 5e-2 0 0 2 0
#exit $?
