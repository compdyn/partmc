#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_aero_total out/emission_part_0001_ out/emission_part_total.txt"
../../extract_aero_total out/emission_part_0001_ out/emission_part_total.txt
echo "../../extract_sectional_aero_total out/emission_exact_ out/emission_exact_total.txt"
../../extract_sectional_aero_total out/emission_exact_ out/emission_exact_total.txt

echo "../../numeric_diff out/emission_part_total.txt out/emission_exact_total.txt 0 5e-2 0 0 2 2"
../../numeric_diff out/emission_part_total.txt out/emission_exact_total.txt 0 5e-2 0 0 2 2
exit $?
