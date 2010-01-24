#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_gas out/nucleate_part_0001_ out/gas.txt"
../../extract_gas out/nucleate_part_0001_ out/gas.txt
echo "../../numeric_diff out/gas.txt out/nucleate_ode_gas.txt 0 5e-2 0 0 2 0"
../../numeric_diff out/gas.txt out/nucleate_ode_gas.txt 0 5e-2 0 0 2 0
exit $?
