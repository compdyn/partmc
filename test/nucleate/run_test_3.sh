#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_gas out/nucleate_part_0001_ out/gas.txt
../../numeric_diff out/gas.txt out/nucleate_ode_gas.txt 0 5e-2 0 0 2 0
exit $?
