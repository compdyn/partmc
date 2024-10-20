#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part_chapman.spec
../../extract_gas out/tchem_chapman_0001
../../numeric_diff --by col --rel-tol 0.4 out/tchem_chapman_0001_gas.txt tchem_chapman_gas_saved.txt
