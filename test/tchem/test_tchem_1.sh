#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part_cb05cl_ae5_with_SIMPOL.spec
../../extract_gas out/tchem_cb05cl_ae5_0001
#../../numeric_diff --by col --rel-tol 0.4 out/tchem_cb05cl_ae5_0001_gas.txt tchem_cb05cl_ae5_gas_saved.txt
