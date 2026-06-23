#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

# sectional run with TChem; the particle run is done in test_tchem_1.sh
../../partmc run_sect_cb05cl_ae5_with_SIMPOL.spec

# compare the sectional gas evolution against the particle-resolved run
../../extract_gas out/tchem_cb05cl_ae5_sect
../../extract_gas out/tchem_cb05cl_ae5_0001
../../numeric_diff --by col --rel-tol 0.4 out/tchem_cb05cl_ae5_0001_gas.txt out/tchem_cb05cl_ae5_sect_gas.txt
