#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

# compare the sectional total aerosol mass evolution against the
# particle-resolved run (sectional run is done in test_tchem_3.sh, particle
# run in test_tchem_1.sh)
../../extract_sectional_aero_time out/tchem_cb05cl_ae5_sect
../../extract_aero_time out/tchem_cb05cl_ae5_0001
../../numeric_diff --by col --rel-tol 0.4 out/tchem_cb05cl_ae5_0001_aero_time.txt out/tchem_cb05cl_ae5_sect_aero_time.txt
