#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../extract_aero_time out/tchem_cb05cl_ae5_0001
../../numeric_diff --rel-tol 0.4 out/tchem_cb05cl_ae5_0001_aero_time.txt tchem_cb05cl_ae5_aero_time_saved.txt
