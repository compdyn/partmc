#!/bin/bash

# exit on error
# set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../test_aq_chem sundial_test_mech.txt sundial_test_init.dat sundial_test_spec_data.dat 1e10 > output.txt
