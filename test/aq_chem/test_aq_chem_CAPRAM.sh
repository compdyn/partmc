#!/bin/bash

# exit on error
# set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../test_aq_chem capram24_red+.txt aq_spec_init_from_pmc.dat aq_spec_data.dat 100.0 > output_CAPRAM.txt
