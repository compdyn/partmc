#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_time out/condense_0001

../../numeric_diff ref_condense_0001_aero_time.txt out/condense_0001_aero_time.txt 0 1e-6
