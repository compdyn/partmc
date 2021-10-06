#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc camp.spec

../../extract_aero_time out/camp_0001
../../numeric_diff --by col --rel-tol 0.5 ref_camp_0001_aero_time.txt out/camp_0001_aero_time.txt
