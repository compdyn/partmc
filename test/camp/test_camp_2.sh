#!/bin/bash
  
# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_gas out/camp_0001
../../numeric_diff --by row --rel-tol 1e-4 ref_camp_0001_gas.txt out/camp_0001_gas.txt
