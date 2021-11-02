#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../extract_env out/camp_0001
../../numeric_diff --by elem --rel-tol 1e-12 ref_camp_0001_env.txt out/camp_0001_env.txt
