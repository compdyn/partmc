#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../build/partmc run_part.spec

../../build/extract_env out/condense_0001
../../build/numeric_diff --by col --rel-tol 1e-4 ref_condense_0001_env.txt out/condense_0001_env.txt
