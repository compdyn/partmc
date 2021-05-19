#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../numeric_diff --by elem --abs-tol 1e-3 ref_groups.txt out/test_groups.txt
