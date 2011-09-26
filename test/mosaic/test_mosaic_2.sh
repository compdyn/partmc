#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_env out/mosaic_0001
../../numeric_diff --by elem --rel-tol 1e-12 ref_mosaic_0001_env.txt out/mosaic_0001_env.txt
