#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../test_poisson_sample 1 50 10000000 > out/poisson_1_approx.dat
../../test_poisson_sample 1 50 0        > out/poisson_1_exact.dat
../../numeric_diff --by col --rel-tol 1e-3 out/poisson_1_exact.dat out/poisson_1_approx.dat
