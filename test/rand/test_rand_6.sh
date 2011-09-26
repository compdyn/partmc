#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../test_binomial_sample 10 0.3 10000000 > out/binomial_1_approx.dat
../../test_binomial_sample 10 0.3 0        > out/binomial_1_exact.dat
../../numeric_diff --by col --rel-tol 1e-3 out/binomial_1_exact.dat out/binomial_1_approx.dat
