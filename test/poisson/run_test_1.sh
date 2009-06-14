#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../test_poisson_sample 1 50 10000000 > out/poisson_1_approx.dat
../../test_poisson_sample 1 50 0        > out/poisson_1_exact.dat
../../numeric_diff out/poisson_1_approx.dat out/poisson_1_exact.dat 0 1e-3
exit $?
