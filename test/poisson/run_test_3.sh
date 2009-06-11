#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../poisson_sample 10 50 10000000 > out/poisson_10_approx.dat
../../poisson_sample 10 50 0        > out/poisson_10_exact.dat
../../numeric_diff out/poisson_10_approx.dat out/poisson_10_exact.dat 0 1e-3
exit $?
