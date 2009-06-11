#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../poisson_sample 30 50 10000000 > out/poisson_30_approx.dat
../../poisson_sample 30 50 0        > out/poisson_30_exact.dat
../../numeric_diff out/poisson_30_approx.dat out/poisson_30_exact.dat 0 1e-3
exit $?
