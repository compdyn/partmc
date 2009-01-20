#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../poisson_sample 20 50 10000000 > poisson_20_approx.dat
../../poisson_sample 20 50 0        > poisson_20_sampled.dat
../../numeric_diff poisson_20_approx.dat poisson_20_sampled.dat 0 1e-3
exit $?
