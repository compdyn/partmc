#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../poisson_sample 30 50 10000000 > poisson_30_approx.dat
../../poisson_sample 30 50 0        > poisson_30_sampled.dat
../../numeric_diff poisson_30_approx.dat poisson_30_sampled.dat 0 1e-3
exit $?
