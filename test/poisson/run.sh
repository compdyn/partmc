#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

cat README

echo ../../poisson_sample 1 50 10000000
../../poisson_sample 1 50 10000000 > poisson_1.d
echo ../../poisson_sample 4 50 10000000
../../poisson_sample 4 50 10000000 > poisson_4.d
echo ../../poisson_sample 10 50 10000000
../../poisson_sample 10 50 10000000 > poisson_10.d
echo ../../poisson_sample 20 50 10000000
../../poisson_sample 20 50 10000000 > poisson_20.d
echo ../../poisson_sample 30 50 10000000
../../poisson_sample 30 50 10000000 > poisson_30.d

echo 'To plot the data run ./plot.pdf and view the generated PDFs'
