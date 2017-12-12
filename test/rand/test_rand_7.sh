#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

((counter = 1))
while [ true ]
do
  echo Attempt $counter

if ! ../../test_binomial_sample 10 0.7 10000000 > out/binomial_2_approx.dat || \
   ! ../../test_binomial_sample 10 0.7 0        > out/binomial_2_exact.dat || \
   ! ../../numeric_diff --by col --rel-tol 2e-3 out/binomial_2_exact.dat out/binomial_2_approx.dat; then
	  echo Failure "$counter"
	  if [ "$counter" -gt 10 ]
	  then
		  echo FAIL
		  exit 1
	  fi
	  echo retrying...
  else
	  echo PASS
	  exit 0
  fi
  ((counter++))
done
