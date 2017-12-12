#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

((counter = 1))
while [ true ]
do
  echo Attempt $counter

if ! ../../test_poisson_sample 10 50 10000000 > out/poisson_10_approx.dat || \
   ! ../../test_poisson_sample 10 50 0        > out/poisson_10_exact.dat || \
   ! ../../numeric_diff --by col --rel-tol 2e-3 out/poisson_10_exact.dat out/poisson_10_approx.dat; then
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
