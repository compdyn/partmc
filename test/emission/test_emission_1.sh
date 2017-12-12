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

if ! ../../partmc run_part.spec || \
   ! ../../partmc run_exact.spec || \
   ! ../../extract_aero_size --num --dmin 1e-8 --dmax 1e-3 --nbin 160 out/emission_part_0001 || \
   ! ../../extract_sectional_aero_size --num out/emission_exact || \
   ! ../../numeric_diff --by col --rel-tol 0.05 out/emission_exact_aero_size_num.txt out/emission_part_0001_aero_size_num.txt; then
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
