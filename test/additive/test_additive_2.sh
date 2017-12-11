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

../../extract_aero_size --mass --dmin 1e-8 --dmax 1e-3 --nbin 160 out/additive_part_0001
../../extract_sectional_aero_size --mass out/additive_exact

if ! ../../numeric_diff --by col --rel-tol 0.2 out/additive_exact_aero_size_mass.txt out/additive_part_0001_aero_size_mass.txt &> /dev/null; then
	  echo Failure "$counter"
	  if [ "$counter" -gt 10 ]
	  then
		  echo FAIL
		  exit 1
	  fi
	  echo retrying...
	  ../../partmc run_part.spec
	  ../../partmc run_exact.spec
  else
	  echo PASS
	  exit 0
  fi
  ((counter++))
done
