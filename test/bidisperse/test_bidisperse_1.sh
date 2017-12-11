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

../../partmc run_part.spec
../../test_bidisperse_extract

../../test_bidisperse_ode

# extract size distributions for plotting
../../extract_aero_size --num --dmin 1e-5 --dmax 1e-3 --nbin 255 out/bidisperse_part_0001

if ! ../../numeric_diff --by col --rel-tol 0.3 out/bidisperse_ode_data.txt out/bidisperse_part_data.txt &> /dev/null; then
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
