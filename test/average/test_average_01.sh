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
   ! ../../bin_average_comp 1e-10 1e-4 24 wet out/average_0001_00000001.nc out/average_comp || \

   ! ../../extract_aero_size --num --dmin 1e-10 --dmax 1e-4 --nbin 24 out/average_0001 || \
   ! ../../extract_aero_size --num --dmin 1e-10 --dmax 1e-4 --nbin 24 out/average_comp_0001 || \
   ! ../../numeric_diff --by col --rel-tol 1e-12 out/average_0001_aero_size_num.txt out/average_comp_0001_aero_size_num.txt; then
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
