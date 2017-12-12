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

if ! ../../extract_aero_time out/average_sizevol_0001 || \
   ! ../../numeric_diff --rel-tol 1e-10 out/average_0001_aero_time.txt out/average_sizevol_0001_aero_time.txt; then
	  echo Failure "$counter"
	  if [ "$counter" -gt 10 ]
	  then
		  echo FAIL
		  exit 1
	  fi
	  echo retrying...
	  if ! ../../partmc run_part.spec; then continue; fi
	  if ! ../../bin_average_comp 1e-10 1e-4 24 wet out/average_0001_00000001.nc out/average_comp; then continue; fi
  else
	  echo PASS
	  exit 0
  fi
  ((counter++))
done
