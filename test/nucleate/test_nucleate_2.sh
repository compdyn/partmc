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

if ! ../../extract_gas out/nucleate_part_0001 || \
   ! ../../numeric_diff --by col --rel-tol 0.05  out/nucleate_ode_gas.txt out/nucleate_part_0001_gas.txt; then
	  echo Failure "$counter"
	  if [ "$counter" -gt 10 ]
	  then
		  echo FAIL
		  exit 1
	  fi
	  echo retrying...
	  if ! ../../partmc run_part.spec; then continue; fi
  else
	  echo PASS
	  exit 0
  fi
  ((counter++))
done
