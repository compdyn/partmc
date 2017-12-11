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

../../extract_gas out/mosaic_0001
if ! ../../numeric_diff --by row --rel-tol 1e-4 ref_mosaic_0001_gas.txt out/mosaic_0001_gas.txt &> /dev/null; then
	  echo Failure "$counter"
	  if [ "$counter" -gt 10 ]
	  then
		  echo FAIL
		  exit 1
	  fi
	  echo retrying...
	  ../../partmc run_part.spec
  else
	  echo PASS
	  exit 0
  fi
  ((counter++))
done
