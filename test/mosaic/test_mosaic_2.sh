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

if ! ../../extract_env out/mosaic_0001 || \
   ! ../../numeric_diff --by elem --rel-tol 1e-12 ref_mosaic_0001_env.txt out/mosaic_0001_env.txt &> /dev/null; then
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
