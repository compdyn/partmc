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

../../extract_aero_size --mass --dmin 1e-7 --dmax 1 --nbin 100 out/sedi_part_0001
../../extract_sectional_aero_size --mass out/sedi_sect

if ! ../../numeric_diff --by col --rel-tol 0.7 out/sedi_sect_aero_size_mass.txt out/sedi_part_0001_aero_size_mass.txt &> /dev/null; then
	  echo Failure "$counter"
	  if [ "$counter" -gt 10 ]
	  then
		  echo FAIL
		  exit 1
	  fi
	  echo retrying...
	  ../../partmc run_part.spec
	  ../../partmc run_sect.spec
  else
	  echo PASS
	  exit 0
  fi
  ((counter++))
done
