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

if ! ../../extract_aero_size --mass --dmin 1e-9 --dmax 1e-3 --nbin 160 out/loss_part_drydep_0001 || \
   ! ../../extract_sectional_aero_size --mass out/loss_exact_drydep || \
   ! ../../numeric_diff --by col --rel-tol 0.15 out/loss_exact_drydep_aero_size_mass.txt out/loss_part_drydep_0001_aero_size_mass.txt; then
	  echo Failure "$counter"
	  if [ "$counter" -gt 10 ]
	  then
		  echo FAIL
		  exit 1
	  fi
	  echo retrying...
	  if ! ../../partmc run_constant_part.spec; then continue; fi
	  if ! ../../partmc run_constant_exact.spec; then continue; fi
          if ! ../../partmc run_volume_part.spec; then continue; fi
          if ! ../../partmc run_volume_exact.spec; then continue; fi
          if ! ../../partmc run_drydep_part.spec; then continue; fi
          if ! ../../partmc run_drydep_exact.spec; then continue; fi
  else
	  echo PASS
	  exit 0
  fi
  ((counter++))
done
