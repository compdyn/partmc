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

if ! ../../extract_aero_particles out/average_compsizevol_0001_00000001.nc || \
   ! sort out/average_compsizevol_0001_00000001_aero_particles.txt | grep -v '^[[:blank:]]*$' > out/average_compsizevol_aero_particles_sorted.txt || \
   ! ../../numeric_diff --by elem --min-col 3 --max-col 4 --rel-tol 3 out/average_aero_particles_sorted.txt out/average_compsizevol_aero_particles_sorted.txt; then
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
