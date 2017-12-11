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

mpirun -v -np 1 ../../partmc run_sect.spec
../../extract_sectional_aero_size --num out/sect
../../extract_sectional_aero_size --mass out/sect
../../extract_sectional_aero_time out/sect

mpirun -v -np 1 ../../partmc run_part_serial.spec
../../extract_aero_size --num --dmin 1e-10 --dmax 1e-4 --nbin 220 out/serial_0001
../../extract_aero_size --mass --dmin 1e-10 --dmax 1e-4 --nbin 220 out/serial_0001
../../extract_aero_time out/serial_0001

if [ ! ../../numeric_diff --by col --rel-tol 0.2 out/sect_aero_size_num.txt out/serial_0001_aero_size_num.txt &> /dev/null ] ||
[ ! ../../numeric_diff --by col --rel-tol 0.2 out/sect_aero_size_mass.txt out/serial_0001_aero_size_mass.txt  &> /dev/null ] ||
[ ! ../../numeric_diff --by col --rel-tol 0.1 out/sect_aero_time.txt out/serial_0001_aero_time.txt &> /dev/null ] then
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
