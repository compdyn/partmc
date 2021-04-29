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

if [[ $1 = "MPI" ]]; then
  exec_str="mpirun -v -np 2 ../../test_rxn_wennberg_no_ro2"
else
  exec_str="../../test_rxn_wennberg_no_ro2"
fi

if ! $exec_str; then
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
