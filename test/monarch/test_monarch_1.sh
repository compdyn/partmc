#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out
# copy the compare file to the output directory
cp simple_comp.txt out/simple_comp.txt

((counter = 1))
while [ true ]
do
  echo Attempt $counter

if [[ $1 == "MPI" ]]; then
  exec_str="mpirun -v -np 2 ../../mock_monarch config_simple.json interface_simple.json out/simple"
else
  exec_str="../../mock_monarch config_simple.json interface_simple.json out/simple"
  #exec_str="nvprof ../../mock_monarch config_simple.json interface_simple.json out/simple"
  #exec_str="../../mock_monarch config_simple_cb05.json ../chemistry/cb05cl_ae5/cb05cl_ae5_init.json out/simple"
fi

  if ! $exec_str; then 
	  echo Failure "$counter"
	  if [ "$counter" -gt 1 ]
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

