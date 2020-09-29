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
cp monarch_mod37_comp.txt out/monarch_mod37_comp.txt

((counter = 1))
while [ true ]
do
  echo Attempt $counter

if [[ $1 == "MPI" ]]; then
  #exec_str="mpirun -v -np 2 ../../mock_monarch config_monarch_mod37.json interface_monarch_mod37.json out/monarch_mod37"
  exec_str="mpirun -v -np 2 ../../mock_monarch config_monarch_cb05_soa.json interface_monarch_cb05_soa.json out/monarch_cb05_soa"
else
  #exec_str="../../mock_monarch config_monarch_mod37.json interface_monarch_mod37.json out/monarch_mod37"
  exec_str="../../mock_monarch config_monarch_cb05_soa.json interface_monarch_cb05_soa.json out/monarch_cb05_soa"
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

