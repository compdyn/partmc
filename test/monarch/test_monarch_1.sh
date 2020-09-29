#!/bin/bash
#SBATCH --job-name=cb05_test
#SBATCH --output=out_cb05.txt
#SBATCH --error=err_cb05_.txt
#SBATCH --ntasks=1
#   #SBATCH --cpus-per-task=4
#   #SBATCH --tasks-per-node=1
##mpi_%j.out ##mpi_%j.err

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
  exec_str="mpirun -v -np 2 ../../mock_monarch config_simple.json interface_simple.json out/simple" #local
  #exec_str="srun ../../mock_monarch config_simple.json interface_simple.json out/simple" #clusterelse
else
  exec_str="../../mock_monarch config_simple.json interface_simple.json out/simple"
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

