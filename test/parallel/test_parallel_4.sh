#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

a=4
b=1
if [ "$(uname)" == "Darwin" ]; then
    b=$(getconf _NPROCESSORS_ONLN)
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    b=$(lscpu --parse=Core,Socket | grep --invert-match '^#' | sort -u | wc -l)
fi
n_proc=$(( a < b ? a : b ))
echo $n_proc

mpirun -v -np $n_proc ../../partmc run_part_parallel_dist_single.spec
../../extract_aero_size --num --dmin 1e-10 --dmax 1e-4 --nbin 220 out/parallel_dist_single_0001
../../extract_aero_size --mass --dmin 1e-10 --dmax 1e-4 --nbin 220 out/parallel_dist_single_0001
../../extract_aero_time out/parallel_dist_single_0001

../../numeric_diff --by col --rel-tol 0.2 out/sect_aero_size_num.txt out/parallel_dist_single_0001_aero_size_num.txt
../../numeric_diff --by col --rel-tol 0.2 out/sect_aero_size_mass.txt out/parallel_dist_single_0001_aero_size_mass.txt
../../numeric_diff --by col --rel-tol 0.1 out/sect_aero_time.txt out/parallel_dist_single_0001_aero_time.txt
