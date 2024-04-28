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
echo $(uname)
echo $a $b
n_proc=$(( a < b ? a : b ))
echo $n_proc

parallel_type=dist

mpirun -np $n_proc ../../partmc run_part_parallel_${parallel_type}.spec
for f in out/parallel_${parallel_type}_0001_????_00000001.nc ; do
    echo "####################################################################"
    echo "####################################################################"
    echo $f
    prefix=${f/_00000001.nc/}
    ../../extract_aero_size --num --dmin 1e-10 --dmax 1e-4 --nbin 220 ${prefix}
    ../../extract_aero_size --mass --dmin 1e-10 --dmax 1e-4 --nbin 220 ${prefix}
    ../../extract_aero_time ${prefix}
    
    ../../numeric_diff --by col --rel-tol 0.4 out/sect_aero_size_num.txt ${prefix}_aero_size_num.txt
    ../../numeric_diff --by col --rel-tol 0.4 out/sect_aero_size_mass.txt ${prefix}_aero_size_mass.txt
    ../../numeric_diff --by col --rel-tol 0.2 out/sect_aero_time.txt ${prefix}_aero_time.txt
done

# #######################################################################
# #######################################################################
# Averaging

../../numeric_average out/parallel_${parallel_type}_aero_size_num.txt out/parallel_${parallel_type}_0001_????_aero_size_num.txt
../../numeric_average out/parallel_${parallel_type}_aero_size_mass.txt out/parallel_${parallel_type}_0001_????_aero_size_mass.txt
../../numeric_average out/parallel_${parallel_type}_aero_time.txt out/parallel_${parallel_type}_0001_????_aero_time.txt

../../numeric_diff --by col --rel-tol 0.2 out/sect_aero_size_num.txt out/parallel_${parallel_type}_aero_size_num.txt
../../numeric_diff --by col --rel-tol 0.2 out/sect_aero_size_mass.txt out/parallel_${parallel_type}_aero_size_mass.txt
../../numeric_diff --by col --rel-tol 0.1 out/sect_aero_time.txt out/parallel_${parallel_type}_aero_time.txt
