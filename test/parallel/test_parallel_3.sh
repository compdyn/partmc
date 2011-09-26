#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

parallel_type=mix

mpirun -v -np 4 ../../partmc run_part_parallel_${parallel_type}.spec
for f in out/parallel_${parallel_type}_0001_????_00000001.nc ; do
    echo "####################################################################"
    echo "####################################################################"
    echo $f
    prefix=${f/_00000001.nc/}
    ../../extract_aero_size_num 1e-10 1e-4 220 ${prefix}
    ../../extract_aero_size_mass 1e-10 1e-4 220 ${prefix}
    ../../extract_aero_time ${prefix}
    
    ../../numeric_diff out/sect_aero_size_num.txt ${prefix}_aero_size_num.txt 0 0.3 0 0 2 0
    ../../numeric_diff out/sect_aero_size_mass.txt ${prefix}_aero_size_mass.txt 0 0.5 0 0 2 0
    ../../numeric_diff out/sect_aero_time.txt ${prefix}_aero_time.txt 0 0.1 0 0 2 2
    ../../numeric_diff out/sect_aero_time.txt ${prefix}_aero_time.txt 0 0.3 0 0 3 3
done

# #######################################################################
# #######################################################################
# Averaging

../../numeric_average out/parallel_${parallel_type}_aero_size_num.txt out/parallel_${parallel_type}_0001_????_aero_size_num.txt
../../numeric_average out/parallel_${parallel_type}_aero_size_mass.txt out/parallel_${parallel_type}_0001_????_aero_size_mass.txt
../../numeric_average out/parallel_${parallel_type}_aero_time.txt out/parallel_${parallel_type}_0001_????_aero_time.txt

../../numeric_diff out/sect_aero_size_num.txt out/parallel_${parallel_type}_aero_size_num.txt 0 0.1 0 0 2 0
../../numeric_diff out/sect_aero_size_mass.txt out/parallel_${parallel_type}_aero_size_mass.txt 0 0.3 0 0 2 0
../../numeric_diff out/sect_aero_time.txt out/parallel_${parallel_type}_aero_time.txt 0 0.05 0 0 2 2
../../numeric_diff out/sect_aero_time.txt out/parallel_${parallel_type}_aero_time.txt 0 0.1 0 0 3 3
