#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

mpirun -v -np 10 ../../partmc run_part_parallel_dist_single.spec
../../extract_aero_size_num 1e-10 1e-4 220 out/parallel_dist_single_0001_ out/parallel_dist_single_aero_size_num.txt
../../extract_aero_size_mass 1e-10 1e-4 220 out/parallel_dist_single_0001_ out/parallel_dist_single_aero_size_mass.txt
../../extract_aero_total out/parallel_dist_single_0001_ out/parallel_dist_single_aero_total.txt

../../numeric_diff out/sect_aero_size_num.txt out/parallel_dist_single_aero_size_num.txt 0 0.3 0 0 2 0
R1=$?
../../numeric_diff out/sect_aero_size_mass.txt out/parallel_dist_single_aero_size_mass.txt 0 0.3 0 0 2 0
R2=$?
../../numeric_diff out/sect_aero_total.txt out/parallel_dist_single_aero_total.txt 0 0.3 0 0 2 2
R3=$?
../../numeric_diff out/sect_aero_total.txt out/parallel_dist_single_aero_total.txt 0 0.3 0 0 3 3
R4=$?
exit $(( $R1 || $R2 || $R3 || $R4 ))
