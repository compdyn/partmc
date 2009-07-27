#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "om-mpirun -v -np 10 ../../partmc run_part_parallel_dist_single.spec"
om-mpirun -v -np 10 ../../partmc run_part_parallel_dist_single.spec
echo "../../extract_aero_size_num 1e-10 1e-4 220 out/parallel_dist_single_0001_ out/parallel_dist_single_aero_size_num.txt"
../../extract_aero_size_num 1e-10 1e-4 220 out/parallel_dist_single_0001_ out/parallel_dist_single_aero_size_num.txt
echo "../../extract_aero_size_mass 1e-10 1e-4 220 out/parallel_dist_single_0001_ out/parallel_dist_single_aero_size_mass.txt"
../../extract_aero_size_mass 1e-10 1e-4 220 out/parallel_dist_single_0001_ out/parallel_dist_single_aero_size_mass.txt
echo "../../extract_aero_total out/parallel_dist_single_0001_ out/parallel_dist_single_aero_total.txt"
../../extract_aero_total out/parallel_dist_single_0001_ out/parallel_dist_single_aero_total.txt

echo "../../numeric_diff out/sect_aero_size_num.txt out/parallel_dist_single_aero_size_num.txt 0 0.3 0 0 2 0"
../../numeric_diff out/sect_aero_size_num.txt out/parallel_dist_single_aero_size_num.txt 0 0.3 0 0 2 0
R1=$?
echo "../../numeric_diff out/sect_aero_size_mass.txt out/parallel_dist_single_aero_size_mass.txt 0 0.3 0 0 2 0"
../../numeric_diff out/sect_aero_size_mass.txt out/parallel_dist_single_aero_size_mass.txt 0 0.3 0 0 2 0
R2=$?
echo "../../numeric_diff out/sect_aero_total.txt out/parallel_dist_single_aero_total.txt 0 0.3 0 0 2 2"
../../numeric_diff out/sect_aero_total.txt out/parallel_dist_single_aero_total.txt 0 0.3 0 0 2 2
R3=$?
echo "../../numeric_diff out/sect_aero_total.txt out/parallel_dist_single_aero_total.txt 0 0.3 0 0 3 3"
../../numeric_diff out/sect_aero_total.txt out/parallel_dist_single_aero_total.txt 0 0.3 0 0 3 3
R4=$?
exit $(( $R1 || $R2 || $R3 || $R4 ))
