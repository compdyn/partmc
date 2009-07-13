#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

parallel_type=central

RT=0

echo "om-mpirun -v -np 10 ../../partmc run_part_parallel_${parallel_type}.spec"
om-mpirun -v -np 10 ../../partmc run_part_parallel_${parallel_type}.spec
for f in out/parallel_${parallel_type}_0001_????_00000001.nc ; do
    prefix=${f/00000001.nc/}
    echo "../../extract_aero_size_num 1e-10 1e-4 220 ${prefix} ${prefix}aero_size_num.txt"
    ../../extract_aero_size_num 1e-10 1e-4 220 ${prefix} ${prefix}aero_size_num.txt
    echo "../../extract_aero_size_mass 1e-10 1e-4 220 ${prefix} ${prefix}aero_size_mass.txt"
    ../../extract_aero_size_mass 1e-10 1e-4 220 ${prefix} ${prefix}aero_size_mass.txt
    echo "../../extract_aero_total ${prefix} ${prefix}aero_total.txt"
    ../../extract_aero_total ${prefix} ${prefix}aero_total.txt
    
    echo "../../numeric_diff out/sect_aero_size_num.txt ${prefix}aero_size_num.txt 0 0.3 0 0 2 0"
    ../../numeric_diff out/sect_aero_size_num.txt ${prefix}aero_size_num.txt 0 0.3 0 0 2 0
    R1=$?
    echo "../../numeric_diff out/sect_aero_size_mass.txt ${prefix}aero_size_mass.txt 0 0.3 0 0 2 0"
    ../../numeric_diff out/sect_aero_size_mass.txt ${prefix}aero_size_mass.txt 0 0.3 0 0 2 0
    R2=$?
    echo "../../numeric_diff out/sect_aero_total.txt ${prefix}aero_total.txt 0 0.3 0 0 2 2"
    ../../numeric_diff out/sect_aero_total.txt ${prefix}aero_total.txt 0 0.3 0 0 2 2
    R3=$?
    echo "../../numeric_diff out/sect_aero_total.txt ${prefix}aero_total.txt 0 0.3 0 0 3 3"
    ../../numeric_diff out/sect_aero_total.txt ${prefix}aero_total.txt 0 0.3 0 0 3 3
    R4=$?
    RT=$(( $RT || $R1 || $R2 || $R3 || $R4 ))
done

echo "../../numeric_average out/parallel_${parallel_type}_aero_size_num.txt out/parallel_${parallel_type}_0001_????_aero_size_num.txt"
../../numeric_average out/parallel_${parallel_type}_aero_size_num.txt out/parallel_${parallel_type}_0001_????_aero_size_num.txt
echo "../../numeric_average out/parallel_${parallel_type}_aero_size_mass.txt out/parallel_${parallel_type}_0001_????_aero_size_mass.txt"
../../numeric_average out/parallel_${parallel_type}_aero_size_mass.txt out/parallel_${parallel_type}_0001_????_aero_size_mass.txt
echo "../../numeric_average out/parallel_${parallel_type}_aero_total.txt out/parallel_${parallel_type}_0001_????_aero_total.txt"
../../numeric_average out/parallel_${parallel_type}_aero_total.txt out/parallel_${parallel_type}_0001_????_aero_total.txt

echo "../../numeric_diff out/sect_aero_size_num.txt out/parallel_${parallel_type}_aero_size_num.txt 0 0.3 0 0 2 0"
../../numeric_diff out/sect_aero_size_num.txt out/parallel_${parallel_type}_aero_size_num.txt 0 0.3 0 0 2 0
R1=$?
echo "../../numeric_diff out/sect_aero_size_mass.txt out/parallel_${parallel_type}_aero_size_mass.txt 0 0.3 0 0 2 0"
../../numeric_diff out/sect_aero_size_mass.txt out/parallel_${parallel_type}_aero_size_mass.txt 0 0.3 0 0 2 0
R2=$?
echo "../../numeric_diff out/sect_aero_total.txt out/parallel_${parallel_type}_aero_total.txt 0 0.3 0 0 2 2"
../../numeric_diff out/sect_aero_total.txt out/parallel_${parallel_type}_aero_total.txt 0 0.3 0 0 2 2
R3=$?
echo "../../numeric_diff out/sect_aero_total.txt out/parallel_${parallel_type}_aero_total.txt 0 0.3 0 0 3 3"
../../numeric_diff out/sect_aero_total.txt out/parallel_${parallel_type}_aero_total.txt 0 0.3 0 0 3 3
R4=$?
exit $(( $RT || $R1 || $R2 || $R3 || $R4 ))
