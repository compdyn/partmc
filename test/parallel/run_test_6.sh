#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

time om-mpirun -v -np 1 ../../partmc run_sect.spec
../../extract_sectional_aero_size_num out/sect_ out/sect_aero_size_num.txt
../../extract_sectional_aero_total out/sect_ out/sect_aero_total.txt

time om-mpirun -v -np 2 ../../partmc run_part_parallel.spec
../../extract_aero_size_num 1e-10 1e-4 220 out/parallel_0001_0001_ out/parallel_0001_aero_size_num.txt
../../extract_aero_size_num 1e-10 1e-4 220 out/parallel_0001_0002_ out/parallel_0002_aero_size_num.txt
../../extract_aero_total out/parallel_0001_0001_ out/parallel_0001_aero_total.txt
../../extract_aero_total out/parallel_0001_0002_ out/parallel_0002_aero_total.txt

time om-mpirun -v -np 1 ../../partmc run_part_serial.spec
../../extract_aero_size_num 1e-10 1e-4 220 out/serial_0001_ out/serial_aero_size_num.txt
../../extract_aero_total out/serial_0001_ out/serial_aero_total.txt
