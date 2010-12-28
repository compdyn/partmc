#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

time mpirun -v -np 1 ../../partmc run_part_serial.spec
time mpirun -v -np 2 ../../partmc run_part_parallel.spec

../../extract_aero_size_num 1e-8 1e-3 160 out/serial_0001_ out/serial_aero_size_num.txt
../../extract_aero_size_mass 1e-8 1e-3 160 out/serial_0001_ out/serial_aero_size_mass.txt
../../extract_aero_species out/serial_0001_ out/serial_aero_species.txt
../../extract_aero_total out/serial_0001_ out/serial_aero_total.txt
../../extract_env out/serial_0001_ out/serial_env.txt
../../extract_gas out/serial_0001_ out/serial_gas.txt

../../extract_aero_size_num 1e-8 1e-3 160 out/parallel_0001_ out/parallel_aero_size_num.txt
../../extract_aero_size_mass 1e-8 1e-3 160 out/parallel_0001_ out/parallel_aero_size_mass.txt
../../extract_aero_species out/parallel_0001_ out/parallel_aero_species.txt
../../extract_aero_total out/parallel_0001_ out/parallel_aero_total.txt
../../extract_env out/parallel_0001_ out/parallel_env.txt
../../extract_gas out/parallel_0001_ out/parallel_gas.txt
