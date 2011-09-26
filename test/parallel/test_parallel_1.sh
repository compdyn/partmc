#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

mpirun -v -np 1 ../../partmc run_sect.spec
../../extract_sectional_aero_size --num out/sect
../../extract_sectional_aero_size --mass out/sect
../../extract_sectional_aero_total out/sect_ out/sect_aero_total.txt

mpirun -v -np 1 ../../partmc run_part_serial.spec
../../extract_aero_size --num --dmin 1e-10 --dmax 1e-4 --nbin 220 out/serial_0001
../../extract_aero_size --mass --dmin 1e-10 --dmax 1e-4 --nbin 220 out/serial_0001
../../extract_aero_total out/serial_0001_ out/serial_aero_total.txt

../../numeric_diff out/sect_aero_size_num.txt out/serial_0001_size_num.txt 0 0.3 0 0 2 0
../../numeric_diff out/sect_aero_size_mass.txt out/serial_0001_size_mass.txt 0 0.3 0 0 2 0
../../numeric_diff out/sect_aero_total.txt out/serial_aero_total.txt 0 0.3 0 0 2 2
../../numeric_diff out/sect_aero_total.txt out/serial_aero_total.txt 0 0.3 0 0 3 3
