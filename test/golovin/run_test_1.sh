#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_mc.spec
../../partmc run_exact.spec

../../extract_num_den out/golovin_mc_0001.nc out/golovin_mc_num.txt
../../extract_num_den out/golovin_exact_0001.nc out/golovin_exact_num.txt
../../extract_mass_den out/golovin_mc_0001.nc out/golovin_mc_mass.txt
../../extract_mass_den out/golovin_exact_0001.nc out/golovin_exact_mass.txt

../../numeric_diff out/golovin_mc_num.txt out/golovin_exact_num.txt 0 1e-2
num_result=$?
../../numeric_diff out/golovin_mc_mass.txt out/golovin_exact_mass.txt 0 1e-2
mass_result=$?
exit $(( num_result || mass_result ))
