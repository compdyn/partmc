#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_aero_size_mass out/brownian_mc_0001.nc out/brownian_mc_size_mass.txt
../../extract_summary_aero_size_mass out/brownian_sect_0001.nc out/brownian_sect_size_mass.txt

../../numeric_diff out/brownian_mc_size_mass.txt out/brownian_sect_size_mass.txt 0 1.0 0 0 2 0
exit $?
