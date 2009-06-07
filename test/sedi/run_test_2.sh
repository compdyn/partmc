#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_state_aero_size_mass 1e-8 1e-2 220 out/sedi_mc_0001_ out/sedi_mc_size_mass.txt"
../../extract_state_aero_size_mass 1e-8 1e-2 220 out/sedi_mc_0001_ out/sedi_mc_size_mass.txt
echo "../../extract_summary_aero_size_mass out/sedi_sect_0001.nc out/sedi_sect_size_mass.txt"
../../extract_summary_aero_size_mass out/sedi_sect_0001.nc out/sedi_sect_size_mass.txt

echo "../../numeric_diff out/sedi_mc_size_mass.txt out/sedi_sect_size_mass.txt 0 0.5 0 0 2 0"
../../numeric_diff out/sedi_mc_size_mass.txt out/sedi_sect_size_mass.txt 0 0.5 0 0 2 0
exit $?
