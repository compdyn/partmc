#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_aero_size_mass out/sedi_mc_0001.nc out/sedi_mc_size_mass.txt
../../extract_summary_aero_size_mass out/sedi_sect_0001.nc out/sedi_sect_size_mass.txt

../../numeric_diff out/sedi_mc_size_mass.txt out/sedi_sect_size_mass.txt 0 0.5 0 0 2 0
exit $?
