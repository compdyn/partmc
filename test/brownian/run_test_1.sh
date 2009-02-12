#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_mc.spec
../../partmc run_sect.spec

../../extract_summary_aero_size_num out/brownian_mc_0001.nc out/brownian_mc_size_num.txt
../../extract_summary_aero_size_num out/brownian_sect_0001.nc out/brownian_sect_size_num.txt

../../numeric_diff out/brownian_mc_size_num.txt out/brownian_sect_size_num.txt 0 0.5 0 0 2 0
exit $?
