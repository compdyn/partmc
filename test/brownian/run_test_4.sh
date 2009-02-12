#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

for f in out/brownian_mc_????.nc ; do
    f1=${f/_mc_/_mc_size_mass_}
    ../../extract_summary_aero_size_mass $f ${f1/.nc/.txt};
done
../../numeric_average out/brownian_mc_size_mass_average.txt out/brownian_mc_size_mass_????.txt

../../extract_summary_aero_size_mass out/brownian_sect_0001.nc out/brownian_sect_size_mass.txt

../../numeric_diff out/brownian_mc_size_mass_average.txt out/brownian_sect_size_mass.txt 0 0.7 0 0 2 0
exit $?
