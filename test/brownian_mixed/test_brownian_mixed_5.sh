#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

# per-species (C) mass distribution, averaged over the Monte Carlo repeats
for f in out/brownian_mixed_part_????_00000001.nc ; do
    ../../extract_aero_size --mass --species C --dmin 1e-10 --dmax 1e-4 \
        --nbin 220 ${f/_00000001.nc/}
done
../../numeric_average out/brownian_mixed_part_aero_size_mass_C_average.txt \
    out/brownian_mixed_part_????_aero_size_mass_C.txt
../../extract_sectional_aero_size --mass --species C out/brownian_mixed_sect
../../numeric_diff --by col --rel-tol 0.2 \
    out/brownian_mixed_sect_aero_size_mass_C.txt \
    out/brownian_mixed_part_aero_size_mass_C_average.txt
