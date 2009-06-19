#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../extract_aero_particle_mass out/mosaic_0001_00000001.nc out/mosaic_initial_particle_masses.txt"
../../extract_aero_particle_mass out/mosaic_0001_00000001.nc out/mosaic_initial_particle_masses.txt
echo "../../numeric_diff true_initial_particle_masses.txt out/mosaic_initial_particle_masses.txt 0 1e-8 0 0 3 0"
../../numeric_diff true_initial_particle_masses.txt out/mosaic_initial_particle_masses.txt 0 1e-8 0 0 3 0
exit $?
