#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_particle_mass out/average_0001_00000001.nc out/average_particles_unsorted.txt
sort out/average_particles_unsorted.txt > out/average_particles.txt
../../extract_aero_particle_mass out/average_comp_0001_00000001.nc out/average_comp_particles_unsorted.txt
sort out/average_comp_particles_unsorted.txt > out/average_comp_particles.txt
../../numeric_diff out/average_particles.txt out/average_comp_particles.txt 0 1e-12 0 0 3 3
