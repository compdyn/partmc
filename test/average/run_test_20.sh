#!/bin/bash

# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_aero_particle_mass out/average_compsizevol_0001_00000001.nc out/average_compsizevol_particles_unsorted.txt
sort out/average_compsizevol_particles_unsorted.txt > out/average_compsizevol_particles.txt
../../numeric_diff out/average_particles.txt out/average_compsizevol_particles.txt 0 0.1 0 0 3 3
