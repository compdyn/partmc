#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../build/partmc run_sedi_part.spec
../../build/partmc run_sedi_sect.spec

../../build/extract_aero_size --num --dmin 1e-7 --dmax 1 --nbin 100 out/sedi_part_0001
../../build/extract_sectional_aero_size --num out/sedi_sect

../../build/extract_aero_size --mass --dmin 1e-7 --dmax 1 --nbin 100 out/sedi_part_0001
../../build/extract_sectional_aero_size --mass out/sedi_sect
