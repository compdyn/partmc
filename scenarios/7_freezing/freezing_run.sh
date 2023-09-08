#!/bin/bash

caseName=git_trans_t1
# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

sed -i "/output_prefix /coutput_prefix output/${caseName}/freezing_part # prefix of output files" run_part.spec

# make the output directory if it doesn't exist
mkdir -p output/$caseName
cp -p run_part.spec output/$caseName
cp -p *.dat output/$caseName

../../build/partmc run_part.spec
#for i in out/freezing_part_*.nc; do ../../build/extract_aero_particles $i; done;
