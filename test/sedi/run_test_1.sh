#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_part.spec
../../partmc run_sect.spec

../../extract_aero_size_num 1e-8 1e-2 220 out/sedi_part_0001_ out/sedi_part_size_num.txt
../../extract_sectional_aero_size_num out/sedi_sect_ out/sedi_sect_size_num.txt

../../numeric_diff out/sedi_part_size_num.txt out/sedi_sect_size_num.txt 0 0.2 0 0 2 0
exit $?
