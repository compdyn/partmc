#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_sect.spec
../../partmc run_exact.spec
../../extract_summary num out/emission_sect_0001.nc out/emission_sect.txt
../../extract_summary num out/emission_exact_0001.nc out/emission_exact.txt
../../numeric_diff out/emission_sect.txt out/emission_exact.txt 0 1e-3
exit $?
