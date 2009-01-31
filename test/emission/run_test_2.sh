#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_sect.spec
../../partmc run_exact.spec
../../test_emission_process out/emission_sect_0001.nc out/emission_sect.txt
../../test_emission_process out/emission_exact_0001.nc out/emission_exact.txt
../../numeric_diff out/emission_sect.txt out/emission_exact.txt 0 1e-3
exit $?
