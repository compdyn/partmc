#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../partmc run_part.spec

../../extract_aero_total out/nucleate_part_0001_ out/aero_total.txt

../../test_nucleate_ode

../../numeric_diff out/aero_total.txt out/nucleate_ode_aero_number.txt 0 5e-2 0 0 2 0
exit $?
