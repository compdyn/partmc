#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_part.spec
../../extract_aero_time out/nucleate_part_0001
../../test_nucleate_ode
../../numeric_diff --by col --min-col 1 --max-col 3 --rel-tol 0.15 out/nucleate_ode_aero_time.txt out/nucleate_part_0001_aero_time.txt
