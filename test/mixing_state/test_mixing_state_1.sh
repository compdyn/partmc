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
../../test_mixing_state_process

../../numeric_diff --by elem --abs-tol 1e-3 ref_all_species.txt out/test_all_species.txt
