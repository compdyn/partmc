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

../../extract_env out/condense_0001_ out/condense_env.txt
../../numeric_diff true_env.txt out/condense_env.txt 0 1e-8
exit $?
