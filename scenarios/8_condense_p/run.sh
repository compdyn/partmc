#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../build/partmc run_part.spec

../../build/extract_env out_we/out_20/condense_0001
../../build/extract_aero_time out_we/out_20/condense_0001
