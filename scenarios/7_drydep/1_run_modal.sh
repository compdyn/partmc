#!/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

mkdir -p out

../../build/partmc drydep_modal.spec

# Now run ./2_process_drydep_modal.sh to process the data