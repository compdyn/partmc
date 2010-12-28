#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

../../numeric_diff out/bidisperse_part_data.txt out/bidisperse_ode_data.txt 0 5e-2 0 0 2 2
exit $?
