#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../numeric_diff out/bidisperse_mc_data.txt out/bidisperse_ode_data.txt 0 3e-2 0 0 3 3
exit $?
