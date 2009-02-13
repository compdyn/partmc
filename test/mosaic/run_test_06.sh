#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_gas out/mosaic_0001.nc out/mosaic_gas_summary.txt
../../extract_state_gas out/mosaic_state_0001_ out/mosaic_gas_state.txt
../../numeric_diff out/mosaic_gas_summary.txt out/mosaic_gas_state.txt 0 1e-8 0 0 0 0
exit $?
