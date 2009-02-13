#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

../../extract_summary_env out/mosaic_0001.nc out/mosaic_env_summary.txt
../../extract_state_env out/mosaic_state_0001_ out/mosaic_env_state.txt
../../numeric_diff out/mosaic_env_summary.txt out/mosaic_env_state.txt 0 1e-8 0 0 0 0
exit $?
