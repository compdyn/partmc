#!/bin/bash

# exit on error
set -e
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory is it doesn't exist
mkdir -p out

# Run the simulation
../../../../camp_box_model config.json out/simple_mech > box_model.log
