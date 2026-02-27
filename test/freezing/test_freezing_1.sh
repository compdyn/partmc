#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

# make the output directory if it doesn't exist
mkdir -p out/ABIFM

# run PartMC freezing test
../../partmc run_part.ABIFM.spec

# extract frozen fraction results to a txt file
../../test_freezing_extract out/ABIFM/freezing_part

# calculate the freezing theoretical results
../../test_freezing_theoretical ABIFM

# compare the PartMC versus theoretical results
../../numeric_diff --rel-tol 0.02 out/ABIFM/freezing_part_frozen_fraction_ensemble_mean.txt out/freezing_theoretical_ABIFM_data.txt
