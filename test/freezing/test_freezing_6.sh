#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}

# make the output directory if it doesn't exist
mkdir -p out/const_melt

# run PartMC freezing test
../../partmc run_part.const.melt.spec

# extract frozen fraction results to a txt file
../../test_freezing_extract out/const_melt/freezing_part

# compare the PartMC versus theoretical results
lines=($(<"out/const_melt/freezing_part_frozen_fraction_ensemble_mean.txt"))
zero="0.000000000000000E+000"

# Check the first and last rows are exactly zero
if [[ "${lines[0]}" != "$zero" || "${lines[-1]}" != "$zero" ]]; then
  echo "Assertion failed: The first or last row is not zero!"
  exit 1
fi

# Check that all other lines are non-zero
for i in "${lines[@]:1:${#lines[@]}-2}"; do
  if [[ "$i" == "$zero" ]]; then
    echo "Assertion failed: A middle row is zero!"
    exit 1
  fi
done
