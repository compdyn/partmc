#!/bin/bash

# Sectional-vs-particle MOSAIC consistency test.
#
# Runs the 1D sectional solver (run_sect) with do_mosaic from the same
# monodisperse initial condition, scenario, and timing as the particle-resolved
# MOSAIC run (run_part.spec, run in test_mosaic_1.sh), then compares the gas
# evolution of the two. This exercises the bin <-> MOSAIC mappers
# (mosaic_from_partmc_sect / mosaic_to_partmc_sect) and the moving-center bin
# remap inside mosaic_timestep_sect against the particle-resolved reference.
#
# The particle run is done in test_mosaic_1.sh (produces out/mosaic_0001).

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_sect_mosaic.spec

# compare the sectional gas evolution against the particle-resolved run
../../extract_gas out/mosaic_0001
../../extract_gas out/mosaic_sect

# Skip column 1 (time) and column 2 (H2SO4): H2SO4 is driven to ~1e-11 ppb by
# condensation in both runs, so its residual differs at the noise level even
# though the absolute values are negligible. Most remaining gas species match
# to ~1e-7 (machine level).
#
# The loose 0.4 tolerance is set by the semivolatile SOA precursors (ALK1/OLE1),
# which differ by up to ~8% early in the run. This is NOT a sectional-vs-particle
# physics difference: it is Monte Carlo sampling noise in the n_part=40 reference,
# whose seed number is 9.75e7 vs the nominal 1e8 (a coarse-sampling artifact).
# That ~2.5% seed error is amplified by the nonlinear semivolatile partitioning
# early on, then decays as condensed mass grows. Raising n_part to a few thousand
# tightens the agreement to ~1%; the moving-center bin remap tracks the
# monodisperse growth essentially exactly. Tighten the tolerance here (or compare
# against a well-resolved reference file) if a stricter check is wanted.
../../numeric_diff --by col --min-col 3 --rel-tol 0.4 \
    out/mosaic_0001_gas.txt out/mosaic_sect_gas.txt
