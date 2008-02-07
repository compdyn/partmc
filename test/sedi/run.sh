#!/bin/sh

cat <<ENDINFO

Sedimentation Exponential Test-case
-----------------------------------

The initial condition is a single species of particle with an
exponential distribution and the simulation using the sedimentation
kernel. No analytical solution is known for this case, so a sectional
code is use as verification.

This test-case demonstrates the use of the sectional model and
comparison of sectional runs with Monte Carlo runs.

ENDINFO
sleep 1

echo ../../src/partmc run_mc.spec
../../src/partmc run_mc.spec
echo ../../src/partmc run_sect.spec
../../src/partmc run_sect.spec

echo ./plot.py
./plot.py
echo "Now view out/sedi_*.pdf"
