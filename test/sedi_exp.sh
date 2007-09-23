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

echo ../src/partmc sedi_exp_mc.spec
../src/partmc sedi_exp_mc.spec
echo ../src/partmc sedi_exp_sect.spec
../src/partmc sedi_exp_sect.spec

echo ./sedi_exp_plot.py
./sedi_exp_plot.py
echo "Now view out/sedi_exp_{num,vol}_{lin,log}.pdf"
