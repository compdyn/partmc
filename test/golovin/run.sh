#!/bin/bash

cat <<ENDINFO

Golovin Test-case
-----------------

This tests the Monte Carlo coagulation code for the Golovin kernel
with an exponential initial condition for which an analytical solution
is known. The Monte Carlo code is looped and average to reduce
error. To improve the quality of the results either the number of
loops or the number of particles can be increased.

This test-case demonstrates generating exact solutions and using
repeated loops of the Monte Carlo code to reduce the error, and serves
as model verification against an analytical solution.

ENDINFO
sleep 1

echo ../../src/partmc run_mc.spec
../../src/partmc run_mc.spec
echo ../../src/partmc run_exact.spec
../../src/partmc run_exact.spec

echo ./plot.py
./plot.py
echo Now view out/golovin_num.pdf and out/golovin_vol.pdf
