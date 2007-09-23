#!/bin/sh

cat <<ENDINFO

Sedimentation Bidisperse Test-case
----------------------------------

Starting from many small particles and one large particle, we simulate
with the Monte Carlo code. For this case, we can write an ODE in one
variable (either number of small particles or volume of large
particle), which we solve as a comparison.

The output plots show both small-particle number and large-particle
volume. We disable doubling to avoid generating multiple large
particles, which would result in large-large coagulations and so the
ODE would no longer be valid. While this means that the simulation is
not all that interesting physically, it is an interesting test-case
that sectional codes cannot do well. The lack of resolution towards
the end of the Monte Carlo simulation is expected.

Because the plots made for this test-case are not number or volume
distributions, but rather numbers and volumes in specific bins, we
need to use a bit of a hack (the test_sedi_bidisperse_state_to_count
program) to process the data. This is not supposed to be a
general-purpose solution.

ENDINFO
sleep 1

echo ../src/partmc sedi_bidisperse_mc.spec
../src/partmc sedi_bidisperse_mc.spec
echo ./sedi_bidisperse_state_to_count
./sedi_bidisperse_state_to_count
echo ./sedi_bidisperse_ode
./sedi_bidisperse_ode

echo ./sedi_bidisperse_plot.py
./sedi_bidisperse_plot.py
echo "Now view out/sedi_bidisperse_*.pdf"
