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
the of the Monte Carlo simulation is entirely expected.

Because the plots made for this test-case are not number or volume
distributions, but rather numbers and volumes in specific bins, we
need to use a bit of a hack (the test_sedi_bidisperse_state_to_count
program) to process the data. This is not supposed to be a
general-purpose solution.

ENDINFO

echo ./partbox test_sedi_bidisperse_mc.spec
./partbox test_sedi_bidisperse_mc.spec
echo ./test_sedi_bidisperse_state_to_count
./test_sedi_bidisperse_state_to_count

echo ./test_sedi_bidisperse_ode
./test_sedi_bidisperse_ode

echo "Plotting number density of small particles"
gnuplot -persist <<ENDNUM
set xlabel "time (s)"
set ylabel "number density of small particles (#/m^3)"
set title "Sedimentation kernel, bidisperse initial condition"
plot "counts_sedi_bidisperse_mc.d" using 1:2 title "Particle resolved"
replot "counts_sedi_bidisperse_ode.d" using 1:2 w l title "ODE"
set terminal postscript eps
set output "plot_sedi_bidisperse_num_small.eps"
replot
ENDNUM
epstopdf plot_sedi_bidisperse_num_small.eps

echo "Plotting number density of small particles (logarithmic)"
gnuplot -persist <<ENDNUM
set logscale y
set xlabel "time (s)"
set ylabel "number density of small particles (#/m^3)"
set title "Sedimentation kernel, bidisperse initial condition"
plot "counts_sedi_bidisperse_mc.d" using 1:2 title "Particle resolved"
replot "counts_sedi_bidisperse_ode.d" using 1:2 w l title "ODE"
set terminal postscript eps
set output "plot_sedi_bidisperse_num_small_log.eps"
replot
ENDNUM
epstopdf plot_sedi_bidisperse_num_small_log.eps

echo "Plotting volume density of big particle"
gnuplot -persist <<ENDVOL
set xlabel "time (s)"
set ylabel "volume density of big particles (m^3/m^3)"
set title "Sedimentation kernel, bidisperse initial condition"
set key right bottom
plot "counts_sedi_bidisperse_mc.d" using 1:3 title "Particle resolved"
replot "counts_sedi_bidisperse_ode.d" using 1:3 w l title "ODE"
set terminal postscript eps
set output "plot_sedi_bidisperse_vol_big.eps"
replot
ENDVOL
epstopdf plot_sedi_bidisperse_vol_big.eps
