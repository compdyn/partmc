#!/bin/sh

cat <<ENDINFO

Dust/Sea-salt Test-case
-----------------------

The initial condition is two log-normal distributions, one of small
dust particles and one of larger salt particles. The run in done in
two phases:

1. Until 800s condensation is active and dominates the evolution. This
   is very expensive, so only a few particles are used. The internal
   state is output at the end of this run.

2. The simulation is restarted from the saved state but with
   condensation turned off, as coagulation now dominates. This is run
   with many particles as it is compartively cheap and needs lots of
   particles to obtain acceptable resolution.

There is no exact solution for this simulation, and 1D sectional codes
cannot determine the particle mixing state. The results are thus
simply plotted on their own, without a comparison.

This test-case demonstates condensation, restarting from saved state,
and processing state for multi-species simulations.

ENDINFO

echo ./partbox test_dust_salt_part1.spec
#./partbox test_dust_salt_part1.spec
echo ./process_out out_dust_salt_part1.d
./process_out out_dust_salt_part1.d

echo ./partbox test_dust_salt_part2.spec
./partbox test_dust_salt_part2.spec
echo ./process_out out_dust_salt_part2.d
./process_out out_dust_salt_part2.d

echo Plotting number density
gnuplot -persist <<ENDNUM
set logscale
set xlabel "radius (m)"
set ylabel "number density (#/m^3)
set title "Testcase dust and seasalt"
plot [1e-7:1e-3] [1e4:1e10] "out_dust_salt_part1_num_avg.d" using 2:3 title "Particle resolved (0 s)"
replot "out_dust_salt_part1_num_avg.d" using 2:7 title "Particle resolved (400 s)"
replot "out_dust_salt_part1_num_avg.d" using 2:11 title "Particle resolved (800 s)"
replot "out_dust_salt_part2_num_avg.d" using 2:13 title "Particle resolved (1000 s)"
set terminal postscript eps
set output "plot_dust_salt_num.eps"
replot
ENDNUM
epstopdf plot_dust_salt_num.eps

echo Plotting volume density
gnuplot -persist <<ENDVOL
set logscale
set xlabel "radius (m)"
set ylabel "volume density (m^3/m^3)"
set title "Testcase dust and seasalt"
set key left top
plot [1e-7:1e-3] [1e-14:1e-4] "out_dust_salt_part1_vol_avg.d" using 2:3 title "Particle resolved (0 s)"
replot "out_dust_salt_part1_vol_avg.d" using 2:7 title "Particle resolved (400 s)"
replot "out_dust_salt_part1_vol_avg.d" using 2:11 title "Particle resolved (800 s)"
replot "out_dust_salt_part2_vol_avg.d" using 2:13 title "Particle resolved (1000 s)"
set terminal postscript eps
set output "plot_dust_salt_vol.eps"
replot
ENDVOL
epstopdf plot_dust_salt_vol.eps
