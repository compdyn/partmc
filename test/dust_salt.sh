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
sleep 1

echo ./partmc test_dust_salt_part1.spec
./partmc test_dust_salt_part1.spec
echo ./process_out out_dust_salt_part1.d
./process_out out_dust_salt_part1.d

echo ./partmc test_dust_salt_part2.spec
./partmc test_dust_salt_part2.spec
echo ./process_out out_dust_salt_part2.d
./process_out out_dust_salt_part2.d

echo ./process_state state_dust_salt_part1_0001_00000000.d
./process_state state_dust_salt_part1_0001_00000000.d
echo ./process_state state_dust_salt_part2_0001_00000800.d
./process_state state_dust_salt_part2_0001_00000800.d
echo ./process_state state_dust_salt_part2_0001_00001000.d
./process_state state_dust_salt_part2_0001_00001000.d

echo Plotting number density
gnuplot -persist <<ENDNUM
set logscale
set xlabel "radius (m)"
set ylabel "number density (#/m^3)
set title "Testcase dust and seasalt"
plot [1e-8:1e-3] [1e3:1e10] "out_dust_salt_part1_num_avg.d" using 2:3 w l title "Monte Carlo (0 s)"
replot "out_dust_salt_part1_num_avg.d" using 2:7 w l title "Monte Carlo (400 s)"
replot "out_dust_salt_part1_num_avg.d" using 2:11 w l title "Monte Carlo (800 s)"
replot "out_dust_salt_part2_num_avg.d" using 2:13 w l title "Monte Carlo (1000 s)"
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
set key top left
plot [1e-8:1e-3] [1e-14:1e-4] "out_dust_salt_part1_vol_avg.d" index 0 using 2:3 w l title "Monte Carlo (0 s)"
replot "out_dust_salt_part1_vol_avg.d" index 0 using 2:7 w l title "Monte Carlo (400 s)"
replot "out_dust_salt_part1_vol_avg.d" index 0 using 2:11 w l title "Monte Carlo (800 s)"
replot "out_dust_salt_part2_vol_avg.d" index 0 using 2:13 w l title "Monte Carlo (1000 s)"
set terminal postscript eps
set output "plot_dust_salt_vol.eps"
replot
ENDVOL
epstopdf plot_dust_salt_vol.eps

echo Plotting composition
gnuplot -persist <<ENDCOM
set logscale
set xlabel "radius (m)"
set ylabel "number density (#/m^3)"
set title "Testcase dust and seasalt"
plot [1e-8:1e-3] [1e3:1e10] "state_dust_salt_part1_0001_00000000_moments_comp.d" using 2:3 w l title "species 1 (0 s)"
replot "state_dust_salt_part1_0001_00000000_moments_comp.d" using 2:4 w l title "species 2 (0 s)"
replot "state_dust_salt_part2_0001_00000800_moments_comp.d" using 2:3 w l title "species 1 (800 s)"
replot "state_dust_salt_part2_0001_00000800_moments_comp.d" using 2:4 w l title "species 2 (800 s)"
replot "state_dust_salt_part2_0001_00000800_moments_comp.d" using 2:5 w l title "mixed (800 s)"
replot "state_dust_salt_part2_0001_00001000_moments_comp.d" using 2:3 w l title "species 1 (1000 s)"
replot "state_dust_salt_part2_0001_00001000_moments_comp.d" using 2:4 w l title "species 2 (1000 s)"
replot "state_dust_salt_part2_0001_00001000_moments_comp.d" using 2:5 w l title "mixed (1000 s)"
set terminal postscript eps
set output "plot_dust_salt_comp.eps"
replot
ENDCOM
epstopdf plot_dust_salt_comp.eps
