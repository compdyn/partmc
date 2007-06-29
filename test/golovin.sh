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

if [[ ! -d out ]]; then mkdir out; fi

echo ../src/partmc golovin_mc.spec
../src/partmc golovin_mc.spec
echo ../src/process_summary out/golovin_mc_summary.d
../src/process_summary out/golovin_mc_summary.d

echo ../src/partmc golovin_exact.spec
../src/partmc golovin_exact.spec
echo ../src/process_summary out/golovin_exact_summary.d
../src/process_summary out/golovin_exact_summary.d

echo Plotting number density
gnuplot -persist <<ENDNUM
set logscale
set xlabel "radius (m)"
set ylabel "number density (#/m^3)
set title "Golovin kernel, exponential initial condition"
plot [1e-7:1e-3] [1e4:1e10] "out/golovin_mc_summary_num_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out/golovin_mc_summary_num_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out/golovin_mc_summary_num_avg.d" using 2:13 title "Monte Carlo (10 mins)"
replot "out/golovin_exact_summary_num_avg.d" using 2:3 w l title "Analytical (0 mins)"
replot "out/golovin_exact_summary_num_avg.d" using 2:8 w l title "Analytical (5 mins)"
replot "out/golovin_exact_summary_num_avg.d" using 2:13 w l title "Analytical (10 mins)"
set terminal postscript eps
set output "out/golovin_plot_num.eps"
replot
ENDNUM
epstopdf out/golovin_plot_num.eps

echo Plotting volume density
gnuplot -persist <<ENDVOL
set logscale
set xlabel "radius (m)"
set ylabel "volume density (m^3/m^3)"
set title "Golovin kernel, exponential initial condition"
set key left top
plot [1e-7:1e-3] [1e-14:1e-4] "out/golovin_mc_summary_vol_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out/golovin_mc_summary_vol_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out/golovin_mc_summary_vol_avg.d" using 2:13 title "Monte Carlo (10 mins)"
replot "out/golovin_exact_summary_vol_avg.d" using 2:3 w l title "Analytical (0 mins)"
replot "out/golovin_exact_summary_vol_avg.d" using 2:8 w l title "Analytical (5 mins)"
replot "out/golovin_exact_summary_vol_avg.d" using 2:13 w l title "Analytical (10 mins)"
set terminal postscript eps
set output "out/golovin_plot_vol.eps"
replot
ENDVOL
epstopdf out/golovin_plot_vol.eps
