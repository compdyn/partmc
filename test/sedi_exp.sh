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
echo ../src/process_summary out/sedi_exp_mc_summary.d
../src/process_summary out/sedi_exp_mc_summary.d

echo ../src/partmc sedi_exp_sect.spec
../src/partmc sedi_exp_sect.spec
echo ../src/process_summary out/sedi_exp_sect_summary.d
../src/process_summary out/sedi_exp_sect_summary.d

echo Plotting number density
gnuplot -persist <<ENDNUM
set logscale x
set xlabel "radius (m)"
set ylabel "number density (#/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] "out/sedi_exp_mc_summary_aero.d" index 0 using 1:2 title "Monte Carlo (0 mins)"
replot "out/sedi_exp_mc_summary_aero.d" index 5 using 1:2 title "Monte Carlo (5 mins)"
replot "out/sedi_exp_mc_summary_aero.d" index 10 using 1:2 title "Monte Carlo (10 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 0 using 1:2 w l title "Sectional (0 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 5 using 1:2 w l title "Sectional (5 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 10 using 1:2 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "out/sedi_exp_plot_num.eps"
replot
ENDNUM
epstopdf out/sedi_exp_plot_num.eps

echo Plotting logarithmic number density
gnuplot -persist <<ENDNUMLOG
set logscale
set xlabel "radius (m)"
set ylabel "number density (#/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] [1e2:1e10] "out/sedi_exp_mc_summary_aero.d" index 0 using 1:2 title "Monte Carlo (0 mins)"
replot "out/sedi_exp_mc_summary_aero.d" index 5 using 1:2 title "Monte Carlo (5 mins)"
replot "out/sedi_exp_mc_summary_aero.d" index 10 using 1:2 title "Monte Carlo (10 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 0 using 1:2 w l title "Sectional (0 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 5 using 1:2 w l title "Sectional (5 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 10 using 1:2 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "out/sedi_exp_plot_num_log.eps"
replot
ENDNUMLOG
epstopdf out/sedi_exp_plot_num_log.eps

echo Plotting volume density
gnuplot -persist <<ENDVOL
set logscale x
set xlabel "radius (m)"
set ylabel "volume density (m^3/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] [0:8e-6]  "out/sedi_exp_mc_summary_aero.d" index 0 using 1:3 title "Monte Carlo (0 mins)"
replot "out/sedi_exp_mc_summary_aero.d" index 5 using 1:3 title "Monte Carlo (5 mins)"
replot "out/sedi_exp_mc_summary_aero.d" index 10 using 1:3 title "Monte Carlo (10 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 0 using 1:3 w l title "Sectional (0 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 5 using 1:3 w l title "Sectional (5 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 10 using 1:3 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "out/sedi_exp_plot_vol.eps"
replot
ENDVOL
epstopdf out/sedi_exp_plot_vol.eps

echo Plotting logarithmic volume density
gnuplot -persist <<ENDVOLLOG
set logscale
set xlabel "radius (m)"
set ylabel "volume density (m^3/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] [1e-20:1e3]  "out/sedi_exp_mc_summary_aero.d" index 0 using 1:3 title "Monte Carlo (0 mins)"
replot "out/sedi_exp_mc_summary_aero.d" index 5 using 1:3 title "Monte Carlo (5 mins)"
replot "out/sedi_exp_mc_summary_aero.d" index 10 using 1:3 title "Monte Carlo (10 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 0 using 1:3 w l title "Sectional (0 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 5 using 1:3 w l title "Sectional (5 mins)"
replot "out/sedi_exp_sect_summary_aero.d" index 10 using 1:3 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "out/sedi_exp_plot_vol_log.eps"
replot
ENDVOLLOG
epstopdf out/sedi_exp_plot_vol_log.eps
