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
echo ../src/process_out out/out_sedi_exp_mc.d
../src/process_out out/out_sedi_exp_mc.d

echo ../src/partmc sedi_exp_sect.spec
../src/partmc sedi_exp_sect.spec
echo ../src/process_out out/out_sedi_exp_sect.d
../src/process_out out/out_sedi_exp_sect.d

echo Plotting number density
gnuplot -persist <<ENDNUM
set logscale x
set xlabel "radius (m)"
set ylabel "number density (#/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] "out/out_sedi_exp_mc_num_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out/out_sedi_exp_mc_num_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out/out_sedi_exp_mc_num_avg.d" using 2:13 title "Monte Carlo (10 mins)"
replot "out/out_sedi_exp_sect_num_avg.d" using 2:3 w l title "Sectional (0 mins)"
replot "out/out_sedi_exp_sect_num_avg.d" using 2:8 w l title "Sectional (5 mins)"
replot "out/out_sedi_exp_sect_num_avg.d" using 2:13 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "plot_sedi_exp_num.eps"
replot
ENDNUM
epstopdf plot_sedi_exp_num.eps

echo Plotting logarithmic number density
gnuplot -persist <<ENDNUMLOG
set logscale
set xlabel "radius (m)"
set ylabel "number density (#/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] [1e2:1e10] "out/out_sedi_exp_mc_num_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out/out_sedi_exp_mc_num_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out/out_sedi_exp_mc_num_avg.d" using 2:13 title "Monte Carlo (10 mins)"
replot "out/out_sedi_exp_sect_num_avg.d" using 2:3 w l title "Sectional (0 mins)"
replot "out/out_sedi_exp_sect_num_avg.d" using 2:8 w l title "Sectional (5 mins)"
replot "out/out_sedi_exp_sect_num_avg.d" using 2:13 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "plot_sedi_exp_num_log.eps"
replot
ENDNUMLOG
epstopdf plot_sedi_exp_num_log.eps

echo Plotting volume density
gnuplot -persist <<ENDVOL
set logscale x
set xlabel "radius (m)"
set ylabel "volume density (m^3/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] [0:8e-6]  "out/out_sedi_exp_mc_vol_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out/out_sedi_exp_mc_vol_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out/out_sedi_exp_mc_vol_avg.d" using 2:13 title "Monte Carlo (10 mins)"
replot "out/out_sedi_exp_sect_vol_avg.d" using 2:3 w l title "Sectional (0 mins)"
replot "out/out_sedi_exp_sect_vol_avg.d" using 2:8 w l title "Sectional (5 mins)"
replot "out/out_sedi_exp_sect_vol_avg.d" using 2:13 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "plot_sedi_exp_vol.eps"
replot
ENDVOL
epstopdf plot_sedi_exp_vol.eps

echo Plotting logarithmic volume density
gnuplot -persist <<ENDVOLLOG
set logscale
set xlabel "radius (m)"
set ylabel "volume density (m^3/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] [1e-20:1e3]  "out/out_sedi_exp_mc_vol_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out/out_sedi_exp_mc_vol_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out/out_sedi_exp_mc_vol_avg.d" using 2:13 title "Monte Carlo (10 mins)"
replot "out/out_sedi_exp_sect_vol_avg.d" using 2:3 w l title "Sectional (0 mins)"
replot "out/out_sedi_exp_sect_vol_avg.d" using 2:8 w l title "Sectional (5 mins)"
replot "out/out_sedi_exp_sect_vol_avg.d" using 2:13 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "plot_sedi_exp_vol_log.eps"
replot
ENDVOLLOG
epstopdf plot_sedi_exp_vol_log.eps
