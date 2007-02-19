#!/bin/sh

echo ./partbox test_sedi_exp_mc.spec
./partbox test_sedi_exp_mc.spec
echo ./process_out out_sedi_exp_mc.d
./process_out out_sedi_exp_mc.d

echo ./partbox test_sedi_exp_sect.spec
./partbox test_sedi_exp_sect.spec
echo ./process_out out_sedi_exp_sect.d
./process_out out_sedi_exp_sect.d

echo Plotting number density
gnuplot -persist <<ENDNUM
set logscale x
set xlabel "radius (m)"
set ylabel "number density (#/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] "out_sedi_exp_mc_num_avg.d" using 2:3 title "Particle resolved (0 mins)"
replot "out_sedi_exp_mc_num_avg.d" using 2:8 title "Particle resolved (5 mins)"
replot "out_sedi_exp_mc_num_avg.d" using 2:13 title "Particle resolved (10 mins)"
replot "out_sedi_exp_sect_num_avg.d" using 2:3 w l title "Sectional (0 mins)"
replot "out_sedi_exp_sect_num_avg.d" using 2:8 w l title "Sectional (5 mins)"
replot "out_sedi_exp_sect_num_avg.d" using 2:13 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "plot_sedi_exp_num.eps"
replot
ENDNUM
epstopdf plot_sedi_exp_num.eps

echo Plotting number density
gnuplot -persist <<ENDNUMLOG
set logscale
set xlabel "radius (m)"
set ylabel "number density (#/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] [1e2:1e10] "out_sedi_exp_mc_num_avg.d" using 2:3 title "Particle resolved (0 mins)"
replot "out_sedi_exp_mc_num_avg.d" using 2:8 title "Particle resolved (5 mins)"
replot "out_sedi_exp_mc_num_avg.d" using 2:13 title "Particle resolved (10 mins)"
replot "out_sedi_exp_sect_num_avg.d" using 2:3 w l title "Sectional (0 mins)"
replot "out_sedi_exp_sect_num_avg.d" using 2:8 w l title "Sectional (5 mins)"
replot "out_sedi_exp_sect_num_avg.d" using 2:13 w l title "Sectional (10 mins)"
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
plot [1e-7:1e-2] [0:8e-6]  "out_sedi_exp_mc_vol_avg.d" using 2:3 title "Particle resolved (0 mins)"
replot "out_sedi_exp_mc_vol_avg.d" using 2:8 title "Particle resolved (5 mins)"
replot "out_sedi_exp_mc_vol_avg.d" using 2:13 title "Particle resolved (10 mins)"
replot "out_sedi_exp_sect_vol_avg.d" using 2:3 w l title "Sectional (0 mins)"
replot "out_sedi_exp_sect_vol_avg.d" using 2:8 w l title "Sectional (5 mins)"
replot "out_sedi_exp_sect_vol_avg.d" using 2:13 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "plot_sedi_exp_vol.eps"
replot
ENDVOL
epstopdf plot_sedi_exp_vol.eps

echo Plotting volume density
gnuplot -persist <<ENDVOLLOG
set logscale
set xlabel "radius (m)"
set ylabel "volume density (m^3/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2] [1e-20:1e3]  "out_sedi_exp_mc_vol_avg.d" using 2:3 title "Particle resolved (0 mins)"
replot "out_sedi_exp_mc_vol_avg.d" using 2:8 title "Particle resolved (5 mins)"
replot "out_sedi_exp_mc_vol_avg.d" using 2:13 title "Particle resolved (10 mins)"
replot "out_sedi_exp_sect_vol_avg.d" using 2:3 w l title "Sectional (0 mins)"
replot "out_sedi_exp_sect_vol_avg.d" using 2:8 w l title "Sectional (5 mins)"
replot "out_sedi_exp_sect_vol_avg.d" using 2:13 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "plot_sedi_exp_vol_log.eps"
replot
ENDVOLLOG
epstopdf plot_sedi_exp_vol_log.eps
