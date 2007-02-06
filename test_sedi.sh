#!/bin/sh

echo ./partbox sedi_mc.spec
./partbox sedi_mc.spec
echo ./process_out out_sedi_mc.d
./process_out out_sedi_mc.d

echo ./partbox sedi_sect.spec
./partbox sedi_sect.spec
echo ./process_out out_sedi_sect.d
./process_out out_sedi_sect.d

echo Plotting number density
gnuplot -persist <<ENDNUM
set logscale x
set xlabel "radius (m)"
set ylabel "number density (#/m^3)"
set title "Sedimentation kernel, exponential initial condition"
plot [1e-7:1e-2]  "out_sedi_mc_num_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out_sedi_mc_num_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out_sedi_mc_num_avg.d" using 2:13 title "Monte Carlo (10 mins)"
replot "out_sedi_sect_num_avg.d" using 2:3 w l title "Sectional (0 mins)"
replot "out_sedi_sect_num_avg.d" using 2:8 w l title "Sectional (5 mins)"
replot "out_sedi_sect_num_avg.d" using 2:13 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "plot_sedi_num.eps"
replot
ENDNUM
epstopdf plot_sedi_num.eps

echo Plotting volume density
gnuplot -persist <<ENDVOL
set logscale x
set xlabel "radius (m)"
set ylabel "volume density (m^3/m^3)"
set title "Sedimentation kernel, exponential initial condition"
set key left top
plot [1e-7:1e-2] [0:1.e-5]  "out_sedi_mc_vol_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out_sedi_mc_vol_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out_sedi_mc_vol_avg.d" using 2:13 title "Monte Carlo (10 mins)"
replot "out_sedi_sect_vol_avg.d" using 2:3 w l title "Sectional (0 mins)"
replot "out_sedi_sect_vol_avg.d" using 2:8 w l title "Sectional (5 mins)"
replot "out_sedi_sect_vol_avg.d" using 2:13 w l title "Sectional (10 mins)"
set terminal postscript eps
set output "plot_sedi_vol.eps"
replot
ENDVOL
epstopdf plot_sedi_vol.eps
