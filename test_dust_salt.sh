#!/bin/sh

echo ./partbox test_dust_salt.spec
./partbox test_dust_salt.spec
echo ./process_out out_dust_salt.d
./process_out out_dust_salt.d

echo Plotting number density
gnuplot -persist <<ENDNUM
set logscale
set xlabel "radius (m)"
set ylabel "number density (#/m^3)
set title "Testcase dust and seasalt"
plot [1e-7:1e-3] [1e4:1e10] "out_dust_salt_num_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out_dust_salt_num_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out_dust_salt_num_avg.d" using 2:13 title "Monte Carlo (10 mins)"
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
plot [1e-7:1e-3] [1e-14:1e-4] "out_dust_salt_vol_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out_dust_salt_vol_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out_dust_salt_vol_avg.d" using 2:13 title "Monte Carlo (10 mins)"
set terminal postscript eps
set output "plot_dust_salt_vol.eps"
replot
ENDVOL
epstopdf plot_dust_salt_vol.eps
