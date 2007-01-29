#!/bin/sh

echo ./run_golovin_fix_hybrid
./run_golovin_fix_hybrid
echo ./process_out out_golovin_fix_hybrid.d
./process_out out_golovin_fix_hybrid.d

echo ./run_golovin_exact
./run_golovin_exact
echo ./process_out out_golovin_exact.d
./process_out out_golovin_exact.d

echo Plotting number density
gnuplot -persist <<ENDNUM
set logscale
set xlabel "radius (m)"
set ylabel "number density (#/m/m^3)"
set title "Golovin kernel, exponential initial condition"
plot [1e-7:1e-3] [1e4:1e10] "out_golovin_fix_hybrid_num_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out_golovin_fix_hybrid_num_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out_golovin_fix_hybrid_num_avg.d" using 2:13 title "Monte Carlo (10 mins)"
replot "out_golovin_exact_num_avg.d" using 2:3 w l title "Analytical (0 mins)"
replot "out_golovin_exact_num_avg.d" using 2:8 w l title "Analytical (5 mins)"
replot "out_golovin_exact_num_avg.d" using 2:13 w l title "Analytical (10 mins)"
set terminal postscript eps
set output "plot_golovin_num.eps"
replot
ENDNUM
epstopdf plot_golovin_num.eps

echo Plotting volume density
gnuplot -persist <<ENDVOL
set logscale
set xlabel "radius (m)"
set ylabel "volume density (m^3/m/m^3)"
set title "Golovin kernel, exponential initial condition"
set key left top
plot [1e-7:1e-3] [1e-14:1e-4] "out_golovin_fix_hybrid_mass_avg.d" using 2:3 title "Monte Carlo (0 mins)"
replot "out_golovin_fix_hybrid_mass_avg.d" using 2:8 title "Monte Carlo (5 mins)"
replot "out_golovin_fix_hybrid_mass_avg.d" using 2:13 title "Monte Carlo (10 mins)"
replot "out_golovin_exact_mass_avg.d" using 2:3 w l title "Analytical (0 mins)"
replot "out_golovin_exact_mass_avg.d" using 2:8 w l title "Analytical (5 mins)"
replot "out_golovin_exact_mass_avg.d" using 2:13 w l title "Analytical (10 mins)"
set terminal postscript eps
set output "plot_golovin_vol.eps"
replot
ENDVOL
epstopdf plot_golovin_vol.eps
