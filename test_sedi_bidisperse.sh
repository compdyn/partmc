#!/bin/sh

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
plot "counts_sedi_bidisperse_mc.d" using 1:2 title "Monte Carlo"
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
plot "counts_sedi_bidisperse_mc.d" using 1:2 title "Monte Carlo"
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
set key left top
plot "counts_sedi_bidisperse_mc.d" using 1:3 title "Monte Carlo"
replot "counts_sedi_bidisperse_ode.d" using 1:3 w l title "ODE"
set terminal postscript eps
set output "plot_sedi_bidisperse_vol_big.eps"
replot
ENDVOL
epstopdf plot_sedi_bidisperse_vol_big.eps
