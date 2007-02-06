#!/bin/sh

#echo ./run_sedi_fix_hybrid_testbi
#./run_sedi_fix_hybrid_testbi
#echo ./process_out out_sedi_fix_hybrid.d
#./process_out out_sedi_fix_hybrid.d

echo ./run_sedi_ode
./run_sedi_ode
echo ./process_out out_sedi_ode.d
./process_out out_sedi_ode.d

#echo Plotting number density
#gnuplot -persist <<ENDNUM
#set logscale x
#set xlabel "radius (m)"
#set ylabel "number density (#/m^3)"
#set title "Sedimentation kernel, bidisperse initial condition"
#plot [1e-7:1e-3]  "out_sedi_fix_hybrid_num_avg.d" using 2:3 title "Monte Carlo (0 mins)"
#replot "out_sedi_fix_hybrid_num_avg.d" using 2:8 title "Monte Carlo (5 mins)"
#replot "out_sedi_fix_hybrid_num_avg.d" using 2:13 title "Monte Carlo (10 mins)"
#plot "out_sedi_ode_num_avg.d" using 2:3 w l title "Ode (0 mins)"
#replot "out_sedi_ode_num_avg.d" using 2:8 w l title "Ode (5 mins)"
#replot "out_sedi_ode_num_avg.d" using 2:13 w l title "Ode (10 mins)"
#set terminal postscript eps
#set output "plot_sedi_num.eps"
#replot
#ENDNUM
#epstopdf plot_sedi_num.eps

#echo Plotting volume density
#gnuplot -persist <<ENDVOL
#set logscale x
#set xlabel "radius (m)"
#set ylabel "volume density (m^3/m^3)"
#set title "Sedimentation kernel, bidisperse initial condition"
#set key left top
#plot [1e-7:1e-3] [0:1.e-5]  "out_sedi_fix_hybrid_mass_avg.d" using 2:3 title "Monte Carlo (0 mins)"
#replot "out_sedi_fix_hybrid_mass_avg.d" using 2:8 title "Monte Carlo (5 mins)"
#replot "out_sedi_fix_hybrid_mass_avg.d" using 2:13 title "Monte Carlo (10 mins)"
#replot "out_sedi_ode_mass_avg.d" using 2:3 w l title "Ode (0 mins)"
#replot "out_sedi_ode_mass_avg.d" using 2:8 w l title "Ode (5 mins)"
#replot "out_sedi_ode_mass_avg.d" using 2:13 w l title "Ode (10 mins)"
#set terminal postscript eps
#set output "plot_sedi_vol.eps"
#replot
#ENDVOL
#epstopdf plot_sedi_vol.eps
