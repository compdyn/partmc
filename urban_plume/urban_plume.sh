#!/bin/sh

cat <<ENDINFO

Urban Plume Test-case
---------------------

This simulates an idealized urban plume with gas and aerosol
chemistry.

ENDINFO
sleep 1

echo ../src/partmc urban_plume.spec
../src/partmc urban_plume.spec
echo ../src/process_summary out/urban_plume_summary.d
../src/process_summary out/urban_plume_summary.d

echo Plotting time history of gas and aerosol
gnuplot -persist <<ENDHISTORY
set xlabel "time (hr)"
set ylabel "gas concentration (ppb)
set y2label "aerosol volume density (#/m^3)
set y2tics
set title "URBAN_PLUME test case"
plot "out/urban_plume_summary_gas.d" using (\$1/3600):7 with lines axes x1y1 title "gas NO2"
replot "out/urban_plume_summary_gas.d" using (\$1/3600):12 with lines axes x1y1 title "gas O3"
replot "out/urban_plume_summary_gas.d" using (\$1/3600):5 with lines axes x1y1 title "gas NH3"
replot "out/urban_plume_summary_aero_total.d" using (\$1/3600):4 with lines axes x1y2 title "aerosol SO4"
replot "out/urban_plume_summary_aero_total.d" using (\$1/3600):5 with lines axes x1y2 title "aerosol NO3"
replot "out/urban_plume_summary_aero_total.d" using (\$1/3600):7 with lines axes x1y2 title "aerosol NH4"
set terminal postscript eps
set output "out/urban_plume_plot_history.eps"
replot
ENDHISTORY
epstopdf out/urban_plume_plot_history.eps
