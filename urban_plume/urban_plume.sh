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
set ylabel "gas concentration (ppb)"
set y2label "aerosol volume density (#/m^3)"
set y2tics
set title "URBAN_PLUME test case"
#plot "out/urban_plume_summary_aero_total.d" using (\$1/3600):2 axes x1y1 title "number dens." with lines
plot [0:8] "out/urban_plume_summary_aero_total.d" using (\$1/3600):3 axes x1y1 title "aerosol SO4 " with lines
#replot "out/urban_plume_summary_gas.d" using (\$1/3600):5 axes x1y1 title "gas NH3" with lines
#replot "out/urban_plume_summary_gas.d" using (\$1/3600):2 axes x1y1 title "gas H2SO4" with lines
#plot "out/urban_plume_summary_gas.d" using (\$1/3600):34 axes x1y1 title "gas ARO1" with lines
#plot "out/urban_plume_summary_gas.d" using (\$1/3600):12 axes x1y1 title "gas O3" with lines
#plot "out/urban_plume_summary_gas.d" using (\$1/3600):6 axes x1y1 title "gas NO" with lines
#replot "out/urban_plume_summary_gas.d" using (\$1/3600):7 axes x1y1 title "gas NO2" with lines
#replot "out/urban_plume_summary_gas.d" using (\$1/3600):24 axes x1y1 title "gas HCHO" with lines
#replot "out/urban_plume_summary_aero_total.d" using (\$1/3600):46 axes x1y2 title "aerosol NH4" with lines
#replot "out/urban_plume_summary_aero_total.d" using (\$1/3600):44 axes x1y2 title "aerosol NO3" with lines
#plot "out/urban_plume_summary_aero_total.d" using (\$1/3600):20 axes x1y1 title "aerosol OC" with lines
#replot "out/urban_plume_summary_aero_total.d" using (\$1/3600):21 axes x1y2 title "aerosol BC" with lines
#plot "out/urban_plume_summary_aero_total.d" using (\$1/3600):28 axes x1y2 title "aerosol ARO1" with lines
#replot "out/urban_plume_summary_aero_total.d" using (\$1/3600):29 axes x1y2 title "aerosol ARO2" with lines
set terminal postscript eps
set output "out/urban_plume_plot_history.eps"
replot
ENDHISTORY
epstopdf out/urban_plume_plot_history.eps
