#!/bin/sh

cat <<ENDINFO

MOSAIC Test-case
----------------

This tests the interface to the MOSAIC chemistry code.

ENDINFO
sleep 1

echo ../src/partmc mosaic.spec
pgdbg ../src/partmc mosaic.spec
echo ../src/process_summary out/mosaic_summary.d
../src/process_summary out/mosaic_summary.d

echo Plotting time history of gas and aerosol
gnuplot -persist <<ENDHISTORY
set xlabel "time (s)"
set ylabel "gas concentration (ppb)
set y2label "aerosol number density (#/m^3)
set title "MOSAIC test case"
plot "out/mosaic_summary_gas.d" using 1:7 with lines axes x1y1 title "gas NO2"
replot "out/mosaic_summary_gas.d" using 1:12 with lines axes x1y1 title "gas O3"
replot "out/mosaic_summary_gas.d" using 1:5 with lines axes x1y1 title "gas NH3"
replot "out/mosaic_summary_aero_total.d" using 1:4 with lines axes x1y2 title "aerosol SO4"
set terminal postscript eps
set output "out/mosaic_plot_history.eps"
replot
ENDHISTORY
epstopdf out/mosaic_plot_history.eps
