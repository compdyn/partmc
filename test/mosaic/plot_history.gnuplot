# run from inside gnuplot with:
# load "weight_plot.gnuplot"
# or from the commandline with:
# gnuplot weight_plot.gnuplot -

set xlabel "time (s)"
set ylabel "gas concentration (ppb)"
set y2label "aerosol mass concentration (kg/m^3)"

set key top left

set ytics nomirror
set y2tics

plot "out/mosaic_gas_summary.txt" using 1:5 axes x1y1 with lines title "gas NH3 summary"
replot "out/mosaic_gas_state.txt" using 1:5 axes x1y1 with points title "gas NH3 state"

replot "out/mosaic_gas_summary.txt" using 1:7 axes x1y1 with lines title "gas NO2 summary"
replot "out/mosaic_gas_state.txt" using 1:7 axes x1y1 with points title "gas NO2 state"

replot "out/mosaic_aero_species_summary.txt" using 1:3 axes x1y2 with lines title "aerosol SO4 summary"
replot "out/mosaic_aero_species_state.txt" using 1:3 axes x1y2 with points title "aerosol SO4 state"

replot "out/mosaic_aero_species_summary.txt" using 1:4 axes x1y2 with lines title "aerosol NO3 summary"
replot "out/mosaic_aero_species_state.txt" using 1:4 axes x1y2 with points title "aerosol NO3 state"

replot "out/mosaic_aero_species_summary.txt" using 1:6 axes x1y2 with lines title "aerosol NH4 summary"
replot "out/mosaic_aero_species_state.txt" using 1:6 axes x1y2 with points title "aerosol NH4 state"
