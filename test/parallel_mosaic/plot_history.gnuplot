# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time / s"
set ylabel "gas mixing ratio / ppb"
set y2label "aerosol mass concentration / (kg/m^3)"

set key top left

set ytics nomirror
set y2tics

plot "out/serial_gas.txt" using 1:5 axes x1y1 with lines title "serial gas NH3", \
     "out/serial_gas.txt" using 1:7 axes x1y1 with lines title "serial gas NO2", \
     "out/serial_aero_species.txt" using 1:3 axes x1y2 with lines title "serial aerosol SO4", \
     "out/serial_aero_species.txt" using 1:4 axes x1y2 with lines title "serial aerosol NO3", \
     "out/serial_aero_species.txt" using 1:6 axes x1y2 with lines title "serial aerosol NH4", \
     "out/parallel_gas.txt" using 1:5 axes x1y1 with points title "parallel gas NH3", \
     "out/parallel_gas.txt" using 1:7 axes x1y1 with points title "parallel gas NO2", \
     "out/parallel_aero_species.txt" using 1:3 axes x1y2 with points title "parallel aerosol SO4", \
     "out/parallel_aero_species.txt" using 1:4 axes x1y2 with points title "parallel aerosol NO3", \
     "out/parallel_aero_species.txt" using 1:6 axes x1y2 with points title "parallel aerosol NH4"
