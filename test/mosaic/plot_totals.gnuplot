# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time (s)"
set ylabel "aerosol number concentration (#/m^3)"
set y2label "aerosol mass concentration (kg/m^3)"

set key top left

set ytics nomirror
set y2tics

set yrange [0:*]
set y2range [0:*]

plot "out/mosaic_aero_total.txt" using 1:2 axes x1y1 with lines title "number", \
     "out/mosaic_aero_total.txt" using 1:3 axes x1y2 with lines title "mass", \
     "true_aero_total.txt" using 1:2 axes x1y1 with points title "true number", \
     "true_aero_total.txt" using 1:3 axes x1y2 with points title "true mass"
