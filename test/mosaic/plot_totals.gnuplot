# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time / s"
set ylabel "aerosol number concentration / (1/m^3)"
set y2label "aerosol mass concentration / (kg/m^3)"

set key top left

set ytics nomirror
set y2tics

set yrange [0:*]
set y2range [0:*]

plot "out/mosaic_0001_aero_time.txt" using 1:2 axes x1y1 with lines title "number", \
     "out/mosaic_0001_aero_time.txt" using 1:3 axes x1y2 with lines title "mass", \
     "ref_mosaic_0001_aero_time.txt" using 1:2 axes x1y1 with points title "ref number", \
     "ref_mosaic_0001_aero_time.txt" using 1:3 axes x1y2 with points title "ref mass"
