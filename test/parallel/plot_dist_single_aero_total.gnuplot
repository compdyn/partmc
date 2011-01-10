# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time / h"
set ylabel "number concentration / (1/m^3)"
set y2label "mass concentration / (kg/m^3)"

set ytics nomirror
set y2tics

set yrange [0:*]
set y2range [0:*]

set key right center

plot "out/sect_aero_total.txt" using ($1/3600):2 axes x1y1 with lines title "num sect", \
     "out/parallel_dist_single_aero_total.txt" using ($1/3600):2 axes x1y1 with lines title "num parallel", \
     "out/serial_aero_total.txt" using ($1/3600):2 axes x1y1 with lines title "num serial", \
     "out/sect_aero_total.txt" using ($1/3600):3 axes x1y2 with lines title "mass sect", \
     "out/parallel_dist_single_aero_total.txt" using ($1/3600):3 axes x1y2 with lines title "mass parallel", \
     "out/serial_aero_total.txt" using ($1/3600):3 axes x1y2 with lines title "mass serial"
