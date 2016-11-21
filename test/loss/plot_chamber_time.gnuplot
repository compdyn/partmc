# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set title "chamber time"

set logscale y
set logscale y2
set xlabel "time / h"
set ylabel "number concentration / (1/m^3)"
set y2label "mass concentration / (kg/m^3)"

set ytics nomirror
set y2tics

set key right top

set xrange [0:12]
set yrange [4e10:1.1e11]
set y2range [1e-5:1e-3]

plot "out/loss_part_chamber_0001_aero_time.txt" using ($1/3600):2 axes x1y1 title "particle number", \
     "out/loss_part_chamber_0001_aero_time.txt" using ($1/3600):3 axes x1y2 title "particle mass", \
     "out/loss_exact_chamber_aero_time.txt" using ($1/3600):2 with lines axes x1y1 title "exact number", \
     "out/loss_exact_chamber_aero_time.txt" using ($1/3600):3 with lines axes x1y2 title "exact mass"
