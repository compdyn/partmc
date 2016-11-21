# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale y
set logscale y2
set xlabel "time / min"
set ylabel "number concentration / (1/m^3)"
set y2label "mass concentration / (kg/m^3)"

set ytics nomirror
set y2tics

set key right top

set xrange [0:10]
set yrange [1e9:1e12]
set y2range [1e-5:1e-2]

plot "out/loss_part_volume_0001_aero_time.txt" using ($1/60):2 axes x1y1 title "particle number", \
     "out/loss_part_volume_0001_aero_time.txt" using ($1/60):3 axes x1y2 title "particle mass", \
     "out/loss_exact_volume_aero_time.txt" using ($1/60):2 with lines axes x1y1 title "exact number", \
     "out/loss_exact_volume_aero_time.txt" using ($1/60):3 with lines axes x1y2 title "exact mass"
