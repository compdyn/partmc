# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale
set xlabel "diameter / m"
set ylabel "number concentration / (1/m^3)"

set key right top

set xrange [1e-9:1e-3]
set yrange [1e0:1e11]

plot "out/loss_part_drydep_0001_aero_size_num.txt" using 1:2 title "particle t = 0 hours", \
     "out/loss_part_drydep_0001_aero_size_num.txt" using 1:26 title "particle t = 24 hours", \
     "out/loss_part_drydep_0001_aero_size_num.txt" using 1:50 title "particle t = 48 hours", \
     "out/loss_exact_drydep_aero_size_num.txt" using 1:2 with lines title "exact t = 0 hours", \
     "out/loss_exact_drydep_aero_size_num.txt" using 1:26 with lines title "exact t = 24 hours", \
     "out/loss_exact_drydep_aero_size_num.txt" using 1:50 with lines title "exact t = 48 hours"
