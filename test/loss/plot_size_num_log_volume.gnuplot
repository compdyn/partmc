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

plot "out/loss_part_volume_0001_aero_size_num.txt" using 1:2 title "particle t = 0 minutes", \
     "out/loss_part_volume_0001_aero_size_num.txt" using 1:7 title "particle t = 5 minutes", \
     "out/loss_part_volume_0001_aero_size_num.txt" using 1:12 title "particle t = 10 minutes", \
     "out/loss_exact_volume_aero_size_num.txt" using 1:2 with lines title "exact t = 0 minutes", \
     "out/loss_exact_volume_aero_size_num.txt" using 1:7 with lines title "exact t = 5 minutes", \
     "out/loss_exact_volume_aero_size_num.txt" using 1:12 with lines title "exact t = 10 minutes"
