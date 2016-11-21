# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set title "chamber size num"

set logscale
set xlabel "diameter / m"
set ylabel "number concentration / (1/m^3)"

set key right top

set xrange [1e-9:1e-3]
set yrange [1e0:1e11]

plot "out/loss_part_chamber_0001_aero_size_num.txt" using 1:2 title "particle t = 0 hours", \
     "out/loss_part_chamber_0001_aero_size_num.txt" using 1:8 title "particle t = 6 hours", \
     "out/loss_part_chamber_0001_aero_size_num.txt" using 1:14 title "particle t = 12 hours", \
     "out/loss_exact_chamber_aero_size_num.txt" using 1:2 with lines title "exact t = 0 hours", \
     "out/loss_exact_chamber_aero_size_num.txt" using 1:8 with lines title "exact t = 6 hours", \
     "out/loss_exact_chamber_aero_size_num.txt" using 1:14 with lines title "exact t = 12 hours"
