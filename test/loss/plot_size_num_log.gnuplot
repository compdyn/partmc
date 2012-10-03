# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale
set xlabel "diameter / m"
set ylabel "number concentration / (1/m^3)"

set key left top

set xrange [1e-7:1e-3]
set yrange [1e5:2e9]

plot "out/additive_part_0001_aero_size_num.txt" using 1:2 title "particle t = 0 hours", \
     "out/additive_part_0001_aero_size_num.txt" using 1:7 title "particle t = 5 minutes", \
     "out/additive_part_0001_aero_size_num.txt" using 1:12 title "particle t = 10 minutes", \
     "out/additive_exact_aero_size_num.txt" using 1:2 with lines title "exact t = 0 hours", \
     "out/additive_exact_aero_size_num.txt" using 1:7 with lines title "exact t = 5 minutes", \
     "out/additive_exact_aero_size_num.txt" using 1:12 with lines title "exact t = 10 minutes"
