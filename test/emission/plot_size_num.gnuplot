# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "diameter / m"
set ylabel "number concentration / (1/m^3)"
set xrange [1e-5:1e-3]

plot "out/emission_part_0001_aero_size_num.txt" using 1:2 title "particle t = 0 hours", \
     "out/emission_part_0001_aero_size_num.txt" using 1:38 title "particle t = 12 hours", \
     "out/emission_part_0001_aero_size_num.txt" using 1:74 title "particle t = 24 hours", \
     "out/emission_exact_aero_size_num.txt" using 1:2 with lines title "exact t = 0 hours", \
     "out/emission_exact_aero_size_num.txt" using 1:38 with lines title "exact t = 12 hours", \
     "out/emission_exact_aero_size_num.txt" using 1:74 with lines title "exact t = 24 hours"
