# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time / s"
set ylabel "mass concentration / (kg/m^3)"

plot "out/deposition_0001_aero_time.txt" using 1:4 title "particle", \
     "out/deposition_data.txt" using 1:4 with lines title "exact"
