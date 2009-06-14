# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale
set xlabel "radius (m)"
set ylabel "number concentration (#/m^3)"

plot "out/bidisperse_part_aero_size_num.txt" using 1:2 with linespoints title "MC initial condition"
