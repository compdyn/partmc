# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "radius (m)"
set ylabel "BC ratio"

plot "out/average_sizenum_particles.txt" using 3:($15/($15+$14)) with points
