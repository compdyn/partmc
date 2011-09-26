# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "diameter / m"
set ylabel "BC ratio"

plot "out/average_compsizenum_0001_00000001_aero_particles.txt" using 3:($15/($15+$14)) with points
