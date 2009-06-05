# run from inside gnuplot with:
# load "weight_plot.gnuplot"
# or from the commandline with:
# gnuplot weight_plot.gnuplot -

set logscale
set xlabel "radius (m)"
set ylabel "number concentration (#/m^3)"

plot "out/bidisperse_mc_aero_size_num.txt" using 1:2 with linespoints title "MC initial condition"
