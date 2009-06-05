# run from inside gnuplot with:
# load "weight_plot.gnuplot"
# or from the commandline with:
# gnuplot weight_plot.gnuplot -

set logscale y
set xlabel "time (s)"
set ylabel "number concentration (#/m^3)"

#set xrange [1e-9:1e-6]
#set yrange [1e7:1e11]

plot "out/bidisperse_ode_data.txt" using 1:2 with lines title "ODE benchmark"
replot "out/bidisperse_mc_data.txt" using 1:2 with points title "particle method"
