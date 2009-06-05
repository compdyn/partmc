# run from inside gnuplot with:
# load "weight_plot.gnuplot"
# or from the commandline with:
# gnuplot weight_plot.gnuplot -

set logscale y
set xlabel "time (s)"
set ylabel "mass concentration (kg/m^3)"

plot "out/bidisperse_ode_data.txt" using 1:3 with lines title "ODE benchmark"
replot "out/bidisperse_mc_data.txt" using 1:3 with points title "particle method"
