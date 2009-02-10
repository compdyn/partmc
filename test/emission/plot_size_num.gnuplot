# run from inside gnuplot with:
# load "weight_plot.gnuplot"
# or from the commandline with:
# gnuplot weight_plot.gnuplot -

set logscale x
set xlabel "radius (m)"
set ylabel "number concentration (#/m^3)"
set xrange [1e-5:1e-4]

plot "out/emission_mc_size_num_summary.txt" using 1:2 title "MC summary t = 0 hours"
replot "out/emission_mc_size_num_summary.txt" using 1:38 title "MC summary t = 12 hours"
replot "out/emission_mc_size_num_summary.txt" using 1:74 title "MC summary t = 24 hours"

replot "out/emission_mc_size_num_state.txt" using 1:2 with lines title "MC state t = 0 hours"
replot "out/emission_mc_size_num_state.txt" using 1:38 with lines title "MC state t = 12 hours"
replot "out/emission_mc_size_num_state.txt" using 1:74 with lines title "MC state t = 24 hours"

replot "out/emission_exact_size_num.txt" using 1:2 with lines title "exact t = 0 hours"
replot "out/emission_exact_size_num.txt" using 1:38 with lines title "exact t = 12 hours"
replot "out/emission_exact_size_num.txt" using 1:74 with lines title "exact t = 24 hours"
