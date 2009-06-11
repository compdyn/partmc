# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot <filename>.gnuplot

set logscale x
set xlabel "radius (m)"
set ylabel "mass concentration (kg/m^3)"

set key left top

set xrange [1e-7:1e-3]
set yrange [0:8e-3]

plot "out/golovin_mc_size_mass.txt" using 1:2 title "MC single t = 0 hours"
replot "out/golovin_mc_size_mass.txt" using 1:7 title "MC single t = 5 minutes"
replot "out/golovin_mc_size_mass.txt" using 1:12 title "MC single t = 10 minutes"

replot "out/golovin_exact_size_mass.txt" using 1:2 with lines title "exact t = 0 hours"
replot "out/golovin_exact_size_mass.txt" using 1:7 with lines title "exact t = 5 minutes"
replot "out/golovin_exact_size_mass.txt" using 1:12 with lines title "exact t = 10 minutes"
