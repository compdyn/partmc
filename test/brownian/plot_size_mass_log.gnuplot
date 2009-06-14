# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale
set xlabel "radius (m)"
set ylabel "mass concentration (kg/m^3)"

set key left top

set xrange [1e-9:1e-6]
set yrange [1e-13:1e-7]

plot "out/brownian_part_size_mass.txt" using 1:2 title "MC single t = 0 hours"
replot "out/brownian_part_size_mass.txt" using 1:14 title "MC single t = 12 hours"
replot "out/brownian_part_size_mass.txt" using 1:26 title "MC single t = 24 hours"

replot "out/brownian_part_size_mass_average.txt" using 1:2 with lines title "MC average t = 0 hours"
replot "out/brownian_part_size_mass_average.txt" using 1:14 with lines title "MC average t = 12 hours"
replot "out/brownian_part_size_mass_average.txt" using 1:26 with lines title "MC average t = 24 hours"

replot "out/brownian_sect_size_mass.txt" using 1:2 with lines title "sectional t = 0 hours"
replot "out/brownian_sect_size_mass.txt" using 1:14 with lines title "sectional t = 12 hours"
replot "out/brownian_sect_size_mass.txt" using 1:26 with lines title "sectional t = 24 hours"
