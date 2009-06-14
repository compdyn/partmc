# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot <filename>.gnuplot

set logscale x
set xlabel "radius (m)"
set ylabel "mass concentration (kg/m^3)"

set key left top

set xrange [1e-7:1e-2]
set yrange [0:1e-2]

plot "out/sedi_mc_size_mass.txt" using 1:2 title "MC single t = 0 hours"
replot "out/sedi_mc_size_mass.txt" using 1:3 title "MC single t = 5 minutes"
replot "out/sedi_mc_size_mass.txt" using 1:4 title "MC single t = 10 minutes"

replot "out/sedi_sect_size_mass.txt" using 1:2 with lines title "sectional t = 0 hours"
replot "out/sedi_sect_size_mass.txt" using 1:3 with lines title "sectional t = 5 minutes"
replot "out/sedi_sect_size_mass.txt" using 1:4 with lines title "sectional t = 10 minutes"
