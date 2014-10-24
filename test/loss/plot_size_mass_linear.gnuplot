# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "diameter / m"
set ylabel "mass concentration / (kg/m^3)"

set key left top

set xrange [1e-7:1e-3]
set yrange [0:5e-4]

plot "out/loss_part_0001_aero_size_mass.txt" using 1:2 title "particle t = 0 hours", \
     "out/loss_part_0001_aero_size_mass.txt" using 1:7 title "particle t = 5 minutes", \
     "out/loss_part_0001_aero_size_mass.txt" using 1:12 title "particle t = 10 minutes", \
     "out/loss_exact_aero_size_mass.txt" using 1:2 with lines title "exact t = 0 hours", \
     "out/loss_exact_aero_size_mass.txt" using 1:7 with lines title "exact t = 5 minutes", \
     "out/loss_exact_aero_size_mass.txt" using 1:12 with lines title "exact t = 10 minutes"
