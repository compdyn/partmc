# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale
set xlabel "diameter / m"
set ylabel "mass concentration / (kg/m^3)"

set key left top

set xrange [1e-9:1e-3]
set yrange [1e-15:1e-3]

plot "out/loss_part_drydep_0001_aero_size_mass.txt" using 1:2 title "particle t = 0 hours", \
     "out/loss_part_drydep_0001_aero_size_mass.txt" using 1:6 title "particle t = 24 hours", \
     "out/loss_part_drydep_0001_aero_size_mass.txt" using 1:10 title "particle t = 48 hours", \
     "out/loss_exact_drydep_aero_size_mass.txt" using 1:2 with lines title "exact t = 0 hours", \
     "out/loss_exact_drydep_aero_size_mass.txt" using 1:6 with lines title "exact t = 24 hours", \
     "out/loss_exact_drydep_aero_size_mass.txt" using 1:10 with lines title "exact t = 48 hours"
