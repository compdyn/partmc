# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time / s"
set ylabel "mass concentration / (kg/m^3)"

plot "out/emission_part_0001_aero_time.txt" using 1:4 title "particle initial", \
     "out/emission_part_0001_aero_time.txt" using 1:5 title "particle background", \
     "out/emission_part_0001_aero_time.txt" using 1:6 title "particle emission", \
     "out/emission_exact_aero_time.txt" using 1:4 with lines title "exact initial", \
     "out/emission_exact_aero_time.txt" using 1:5 with lines title "exact background", \
     "out/emission_exact_aero_time.txt" using 1:6 with lines title "exact emission"
