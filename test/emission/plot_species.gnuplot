# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot <filename>.gnuplot

set xlabel "time (s)"
set ylabel "mass concentration (kg/m^3)"

plot "out/emission_mc_species.txt" using 1:2 title "MC initial"
replot "out/emission_mc_species.txt" using 1:3 title "MC background"
replot "out/emission_mc_species.txt" using 1:4 title "MC emission"

replot "out/emission_exact_species.txt" using 1:2 with lines title "exact initial"
replot "out/emission_exact_species.txt" using 1:3 with lines title "exact background"
replot "out/emission_exact_species.txt" using 1:4 with lines title "exact emission"
