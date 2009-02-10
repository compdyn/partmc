# run from inside gnuplot with:
# load "weight_plot.gnuplot"
# or from the commandline with:
# gnuplot weight_plot.gnuplot -

set xlabel "time (s)"
set ylabel "mass concentration (kg/m^3)"

plot "out/emission_mc_species_summary.txt" using 1:2 title "MC summary initial"
replot "out/emission_mc_species_summary.txt" using 1:3 title "MC summary background"
replot "out/emission_mc_species_summary.txt" using 1:4 title "MC summary emission"

replot "out/emission_mc_species_state.txt" using 1:2 with lines title "MC state initial"
replot "out/emission_mc_species_state.txt" using 1:3 with lines title "MC state background"
replot "out/emission_mc_species_state.txt" using 1:4 with lines title "MC state emission"

replot "out/emission_exact_species.txt" using 1:2 with lines title "exact initial"
replot "out/emission_exact_species.txt" using 1:3 with lines title "exact background"
replot "out/emission_exact_species.txt" using 1:4 with lines title "exact emission"
