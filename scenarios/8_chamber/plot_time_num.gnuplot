# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set title "chamber total number concentration"

set xlabel "time / min"
set ylabel "number concentration / (1/cm^3)"

unset key

set xrange [0:350]
set yrange [0:1.4e5]

plot "out/chamber_tot_num_conc.txt" using ($1/60):($2/1e6):($3/1e6) with errorbars
