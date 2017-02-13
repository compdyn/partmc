# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set title "chamber number size distribution"

set logscale x
set xlabel "mobility diameter / nm"
set ylabel "number concentration / (1/cm^3)"

set key right top
set grid

set xrange [10:1000]
#set yrange [1e0:1e11]

plot "out/chamber_00000001_num_dist.txt" using ($1*1e9):($2/1e6):($3/1e6) with errorbars title "PartMC simulation, t = 0 min", \
     "out/chamber_00000002_num_dist.txt" using ($1*1e9):($2/1e6):($3/1e6) with errorbars title "PartMC simulation, t = 7 min", \
     "out/chamber_00000011_num_dist.txt" using ($1*1e9):($2/1e6):($3/1e6) with errorbars title "PartMC simulation, t = 70 min", \
     "out/chamber_00000031_num_dist.txt" using ($1*1e9):($2/1e6):($3/1e6) with errorbars title "PartMC simulation, t = 210 min", \
     "exp_aero_size_num.txt" using ($1*1e9):($2/1e6) with lines title "barrel experiment 1, t = 0 min", \
     "exp_aero_size_num.txt" using ($1*1e9):($3/1e6) with lines title "barrel experiment 1, t = 7 min", \
     "exp_aero_size_num.txt" using ($1*1e9):($12/1e6) with lines title "barrel experiment 1, t = 70 min", \
     "exp_aero_size_num.txt" using ($1*1e9):($32/1e6) with lines title "barrel experiment 1, t = 210 min"
