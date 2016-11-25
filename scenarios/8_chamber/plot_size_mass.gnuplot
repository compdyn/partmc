# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set title "chamber number size distribution"

set logscale x
set xlabel "mobility diameter / nm"
set ylabel "mass concentration / (ug/m^3)"

set key right top
set grid

set xrange [10:1000]
#set yrange [1e0:1e11]

plot "out/chamber_00000001_mass_dist.txt" using ($1*1e9):($2*1e9):($3*1e9) with errorbars title "t = 0 hours", \
     "out/chamber_00000002_mass_dist.txt" using ($1*1e9):($2*1e9):($3*1e9) with errorbars title "t = 7 minutes", \
     "out/chamber_00000011_mass_dist.txt" using ($1*1e9):($2*1e9):($3*1e9) with errorbars title "t = 70 minutes", \
     "out/chamber_00000031_mass_dist.txt" using ($1*1e9):($2*1e9):($3*1e9) with errorbars title "t = 210 minutes", \
     "exp_aero_size_mass.txt" using ($1*1e9):($2*1e9) with lines title "barrel experiment 1, t = 0 min", \
     "exp_aero_size_mass.txt" using ($1*1e9):($3*1e9) with lines title "barrel experiment 1, t = 7 min", \
     "exp_aero_size_mass.txt" using ($1*1e9):($12*1e9) with lines title "barrel experiment 1, t = 70 min", \
     "exp_aero_size_mass.txt" using ($1*1e9):($32*1e9) with lines title "barrel experiment 1, t = 210 min"
