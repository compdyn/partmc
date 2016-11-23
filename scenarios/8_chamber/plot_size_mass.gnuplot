# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set title "chamber number size distribution"

set logscale x
set xlabel "diameter / nm"
set ylabel "mass concentration / (ug/m^3)"

set key right top

set xrange [10:1000]
#set yrange [1e0:1e11]

plot "out/chamber_00000001_mass_dist.txt" using ($1*1e9):($2*1e9):($3*1e9) with errorbars title "t = 0 hours", \
     "out/chamber_00000002_mass_dist.txt" using ($1*1e9):($2*1e9):($3*1e9) with errorbars title "t = 7 minutes", \
     "out/chamber_00000011_mass_dist.txt" using ($1*1e9):($2*1e9):($3*1e9) with errorbars title "t = 70 minutes", \
     "out/chamber_00000031_mass_dist.txt" using ($1*1e9):($2*1e9):($3*1e9) with errorbars title "t = 210 minutes"
